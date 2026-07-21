# k-blocking a `kji`-ordered loop (no k-recurrence)

This is the counterpart to `references/loop-blocking-guide.md` (the "jki"
pattern), for the opposite structural case. Read that doc first if you
haven't — the two patterns share vocabulary (`do concurrent`, `DO_LOCALITY`,
data mapping) and several correctness/process patterns (commit discipline,
runtime-configurable block sizes, diagnosing a hang) covered once in
`references/loop-blocking-shared-patterns.md`, which this doc assumes rather
than re-explains.

**jki** (loop-blocking-guide.md): `do j; do k; do i`, and `k` carries a real
top-down/bottom-up recurrence (e.g. an isopycnal-slope sweep). `k` must stay
strictly sequential; only `j`/`i` can be blocked or parallelized.

**kji** (this doc): `do k; do j; do i`, and each layer's calculation is
*independent of every other layer* — nothing carries over from `k` to `k+1`.
Because there's no recurrence, `k` itself is a legitimate axis to block *and*
collapse into the parallel index set alongside `j`/`i`, not just something to
loop over sequentially outside a `j`/`i` `do concurrent`. This matters
because a naive per-layer port (`do k=1,nz` sequential on the host, with only
the inner `j`/`i` loop as `do concurrent`) launches one (or several) small
kernels per layer — `nz` times whatever the loop body needs. Grouping
`nkblock` layers together and collapsing `k` into the same `do concurrent`
as `j`/`i` turns that into far fewer, larger kernel launches.

Validated on `MOM_hor_visc.F90`'s `horizontal_viscosity` (branch
`kblock-horvisc-fromgfdl`), which computes the horizontal viscosity
tendencies layer-by-layer with no vertical dependency.

## When this applies

- The loop is `do k=1,nz` (or an equivalent per-layer loop) with no
  reference to `k-1`/`k+1` state carried across iterations — check this
  carefully; a loop that merely *looks* like `kji` but reads/writes
  `array(...,k-1)` or accumulates something from a previous `k` has a real
  recurrence and is a `jki` case instead (see that guide).
- The loop body doesn't call a MOM_EOS routine in a way that would apply
  (see "Before porting any loop: check for MOM_EOS calls" in `SKILL.md`) —
  if it does, resolve that first, independent of k-blocking.

## The target shape

`marshall/cont-blocks-rebase` (Marshall Ward's upstream rebase of this
blocking work onto `dev/gfdl`) names the k-block bounds `k_start`/`k_end`
(not `ks`/`ke`) — use that naming, since it's the convention the rest of
this blocking work is converging on. That branch's `PPM_reconstruction_x`/
`_y` don't yet carry `do concurrent` (it's the CPU-only blocked form — see
the worked example below); adding the offload collapse on top of the same
bounds looks like:

```fortran
do k_start=1,nz,nkblock
  k_end = min(nz, k_start+nkblock-1)
  do concurrent (k=k_start:k_end, j=js:je, i=is:ie) DO_LOCALITY(local(kk))
    kk = k - k_start + 1
    ! original per-k body: arrays already indexed by the real k (h_in(i,j,k),
    ! etc.) need no change; only the new block-relative scratch arrays use kk
  enddo
enddo
```

The `do concurrent` header carries the **real, absolute `k`**, and the
block-relative `kk` is derived inside the loop body only where it's needed —
to index a scratch array sized `(...,nkblock)` rather than `(...,nz)`. This
reads more naturally than the reverse (`kk` in the header, `k =
k_start+kk-1` derived inside) whenever most of the loop body's arrays are
already indexed by the real `k` anyway, which is the common case here.

Contrast with the `jki` shape: there, `k` remains a literal sequential `do`
wrapping `j`/`i` tiles. Here, `k` sits *inside* the same `do concurrent` as
`j`/`i` — the whole point is that the GPU can parallelize across layers too.

The rest of this doc (step-by-step, gotchas) still refers to the k-block
bounds as `ks`/`ke`, matching `MOM_hor_visc.F90`'s current form — read those
as shorthand for `k_start`/`k_end` above.

## Step-by-step (bisectable commits)

Follow the commit discipline in `references/loop-blocking-shared-patterns.md`.
Don't take a real branch's commit boundaries
as license to go coarser than this — a prior k-blocking branch used a
coarser split, and that's not something to emulate for granularity, only for
the technical lessons below.

1. **Add the `nkblock` trailing dimension to every local scratch array that
   needs to carry a value across the loop body**, at compile-time
   `nkblock=1`, indexed with a literal `1` everywhere (bit-identical infra
   commit — `kk` doesn't exist yet, every reference is still to index `1`).
   Arrays already indexed by the real `k` (not scratch, e.g. subroutine
   arguments) don't need touching. If the loop already has some pre-existing
   `do concurrent (j=..., i=...)` sub-loops from an earlier, simpler per-k
   port, it can be easier to first convert them back to plain nested
   `do`/`do` for this commit, so the array-shape diff isn't tangled with
   concurrency-syntax changes — reintroduce `do concurrent` in a later
   commit once the shapes are settled.
2. **Introduce the `ks`/`ke` block loop**, `nkblock` still `1`, `kk` always
   `1`, `kk = k - ks + 1`. Still bit-identical.
3. **Collapse `k` into the `j`/`i` loops as a single `do concurrent`**,
   *except* where it can't be:
   - any loop body that calls a non-`pure` subroutine — `do concurrent`
     cannot contain one, so keep that loop an ordinary sequential `do`.
   - any loop that accumulates a running value *across* `k` within the block
     (not just within `j`/`i`) — collapsing `k` into `do concurrent` races on
     the accumulator unless it can be expressed as a
     `DO_LOCALITY(reduce(...))` clause. If it can't, keep `k` sequential for
     that loop.
4. **Add tailored `DO_LOCALITY` clauses** per loop — `local()` for any
   scalar (like `k` itself, or a temporary) that's assigned inside the
   collapsed loop, `reduce()` for any scalar accumulated across iterations.
   Don't blanket-copy one loop's clause list onto another; the compiler will
   tell you what's missing (`must appear in a SHARED or PRIVATE clause`
   equivalent for `do concurrent` is a locality-clause omission).
5. **Add data-mapping directives** for the newly `nkblock`-shaped locals
   (`!$omp target enter/exit data map(alloc: ...)`), per
   `references/data-mapping-conventions.md`.
6. **Promote `nkblock` to a runtime CS parameter** — see below.
7. **Test with `nkblock` set to something greater than 1** — not just left
   at the compile-time-equivalent default. This step is not optional; see
   the first gotcha below for why testing only ever at `nkblock=1` is not
   sufficient to call the port correct.

## Correctness gotchas

**A hardcoded index `1` left behind after promoting an array to carry the
`nkblock` dimension is invisible at `nkblock=1` and wrong for `nkblock>1`.**
This is the stale-index failure mode described in
`references/loop-blocking-shared-patterns.md` — in `MOM_hor_visc.F90` it
showed up twice independently: an anisotropic-viscosity correction that only updated
`str_xx(:,:,1)` instead of looping `kk=1,kmax` (silently losing the
contribution from every layer but the block's first once `HORVISC_NKBLOCK
> 1`), and a Leith+E initialization (`m_leithy(:,:,1) = 0.0`) that should
have zeroed the whole block.

**`do concurrent` can't contain a call to a non-`pure` procedure, and can't
safely collapse a dimension carrying a genuine cross-iteration dependency**
— see `SKILL.md`'s `do concurrent` notes. `MOM_hor_visc.F90`'s MEKE-reduction
loops kept a serial `kk` step for exactly this reason.

**A one-time initialization guarded by `if (k==1) then ... endif` inside a
loop you're about to k-block is a race, not a safe "runs once" idiom** — see
`references/loop-blocking-shared-patterns.md`. `MOM_hor_visc.F90`
zeroed `MEKE%mom_src`/`mom_src_bh`/`GME_snk` this way; the fix hoisted the
zeroing entirely out of the k-loop into a plain sequential block that runs
once before the block loop starts.

**Guard conditions can get silently dropped when code is relocated during
the refactor.** Moving an accumulation to sit next to a different
conditional block is easy to get subtly wrong (dropping a guard that used to
apply). Diff carefully against the pre-refactor logic when relocating a
conditional accumulation, not just against what compiles.

## Making the block size runtime-configurable

Same pattern as `references/loop-blocking-shared-patterns.md`, applied to
`nkblock`: a public field on the `*_CS` type defaulting
to `0` ("not configured — use the full column as one block"), a resolver
that substitutes `GV%ke` when `CS%nkblock==0`, and the resolved value passed
explicitly into the ported subroutine as an `intent(in)` argument. Use
`get_param(..., layoutParam=.true.)` for the namelist override (e.g.
`HORVISC_NKBLOCK`). The platform-conditional compile-time default is the
same idea as before: `0` (⇒ full column) under `__NVCOMPILER_OPENMP_GPU`, a
small literal (e.g. `1`) otherwise.

## Worked examples

- **`MOM_continuity_PPM.F90`'s `PPM_reconstruction_x`/`PPM_reconstruction_y`**
  (already in this file — no branch needed): the clean, settled final state
  of this exact pattern, worth reading before diving into a WIP branch. A
  few things it does slightly differently from the step-by-step above, worth
  matching:
  - The `do concurrent` header uses the **real, absolute `k`** as its index
    (`do concurrent (k=ks:ke, j=jsl:jel, i=isl:iel) ...`), not a
    block-relative `kk`. The block-relative index is derived *inside* the
    loop body only where it's actually needed, to index the block-sized
    scratch array: `kk = k - ks + 1`. This reads more naturally than
    introducing `kk` in the `do concurrent` header and computing
    `k = kstart+kk-1` — prefer this form (`k` in the header, `kk` derived
    inside) when the loop body's non-scratch arrays are naturally indexed by
    the real `k` anyway.
  - `nkblock` is resolved once at the top of the *caller* with a plain `if`
    (`nkblock = CS%nkblock ; if (nkblock == 0) nkblock = nz`) and passed
    into the subroutine as an ordinary `intent(in)` argument — not a `merge()`
    inlined into a declaration. This is the form to follow.
  - The block-relative scratch array (`slp`) is declared `max(1,nkblock)`-sized
    and mapped with a single `!$omp target enter data`/`exit data` pair
    bracketing the *whole subroutine* (outside the `do ks=1,nz,nkblock`
    loop), not re-mapped per block.
  - A helper subroutine that does the final limiting (`PPM_limit_CW84`/
    `PPM_limit_pos`) is called once per k-block, taking the block's `ks`/`ke`
    as plain integer arguments — a simple alternative to constructing a
    `dom`-style tuple when the callee just needs a sub-range, not a
    re-based index.
- **`MOM_hor_visc.F90`'s `horizontal_viscosity`** (branch
  `kblock-horvisc-fromgfdl`, merge-base `b8c471cfa`): the in-progress branch
  the gotchas above came from — useful for seeing the mistakes and their
  fixes, not just the destination. The branch also contains several commits
  that extract sub-blocks of the subroutine into separate helper subroutines
  (`hor_visc_Leithy_Ah`, `hor_visc_GME_setup`, and others) — those are pure
  code-organization moves, not part of the k-blocking technique, and aren't
  the source of anything in this doc. Use `git log`/`git show` on the
  k-blocking, do-concurrent-conversion, race-fix, and runtime-parameter
  commits in that range for concrete diffs of the techniques above.
