# Loop Blocking Guide

How to convert a MOM6 `jki`-ordered loop nest (`do j; do k; do i`, with `k`
carrying a top-down or bottom-up recurrence that can't be parallelized) into
the tiled/blocked form used by `MOM_continuity_PPM.F90`, as **prep work
before** adding `do concurrent`/`!$omp target` (see the main `SKILL.md`).
Validated end-to-end on `MOM_isopycnal_slopes.F90`'s `calc_isoneutral_slopes`
(both its zonal/u-point and meridional/v-point slope loops), which produced
the pattern and gotchas below.

## When this applies

A loop is a candidate if:

- It's ordered `do j=... ; do k=... ; do i=...` (or `I`/`J` at u/v-points) —
  `k` must stay sequential (a recurrence), so only `j` and `i` can be
  blocked/parallelized.
- The `i`-loop body calls a `MOM_EOS` routine (e.g. `calculate_density_derivs`,
  `calculate_density_second_derivs`) once per row via its 1-D interface. This
  is the main payoff of blocking: promoting that call to the 2-D interface
  (`calculate_density_derivs_2d`-style) lets one call process a whole tile
  instead of one row, which is what makes a later GPU port worthwhile (see
  "Before porting any loop: check for MOM_EOS calls" in `SKILL.md`).

If the loop has no EOS call and no `k`-recurrence, plain `do concurrent`
collapse (main `SKILL.md`, "Choosing a loop construct") is probably enough —
this heavier blocking refactor is for the specific `jki` + per-row-EOS shape.

## The target shape

Use the block-indexing convention introduced in `marshall/cont-blocks-rebase`
(Marshall Ward's upstream rebase of this blocking work onto `dev/gfdl`) —
this supersedes any `jstart`/`istart`/`imax`-style ad hoc naming from earlier
GPU-branch-only work:

- `jsb:jeb` / `isb:ieb` — **center**-indexed block bounds in the j/i
  directions.
- `JsbB:JebB` / `IsbB:IebB` — **face**-indexed block bounds, for the
  meridional/J and zonal/I directions respectively.
- `jj`, `ii` — block-local center-index loop variables.
- `JJ`, `II` — block-local face-index loop variables. These are the *same
  Fortran variables* as `jj`/`ii` (identifiers are case-insensitive) —
  capitalization is mnemonic only, exactly like the existing `i`/`I`, `j`/`J`
  staggering convention (`code-style-guide.md`) extended down to the
  block-local index. A given subroutine declares whichever pair its point
  type actually needs, not both.
- `jje`, `iie` / `JJe`, `IIe` — block-local valid end bounds (center / face).

For a `jki` loop at u-points (center `j`, face `I` — the shape of
`MOM_isopycnal_slopes.F90`'s zonal loop):

```fortran
do jsb=js,je,njblock ; do IsbB=is,ie,niblock
  jeb = min(je, jsb+njblock-1) ; jje = jeb-jsb+1
  IebB = min(ie, IsbB+niblock-1) ; IIe = IebB-IsbB+1
  do K=nz,2,-1              ! stays sequential — this is the recurrence
    ! fill phase: do jj=1,jje ; do II=1,IIe ; I=IsbB+II-1 ; j=jsb+jj-1
    !   fill the batched-EOS input arrays (T/S/pres) for the whole tile
    ! one batched calculate_density_derivs(_2d) call over the whole tile
    ! main phase: do jj=1,jje ; do II=1,IIe ; I=IsbB+II-1 ; j=jsb+jj-1
    !   the original per-(j,K,i) body, reading the batched EOS outputs
  enddo
enddo ; enddo
```

**Style convention: chain the outer/inner block-loop headers onto one
line.** The outer (`jsb`) and inner (`IsbB`) loop openers share a line, with
the outer block's bound/count pair (`jeb`/`jje`) computed as the body's first
statement and the inner block's pair (`IebB`/`IIe`) as its second — same
`j`-then-`i` ordering as the loop nest itself — followed by a blank line and
then the loop's actual work. Close both loops the same way, `enddo ; enddo`
on one line. This matches the `do j=... ; do i=...` single-line nested-loop
idiom already used for short loops elsewhere in this codebase, extended here
to a loop nest whose body spans many lines.

At v-points (the meridional loop), the center/face pairing swaps: the outer
block becomes face-indexed (`JsbB:JebB`, block-local `JJ`/`JJe`) and the
inner block becomes center-indexed (`isb:ieb`, block-local `ii`/`iie`) — same
structure, different letters, per which coordinate is staggered at that
point type.

`njblock`/`niblock` are compile-time `integer, parameter`s during this pass
(CPU-only restructuring, no offload directives yet) — see "Preserving CPU
and GPU performance together" in `SKILL.md` for why per-platform block sizes
come later, once offload is actually being added. The rest of this guide
(step-by-step, gotchas) still refers to the block bounds generically as
`jstart`/`istart`/`jmax`/`imax`/`ii`/`jj` — read those as shorthand for
whichever of the `jsb`/`JsbB`-style pairs above actually applies to your
loop's point type.

## Step-by-step (bisectable commits)

Follow the commit discipline in `references/loop-blocking-shared-patterns.md`
— small, individually-verified commits, commit only once a build+run
reproduces the baseline `ocean.stats` bit-for-bit. Every commit below is
provably bit-identical by construction — block size never changes per-point
arithmetic or ordering, only how work is chunked — until offload directives
are introduced later.

1. **j-blocking skeleton.** Introduce `jstart`/`jend`/`jmax`/`jj` with
   `njblock` fixed at `1`. Convert `do j=js,je` into the `jstart` loop plus
   an inner `do jj=1,jmax ; j=jstart+jj-1`. No array shape changes yet —
   this is pure infrastructure and should be its own commit (bisect marks it
   trivially good).
2. **Promote arrays, split fill/main phases.** Give the arrays that carry
   data from the fill phase to the main phase (the EOS call's inputs and
   outputs, plus anything else read across the `jj` boundary) an added
   `njblock`-sized dimension, indexed by `jj`. Split the loop body into a
   fill pass and a main pass over the same `jj` range. Scratch arrays that
   are filled and consumed *within* the same `jj` iteration don't need
   promoting — see the shared-`dom` gotcha below before deciding which
   arrays need it.
3. **Batch the EOS call.** Replace the per-`jj` 1-D `calculate_density_derivs`
   call with a single 2-D-interface call spanning the whole `jj` block,
   using a `dom` argument shaped `(2,2)`: one row for the existing `i`-range,
   one row for `[1, jmax]`. A simpler intermediate is also valid if you want
   a smaller diff to verify at this step: keep the call on its 1-D interface,
   called once per row inside a `do j=jstart,jend` loop, but slice the array
   argument to the tile's `i`-range directly (`T_EOS(istart:iend,j)` rather
   than the full-domain `T_EOS(:,j)`) instead of passing an explicit `dom`.
   Fortran re-bases an array section to start at index 1, so this needs no
   manually-constructed domain at all — it's a legitimate stepping stone
   before committing to the full 2-D-interface-plus-`dom` batching in a
   later commit (confirmed in `MOM_set_viscosity.F90`'s `set_viscous_ML`,
   which used exactly this slicing form for several commits before
   introducing an explicit `EOSdom(2,2)`).
4. **i-blocking skeleton.** Add `istart`/`iend`/`niblock` (parameter, not yet
   promoted) alongside `jstart`, with `niblock` sized to actually amortize
   the batched EOS call's overhead — see the gotcha below, **do not** start
   this at `niblock=1` the way step 1 safely could for `j`. `MOM_continuity_PPM.F90`'s
   CPU `default_niblock = 32` is the established precedent. Keep every
   promoted array's declared shape full-width in this commit — only the
   loops narrow to `istart:iend`.
5. **Shrink arrays to the i-block.** Only now replace the full-width leading
   dimension on the arrays promoted in step 2 with `niblock`, indexed by
   `ii`. This must come strictly after step 4, never combined with or ahead
   of it — see the ordering gotcha below.

### Applying this to a second, structurally-similar loop

MOM6 subroutines often have a zonal (u-point) and meridional (v-point)
version of the same slope/flux calculation back to back. Once steps 1-5 are
validated on the first one:

- **Reuse the same subroutine-local blocking variables** (`jstart`, `istart`,
  `niblock`, `ii`, etc.) for the second loop instead of declaring a second,
  suffixed set (`jstart_v`). They're plain subroutine-scope locals and the
  two loops run sequentially — no aliasing risk, and duplicating them is
  pure churn.
- **Apply the validated design directly, without re-interviewing** — the
  design questions (two-phase split vs. single-phase, which arrays to
  promote, block size) were already settled by the first loop. Re-asking
  them for a structurally-identical second loop just re-litigates a decided
  question.
- **Still flag genuine structural differences as you find them.** They won't
  always be identical: e.g. a zonal loop's `i`-stencil may need a "+1"
  widened array (h-point values at both `I` and `I+1`) that a meridional
  loop doesn't need if its analogous stencil is already two separate arrays
  (row `j` and row `j+1`) rather than one array read at two adjacent
  indices. Mirroring the pattern is not the same as assuming the two loops
  are identical.
- **Temp/scratch data arrays consumed by the EOS call can also be shared
  between the two loops**, not just the loop-control integers — this is a
  design choice, not a requirement. `set_viscous_ML` reuses the same
  `T_EOS`/`S_EOS`/`press`/`dR_dT`/`dR_dS` arrays for both its u-point and
  v-point loops (rather than `MOM_isopycnal_slopes.F90`'s choice of separate
  `_u`-suffixed copies), which is safe as long as the two loops run
  sequentially and each tile-fill is fully consumed before the other loop's
  fill overwrites it. Sharing saves the array-duplication churn; separate
  copies make each loop's diff independently readable and avoid having to
  reason about the sharing invariant. Either is fine — pick based on whether
  the surrounding code already gives the two loops distinct array names.

If the two loops share a promoted array (rather than using separate `_u`/`_v`
copies), the shrink step (5 above) can't happen until **both** loops' `istart`
skeletons exist — shrinking after only one loop's skeleton is in place
reintroduces the out-of-bounds-write ordering bug below for whichever loop
hasn't been narrowed yet.

## Correctness gotchas

**i-blocking an EOS call needs a moderate block size, unlike j-blocking.**
If the original code already called the EOS routine once per row (`j`/`k` as
the only outer indices, full `i`-range in one call), then `njblock=1`
exactly reproduces that granularity — safe. But moving the same call inside
an `istart` loop, even at `niblock=1`, turns "one call per row" into "one
call per point": call count scales with row width while useful work per call
shrinks by the same factor. Each `calculate_density_derivs_2d`-style call
has real fixed overhead (polymorphic EOS dispatch, plus automatic-array
temporaries sized to the full array passed in, not the `dom` sub-range) —
enough that `niblock=1` turned a ~100s `benchmark_ALE` run into one that
didn't finish in 20+ minutes. Use a block size like `niblock=32` instead;
it stays bit-identical (block size doesn't change arithmetic) and amortizes
the per-call overhead.

**Every array passed to one EOS call must share the same index base.** If
one argument to a `calculate_density_derivs`/`calculate_density_second_derivs`
call is block-relative (promoted to `niblock`-sized, indexed from 1) while
another argument to the *same call* is still full-domain (un-promoted,
globally indexed), the shared `dom` argument indexes them inconsistently
once the block advances past the first tile — their "index 1" no longer
refers to the same global point. Either promote every argument to a given
EOS call consistently, or leave all of them un-promoted; don't mix.

**Don't shrink an array's declared shape before the loop that narrows its
fill range actually exists (step 5 before step 4).** If an array is declared
`(niblock, njblock)` but its fill loop still runs the original, unnarrowed
`do i=is,ie` range because the `istart`/`iend` skeleton commit hasn't landed
yet, every fill beyond index `niblock` is an out-of-bounds write. This
doesn't necessarily crash — it can manifest as the model hanging near
timestep 0 with climbing CPU time, which reads like a performance problem,
not memory corruption. Keep arrays at their full-width type through the
skeleton commit; shrink only in the next commit once every loop touching
that array is already bounded to `istart:iend`.

**Once a boolean mask array (e.g. `do_i`) is shrunk to tile-relative
indices, compute `ii`/`jj` before testing it, not inline in the loop
header.** This is the stale-index-after-promoting-an-array failure mode
described in `references/loop-blocking-shared-patterns.md`, wearing a
mask-specific costume: `do I=istart,iend ; if (do_i(I,j)) then`
works before the shrink (mask is still indexed globally); after shrinking
`do_i` to `(niblock,njblock)` the same line needs `ii`/`jj` computed first —
`do I=istart,iend ; ii=I-istart+1 ; jj=j-jstart+1 ; if (do_i(ii,jj)) then` —
since the mask can no longer be read with the loop's raw index.

## Distinguishing a real hang from transient slowness

See `references/loop-blocking-shared-patterns.md` — the CPU-time thresholds
and diagnosis approach there apply here unchanged. Both
real hangs found in this loop (the EOS-call-count blowup and the
out-of-bounds write above) showed the same 10+-minutes-still-climbing
signature described there.

## After blocking: adding offload directives

The steps above are deliberately scoped to pure CPU-side restructuring —
adding `do concurrent`/`!$omp target` to the now-blocked tile loops is a
separate, later pass (see `SKILL.md`, "Choosing a loop construct"). Once you
get there, on this specific `jki` shape a few more things come up:

**Wrap every `local()`/`reduce()` locality clause in `DO_LOCALITY()`, not
raw** — see `SKILL.md`'s `do concurrent` notes for why. `reduce()` is what
handles the "did any column still need work" style flag that a sequential
`jki` loop would otherwise set with a plain scalar inside the loop body.

**A `real, dimension(niblock,njblock) :: ...` local, sized from dummy
arguments, is not the same "automatic array" the `do concurrent` compiler-bug
caveat in `SKILL.md` warns about.** That caveat is specifically about arrays
privatized *inside* a `local()` clause (one instance per iteration); a
subroutine-scope local merely *sized* by a dummy argument and reused across
the whole tile loop doesn't hit it. `set_viscous_ML`'s tile-sized temp arrays
are declared this way and offload cleanly.

**Reduce kernel-launch count when a sequential `k`-loop's body splits into
k-invariant stages.** If a `do k=1,nz` loop's body only differs by whether
`k` is above/below some threshold (e.g. `k<=nkml` vs. `k>nkml`), and an
expensive call (like the batched EOS call) only fires at one specific `k`,
opening/closing a `do concurrent`/`!$omp target` region *inside* every `k`
iteration launches `nz` times (times however many separate regions the body
has) when 3 launches would do. Restructure into three stages instead: one
`!$omp target ... !$omp end target` region wrapping `do k=1,min(nkml,nz)`
(with `!$omp loop collapse(2)` per-`k` work inside, not a fresh `target` per
`k`), the single EOS call/branch on its own outside any `k`-loop, then a
second `!$omp target` region wrapping `do k=nkml+1,nz`. `set_viscous_ML`'s
"refactor inner k loops" commit did exactly this and cut kernel launches from
`nz*3` down to 3. There's real code duplication in splitting a loop body this
way — accept it; it's the point.

**The same restructuring applies to a `domore`/`exit`-style recurrence, even
without a fixed `k`-threshold to split on.** Wrap the *entire*
`do k=nz,1,-1 ... exit` loop in one `!$omp target` region (see SKILL.md's
"Choosing a loop construct", option 4), with `!$omp loop collapse(n)` doing
the per-`k` masked work inside, instead of a `do concurrent` kernel launched
fresh every `k`. Measured directly on one such recurrence with a
data-transfer-tracing log (see SKILL.md's "Verifying a port"): moving from
one `do concurrent` per `k` to one `!$omp target` per whole recurrence cut
the subroutine's total transfer/launch events by roughly 6x (~8300 to ~1300
in one measured case), with the recurrence's own launch count dropping from
one-per-`k`-per-tile down to one-per-call.

**A sequential CPU-only early-exit (`if (.not.do_any) exit`) needs the
`#ifndef __NVCOMPILER_OPENMP_GPU` guard only once the `k`-loop's body lives
inside a persistent `!$omp target` region shared across every `k`, as in the
restructuring just above** — not simply because offload directives exist
nearby. A `do concurrent` kernel launched fresh each `k`, with the sequential
`k`-loop itself remaining ordinary host Fortran in between launches, can
check `if (.not.do_any) exit` exactly as on CPU: the exit runs as genuine
host code between kernel launches, so it needs no guard for correctness
there. Once the whole recurrence moves inside one target region, though, the
early-exit interacts with the region's control flow differently, and
`set_viscous_ML` compiles out the tracking flag and the `exit` for GPU
builds with `#ifndef __NVCOMPILER_OPENMP_GPU` / `#endif`, always running the
full `k` range instead — safe as long as the mask (`do_i`/`do_any`) still
gates each column's work, making the extra iterations no-ops.
Preprocessor-gating a CPU-only sequential optimization like this is the
standard way to keep both paths correct when their control flow has to
genuinely differ, not just their tile size. This guard is also worth adding
purely for *performance*, even in the `do concurrent`-per-`k` form where
it's not required for correctness: checking the flag every `k` means a
host↔device round trip on that scalar every iteration, which can dominate a
subroutine's transfer count (measured as ~94% of all transfer events in one
such subroutine before this guard was added).

**A single intrinsic call (`exp`, and likely other transcendentals) is not
guaranteed bitwise-identical between CPU and GPU** — see
`code-style-guide.md`'s arithmetic-reproducibility notes. `set_viscous_ML`
hit this and worked around it by precomputing the affected value on the host
before entering the offloaded region rather than calling it inside the
kernel.

## Making the block size runtime-configurable

See `references/loop-blocking-shared-patterns.md` for the general pattern
(`CS` zero-sentinel field, resolver, `intent(in)` argument,
`layoutParam=.true.`). `MOM_continuity_PPM.F90` and `MOM_set_viscosity.F90`'s
`set_viscous_ML` are both settled examples of it applied to
`niblock`/`njblock`; the resolver in the latter is
`viscous_ML_block_sizes(CS, G, niblock, njblock)`, substituting the full
local-domain extent (`G%iedB-G%isdB+1` / `G%jedB-G%jsdB+1`) when either
`CS%niblock`/`CS%njblock` is `0`.

## Worked examples

- `MOM_isopycnal_slopes.F90`'s `calc_isoneutral_slopes` (zonal then
  meridional slope loops): blocked in 10 commits (skeleton → promote+split →
  batch EOS, ×2, one per loop, then skeleton → shrink for the i-dimension,
  ×2, one per loop), each verified with a bit-for-bit `ocean.stats` match
  before committing. Stopped short of adding offload directives.
- `MOM_set_viscosity.F90`'s `set_viscous_ML` (branch
  `port-set_viscous_ML-tile`, commits `42094b2a8..17b3cc671`): the same
  `jki` blocking pattern carried all the way through offload directives,
  data mapping, the kernel-launch-count reduction, and the runtime-configurable
  block size described above. Use `git log`/`git show` on the individual
  commits in that range for concrete diffs of any of the techniques on this
  page.
