# Loop blocking (jki): offload, runtime config, worked examples

Continuation of `references/loop-blocking-guide.md` — read that first for
the target shape and the step-by-step blocking recipe. Covers what comes
after a `jki` loop is blocked: adding offload directives, making the block
size runtime-configurable, and the worked examples the whole `jki` recipe
(this file plus `loop-blocking-guide.md` and `loop-blocking-gotchas.md`)
was validated against.

## After blocking: adding offload directives

The blocking recipe in `references/loop-blocking-guide.md` is deliberately
scoped to pure CPU-side restructuring — adding `do concurrent`/`!$omp
target` to the now-blocked tile loops is a separate, later pass (see
`SKILL.md`, "Choosing a loop construct"). Once you get there, on this
specific `jki` shape a few more things come up:

**Wrap every `local()`/`reduce()` locality clause in `DO_LOCALITY()`, not
raw** — see `references/loop-constructs.md`'s `do concurrent` notes for
why. `reduce()` is what handles the "did any column still need work" style
flag that a sequential `jki` loop would otherwise set with a plain scalar
inside the loop body.

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
