---
name: mom6-gpu-loop-porting
description: Use this skill whenever porting MOM6 (Modular Ocean Model 6) Fortran code to run on GPU, or reviewing/debugging existing GPU-offloaded MOM6 code. This covers converting loops to `do concurrent` or OpenMP `!$omp target` offload directives, deciding where `!$omp target enter/exit data` and `target update` directives belong, handling MOM6 control-structure (`*_CS`) fields and derived types like `tv`/`G`/`US`/`GV` on the device, working around the MOM_EOS polymorphic-type GPU limitation, blocking/tiling a `jki`-ordered loop (`njblock`/`niblock`, à la `MOM_continuity_PPM.F90`) or k-blocking a `kji`-ordered loop with no k-recurrence (`nkblock`, `k_start`/`k_end`, collapsing `k` into `do concurrent`, à la `MOM_hor_visc.F90`) as prep before offload, and preserving CPU performance while adding GPU support (e.g. loop tiling / "blocks"). Trigger this for tasks like "port this subroutine to GPU", "add OpenMP offload to this loop", "block/tile this jki loop", "k-block this loop", "port this kji loop to GPU", "why does this only break when NKBLOCK>1", "why is this GPU kernel giving different answers than the CPU", "how should I map this CS field to the device", or any work touching `!$omp target`, `do concurrent`, loop blocking, or GPU offload anywhere in the MOM6 codebase — even if the user doesn't use the word "skill" or spell out these directive names explicitly.
---

# MOM6 GPU Loop Porting

## Status and scope

MOM6's GPU port is a work in progress, validated only on the`double_gyre`,
`benchmark` and `benchmark_ALE` configurations — the codebase is large, and
only the paths those configurations exercise are ported so far. Treat
everything below as current best practice, not a finished spec: for a pattern
not covered here, make a reasonable call consistent with this philosophy, flag
it to the user, and suggest adding it once confirmed to work.

Overall approach: compute on GPU via `do concurrent` or OpenMP `target`
loop constructs, manage data movement explicitly (never rely on
compiler-inferred transfers), and port only the branches the target
configurations exercise.

## Before porting any loop: check for MOM_EOS calls

GPUs can't execute a polymorphic Fortran type's deferred procedures from
inside a kernel — `MOM_EOS` subroutines are typically dispatched this way
(an abstract `EOS` type selecting the equation-of-state routine at
runtime). This doesn't make MOM_EOS-adjacent code unportable: a kernel can
launch inside an EOS subroutine on data already resident on device, as
long as *that subroutine's own kernels* never touch the polymorphic `EOS`
type internally. The constraint is the dispatch happening inside a kernel,
not EOS logic being off-limits.

Before converting a loop:

1. Check whether the loop body calls a MOM_EOS subroutine.
2. If it does, check whether it separates cleanly from the rest of the
   loop (hoisted out, or called ahead into a temporary array) without
   substantially restructuring the surrounding code.
3. If separation is easy, do it, and port the surrounding loop normally.
4. If separation isn't easy — a real refactor is needed — **stop**. Don't
   attempt the port. Explain to the user why the EOS call blocks a
   straightforward port, sketch 1-2 refactor approaches, and ask how
   they'd like to proceed before writing offload code. Guessing risks a
   kernel that won't compile, crashes on device, or an unwanted refactor.

## Choosing a loop construct

Try options 1-3 in order. Each fallback exists to solve a specific, concrete
problem with the option above it, not out of general preference — don't
skip to a fallback unless its trigger actually applies. Option 4 addresses a
different axis entirely (kernel-launch overhead across a sequential loop,
not which construct can execute a given parallel loop) and can apply
regardless of which of 1-3 you're using for the parallel work inside it.

### 1. `do concurrent` (default choice)

`do concurrent` is preferred by default: it's standard Fortran, so code
stays portable and readable off-GPU too, and it currently gets the best
`nvfortran` code generation (the compiler this port targets).

```fortran
do concurrent( j=js:je, i=is:ie ) DO_LOCALITY(local(tmp))
  tmp = some_expression(i, j)
  output(i,j) = tmp * other_array(i,j)
enddo
```

Notes on `do concurrent` syntax, since it differs from an ordinary `do` loop
in easy-to-miss ways:

- Loop bounds are separated by a **colon**, not a comma:
  `do concurrent( i = is:ie, j = js:je )` — not `do i=is,ie; do j=js,je`.
- Collapse multi-dimensional loops into a single `do concurrent` statement
  instead of nesting separate `do concurrent` loops per dimension. Nesting
  one `do concurrent` inside another is fine (e.g. an outer 2-D horizontal
  loop with an inner k-loop that itself needs to be a `do concurrent`) —
  avoid only splitting one multi-index statement into several nested ones.
- Declare loop-private scalar temporaries with `local(tmp)` so GPU
  parallelization doesn't race on them; accumulate a scalar reduction
  (e.g. a logical-or "did any column still need work" flag) with
  `reduce(.or.:some_flag)` the same way.
- **Always wrap the locality clause in the `DO_LOCALITY(...)` macro**
  (`#include`d from `src/framework/do_concurrent_compat.h`), never bare:
  `DO_LOCALITY(local(tmp))`, `DO_LOCALITY(reduce(.or.:do_any))`.
  `MOM_continuity_PPM.F90` — this pattern's reference file — never writes a
  bare clause. The macro expands to the clause only when
  `HAVE_FC_DO_CONCURRENT_LOCAL` is defined, else to nothing, so the same
  source stays portable to compilers without locality-specifier support;
  writing the clause directly silently drops that fallback.
- **`DO_LOCALITY` takes exactly one locality-clause argument** — never
  combine `local(...)` and `reduce(...)` into a single call, e.g.
  `DO_LOCALITY(local(a,b), reduce(.or.:c))`: the comma inside one macro
  invocation breaks the C preprocessor's argument count, since
  `#define DO_LOCALITY(X) X` only takes one parameter. When a loop needs
  both, chain two separate `DO_LOCALITY(...)` calls with `&` continuation
  instead:
  ```fortran
  do concurrent (j=js:je, i=is:ie, do_i(i,j)) &
        DO_LOCALITY(local(tmp1,tmp2)) &
        DO_LOCALITY(reduce(.or.:do_any))
  ```
  matching precedent in `MOM_set_viscosity.F90` and `MOM_coms.F90`, which
  both chain separate calls rather than combining clauses in one.
- Remove any pre-existing `!$OMP parallel do` directive when converting a
  loop to `do concurrent` — the two shouldn't coexist on the same loop.
- `do concurrent` supports a mask condition, which is a natural fit for the
  common MOM6 pattern of `do i; do j; if (mask(i,j) > 0) ...`:

  ```fortran
  do concurrent( j=js:je, i=is:ie, mask(i,j) > 0 )
    output(i,j) = some_expression(i,j)
  enddo
  ```

  Prefer this masked form over an `if` as the loop body's first statement
  when the condition simply gates whether the iteration runs at all.

**The one case `do concurrent` can't be used**: a loop needing to
privatize an *automatic array* per thread (an array-valued temporary sized
from loop or subroutine bounds, unlike the scalar temporaries `local(...)`
is meant for). As of `nvfortran` 26.5 there's a known compiler bug here —
don't fight it, fall through to the next construct. Since it's a live bug,
check whether a newer `nvfortran` release has fixed it before reflexively
falling back, in case this exception has gone stale.

**`do concurrent` also can't contain a call to a non-`pure` procedure, or
safely collapse a dimension with a genuine cross-iteration dependency.**
Check every loop body for (a) calls to subroutines/functions that aren't
`pure`, and (b) anything accumulated *across* iterations of the dimension
being folded in, rather than freshly computed per-iteration. Either one
means that dimension must stay an ordinary sequential `do` instead of
joining the `do concurrent` index set — unless the accumulation can be
expressed as a `DO_LOCALITY(reduce(...))` clause, which puts it back on the
table. This comes up most often when collapsing a layer dimension (`k`)
that used to be a separate sequential loop into the same `do concurrent` as
`j`/`i` — see `references/kji-loop-blocking-guide.md` for a worked example
(`MOM_hor_visc.F90`'s MEKE-reduction loops kept a serial step for this).

### 2. `!$omp target teams loop collapse(n)` (fallback for automatic-array privatization)

When `do concurrent` isn't usable because of the automatic-array
privatization issue above, use `!$omp target teams loop` instead, with
`collapse(n)` for the `n` parallel loop levels:

```fortran
!$omp target teams loop collapse(2)
do j=js,je ; do i=is,ie
  ! loop body, including whatever needed the automatic array
enddo ; enddo
```

This is the preferred OpenMP fallback over the option below, based on how
it performs with `nvfortran` in practice.

### 3. `!$omp target distribute parallel do collapse(n)` (fallback for when `target teams loop` fails)

If `target teams loop` itself fails to compile or run correctly for a
given loop (e.g. a compiler bug unrelated to the automatic-array issue),
`!$omp target distribute parallel do collapse(n)` has worked as a
substitute in some cases:

```fortran
!$omp target distribute parallel do collapse(2)
do j=js,je ; do i=is,ie
  ! loop body
enddo ; enddo
```

Reach for this only after `target teams loop` has actually been tried and
shown to fail — it's a second-line fallback, not an alternative default.

### 4. `!$omp target` + `!$omp loop collapse(n)` (collapsing a sequential recurrence into one kernel launch)

Unlike options 2-3, this isn't a fallback from `do concurrent`'s capability
— it addresses a different problem entirely. A genuine sequential
recurrence (e.g. the `do k=nz,1,-1` loops covered in
`references/loop-blocking-guide.md`, with a real top-down/bottom-up
dependency) would otherwise launch one `do concurrent` kernel *per
iteration* of that outer sequential loop. Each individual kernel can be
fast, but if the outer loop runs many times (once per vertical level, for
example), the per-iteration launch/sync overhead adds up. This can only be used
if there are no routine calls within the loop.

Wrap the *whole* sequential loop in one `!$omp target` region, and use
`!$omp loop collapse(n)` for the per-iteration parallel work inside it,
instead of a fresh `do concurrent` kernel each iteration:

```fortran
!$omp target
do k=nz,1,-1
  !$omp loop collapse(2) private(ii,jj)
  do j=js,je ; do i=is,ie
    ii = i-is+1 ; jj = j-js+1
    if (mask(ii,jj)) then
      ! per-(k,j,i) work
    endif
  enddo ; enddo
enddo
!$omp end target
```

This is the pattern `MOM_continuity_PPM.F90`'s `zonal_flux_adjust`/
`meridional_flux_adjust` use for their Newton-iteration recurrence. A few
differences from `do concurrent` worth knowing:

- **No inline mask condition.** `!$omp loop` doesn't support `do
  concurrent`'s trailing `, mask(i,j) > 0` condition — use an ordinary
  `if (...) then ... endif` wrapping the loop body instead.
- **No `DO_LOCALITY` macro.** `!$omp loop`'s `private(...)` clause is plain
  OpenMP syntax, not `do concurrent`'s locality-specifier mechanism, so it
  doesn't need the compatibility macro at all.

If the recurrence has a `domore`/`exit`-style early-out, see "Reduce
kernel-launch count" in `references/loop-blocking-guide.md` for how that
interacts with this pattern.

## Data mapping: what's already resident, and what you must map yourself

Two large derived types worth knowing about up front:

- **`G` (grid), `US` (unit scaling), `GV` (vertical grid)** and their
  constituent arrays can be assumed already resident on device everywhere
  in the model, unless told otherwise for a specific case. You generally
  don't need mapping directives for these.
- **`tv` (thermodynamic variables)**, by contrast, is *not* yet
  persistently resident. Any subroutine touching `tv` or its constituent
  arrays (e.g. `tv%S`, `tv%T`) needs them explicitly mapped to device from
  outside the subroutine that uses them.

When mapping a derived type like `tv`, order matters as it does for
`*_CS` structs (see `references/data-mapping-conventions.md`, Rule 2): map
the parent type first, then its constituent arrays.

**The gotcha that has caused segfaults in practice**: once a constituent
array is mapped back with its own `target update from(...)` (because a
subroutine modified it), do **not** subsequently re-map the *parent type*
`tv` itself. Doing so afterward can clobber the arrays' host contents, or
leave them pointing at a device address that no longer exists — a
segfault or accelerator failure, not usually a silent wrong-answer bug, so
it's easy to miss until a run crashes. If only the arrays changed, map
only the arrays back.

For the four rules governing *where in the code* a mapping directive belongs
(subroutine-local temporaries, static vs. dynamic `*_CS` fields, and
subroutine argument arrays), see `references/data-mapping-conventions.md`.

## Preserving CPU and GPU performance together

A port is only useful if it (a) doesn't change answers relative to the CPU
path, (b) doesn't change answers relative to `dev/gpu`, our development
branch, and (c) doesn't slow the CPU path down. All three matter equally —
a numerically-correct port that regresses CPU performance is unacceptable,
and so is a fast GPU kernel that quietly reorders a summation (see the
arithmetic-reproducibility notes in `references/code-style-guide.md` — a
common way for (a)/(b) to fail silently).

**Loop tiling ("blocks") is the main tool for satisfying (c) without
sacrificing GPU throughput.** The idea: keep the same tile/block size on
both platforms, but choose different sizes per platform — cache-friendly
on CPU (often the default block size in `i`, size 1 in `j` and `k`), large
on GPU (often the whole default-sized 3-D array as one block, to maximize
parallel work per kernel launch). `MOM_continuity_PPM.F90` and
`MOM_CoriolisAdv.F90` are already ported this way and are good references
for how the tiling is structured.

See `references/loop-blocking-shared-patterns.md` for the commit
discipline, common pitfalls, and runtime-configurable block-size pattern
shared by both directions of blocking; the two paragraphs below cover what's
specific to each loop shape.

If the loop you're porting is `jki`-ordered (a `k`-recurrence, so only
`j`/`i` can be blocked) and isn't tiled yet, do that restructuring *first*,
as its own bisectable pass, before adding any offload directive — see
`references/loop-blocking-guide.md` for the step-by-step recipe
(introducing `njblock`/`niblock`, batching a per-row `MOM_EOS` call to its
2-D interface, and the shape-specific gotchas found).

If instead the loop is `kji`-ordered — each vertical layer's calculation is
independent, with no `k-1`/`k+1` recurrence — `k` itself can be blocked
*and* collapsed into the same `do concurrent` as `j`/`i`, rather than just
looped over sequentially outside it. See
`references/kji-loop-blocking-guide.md` for the k-blocking recipe
(`nkblock`, `k_start`/`k_end`, collapsing `k` into the parallel index set)
and the shape-specific gotchas found blocking `MOM_hor_visc.F90`'s
`horizontal_viscosity` this way.

**If a subroutine or module looks like it needs a substantial refactor to
get CPU and GPU performance to coexist** (not just a loop-construct swap),
don't push through it solo. Explain to the user why a refactor looks
necessary, sketch a couple of approaches, and ask how they'd like to
proceed — the same escalation pattern as the MOM_EOS case above.

## Blocking a loop: shared patterns and gotchas

The `jki` and `kji` blocking recipes (`references/loop-blocking-guide.md`,
`references/kji-loop-blocking-guide.md`) share the same commit discipline
and hit some of the same failure modes, regardless of which dimension gets
blocked. See `references/loop-blocking-shared-patterns.md` for that shared
material — commit discipline, extracting the target routine into a
submodule before blocking starts, the block-size-1-insufficient gotcha, the
`if(k==1)`-guarded-init race, the runtime-configurable-block-size pattern,
and how to diagnose a hang during a blocking pass. Each loop-shape guide
covers only what's specific to its own shape.

## Verifying a port

Before considering a port done:

- **Correctness (a and b above)**: compare the `ocean.stats` file from a
  run with the port against one without it (or from `dev/gpu`). They
  should match. If they don't, suspect a reordered reduction, a
  mapping-order bug (see the `tv` gotcha above), or a race from a missing
  `local(...)` clause before looking elsewhere.
- **Performance (c above)**: check this with the CPU clocks already built
  into `MOM_cpu_clock` — see `references/performance-timing.md` for how to
  find an existing clock for the component you're porting, or add one if
  it doesn't exist. Separately, a scripts directory is intended for
  testing target configurations against target compiler settings for
  broader performance regressions. It may not be populated yet — check
  whether it exists and has usable scripts before assuming it's available,
  and ask the user if you can't find it.
- **Data-transfer volume**: `nvfortran`'s `NV_ACC_NOTIFY` environment
  variable (set to `3`, or `2` for a lighter version) logs every host↔device
  upload, download, and kernel launch to stderr, tagged with source file,
  line, function name (lowercased), and — for transfers — variable name and
  byte count. Capture a run with it enabled, then filter the log for the
  subroutine you ported (its lowercased name) to isolate its activity from
  the rest of the run. Two things worth checking for: (1) a small scalar or
  array transferred far more often than expected (e.g. once per loop
  iteration rather than once per call) usually means a host-visible control
  variable is round-tripping unnecessarily — see the `domore`/`exit`
  discussion in `references/loop-blocking-guide.md`'s kernel-launch-count
  section; (2) a Rule-4 argument deliberately left unmapped showing up with
  *many* transfers, rather than none, means it isn't actually already
  resident the way you assumed, and needs real mapping after all. The log
  can be large (100MB+ for a full run) — don't read it directly, only
  filter/count against it.

## Reference files

- `references/loop-blocking-guide.md` — converting a `jki`-ordered loop
  (`k`-recurrence, only `j`/`i` blockable) into the `MOM_continuity_PPM.F90`
  tiled form before adding offload directives: the `njblock`/`niblock`
  skeleton-then-promote commit sequence, batching a per-row `MOM_EOS` call
  to its 2-D interface, and the shape-specific gotchas (EOS call-count
  blowup, shared-`dom` index consistency, array-shrink-before-loop-
  narrowing ordering) found blocking `MOM_isopycnal_slopes.F90` — see
  `references/loop-blocking-shared-patterns.md` for the shared pitfalls.
- `references/kji-loop-blocking-guide.md` — k-blocking a `kji`-ordered loop
  (no `k`-recurrence, so `k` collapses into the same `do concurrent` as
  `j`/`i`): the `nkblock`/`k_start`/`k_end` commit sequence and the
  shape-specific gotchas found blocking `MOM_hor_visc.F90`'s
  `horizontal_viscosity` — see `references/loop-blocking-shared-patterns.md`
  for the shared pitfalls.
- `references/loop-blocking-shared-patterns.md` — the commit discipline,
  extracting the target routine into a submodule before blocking starts,
  block-size-1-insufficient gotcha, `if(k==1)`-guarded-init race,
  runtime-configurable-block-size pattern, and hang-diagnosis heuristics
  shared by both loop-blocking guides above.
- `references/data-mapping-conventions.md` — the detailed rules for where
  `enter/exit data` and `target update` directives go for subroutine locals,
  static vs. dynamic `*_CS` fields, and subroutine argument arrays, plus a
  couple of `nvfortran`-specific mapping quirks (cheap zero-init `alloc` vs.
  `to`, and why conditionally-unused arrays still need at least an `alloc`).
- `references/performance-timing.md` — how to find or add a `MOM_cpu_clock`
  timer to check that a port hasn't regressed CPU performance.
- `references/code-style-guide.md` — a summary of MOM6's general Fortran
  style conventions (indentation, naming, loop-index letter conventions,
  arithmetic reproducibility) worth following in any code you write or
  restructure while porting. Full source:
  https://github.com/mom-ocean/MOM6/wiki/Code-style-guide
