---
name: mom6-gpu-loop-porting
description: Use this skill whenever porting MOM6 (Modular Ocean Model 6) Fortran code to run on GPU, or reviewing/debugging existing GPU-offloaded MOM6 code. This covers converting loops to `do concurrent` or OpenMP `!$omp target` offload directives, deciding where `!$omp target enter/exit data` and `target update` directives belong, handling MOM6 control-structure (`*_CS`) fields and derived types like `tv`/`G`/`US`/`GV` on the device, working around the MOM_EOS polymorphic-type GPU limitation, and preserving CPU performance while adding GPU support (e.g. loop tiling / "blocks"). Trigger this for tasks like "port this subroutine to GPU", "add OpenMP offload to this loop", "why is this GPU kernel giving different answers than the CPU", "how should I map this CS field to the device", or any work touching `!$omp target`, `do concurrent`, or GPU offload anywhere in the MOM6 codebase — even if the user doesn't use the word "skill" or spell out these directive names explicitly.
---

# MOM6 GPU Loop Porting

## Status and scope

MOM6's GPU port is a work in progress. The conventions here have been
validated on the `benchmark` and `benchmark_ALE` configurations, but the
codebase is large and only the code paths those configurations exercise have
been ported so far. Treat everything below as the current best practice, not
a finished spec — when you hit a pattern this document doesn't cover, make a
reasonable decision consistent with the philosophy here, flag it to the user,
and suggest adding it to this skill once it's confirmed to work.

The overall approach is: compute on the GPU with `do concurrent` or OpenMP
`target` loop constructs, manage data movement explicitly (don't rely on the
compiler to infer transfers), and only port the branches that are actually
exercised by the target configurations.

## Before porting any loop: check for MOM_EOS calls

GPUs can't execute through a polymorphic Fortran type's deferred procedures
from inside a kernel — `MOM_EOS` subroutines are typically dispatched this
way (an abstract `EOS` type selecting the equation-of-state routine at
runtime). This doesn't mean MOM_EOS-adjacent code is unportable: a kernel can
still be launched within an EOS subroutine on data that's already on the device, as
long as *that subroutine's own kernels* never touch the polymorphic `EOS`
type or its deferred procedures internally. The constraint is specifically
about the polymorphic dispatch happening inside a kernel, not about EOS logic
being off-limits.

So, before converting a loop:

1. Check whether the loop body calls a MOM_EOS subroutine.
2. If it does, check whether that call can be cleanly separated from the rest
   of the loop (e.g. hoisted into its own loop, called ahead of time into a
   temporary array) without a substantial restructuring of the surrounding
   code.
3. If separation is easy — do it, and port the surrounding loop normally.
4. If separation is *not* easy — a real refactor would be needed — **stop**.
   Don't attempt the port. Explain to the user why the EOS call blocks a
   straightforward port, sketch 1-2 plausible refactor approaches, and ask
   how they'd like to proceed before writing any offload code. Guessing here
   risks silently producing a kernel that won't compile or will crash on
   device, or a big unwanted refactor the user didn't ask for.

## Choosing a loop construct

Try these in order. Each fallback exists because of a specific, concrete
problem with the option above it — not because of general preference — so
don't skip straight to a fallback unless the described trigger actually
applies.

### 1. `do concurrent` (default choice)

`do concurrent` is preferred by default: it's standard Fortran (so the code
stays portable and readable even off-GPU), and it currently gets the best
code generation from `nvfortran`, the compiler this port targets.

```fortran
do concurrent( j=js:je, i=is:ie ) local(tmp)
  tmp = some_expression(i, j)
  output(i,j) = tmp * other_array(i,j)
enddo
```

Notes on `do concurrent` syntax, since it differs from an ordinary `do` loop
in easy-to-miss ways:

- Loop bounds are separated by a **colon**, not a comma:
  `do concurrent( i = is:ie, j = js:je )` — not `do i=is,ie; do j=js,je`.
- Collapse multi-dimensional loops into a single `do concurrent` statement
  rather than nesting separate `do concurrent` loops for each dimension.
  Nesting a `do concurrent` inside another `do concurrent` is fine (e.g. an
  outer 2-D horizontal loop containing an inner k-loop that itself needs to
  be a `do concurrent`) — it's collapsing what should be one multi-index
  statement into several separate ones that should be avoided.
- Declare loop-private scalar temporaries with `local(tmp)` so the GPU
  parallelization doesn't race on them.
- Remove any pre-existing `!$OMP parallel do` directive when converting a
  loop to `do concurrent` — the two shouldn't coexist on the same loop.
- `do concurrent` supports a mask condition, which is a natural fit for the
  common MOM6 pattern of `do i; do j; if (mask(i,j) > 0) ...`:

  ```fortran
  do concurrent( j=js:je, i=is:ie, mask(i,j) > 0 )
    output(i,j) = some_expression(i,j)
  enddo
  ```

  Prefer this masked form over keeping an `if` as the first statement inside
  the loop body when the condition is simply gating whether the iteration
  runs at all.

**The one case where `do concurrent` cannot be used**: a loop that needs to
privatize an *automatic array* per thread (i.e. an array-valued temporary
sized from loop or subroutine bounds, as opposed to the scalar temporaries
`local(...)` is meant for). As of the current `nvfortran` release (26.5),
there's a known compiler bug in this scenario. If you hit this, don't fight
it — fall through to the next construct. Since this is a live compiler bug,
check whether it's been resolved in a newer `nvfortran` release before
reflexively falling back, in case this exception has become stale.

### 2. `!$omp target teams loop collapse(n)` (fallback for automatic-array privatization)

When `do concurrent` isn't usable because of the automatic-array
privatization issue above, use `!$omp target teams loop` instead, with
`collapse(n)` for the `n` loop levels that are actually parallel:

```fortran
!$omp target teams loop collapse(2)
do j=js,je
  do i=is,ie
    ! loop body, including whatever needed the automatic array
  enddo
enddo
```

This is the preferred OpenMP fallback (over the option below) based on how it
performs with `nvfortran` in practice.

### 3. `!$omp target distribute parallel do collapse(n)` (fallback for when `target teams loop` fails)

If `target teams loop` itself fails to compile or run correctly for a given
loop (e.g. due to a compiler bug rather than the automatic-array issue),
`!$omp target distribute parallel do collapse(n)` has worked as a substitute
in some of these cases:

```fortran
!$omp target distribute parallel do collapse(2)
do j=js,je
  do i=is,ie
    ! loop body
  enddo
enddo
```

Reach for this only after `target teams loop` has actually been tried and
shown to fail — it's a second-line fallback, not an alternative default.

## Data mapping: what's already resident, and what you must map yourself

Two large derived types worth knowing about up front:

- **`G` (grid), `US` (unit scaling), `GV` (vertical grid)** and their
  constituent arrays can be assumed to already be resident on the device
  everywhere in the model, unless a user tells you otherwise for a specific
  case. You generally don't need to add mapping directives for these.
- **`tv` (thermodynamic variables)**, by contrast, is *not* yet persistently
  resident. Any subroutine that touches `tv` or its constituent arrays (e.g.
  `tv%S`, `tv%T`) needs those explicitly mapped to the device, from outside
  the subroutine that uses them.

When mapping a derived type like `tv`, order matters the same way it does for
`*_CS` structs (see `references/data-mapping-conventions.md`, Rule 2): map the
parent type first, then its constituent arrays.

**The gotcha that has caused segfaults in practice**: once you've mapped a
constituent array back with its own `target update from(...)` (because a
subroutine modified it), do **not** subsequently re-map the *parent type*
`tv` itself. Doing so after the arrays have already been updated can clobber
the array contents on the host, or leave them pointing at a device address
that no longer exists, which shows up as a segfault or an accelerator
failure — not usually as a silent wrong-answer bug, so it's easy to miss
until a run crashes. If only the arrays changed, map only the arrays back.

For the four rules governing *where in the code* a mapping directive belongs
(subroutine-local temporaries, static vs. dynamic `*_CS` fields, and
subroutine argument arrays), see `references/data-mapping-conventions.md`.

## Preserving CPU and GPU performance together

A port is only useful if it (a) doesn't change answers relative to the CPU
path, (b) doesn't change answers relative to `dev/gpu`, our development
branch, and (c) doesn't slow the CPU path down. All three matter equally —
a numerically-correct port that regresses CPU performance is not acceptable,
and neither is a fast GPU kernel that quietly reorders a summation (see the
arithmetic-reproducibility notes in `references/code-style-guide.md` — this
is a common way for (a)/(b) to fail silently).

**Loop tiling ("blocks") is the main tool for satisfying (c) without
sacrificing GPU throughput.** The idea: keep the same tile/block size for the
loop on both platforms, but choose different sizes per platform —
cache-friendly on CPU (often the default block size in `i`, and size 1 in `j`
and `k`), and large on GPU (often the whole default-sized 3-D array as one
block, to maximize parallel work per kernel launch). `MOM_continuity_PPM.F90`
and `MOM_CoriolisAdv.F90` have already been ported with this pattern
successfully and are good references for how the tiling is structured.

**If a subroutine or module looks like it'll need a substantial refactor to
get CPU and GPU performance to coexist** (not just a loop-construct swap),
don't push through it solo. Explain to the user why a refactor looks
necessary, sketch out a couple of ways it could be approached, and ask how
they'd like to proceed — the same escalation pattern as the MOM_EOS case
above.

## Verifying a port

Before considering a port done:

- **Correctness (a and b above)**: compare the `ocean.stats` file from a run
  with the port against `ocean.stats` from a run without it (or from
  `dev/gpu`). They should match. If they don't, suspect a reordered
  reduction, a mapping-order bug (see the `tv` gotcha above), or a race from
  a missing `local(...)` clause before looking elsewhere.
- **Performance (c above)**: check this with the CPU clocks already built
  into `MOM_cpu_clock` — see `references/performance-timing.md` for how to
  find an existing clock for the component you're porting, or add one if it
  doesn't exist yet. Separately, there's also a scripts directory intended
  for testing the target configurations against our target compiler settings
  for performance regressions more broadly. It may not be populated yet —
  check whether it exists and has usable scripts before assuming it's
  available, and ask the user if you can't find it.

## Reference files

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
