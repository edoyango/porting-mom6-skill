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

Four constructs are available, tried in order — see
`references/loop-constructs.md` for the full decision tree: exact syntax,
the `DO_LOCALITY` macro rules, the automatic-array-privatization compiler
bug that rules `do concurrent` out for some loops, and code examples for
each.

1. **`do concurrent`** — default choice: standard Fortran, best `nvfortran`
   code generation. Can't privatize an automatic array per thread (a live
   `nvfortran` 26.5 compiler bug) or contain a call to a non-`pure`
   procedure.
2. **`!$omp target teams loop collapse(n)`** — fallback when `do concurrent`
   hits the automatic-array issue above. Preferred OpenMP fallback based on
   `nvfortran` performance in practice.
3. **`!$omp target distribute parallel do collapse(n)`** — second-line
   fallback, only once `target teams loop` has actually been tried and shown
   to fail for a given loop.
4. **`!$omp target` wrapping a sequential loop, with `!$omp loop
   collapse(n)` inside it** — addresses a different axis than 1-3 (kernel-
   launch overhead across a genuine sequential recurrence, e.g. the `do
   k=nz,1,-1` loops in `references/loop-blocking-guide.md`), not which
   construct can execute a parallel loop, so it can apply regardless of
   which of 1-3 handles the parallel work inside it. Only usable if there
   are no routine calls within the loop. If the recurrence has a
   `domore`/`exit`-style early-out, see "Reduce kernel-launch count" in
   `references/loop-blocking-offload.md` for how that interacts with this
   pattern.

Each fallback exists to solve a specific, concrete problem with the option
above it, not out of general preference — don't skip to a fallback unless
its trigger actually applies.

## Data mapping: what's already resident, and what you must map yourself

**`G` (grid), `US` (unit scaling), `GV` (vertical grid)** and their
constituent arrays can be assumed already resident on device everywhere in
the model — you generally don't need mapping directives for these. **`tv`
(thermodynamic variables)**, by contrast, is *not* yet persistently
resident: any subroutine touching `tv` or its constituent arrays (`tv%S`,
`tv%T`) needs them explicitly mapped to device from outside the subroutine
that uses them, parent type first then constituent arrays — same ordering
as `*_CS` structs (`references/data-mapping-conventions.md`, Rule 2).

**Segfault gotcha**: once a constituent array is mapped back with its own
`target update from(...)` (because a subroutine modified it), do **not**
subsequently re-map the *parent type* `tv` itself — doing so can clobber
the arrays' host contents or leave them pointing at a freed device
address. If only the arrays changed, map only the arrays back.

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
sacrificing GPU throughput**: same tile/block size on both platforms, but a
different size per platform — cache-friendly on CPU (often the default
block size in `i`, size 1 in `j`/`k`), large on GPU (often the whole
default-sized 3-D array as one block, to maximize work per kernel launch).
`MOM_continuity_PPM.F90` and `MOM_CoriolisAdv.F90` are already ported this
way and are good references for how the tiling is structured.

See `references/loop-blocking-shared-patterns.md` for the commit
discipline, pitfalls, and runtime-configurable block-size pattern shared by
both directions of blocking below.

If the loop is `jki`-ordered (a `k`-recurrence, so only `j`/`i` can be
blocked) and isn't tiled yet, do that restructuring *first*, as its own
bisectable pass, before adding any offload directive — see
`references/loop-blocking-guide.md` for the step-by-step recipe.

If instead the loop is `kji`-ordered — each vertical layer's calculation is
independent, with no `k-1`/`k+1` recurrence — `k` itself can be blocked
*and* collapsed into the same `do concurrent` as `j`/`i`, rather than
looped over sequentially outside it. See
`references/kji-loop-blocking-guide.md` for the k-blocking recipe.

**If a subroutine or module looks like it needs a substantial refactor to
get CPU and GPU performance to coexist** (not just a loop-construct swap),
don't push through it solo — explain why to the user, sketch a couple of
approaches, and ask how they'd like to proceed, the same escalation
pattern as the MOM_EOS case above.

## Verifying a port

Before considering a port done:

- **Correctness (a and b above)**: compare the `ocean.stats` file from a
  run with the port against one without it (or from `dev/gpu`). They
  should match. If they don't, suspect a reordered reduction, a
  mapping-order bug (see the `tv` gotcha above), or a race from a missing
  `local(...)` clause before looking elsewhere.
- **Performance (c above)**: check this with the CPU clocks already built
  into `MOM_cpu_clock` — see `references/performance-timing.md` for how to
  find an existing clock, or add one if it doesn't exist. Separately, a
  scripts directory is intended for testing target configurations against
  target compiler settings for broader performance regressions. It may not
  be populated yet — check whether it exists before assuming it's
  available, and ask the user if you can't find it.
- **Data-transfer volume**: `nvfortran`'s `NV_ACC_NOTIFY` environment
  variable (set to `3`, or `2` for lighter) logs every host↔device transfer
  and kernel launch to stderr, tagged with source location and — for
  transfers — variable name and byte count. Capture a run with it enabled,
  filter the log for the subroutine you ported. Watch for: (1) a small
  scalar/array transferred far more often than expected (e.g. once per loop
  iteration, not once per call) — a host-visible control variable
  round-tripping unnecessarily, see the `domore`/`exit` discussion in
  `references/loop-blocking-offload.md`'s kernel-launch-count section; (2) a
  Rule-4 argument deliberately left unmapped showing *many* transfers rather
  than none — it isn't actually already resident. The log can be large
  (100MB+) — don't read it directly, only filter/count against it.

## Reference files

- `references/loop-constructs.md` — the four loop-construct options in
  full: syntax, `DO_LOCALITY` rules, and the automatic-array compiler bug.
- `references/loop-blocking-guide.md` — tiling a `jki`-ordered
  (`k`-recurrence) loop into the `MOM_continuity_PPM.F90` form: target
  shape and the `njblock`/`niblock` step-by-step commit sequence.
- `references/loop-blocking-gotchas.md` — applying the `jki` recipe to a
  second, structurally-similar loop, its correctness gotchas (EOS
  call-count blowup, shared-`dom` index consistency, array-shrink
  ordering), and distinguishing a real hang from transient slowness.
- `references/loop-blocking-offload.md` — adding offload directives once a
  `jki` loop is blocked, making the block size runtime-configurable, and
  worked examples (`MOM_isopycnal_slopes.F90`, `MOM_set_viscosity.F90`).
- `references/kji-loop-blocking-guide.md` — k-blocking a `kji`-ordered
  (no-recurrence) loop: the `nkblock`/`k_start`/`k_end` commit sequence,
  gotchas, and worked examples from `MOM_hor_visc.F90`.
- `references/loop-blocking-shared-patterns.md` — commit discipline,
  submodule-extraction, block-size-1, `if(k==1)`-race, runtime-config, and
  hang-diagnosis patterns shared by every blocking guide above.
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
