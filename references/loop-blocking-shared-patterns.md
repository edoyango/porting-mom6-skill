# Loop blocking: shared patterns and gotchas

The `jki` and `kji` blocking recipes (`references/loop-blocking-guide.md`,
`references/kji-loop-blocking-guide.md`) share the same discipline and hit
some of the same failure modes, regardless of which dimension gets
blocked. This doc covers what's shared; each loop-shape guide covers what's
specific to its own shape.

**Commit discipline.** Land blocking work as small, individually-verified
commits: build, run `benchmark_ALE`, sha256-compare `ocean.stats` against
a known-good baseline, and commit only on a match. Don't proceed while the
current step is unverified — every commit should be independently
bisectable. Every step in both recipes is bit-identical to the
pre-blocking baseline by construction (block size changes only how work is
chunked, not per-point arithmetic or ordering) until offload directives
appear, so a mismatch always means a real bug in that commit, not
accumulated drift.

**Extract the target routine into a submodule before starting a blocking
pass.** The commit discipline above means many build/test cycles in a row.
If the routine still lives directly in its parent module, editing its body
still regenerates that module's `.mod` interface file on most builds, which
forces every other module that `use`s it to rebuild too — on a large,
widely-`use`d module that tax lands on every single commit in the
sequence. Splitting the routine out into a `submodule` of its parent module
keeps the parent's public interface (and `.mod` file) untouched by body-only
edits, so each iteration only rebuilds the submodule and relinks. Do this
once, as its own build/test-verified commit, before the blocking sequence
starts — not partway through — so every commit in the bisectable sequence
builds under the same, faster footprint. Once the blocking pass is accepted,
fold the routine back into its parent module as its own final build/test-
verified commit — the submodule split is a build-speed convenience for the
duration of the blocking work, not a structural change worth keeping in the
codebase afterward.

**A block-size-1 test is necessary but not sufficient.** Once an array
gains a block-sized dimension, every read/write site needs the
block-relative index. A missed site still compiles and stays bit-for-bit
correct at block size 1 (index 1 is the only valid index, so a stale
reference is accidentally right) — it's wrong only once the block size
exceeds 1. Since block-size-1 is how every intermediate commit above gets
verified, this bug class is invisible until you test a larger size. Always
verify at least once with a block size greater than 1 — ideally one that
doesn't evenly divide the loop extent, to exercise a partial final block
too — before considering a blocking port done. Found twice independently
in `MOM_hor_visc.F90` (a hardcoded index `1` left in a promoted array) and
once in `MOM_isopycnal_slopes.F90`'s mask handling (a loop header reading a
shrunk mask array with the un-converted global index) — see the reference
files for the concrete diffs.

**A one-time initialization guarded by `if (k==1) then ... endif`** (or
the equivalent for whichever dimension is being blocked) **is a race once
that loop is blocked or collapsed, not a safe "runs once" idiom.** It
relied on the old strictly-sequential loop guaranteeing the guarded
iteration runs alone, before anything else — once blocked, and especially
once inside a `do concurrent`, there's no such guarantee, and the
initialization can interleave with the accumulating writes in any order.
Fix by hoisting it entirely out of the loop into a plain sequential block
that runs once before the block loop starts, decoupling "runs once" from
"happens to be first in a loop about to be parallelized."

**Making a block size runtime-configurable.** Once validated with a
compile-time `integer, parameter`, this is the pattern
`MOM_continuity_PPM.F90` (and later `MOM_set_viscosity.F90`,
`MOM_hor_visc.F90`) converged on for a per-platform runtime choice — small
on CPU, large on GPU, from the *same* binary:

1. Add the block size as a **public integer field on the `*_CS` type**,
   defaulting to `0` — a sentinel meaning "not configured, use the full
   domain/column as one block."
2. Add a small resolver (a subroutine, or a plain `if` at the top of the
   caller) that substitutes the full extent when the field is `0`.
3. Change the array declarations inside the ported subroutine from
   `integer, parameter` to an `intent(in)` dummy argument.
4. At the call site, resolve first, then pass the value explicitly:
   `if (niblock==0) niblock = ie-is+1 ; call foo(..., niblock)`. Pass the
   resolved integer as a plain argument — don't resolve it inline with a
   `merge()` baked into a declaration's bound expression; an explicit,
   named argument is what every settled example in this codebase does.
5. Read the namelist override with `get_param(..., layoutParam=.true.)` —
   this marks it as affecting parallel decomposition/performance, not
   physics, the correct category for a block-size knob.

**Diagnosing a hang during a blocking pass.** MOM6's output only appears once
the run exits, so checking mid-run is useless — check `ps aux` for the MOM6
process's accumulated CPU time instead. A normal `benchmark_ALE` run completes
in roughly 100-103s; 150-250s (1.5-2.5x baseline) can also occur from correct
code. 10+ minutes, still climbing, stuck at timestep 0 is a real signature
and justifies killing the process to look for a defect — usually an accidental
blowup in call count (e.g. i-blocking an EOS call at too small a block size)
or an out-of-bounds write (e.g. an array shrunk to its block size before the
loop narrowing its fill range exists).
