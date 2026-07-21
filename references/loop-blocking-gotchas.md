# Loop blocking (jki): reuse and correctness gotchas

Continuation of `references/loop-blocking-guide.md` — read that first for
the target shape and the step-by-step blocking recipe this assumes. Covers
applying the validated recipe to a second, structurally-similar loop, the
correctness gotchas found blocking `MOM_isopycnal_slopes.F90`, and telling
a real hang from transient slowness during the pass.

## Applying this to a second, structurally-similar loop

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
copies), the shrink step (step 5 in `references/loop-blocking-guide.md`)
can't happen until **both** loops' `istart` skeletons exist — shrinking
after only one loop's skeleton is in place reintroduces the
out-of-bounds-write ordering bug below for whichever loop hasn't been
narrowed yet.

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

**Don't leave the step-1 `jj`-loop wrapping the whole `K`-recurrence once
i-blocking (step 4) is added.** It must relocate to sit directly inside
`K`, paired with the new `ii` loop, as the actual loop-control variables —
see "The target shape" in `references/loop-blocking-guide.md` for the rule
and why. Leaving it wrapping `K` still compiles and gives bit-identical
answers at `njblock=1` (it degenerates to one trivial iteration), so the
bit-for-bit `ocean.stats` check doesn't catch it — this was caught by
inspection, not a failed test, when applying the recipe to
`MOM_thickness_diffuse.F90`'s `thickness_diffuse_full`.

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
