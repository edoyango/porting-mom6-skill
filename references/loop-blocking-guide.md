# Loop Blocking Guide

How to convert a MOM6 `jki`-ordered loop nest (`do j; do k; do i`, with `k`
carrying a top-down or bottom-up recurrence that can't be parallelized) into
the tiled/blocked form used by `MOM_continuity_PPM.F90`, as **prep work
before** adding `do concurrent`/`!$omp target` (see the main `SKILL.md`).
Validated end-to-end on `MOM_isopycnal_slopes.F90`'s `calc_isoneutral_slopes`
(both its zonal/u-point and meridional/v-point slope loops), which produced
the pattern below. See `references/loop-blocking-gotchas.md` for the
gotchas found applying it, and `references/loop-blocking-offload.md` for
adding offload directives afterward.

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
    ! fill phase — paired on the same loop nest as the main phase below:
    do jj=1,jje ; do II=1,IIe
      j = jsb+jj-1 ; I = IsbB+II-1
      ! fill the batched-EOS input arrays (T/S/pres) for the whole tile
    enddo ; enddo
    ! one batched calculate_density_derivs(_2d) call over the whole tile
    do jj=1,jje ; do II=1,IIe
      j = jsb+jj-1 ; I = IsbB+II-1
      ! the original per-(j,K,i) body, reading the batched EOS outputs
    enddo ; enddo
  enddo
enddo ; enddo
```

**`jj`/`II` are the actual loop-control variables, directly inside
whichever loop already ranges over `j`/`i`** — the `K`-recurrence itself
above, or a k-independent fill/solve loop — never hoisted into a
standalone wrapper loop between the outer block loop and that loop. Native
`j`/`I` are computed from them as the loop body's first statement (`j =
jsb+jj-1 ; I = IsbB+II-1`), not the other way around, so the block-local
index stays the one real loop variable everywhere — including the
k-independent `do concurrent (jj=1:jje, II=1:IIe)`-style loops elsewhere in
the routine — with no index-variable-role-swap needed later. This is easy
to get wrong at `njblock=1` (see the gotcha in
`references/loop-blocking-gotchas.md`).

**Style convention: chain the outer/inner block-loop headers onto one
line.** `jsb`/`IsbB` share a line; the outer block's bound/count pair
(`jeb`/`jje`) is computed as the body's first statement, the inner pair
(`IebB`/`IIe`) as its second, followed by a blank line then the loop's
actual work. Close both loops the same way, `enddo ; enddo` on one line —
matching the `do j=... ; do i=...` single-line idiom already used for short
loops in this codebase, extended here to a body spanning many lines.

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
   trivially good). This `jj`-loop is a temporary placeholder here, still
   wrapping the whole `K`-recurrence — step 4 relocates it.
2. **Promote arrays, split fill/main phases.** Give the arrays that carry
   data from the fill phase to the main phase (the EOS call's inputs and
   outputs, plus anything else read across the `jj` boundary) an added
   `njblock`-sized dimension, indexed by `jj`. Split the loop body into a
   fill pass and a main pass over the same `jj` range. Scratch arrays that
   are filled and consumed *within* the same `jj` iteration don't need
   promoting — see the shared-`dom` gotcha in
   `references/loop-blocking-gotchas.md` before deciding which arrays need
   it.
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
   the batched EOS call's overhead — see `references/loop-blocking-gotchas.md`
   for why; **do not** start this at `niblock=1` the way step 1 safely could
   for `j`. `MOM_continuity_PPM.F90`'s CPU `default_niblock = 32` is the
   established precedent. Keep every promoted array's declared shape
   full-width in this commit — only the loops narrow to `istart:iend`.
   **Also relocate step 1's `jj`-loop here**, from wrapping the whole
   `K`-recurrence to sitting directly inside `K`, paired with the new `ii`
   loop on the same line, matching the target shape above — see the gotcha
   in `references/loop-blocking-gotchas.md` for why this is easy to miss.
5. **Shrink arrays to the i-block.** Only now replace the full-width leading
   dimension on the arrays promoted in step 2 with `niblock`, indexed by
   `ii`. This must come strictly after step 4, never combined with or ahead
   of it — see `references/loop-blocking-gotchas.md`.

## Next

Once this recipe is validated on a first loop, see
`references/loop-blocking-gotchas.md` for applying it to a second,
structurally-similar loop and the correctness gotchas found doing so, and
`references/loop-blocking-offload.md` for adding offload directives, making
the block size runtime-configurable, and worked examples.

