# Choosing a GPU loop construct

How to pick which construct (`do concurrent`, `!$omp target teams loop`,
`!$omp target distribute parallel do`, or a wrapping `!$omp target` +
`!$omp loop`) to use for a given parallel loop when porting to GPU — see
the main `SKILL.md` for how this fits into the overall porting workflow
(the MOM_EOS check that must happen first, data mapping, blocking, and
verification).

Try options 1-3 in order. Each fallback exists to solve a specific, concrete
problem with the option above it, not out of general preference — don't
skip to a fallback unless its trigger actually applies. Option 4 addresses a
different axis entirely (kernel-launch overhead across a sequential loop,
not which construct can execute a given parallel loop) and can apply
regardless of which of 1-3 you're using for the parallel work inside it.

## 1. `do concurrent` (default choice)

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

## 2. `!$omp target teams loop collapse(n)` (fallback for automatic-array privatization)

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

## 3. `!$omp target distribute parallel do collapse(n)` (fallback for when `target teams loop` fails)

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

## 4. `!$omp target` + `!$omp loop collapse(n)` (collapsing a sequential recurrence into one kernel launch)

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
kernel-launch count" in `references/loop-blocking-offload.md` for how that
interacts with this pattern.

