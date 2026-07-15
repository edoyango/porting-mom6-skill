# Data mapping conventions

These are the conventions for getting data onto and off of the device. They
govern *where in the code* a `!$omp target ... data` or `!$omp target update`
directive belongs, not which loop construct to use (see `SKILL.md` for the
loop-construct decision tree).

This document was inherited from an earlier porting pass validated on the
`benchmark` and `benchmark_ALE` configurations. Treat it as a strong default,
not an exhaustive rulebook — extend it as more of the model is ported.

## Rule 1 — Subroutine-local variables

Local arrays and scalars that are only needed inside a single subroutine should
be allocated on the device at the start of that subroutine and deleted at the
end.

```fortran
! At the top of the subroutine (after index setup):
!$omp target enter data map(alloc: local_array_a, local_array_b, local_scalar)

! ... body of subroutine, including do concurrent loops ...

! At the bottom of the subroutine (before end subroutine):
!$omp target exit data map(delete: local_array_a, local_array_b, local_scalar)
```

A scalar used only as a `local` inside a `do concurrent` loop does **not** need
to appear in this subroutine-level `map(alloc:)` directive — the `local(...)`
clause on the loop itself is sufficient and preferred.

## Rule 2 — CS variables that are static at runtime

Fields inside a `*_CS` control structure that are set once during
initialization and never modified again (e.g. precomputed grid metrics such as
`f2_dx2_h`, `beta_dx2_h`) should be:

- Mapped to the device during the corresponding `*_init` subroutine, **after**
  the parent struct itself has been mapped.
- Deleted from the device during the corresponding `*_end` subroutine,
  **before** the parent struct is deleted.

Mapping order matters — the parent struct goes to the device before its
allocatable component arrays, and comes off after them:

```fortran
! In *_init (at the very end, after all allocations and data filling):
!$omp target enter data map(to: CS)                  ! parent struct first
if (CS%some_flag) then
  !$omp target enter data map(to: CS%static_array_a, CS%static_array_b)
  !$omp target enter data map(alloc: CS%output_array) ! output: alloc, not 'to'
endif

! In *_end (before any host deallocates):
if (allocated(CS%output_array)) &
  !$omp target exit data map(delete: CS%output_array, CS%static_array_a, CS%static_array_b)
! ... host deallocates ...
!$omp target exit data map(delete: CS)               ! parent struct last
```

Use `map(to:...)` for arrays already filled at init time and only read on the
GPU. Use `map(alloc:...)` for arrays whose device values are written by GPU
loops during the run (the host values at init time aren't needed on the
device).

**Compiler-specific optimization**: with `nvfortran`, an array allocated on
the device via `map(alloc:...)` comes back zero-initialized. If the host
array is already zero at the point you're mapping it — a common case for
freshly-allocated diagnostic or accumulator arrays — you can map it with
`map(alloc:...)` instead of `map(to:...)` and skip the host-to-device copy
entirely, since `alloc` is a cheaper operation than `to`. Only make this
substitution when you've actually confirmed the host value is zero at that
point in the code; if it could be anything else, `map(to:...)` is still
required.

This same parent-struct-first, arrays-second ordering applies to any derived
type that needs to be mapped explicitly — see "Derived-type residency" in
`SKILL.md` for the `tv` example and a mapping-order gotcha that has caused
segfaults in practice.

## Rule 3 — CS variables that change between calls

Some CS fields are recomputed on the **CPU** between successive calls to an
offloaded subroutine (e.g. by a sub-call inside the subroutine itself, or by
other model code that runs on the CPU before the next call). These should be:

- Mapped to the device *before* the call using `!$omp target update to(...)`,
  placed within the CPU clock region surrounding the call in the caller (e.g.
  `MOM.F90`).
- Mapped from the device *after* the call using `!$omp target update
  from(...)`, or released with `!$omp target exit data map(release:...)` if
  the caller no longer needs the host copy — also within the same clock
  region.

```fortran
! In MOM.F90 (example):
call cpu_clock_begin(id_clock_some_routine)
!$omp target update to(CS%some_cs_field_updated_on_cpu)
call some_offloaded_routine(...)
!$omp target update from(CS%some_cs_output_array)
call cpu_clock_end(id_clock_some_routine)
```

If a CS field is instead updated *inside* the offloaded subroutine (e.g. by a
CPU sub-call within it), put the `target update to(...)` inside the
subroutine, immediately after that CPU call — not in the caller.

## Rule 4 — Other subroutine inputs and outputs

Arrays passed as arguments to an offloaded subroutine (rather than living in a
CS struct) should be mapped to/from the device around the call site in the
caller, again within the clock region:

```fortran
call cpu_clock_begin(id_clock_routine)
!$omp target update to(input_array)
call offloaded_routine(input_array, output_array, ...)
!$omp target update from(output_array)
call cpu_clock_end(id_clock_routine)
```

Many of these arrays (e.g. `h`, grid arrays in `G`) are already persistently
mapped on the device for other parts of the model — check `MOM.F90`
init/finalization before adding a new `map` directive, to avoid double-mapping.

Porting of rule-4 inputs has been deferred for many subroutines because these
arrays tend to involve complicated ownership/lifetime questions. When in
doubt, leave rule-4 inputs for a dedicated follow-up pass rather than guessing.

## Summary table

| Variable category | Where to map to device | Where to map from device | Where to delete |
|---|---|---|---|
| Subroutine-local arrays/scalars | Start of subroutine (`enter data map(alloc:)`) | As needed before CPU reads inside subroutine | End of subroutine (`exit data map(delete:)`) |
| Static CS arrays (set at init, read-only at runtime) | End of `*_init` (`enter data map(to:)`, after parent struct) | Not needed (never modified on device) | Start of `*_end`, before host `deallocate` |
| Dynamic CS output arrays (written by GPU loops) | End of `*_init` (`enter data map(alloc:)`, after parent struct) | After GPU loops, unconditionally, before any CPU read | Start of `*_end`, before host `deallocate` |
| CS fields updated by CPU between calls | Before the offloaded call in caller (`target update to`) | After the call in caller (`target update from`) | In `*_end` |
| CS fields updated by CPU *inside* the offloaded subroutine | Inside the subroutine, after the CPU update (`target update to`) | Inside the subroutine, before any subsequent CPU read | In `*_end` |
| Subroutine argument arrays | Before call in caller (`target update to`) | After call in caller (`target update from`) | Managed by whoever owns the array |
| `tv` and its constituent arrays (not yet persistently resident) | Outside the subroutine that calls into it; parent type first, then arrays | Map only the modified arrays back, never the parent type after arrays have been remapped | Wherever `tv` is torn down |

## nvfortran mapping quirks worth knowing

Beyond the zero-initialization behavior noted under Rule 2, there's a second
`nvfortran` quirk that affects arrays referenced only inside conditionally-run
code — e.g. something guarded by `if (CS%some_condition) then` where that
branch isn't exercised in every configuration:

**Even if a kernel doesn't end up touching such an array in a given
configuration, that array still needs to be at least `map(alloc:...)`'d on
the device.** If it's left completely unmapped, `nvfortran` appears to fall
back to managing the transfer itself implicitly, silently mapping the array
back and forth around the kernel — which can be a significant, easy-to-miss
performance regression for large arrays, especially in big kernels with many
conditional branches where several such arrays might be affected at once. The
fix is cheap: map it (`alloc` is enough if nothing needs to be on the device
yet) regardless of whether the current configuration's branches will
reference it.

## Worked examples

For the loop-tiling ("blocks") approach to preserving both CPU and GPU
performance, see `MOM_continuity_PPM.F90` and `MOM_CoriolisAdv.F90`, which have
already been ported successfully with that pattern. As other modules land and
stabilize, add pointers to them here.
