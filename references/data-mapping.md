# Data mapping and host/device transfers

Where `!$omp target enter/exit data` and `target update` go, which map kind to use, and the
mistakes that pass every benchmark. Findings are for nvfortran 26.3 unless marked. Tags as in
`loop-constructs.md`. Reproducers: uwagura/nvfortran-mres repo, `data_mapping/`.

Contents: 1 check what is resident · 2 where each kind of array is mapped · 3 map kinds ·
4 derived types · 5 transfers to and from host-only code · 6 hazards · 7 auditing a routine

## 1. Check what is already resident

Before adding a map for anything, grep for existing ones:

```bash
grep -rn "<name>" src/ --include=*.F90 | grep -iE "omp.*(enter data|exit data|update)"
```

Once an array is resident, code elsewhere may only `update` it (section 3). The usual state:
- `G`: a fixed list of metric and mask arrays is mapped in `initialize_MOM` (`MOM.F90`), not all
  of `G`. A `G%` array missing from that list needs its own map.
- `GV`: only `GV`, `GV%Rlay` and `GV%g_prime`. Other `GV` arrays are not resident.
- `US`: mapped whole in `MOM.F90`.
- `tv`: not persistently resident; it is mapped around the calls that need it. Map `tv` and
  its arrays from outside the routine that uses them.
- `visc` (`vertvisc_type`) members: persistently mapped in `set_visc_init`.

**Consider making a type persistently resident as part of a port.** When most of the code that
reads a control structure or derived type runs on the device, map its arrays once in `_init`
(released in `_end`) instead of around each call. Then bracket the remaining host-only consumers
with `update from`/`update to`, guarded by the same condition that selects the host path.

## 2. Where each kind of array is mapped

| Kind | Map | Unmap |
|---|---|---|
| CS member (`ALLOCABLE_`) | in `<mod>_init`, next to its `ALLOC_`, after the parent CS | in `<mod>_end`, before `DEALLOC_`, in reverse order, `delete` |
| Subroutine scratch | at routine entry, `map(alloc:)`, **below every early `return`** | at return, `delete` or `release`; earlier once its last use has passed |
| Dummy argument | not in the callee: the caller (or the array's owner) maps it | same |
| Pointer member | `map(to:)`, guarded `if (associated(x))` | mirrored, same guard |

- **Every array a kernel references needs an explicit map**, even one used only in a branch
  that never runs in your configuration. Without one, nvfortran inserts an implicit copy and
  transfers the whole array at every launch: an unmapped 8 MB array referenced only under
  `if (.false.)` was uploaded on each of three launches `[run-verified]`. Answers stay right, so
  only `-Minfo` (`Generating implicit copyin(...)`) or `NVCOMPILER_ACC_NOTIFY=2` shows it.
  `map(alloc:)` is enough for an array the branch never reads.
- Inputs and outputs of a ported routine called from host code: `update to` before the call and
  `update from` after it, in the caller, inside the call's `cpu_clock` region.
- If a CS field is changed by host code *inside* the ported routine, put the `update to` right
  after that host code, not in the caller.
- New persistent arrays: `allocatable` (macro form), not `pointer`. Use `pointer` only for
  restart-registry targets or real aliasing.
- Guard maps with the intrinsic that matches the declaration: `allocated()` for allocatables,
  `associated()` for pointers.

## 3. Map kinds

**`map(alloc:)` does not give zeros.** The device copy holds whatever was there before
`[run-verified]`. A fresh allocation happened to read as zero, but after any same-size block had
been freed, `map(alloc:)` returned that block's old values. That included a host-zeroed automatic
mapped `alloc` in a routine called repeatedly. Never rely on `alloc` for a zero device array.

**Preferred: `map(alloc:)` plus device-side initialization.** It needs no host-to-device
transfer and is always correct. For an array that must start at zero, allocate it and zero it in
a `do concurrent`, as `zonal_mass_flux` does with `zeros` (`MOM_continuity_PPM.F90`). Use
`map(to:)` only when the device reads values the host computed.

- `map(to:)` is required whenever the device reads host-set scalars, pointer descriptors, or
  `associated()`/`allocated()` state of a struct. `map(alloc: Reg, Reg%Tr(:))` read garbage
  descriptors and changed answers with the number of GPUs (`a774eb331`).
- **`map(to:)` on an object that is already present copies nothing**: it only raises the
  reference count `[run-verified]`. Likewise `exit data map(from:)` on an object still mapped
  elsewhere does not copy back. Refresh resident data only with `target update to/from`.
- **Neither `delete` nor `release` copies back.** Use `update from` or `map(from:)` first if the
  host needs the values.
- **`delete` forces the reference count to zero**, destroying outer mappings of the same object.
  Use `delete` only in the scope that owns the mapping (`_end`, or the routine that mapped its own
  scratch). Use `release` for anything another scope may also have mapped.
- `target update` of data that is not present does nothing, with no error.

## 4. Derived types

- Map the parent before its members, and unmap in reverse:
  `map(to: CS)` (or `alloc` for a shell whose `_init` fills it), then `map(to: CS%a)`.
  Plain scalar members travel with the parent.
- A pointer component needs the parent and component together, which attaches the device copy:
  `map(to: tv, tv%eqn_of_state)` `[run-verified]`. Without it the kernel faults with
  `CUDA_ERROR_ILLEGAL_ADDRESS`. Put the directive where `tv` becomes resident, not in the leaf.
- **Whole-struct `target update to(CS)` only before any member is attached.** After a member is
  attached it overwrites the device copy's descriptor with the host address, and the next kernel
  touching the member dies with `CUDA_ERROR_ILLEGAL_ADDRESS` `[run-verified]`. The `_init`
  routines that use it (`hor_visc_init`, `barotropic_init`, `VarMix_init`) map their members
  after it. To refresh a scalar later, update the scalar: `update to(CS%dtbt)`.
- **Never copy a parent back to the host** (`update from(CS)` or `exit data map(from: CS)`) once
  a member is attached. Both hung the test program `[run-verified]`. Copy back only the member
  arrays. Re-entering the parent with another `map(to: CS)` is harmless but does nothing
  (section 3).
- Don't map an array of structs element by element. Each element's members are a separate attach,
  which roughly doubled GPU time in tracer horizontal diffusion; flatten them to one array with
  an extra dimension (`1865612de`) `[source-only]`.
- Two shifted sections of one array (`ent_s(:,:,1:nz)` and `ent_s(:,:,2:nz+1)` passed as
  `ea`/`eb`) overlap in the present table and fail as "partially present". Map the whole parent
  once, at the call site.

## 5. Transfers to and from host-only code

Host-only consumers need a dominating `update from` for every array they read that the device
wrote, and an `update to` afterwards for every array they write that the device reads next:
- `post_data` and the diag mediator; `hchksum`/`uvchksum`/... (`MOM_checksums.F90`);
  `save_restart`. A newly device-resident field registered for restarts needs an `update from`
  before restarts are written.
- `pass_var`, `pass_vector`, and `do_group_pass` or `start/complete_group_pass` without
  `omp_offload=.true.`. With the flag, they work on device data.
- Unported modules, `MOM_error`, any I/O.
- `!$OMP parallel do` is **host** code (CPU threads), easily misread as an offload directive.

Placement:
- Put the transfer next to the producer, covering exactly the arrays the consumer reads. A
  transfer of a different array, or under a different condition, does not cover yours.
- **Each branch of an if/elseif chain needs its own transfer.** In `tracer_hordiff` the device
  branch computing `khdt_x`/`Kh_u` had no `update from`, while its two host siblings had
  `update to`, so later host code read unwritten memory (`979be73e6`).
- **A missing `update from` can leave answers bitwise identical.** If the stale values reach only
  control flow (an iteration count, a limiter), the run misbehaves instead: the `khdt_x` bug made
  one timestep run 89,159 iterations. Matching `ocean.stats` does not clear a transfer.
- Diagnostics: one `update from(a, b) if (CS%debug .or. CS%id_a > 0 .or. CS%id_b > 0)` covering
  all consumers, then the individual `if (id > 0) call post_data(...)`.
- At coarse sync points in `step_MOM` a blanket transfer of the state is the norm; match it rather
  than adding per-field ones.
- Passing a derived-type member (`CS%tv%T`) into a halo exchange or kernel without mapping it
  makes the runtime issue many small implicit transfers. Bracket the call with an explicit
  `map(to: CS%tv, CS%tv%T, CS%tv%S)` / `map(from: ...)`, as around `post_diabatic_halo_updates`
  in `step_MOM`.

Halo exchanges on device-resident fields:
- **Preferred:** a group pass. Register the fields once on a `group_pass_type` in the CS
  (`create_group_pass(CS%pass_x, field, G%Domain, halo=...)`, several fields on one handle for one
  message), and exchange them in place with `do_group_pass(CS%pass_x, G%Domain, omp_offload=.true.)`.
  The fields must stay mapped while the pass is used.
- `pass_var` and `pass_vector` have no `omp_offload` argument and are host-only. Bracket them with
  `update from` before and `update to` after, or convert them to an offloaded group pass.
- `start_group_pass` / `complete_group_pass` (the `NONBLOCKING_UPDATES` path) take
  `omp_offload=.true.` like `do_group_pass`. Pass it to **both** calls of a pair (FMS stops with a
  fatal error if they differ) and drop the `update from`/`update to` around them. Before dropping
  them, check that nothing between the start and the complete reads or writes the group's fields on
  the host, and that no device kernel writes them there (the old `update to` would have overwritten
  that). The device copy of the message buffer is created in start and released at the end of
  complete, so a start without its complete leaks it.
- Size `halo=` to the stencil of the kernels that read the field (`dynamics_split_RK2` computes
  `cor_stencil`, `vel_stencil` and `cont_stencil` for this). If a ported kernel widens a stencil,
  widen the matching pass.

## 6. Hazards that pass every benchmark

- **Unbalanced enter/exit.** After editing a map list, grep the routine for the matching exit.
  A leaked device allocation per call (`15ca2a25f`) shows no symptom until memory runs out.
- **`enter data` above an early `return`.** When the return fires, the entry leaks; the freed
  host block is then reused by another routine's automatic, whose own `map(alloc:)` fails with
  `variable in data clause is partially present` in a module you never touched. Read the
  present-table dump for the stale entry and its `line:`.
- **Hoisting a transfer to the caller while the callee still initializes on the host.** A
  `map(alloc:)` at the call site freezes the device copy; a host `cTKE(:,:,:) = 0.0` left inside
  the callee never reaches it, and a kernel accumulating into `cTKE` reads garbage. Initialize on
  the device instead. Also remove any `exit data map(release:)` left in the callee for data the
  caller now owns.
- **Persistently mapped arrays wrapped in a local enter/exit pair.** `exit data map(from:)` only
  decrements the count, so nothing comes back and the host copy stays stale. Symptom: one output
  is all zeros while the kernel's other outputs are right. Use `update from`.

## 7. Auditing a routine

Track, for each array, whether the host copy, the device copy, or both are current, walking the
code in execution order, including both sides of every branch:

| Event | Host | Device |
|---|---|---|
| `enter data map(to:)` (not already present) | current | current |
| `enter data map(alloc:)` | current | **garbage** |
| host write (`do`, `!$OMP parallel do`, host call) | current | **stale** |
| device write (`do concurrent`, `!$omp target`) | **stale** | current |
| `update to` / `update from` | current | current |
| `exit data map(delete:/release:)` | unchanged, possibly stale | gone |

A host read while the host copy is stale needs an `update from` before it; a device read while
the device copy is stale or garbage needs an `update to`, or a device-side initialization. When
the write and the read sit in sibling branches, the condition that reaches the read without the
transfer is the bug report.
