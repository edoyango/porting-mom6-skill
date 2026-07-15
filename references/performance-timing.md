# Assessing performance with MOM_cpu_clock

MOM6 already has a lightweight timing framework in `MOM_cpu_clock`, used
throughout the model to measure how long a given component takes. This is
the primary tool for checking whether a port has held CPU performance
steady (requirement (c) in `SKILL.md`'s "Preserving CPU and GPU performance
together" section) — it's the same clock mechanism, in fact, that shows up
in the Rule 3 / Rule 4 code examples in `data-mapping-conventions.md`, where
`cpu_clock_begin`/`cpu_clock_end` already bracket the mapping directives and
the call site.

## Step 1: check whether a clock already exists for the component you're porting

Two signs indicate one does:

- A `call cpu_clock_begin(id_clock_X)` / `call cpu_clock_end(id_clock_X)`
  pair already wraps the call site of the subroutine you're porting, in its
  caller.
- The end-of-run timing summary (typically near the end of the model's
  stdout/log output) lists a timer whose name matches, or clearly
  corresponds to, the subroutine or component you're porting.

If both are true, you already have a timing baseline. Capture that number
from a run before your change lands, then compare it to the same clock's
value after the port to confirm CPU timing hasn't regressed.

## Step 2: if no clock exists yet, define one

- Declare an integer clock id — module-level or as a `*_CS` field, depending
  on how similar existing clock ids are declared in that module.
- Register it in the module's `_init` routine with `cpu_clock_id(...)`,
  giving it a descriptive name so it's identifiable in the end-of-run
  summary.
- Wrap the call site of the subroutine you're porting (in the caller) with
  `cpu_clock_begin(id_clock_X)` / `cpu_clock_end(id_clock_X)`. This is
  usually the same clock region that already brackets (or should bracket)
  the `target update to/from` directives from Rules 3 and 4 in
  `data-mapping-conventions.md` — one region, one clock, covering both the
  data movement and the call.

## Step 3: use it

Compare the clock's CPU-only timing before and after the port to confirm (c)
hasn't regressed. If you also want a sense of GPU speedup, the same
clock/timer name recorded from a GPU run gives you a like-for-like
comparison against the CPU baseline.

This is a per-component, fine-grained check. It complements, rather than
replaces, whatever the target-configuration regression-testing scripts turn
out to check once that scripts directory is populated (see `SKILL.md`).
