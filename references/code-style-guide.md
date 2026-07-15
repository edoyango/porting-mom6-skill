# MOM6 code style guide (summary)

Full source: https://github.com/mom-ocean/MOM6/wiki/Code-style-guide

Ported code should follow the same conventions as the rest of MOM6. Highlights
most relevant to GPU porting work:

- **Indentation**: 2 spaces, consistent within a block. Continuation lines
  indent at least 4 spaces (extra space is fine to aid alignment).
- **Line length**: aim for 100 characters for code (excluding comments); never
  exceed 120 characters including comments.
- **Loop index letters follow the staggering convention**: lowercase `i`,`j`
  for cell-center references; uppercase `I`,`J` for cell-edge/staggered
  references (north-east staggering, so `I` means `i+1/2`); uppercase `K` for
  interface-level references (`K` means `k-1/2`). Keep this convention when
  writing new `do concurrent` loop bounds — don't introduce new index letters.
- **Local variable declarations** go after dummy-argument declarations,
  separated by a `! Local variables` comment. Prefer descriptive
  multi-character snake_case names; any short/abbreviated name must be
  commented, and real variables should have units in the comment. This applies
  to new locals introduced when restructuring a loop for GPU offload (e.g. a
  tile-loop index or a temporary introduced for tiling).
- **Block termination**: `do`/`enddo` and `if`/`endif` (not `end do`/`end if`);
  other constructs (`subroutine`, `function`, `module`, `type`, ...) use a
  separated `end <token>`.
- **No module-level global data**, no implicit typing (`implicit none ;
  private` in every module), and all module `use` statements must carry `,
  only:`. These constraints shouldn't change when a subroutine is ported, but
  worth double-checking if a refactor moves code between modules.
- **Array syntax**: whole-array assignment and identical-copy syntax
  (`tv%S(:,:,:) = 0.`, `S_tmp(:,:,:) = tv%S(:,:,:)`) are fine, but bare
  scalar-style assignment without colons (`tv%S = 0.`) is not permitted, and
  array-syntax math that touches halos is disallowed because halo values
  aren't guaranteed valid. This matters when rewriting a loop nest as an array
  expression instead of a `do concurrent` — prefer `do concurrent` over
  whole-array math expressions for exactly this reason.
- **Arithmetic reproducibility**: additions should be explicitly paired with
  parentheses (`(a + b) + c`, not `a + b + c`); avoid `sum()`, `prod()`, and
  `matmul()` intrinsics in favor of explicit loops with a defined summation
  order; avoid the `**` exponent operator for integer powers (use repeated
  multiplication) and for cube roots (use the MOM6 `cuberoot` intrinsic).
  **This is directly relevant to GPU porting**: reordering a reduction to fit
  a `do concurrent` or `target teams loop` construct must preserve the
  original summation order, or the ported loop will produce answers that
  differ from the CPU / `dev/gpu` baseline even though the physics is
  unchanged. When collapsing or reassociating a loop for offload, check that
  the reassociation doesn't change floating-point order of operations for any
  sum, product, or exponentiation in the loop body.

For anything not covered above (whitespace details, documentation
conventions, etc.), check the linked wiki page directly rather than guessing.
