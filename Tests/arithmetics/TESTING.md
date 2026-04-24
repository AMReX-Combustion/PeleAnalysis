# arithmetics Test Suite

Tests for `arithmetics.cpp` using synthetic constant-field plotfiles generated
with `generateTestPlt`.  Value correctness is verified with a round-trip
technique: applying inverse operations and reading the result with the included
`check_plt_value.py` helper.

---

## Running the test suite

```bash
cd Tests/arithmetics
./run_tests.sh               # compile + run all tests
./run_tests.sh --no-compile  # skip build, use existing executables in Src/
```

The script compiles 3D serial and 3D MPI executables plus `generateTestPlt`,
generates two synthetic plotfiles, and runs all tests.  All artefacts land in
`testrun/`.  MPI tests are skipped automatically when `mpirun`/`mpiexec` is
absent.  Value-correctness tests require Python 3.

---

## Manual setup

### Compile

```bash
cd Src/
make EBASE=generateTestPlt DIM=3 DEBUG=FALSE USE_MPI=FALSE COMP=gnu   # generateTestPlt3d.gnu.ex
make EBASE=arithmetics     DIM=3 DEBUG=FALSE USE_MPI=FALSE COMP=gnu   # arithmetics3d.gnu.ex
make EBASE=arithmetics     DIM=3 DEBUG=FALSE USE_MPI=TRUE  COMP=gnu   # arithmetics3d.gnu.MPI.ex
```

### Generate test plotfiles

From `Tests/arithmetics/testrun/` (or any working directory):

```bash
generateTestPlt3d.gnu.ex ../generate_plt_arith.inp
generateTestPlt3d.gnu.ex ../generate_plt_arith_divzero.inp
```

| Plotfile | Fields | Values |
|---|---|---|
| `plt_arith` | `varA`, `varB` | 4.0, 2.0 — constant everywhere |
| `plt_arith_divzero` | `varA`, `varB` | 4.0 constant; `varB` = plane_step (0.0 for x<0.5, 2.0 for x≥0.5) |

Both are 8×8×8 grids on [0,1]³ with Cartesian coordinates and no AMR.  The
integer values (4 and 2, both exact in IEEE 754) make round-trip results
bit-exact, so round-trip comparisons use zero tolerance.

---

## Test Matrix

| # | Input file | Feature tested | Pass criterion |
|---|---|---|---|
| 1 | `test1_add.inp` | `operator=add` | exit 0; `plt_t1_add/` exists; Header nComp=3; "sum" in Header |
| 2 | `test2_subtract.inp` | `operator=subtract` | exit 0; `plt_t2_subtract/` exists; Header nComp=3; "diff" in Header |
| 3 | `test3_multiply.inp` | `operator=multiply` | exit 0; `plt_t3_multiply/` exists; Header nComp=3; "prod" in Header |
| 4 | `test4_divide.inp` | `operator=divide` (no zeros) | exit 0; `plt_t4_divide/` exists; Header nComp=3; "quot" in Header |
| 5 | `test5_divzero_check.inp` | divide + zeros, `checkDivByZero=1` | exits **non-zero**; no output dir |
| 6 | `test6_divzero_skip.inp` | divide + zeros, `checkDivByZero=0` | exit 0; output dir exists |
| 7 | inline args | default output naming | `plt_arith_add/` created |
| 7b | `test7_outfile_override.inp` | `outfile` override | `plt_custom_output/` created |
| 8 | `test8a/b_roundtrip_add.inp` | add / subtract inverse | `recovered` = 4.0 ± 0 in all cells |
| 9 | `test9a/b_roundtrip_mul.inp` | multiply / divide inverse | `recovered` = 4.0 ± 0 in all cells |
| 10 | inline args | `inVarAName` not in plotfile | exits non-zero |
| 11 | inline args | `inVarBName` not in plotfile | exits non-zero |
| 12 | inline args | invalid `operator` | exits non-zero |
| MPI-2 | inline args | add with 2 MPI ranks | `sum`=6.0, `varA`=4.0 in all cells |
| MPI-4 | inline args | add with 4 MPI ranks | `sum`=6.0, `varA`=4.0 in all cells |
| MPI-div | `test5_divzero_check.inp` | divide + zeros, 4 MPI ranks | all ranks abort (non-zero exit) |

---

## Detailed Expected Outputs

### Tests 1–4 — Operator structural checks

Output directory: `plt_t{N}_{operator}/`

```
plt_t1_add/
├── Header         ← text: line 2 = "3", lines 3–5 = "varA", "varB", "sum"
└── Level_0/
    ├── Cell_H     ← FAB layout
    └── Cell_D_00000   ← binary: 3 × 512 doubles (8×8×8 grid)
```

The Header line 2 (number of components) must equal `nCompIn + 1 = 3`.
All three variable names must appear on separate lines immediately after.
Original variables (`varA`, `varB`) are copied unchanged; the derived component
is appended last.

### Test 5 — Divide-by-zero detection

`plt_arith_divzero` has `varB = 0` on the left half of the domain (64 of 512
cells).  With `checkDivByZero=1` the `ReduceOps<ReduceOpLogicalOr>` reduction
detects the zeros across all MPI ranks and calls `amrex::Abort()`.  The tool
must exit non-zero and must **not** write an output plotfile.

### Test 6 — Divide-by-zero bypass

With `checkDivByZero=0` the zero check is skipped entirely.  `MultiFab::Divide`
in AMReX does not guard against zeros; cells on the right half (varB=2) produce
`4/2=2`, cells on the left half (varB=0) produce IEEE infinity or NaN depending
on varA.  The tool must complete and write an output plotfile.

### Tests 7 / 7b — Output naming

Without `outfile`, the default name is `<getFileRoot(infile)>_<operator>`.  For
`infile=plt_arith` and `operator=add` this is `plt_arith_add`.  With
`outfile=plt_custom_output` the default is overridden.

### Tests 8–9 — Round-trip value verification

**Add/subtract inverse** (`varA=4, varB=2`):
```
pass 1:  sum  = varA + varB     = 4 + 2 = 6   (stored in plt_t8a_add)
pass 2:  recovered = sum - varB = 6 - 2 = 4   (stored in plt_t8b_sub)
```
`recovered` must equal `varA = 4.0` in every cell.  IEEE 754 guarantees this
is bit-exact for integers representable as `double`.

**Multiply/divide inverse** (`varA=4, varB=2`):
```
pass 1:  prod      = varA * varB    = 4 * 2 = 8   (stored in plt_t9a_mul)
pass 2:  recovered = prod / varB    = 8 / 2 = 4   (stored in plt_t9b_div)
```
Same reasoning: exact for small integer powers of 2.

Value verification uses `check_plt_value.py`, which reads the AMReX binary FAB
files directly (component-major double storage, endian-detected from FAB header)
and checks every cell against the expected value.

### Tests 10–12 — Error handling

Each passes an invalid parameter combination to the serial executable and
expects a non-zero exit code (AMReX `Abort()`):

| Test | Bad input | Expected message |
|---|---|---|
| 10 | `inVarAName=no_such_var` | `Variable no_such_var not found in file plt_arith` |
| 11 | `inVarBName=no_such_var` | `Variable no_such_var not found in file plt_arith` |
| 12 | `operator=modulo` | `operator must be one of: add, subtract, multiply, divide` |

### MPI tests

The MPI divide-by-zero test (MPI-div) verifies that **all** ranks abort when
zeros are detected.  Before the GPU-compatible refactor, `hasZero` was a local
host variable never reduced across ranks; the abort would only fire on ranks
that owned cells with zero values, causing an MPI deadlock.  The test passes
only if `mpirun` itself returns non-zero (i.e. `MPI_Abort` propagated to all
ranks).

---

## `check_plt_value.py` — FAB reader

The script reads an AMReX plotfile and checks that every cell of a named
variable equals a target value.

```bash
python3 check_plt_value.py <plt_dir> <var_name> <expected> [<tolerance>]
```

It handles:
- Multi-level AMR (iterates over all `Level_N/` directories)
- Multiple FABs per level (MPI decomposition; reads all `FabOnDisk` entries from `Cell_H`)
- Little- and big-endian FAB data (detected from the FAB header line)

Exit codes: `0` = all values within tolerance, `1` = value mismatch, `2` = usage error.

---

## Known Pitfalls

### Pitfall 1 — Divide-by-zero MPI correctness (fixed)

Before the `ReduceOps` refactor, `hasZero` was local to each MPI rank and not
reduced.  A zero on rank 2 would cause only rank 2 to abort while other ranks
continued, resulting in an MPI deadlock.  `ParallelDescriptor::ReduceIntMax`
now ensures all ranks agree before the assertion.

### Pitfall 2 — `refRatios` hardcoded (fixed)

The plotfile write originally used a hardcoded `{2,2,2}` refinement ratio for
all levels.  For AMR data with actual refinement ratios other than 2, the output
plotfile would be structurally incorrect.  The tool now reads
`amrData.RefRatio()[lev]` per level.

### Pitfall 3 — `getFileRoot` on empty `infileName`

If `infile` is not set (empty string), `Tokenize` on `/` returns an empty
vector and `tokens[tokens.size()-1]` is undefined behaviour.  AMReX
`ParmParse::get("infile", ...)` will abort before this path is reached, but
the underlying issue remains.

### Pitfall 4 — Divide-by-zero with `checkDivByZero=0`

IEEE 754 division by zero produces `+Inf` (for positive numerator) or `-Inf`
(negative) or `NaN` (0/0).  The tool does not warn when `checkDivByZero=0`
produces these values.  The output plotfile will contain IEEE special values
silently.
