# Testing `analysis_util`

## Overview

This test suite validates the generic utility functions in `PeleAnalysis/Tools/Util/`
(`analysis_util` namespace). Tests are exercised via a C++ executable built from
`test_analysis_util.cpp` and driven by `run_tests.sh`.

## Running

```
cd PeleAnalysis/Tests/analysis_util
./run_tests.sh            # build + run all tests (serial and MPI)
./run_tests.sh --no-compile   # skip build, reuse existing executables
```

Artifacts are written to `testrun/`. Python 3 is required for plotfile value
assertions. MPI tests are skipped if `mpirun`/`mpiexec` is not found.

## Build system

The test executable is built from `test_analysis_util.cpp` in this directory using
the one-line `GNUmakefile` which includes `../../Src/GNUmakefile` verbatim.
`PELE_ANALYSIS_HOME` is overridden from the command line to fix the relative path
in `Src/GNUmakefile`. `Blocs := .` in `Src/GNUmakefile` expands to the current
directory, so the test source is found automatically. No build logic is duplicated.

## Test groups

### Group 1 — String utilities

| ID | What | Input → Expected |
|---|---|---|
| T-STR-1 | `get_file_root` with path | `"path/to/plt"` → `"plt"` |
| T-STR-2 | `get_file_root` bare filename | `"plt"` → `"plt"` |
| T-STR-3 | `get_file_root` root-relative | `"/plt"` → `"plt"` |
| T-STR-4 | `find_var_index` middle | `{"a","b","c"}`, `"b"` → `1` |
| T-STR-5 | `find_var_index` first | `"a"` → `0` |
| T-STR-6 | `find_var_index` last | `"c"` → `2` |
| T-STR-7 | `find_var_index` not found, no abort | `"d"` → `-1` |
| T-STR-8 | `find_var_index` exact match | `{"var","var_extra"}`, `"var"` → `0` |
| T-STR-9 | `parse_title` + `parse_var_names` (space) | `"My Title\nvar1 var2 var3"` |
| T-STR-10 | `parse_var_names` comma+space | `"var1, var2, var3"` → 3 tokens |
| T-STR-A1 | `find_var_index` abort on missing (subprocess) | non-zero exit |

### Group 2 — AMReX plotfile round-trip

All tests create MultiFabs programmatically, write with `write_plotfile`, and
read back with `read_plotfile`.

| ID | What | Verification |
|---|---|---|
| T-PLT-1 | Single-variable constant=3.14 | All cells == 3.14 (`check_plt_value.py`) |
| T-PLT-2 | Two-variable file, load one | Correct component value |
| T-PLT-3 | `PlotfileData.var_names` has all file vars | `size()==2` for 2-var file |
| T-PLT-4 | `finest_level` cap | `n_lev==1` when capped at 0 |
| T-PLT-5 | `n_grow=1` ghost cells | `nGrowVect() == IntVect(1)` |
| T-PLT-6 | `is_per={1,0,0}` periodicity | `isPeriodic(0)==true`, `isPeriodic(1)==false` |
| T-PLT-7 | Default `ref_ratios` (empty) | Single-level write succeeds |
| T-PLT-8 | Header nComp matches var count | `sed -n '2p' Header == "2"` |
| T-PLT-9 | Returned data works with `init()` | `initialized==true` after call |
| T-PLT-A1 | Unknown variable → abort (subprocess) | non-zero exit |
| T-PLT-A2 | Non-existent file → abort (subprocess) | non-zero exit |

### Group 3 — MEF file round-trip

Tests create a `MEFData` struct programmatically (4 nodes, 4 components, 2 triangles),
write with `write_mef`, read back with `read_mef`.

| ID | What | Verification |
|---|---|---|
| T-MEF-1 | Title preserved | `data_out.title == data_in.title` |
| T-MEF-2 | `var_names` preserved | Count and each name equal |
| T-MEF-3 | `n_elts`/`nodes_per_elt` preserved | Both equal after round-trip |
| T-MEF-4 | Node values preserved | All cells match within 1e-12 |
| T-MEF-5 | Connectivity data preserved | All integers equal |
| T-MEF-6 | Single-triangle edge case | `n_elts==1`, `connectivity.size()==3` |
| T-MEF-7 | Multi-component nodes (nComp=6) | `nComp==6` after round-trip |
| T-MEF-8 | Title with spaces | `"My surface title with spaces"` round-trips |
| T-MEF-A1 | Non-existent file → abort (subprocess) | non-zero exit |

### Group 4 — `get_covered_mf`

Tests verify the fine-covered mask after `get_covered_mf()` was relocated from
`analysis_util_integrate.cpp` to `analysis_util.cpp`.

| ID | What | Verification |
|---|---|---|
| T-COV-1 | Single-level: all uncovered | `mask[0].sum() == 8^3 = 512` |
| T-COV-2 | Two-level: fine-covered coarse cells = 0 | Sum == `512 - 4^3 = 448` |
| T-COV-4 | `integrate()` still resolves `get_covered_mf` | No link error; no crash |

### Group 5 — `integrate` regression

Confirms integration behavior is unchanged after the `get_covered_mf` relocation.

| ID | What | Verification |
|---|---|---|
| T-INT-1 | All-axis integral of const=3 field | Result ≈ `3.0 * 1.0 = 3.0` (domain vol = 1) |
| T-INT-2 | X-axis integral only | Result size = `ny*nz = 64`; each ≈ `0.375` |

### MPI tests (Phase 6)

| ID | np | Verification |
|---|---|---|
| T-PLT-MPI-1 | 2, 4 | Plotfile values match serial (via `check_plt_value.py`) |
| T-MEF-MPI-1 | 2, 4 | MEF file written only by IOProcessor; file parseable after run |
| T-INT-MPI-1 | 2, 4 | T-INT-1 passes on MPI run |

## Known pitfalls

1. **`get_file_root` bare filename** — `Tokenize("name", "/")` returns a single-element
   vector; index `[size-1]` == `[0]` == `"name"`. Covered by T-STR-2.

2. **Stream position after `>> n_elts >> nodes_per_elt`** — `>>` leaves the stream
   after the numbers but before the newline. `FArrayBox::readFrom` calls `>>` internally
   which skips leading whitespace. Covered by T-MEF-3.

3. **`write_plotfile` single-level `ref_ratios` size** — `std::max(n_lev-1, 0)` in
   `write_plotfile` ensures the vector is empty (not negative) for single-level data.
   Covered by T-PLT-7.

4. **`find_var_index` exact match** — must not match `"var"` against `"var_extra"`.
   Covered by T-STR-8.

5. **`nlev` initialization** — `int nlev = -1` (not zero) is required for `init()` to
   correctly set `nlev` from the MultiFab vector. Covered by T-COV and T-INT which call
   `init()` and rely on `nlev` being correct.

6. **`write_mef` IOProcessor guard** — without `if (IOProcessor()) return;`, all MPI
   ranks try to open and write the same file simultaneously. Covered by T-MEF-MPI-1.

7. **`get_covered_mf` call from `integrate`** — after moving the implementation to
   `analysis_util.cpp`, the function must still resolve correctly from
   `analysis_util_integrate.cpp` (same namespace). Verified at link time by T-COV-4
   and T-INT-1.
