# jpdf Test Suite

Tests for `jpdf.cpp` using synthetic plotfiles generated with `generateTestPlt`.

---

## Running the test suite

The fastest way to run everything is with the automated script:

```bash
cd Src/TestFiles
./run_tests.sh            # compile then run all tests
./run_tests.sh --no-compile  # skip build, use existing executables
```

The script compiles 2D/3D serial and 3D MPI executables, generates synthetic
plotfiles, runs all tests, and prints a PASS/FAIL summary. All output artefacts
are written to `TestFiles/testrun/` so the source tree stays clean. MPI tests
require `mpirun` or `mpiexec` and are skipped automatically if neither is found.

---

## Manual setup

### 1. Compile tools

Executables are named `<EBASE><DIM>d.<COMP>[.MPI].ex` by AMReX.

```bash
cd Src/
make EBASE=generateTestPlt DIM=2 DEBUG=FALSE USE_MPI=FALSE COMP=gnu   # generateTestPlt2d.gnu.ex
make EBASE=generateTestPlt DIM=3 DEBUG=FALSE USE_MPI=FALSE COMP=gnu   # generateTestPlt3d.gnu.ex
make EBASE=jpdf            DIM=2 DEBUG=FALSE USE_MPI=FALSE COMP=gnu   # jpdf2d.gnu.ex
make EBASE=jpdf            DIM=3 DEBUG=FALSE USE_MPI=FALSE COMP=gnu   # jpdf3d.gnu.ex
make EBASE=jpdf            DIM=3 DEBUG=FALSE USE_MPI=TRUE  COMP=gnu   # jpdf3d.gnu.MPI.ex
```

### 2. Generate test plotfiles

From `Src/`:

```bash
./generateTestPlt3d.gnu.ex TestFiles/generate_plt_A.inp
./generateTestPlt3d.gnu.ex TestFiles/generate_plt_B.inp
```

This creates `plt_jpdf_A/` and `plt_jpdf_B/` in the working directory. Both contain three fields on a 16×16×16 grid in [0,1]³:

| Field | Description | Values |
|-------|-------------|--------|
| `var_plane` | Hard step at x=0.5 | 0.0 (x<0.5) or 1.0 (x≥0.5) |
| `var_sine` | sin(2π x)×0.5 + 0.5 | continuous, range ≈ [0, 1] |
| `var_cond` | Hard step at y=0.5 | 0.0 (y<0.5) or 1.0 (y≥0.5) |

`plt_jpdf_B` is identical to `plt_jpdf_A` and is used for the averaging test.

### 3. Run a test

```bash
./jpdf2d3d.gnu.ex TestFiles/test1_all_outputs.inp
```

---

## Test Matrix

| # | Input file | Feature(s) tested | Pass criterion |
|---|-----------|-------------------|----------------|
| 1 | `test1_all_outputs.inp` | All six output formats | Six output file types present in `plt_jpdf_A_test1/`; PDF non-zero only in columns for var_plane=0 and var_plane=1 |
| 2 | `test2_condmean.inp` | 2D conditional mean | `condMean_var_cond_on_var_plane_var_sine.dat` created; all non-zero entries ≈ 0.5 |
| 3 | `test3_condmean_duplicate.inp` | Duplicate in condMean_vars | Tool runs without error; verbose prints "already on position 2"; no variable loaded twice |
| 4 | `test4_useminmax.inp` | `useminmax` range override + clamping | Verbose reports `v1g > 0`; plotfile Header axis min/max = -0.5 and 0.5 |
| 5 | `test5_conditioning_1.inp` | `do_conditioning=1` | JPDF non-zero only at var_plane=1 column; sum(PDF) = 1.0 |
| 6 | `test6_conditioning_2.inp` | `do_conditioning=2` (c(1-c)) | JPDF zero for var_sine near 0 and 1; sum(PDF) = 1.0 |
| 7 | `test7_normcval.inp` | `norm_cVal=1` | Only back-half (var_cond=1) cells contribute; sum(PDF) = 1.0 |
| 8 | `test8_average.inp` | `do_average=1` | `JPDFAverage_test8/` created; averaged PDF matches per-file PDF |

---

## Detailed Expected Outputs

### Test 1 — All output formats

Output directory: `plt_jpdf_A_test1/`

| File | Description |
|------|-------------|
| `Header` + `Level_0/Cell*` | AMReX plotfile with 2 components (PDF and log-PDF) |
| `Pdf_var_plane_var_sine.gpd` | 3-column gnuplot file: v1 v2 pdf |
| `Pdf_var_plane_var_sine.dat` | 16×16 MATLAB matrix (row=v1i, col=v2i) |
| `Pdf_var_plane_x.dat` | 16 bin-centre values for var_plane axis |
| `Pdf_var_sine_x.dat` | 16 bin-centre values for var_sine axis |
| `PdfX1_var_plane_var_sine.dat` | 16×16 volume-weighted mean of var_plane per bin |
| `PdfX2_var_plane_var_sine.dat` | 16×16 volume-weighted mean of var_sine per bin |
| `Pdf_var_plane_var_sine.tpd` | Tecplot: N=256 nodes, E=225 quads |
| `Pdf_var_plane_var_sine.fab` | Binary FArrayBox (4 components: v1, v2, log-pdf, pdf) |
| `Scatter_var_plane_var_sine.dat` | 2-column x y for non-zero bins only |

Structural checks:
- PDF is non-zero only in bin column 0 (var_plane=0.0) and column 15 (var_plane=1.0 — clamped, see Pitfall 2).
- `sum(PDF matrix) = 1.0`.
- Log-PDF component: `log(pdf + 1e-7)` for each bin.

### Test 2 — 2D conditional mean

Output directory: `plt_jpdf_A_test2/`

- `condMean_var_cond_on_var_plane_var_sine.dat`: 16×16 matrix.
- `var_cond` is independent of x-axis (no correlation with var_plane or var_sine), so each y-slice of the domain contributes equally to every (var_plane, var_sine) bin. The mean of a 50/50 mix of 0.0 and 1.0 is **0.5**.
- Check: all entries in the condMean matrix that correspond to non-zero PDF bins should equal 0.5 ± floating-point noise.
- Empty bins (zero PDF): condMean value set to bin-centre of var_plane or var_sine (not 0), per the empty-bin fallback in the code.

### Test 3 — condMean duplicate detection

Output directory: `plt_jpdf_A_test3/`

- 3 JPDF pairs → 15 MATLAB PDF files + 3 condMean files:
  - `condMean_var_cond_on_var_plane_var_sine.dat`
  - `condMean_var_cond_on_var_plane_var_cond.dat`
  - `condMean_var_cond_on_var_sine_var_cond.dat`
- Tool must complete without error.
- Verbose line: `"condMeanVar 'var_cond' already on position 2"` confirms deduplication.

### Test 4 — useminmax + bin clamping

Output directory: `plt_jpdf_A_test4/`

- var_plane range overridden to [-0.5, 0.5]; auto-range would be [0.0, 1.0].
- Cells with var_plane=1.0 map to bin index `int(16 × 1.5/1.0) = 24`, clamped to 15.
- Verbose reports `v1g = 2048` (16×16×8 cells with var_plane=1.0 clamped above range).
- Plotfile `Header` last two lines (axis ranges) must show `-5.000...e-01 5.000...e-01` for var1.

### Test 5 — do_conditioning=1

Output directory: `plt_jpdf_A_test5/`

- cVar=0 (var_plane), cMin=0.9, cMax=1.1: only right half of domain (var_plane=1.0) contributes.
- JPDF matrix: non-zero only in column 15.
- `sum(PDF matrix) = 1.0` (normalization fix applied).
- The log-PDF component in the plotfile should show a real value (not just `log(1e-7)`) at bin 15.

### Test 6 — do_conditioning=2

Output directory: `plt_jpdf_A_test6/`

- cVar=1 (var_sine). The code normalises first, then applies c(1−c): with cNormMax=1.0, var_sine maps to [0,1]; c(1−c) then gives [0, 0.25]. cMin=0.1 keeps cells where c(1−c) ≥ 0.1, i.e., var_sine ∈ [~0.13, ~0.87].
- JPDF matrix: zero (or near-zero) in the bins corresponding to the extreme tails of var_sine.
- `sum(PDF matrix) = 1.0`.

### Test 7 — norm_cVal=1

Output directory: `plt_jpdf_A_test7/`

- cVar=2 (var_cond), norm_cVal=1, cMin=0.9, cMax=1.1: only back-half cells (var_cond=1.0) contribute.
- 3 JPDF pairs written; the (var_plane, var_sine) PDF covers both halves of x (both 0 and 1 values of var_plane present in back half).
- `sum(PDF)` per pair = 1.0.

### Test 8 — do_average

Output directories: `plt_jpdf_A_test8/`, `plt_jpdf_B_test8/`, `JPDFAverage_test8/`

- Both input plotfiles are identical.
- Averaged PDF = per-file PDF (to floating-point precision).
  - Accumulated: `binAv = bin_A + bin_B = 2 × bin_A`
  - Normalised: `binAv / (domainVol × 2) = bin_A / domainVol`
- Check: element-wise difference between `JPDFAverage_test8` PDF and `plt_jpdf_A_test8` PDF should be ≈ 0.

---

## Known Pitfalls

### Pitfall 1 — Normalization (fixed)

**Bug (before fix):** When `do_conditioning > 0`, only a subset of cells contribute to the PDF, but `bin[i] /= domainVol` divided by the *full* domain volume. This caused `sum(PDF) < 1`.

**Fix applied:** When `do_conditioning > 0`, the normalization factor is replaced by the total accumulated contributed volume (`sum(bin[i])` before division). This makes `sum(PDF) = 1` regardless of how many cells pass the conditioning mask.

### Pitfall 2 — vMax clamping (expected behaviour)

When a cell value equals `vMax` exactly, `v1i = int(nBins × 1.0) = nBins`, which is clamped to `nBins - 1`. The verbose output reports this as `v1g > 0`. This is intentional (values at the exact maximum end up in the last bin) and is not a bug. Tests 1, 4, 5 will all show non-zero `v1g` for `var_plane`.

### Pitfall 3 — Empty bin fallback

When a bin has no cells (`div == 0`), the volume-weighted mean arrays `binX1` and `binX2` are set to the bin-centre coordinates `v1` and `v2` (not 0). Conditional mean arrays are left at their initialised value of 0. This is important when interpreting MATLAB output: a condMean value of 0 in an empty bin does not mean the conditioning variable is zero there — the bin simply received no data.

### Pitfall 4 — log-PDF small offset

The log-PDF is computed as `log(pdf + 1e-7)`. For bins with PDF = 0, the stored log value is `log(1e-7) ≈ -16.1`. This hardcoded floor compresses the log-space dynamic range when the true PDF minimum is very small but non-zero.

### Pitfall 5 — Slash-in-variable-name (ProtectSlashes)

Variables with `/` in their name (e.g. `Y(OH)`) have slashes replaced by `_` in output filenames via `ProtectSlashes`. This cannot be tested with `generateTestPlt` because ParmParse interprets `.` and `/` as key separators. Testing requires a hand-crafted AMReX plotfile where the variable name string contains a slash.

---

## Reading Output Files

**AMReX plotfile:** Use `amrvis2d` or Python `yt`:
```python
import yt
ds = yt.load("plt_jpdf_A_test1")
```
The plotfile has 2 components: `Pdf_var_plane_var_sine` and `Pdf_var_plane_var_sine (log)`. The Header encodes the variable axis ranges in the last `nVars_jpdf` lines.

**MATLAB `.dat` (PDF matrix):** Row index = v1i (first variable), column index = v2i (second variable). The axis ranges come from the corresponding `Pdf_<var>_x.dat` files.

**gnuplot `.gpd`:** Three columns: v1 v2 pdf. Plot with:
```gnuplot
set pm3d; splot 'Pdf_var_plane_var_sine.gpd' with pm3d
```

**Tecplot `.tpd`:** Quad mesh; open directly in Tecplot or ParaView.
