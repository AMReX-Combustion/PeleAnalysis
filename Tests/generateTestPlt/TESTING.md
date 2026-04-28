# generateTestPlt Test Suite

Tests for `Src/generateTestPlt.cpp` — the synthetic AMReX plotfile generator.

## Running

```bash
cd PeleAnalysis/Tests/generateTestPlt
./run_tests.sh              # build + run
./run_tests.sh --no-compile # skip build, reuse existing executable
```

Output artefacts (plotfiles, build logs) are written to `testrun/` and are not tracked by git.

## Test matrix

### Phase 2 — Structural (single-level)

| Test | Input | What is checked |
|------|-------|-----------------|
| T1   | `gen_t1_all_fields.inp` | nComp=10 in Header; all 10 variable names present |

### Phase 3 — Value correctness (single-level)

| Test | Input | What is checked |
|------|-------|-----------------|
| T2 | `gen_t2_const.inp` | All cells of a `constant` field equal 7.0 |
| T3 | `gen_t3_degenerate.inp` | Five degenerate configurations all produce 5.0 everywhere: direct constant, sphere far outside domain, ring with radii outside domain, double-plane with positions outside domain, sine with zero amplitude |
| T4 | `gen_t4_planes.inp` | `plane_step` on x, y, z axes each produces range exactly [0, 1] |

### Phase 4 — AMR structural

| Test | Input | What is checked |
|------|-------|-----------------|
| T5  | `gen_t5_amr_inbox.inp`   | `in_box_lo/hi` criterion → 2 levels |
| T6  | `gen_t6_amr_vgt.inp`     | `value_greater` criterion → 2 levels |
| T7  | `gen_t7_amr_vlt.inp`     | `value_less` criterion → 2 levels |
| T8  | `gen_t8_amr_adj.inp`     | `adjacent_difference_greater` criterion → 2 levels |
| T9  | `gen_t9_amr_3level.inp`  | `max_level=2`, InBox → 3 levels |
| T10 | `gen_t10_amr_refrat4.inp`| `ref_ratio=4`, InBox → 2 levels |
| T11 | `gen_t11_amr_notags.inp` | Threshold above any field value → no tags → 1 level only (no Level_1 directory) |

### Phase 5 — AMR value correctness

| Test | Input | What is checked |
|------|-------|-----------------|
| T12 | `gen_t12_amr_const.inp` | `constant` field value 3.14 is identical on every AMR level |

### Phase 6 — `output_name` (issue #57)

| Test | Input | What is checked |
|------|-------|-----------------|
| T13 | `gen_t13_output_name.inp` | `myvar.output_name = Y(OH)` writes `Y(OH)` (not `myvar`) to the plotfile Header |

This addresses the second gap reported in issue #57: variable names containing `/` or other ParmParse-reserved characters can now be produced by using a safe ParmParse prefix and overriding the plotfile name via `output_name`.

### Phase 7 — Error handling

| Test | Input | What is checked |
|------|-------|-----------------|
| err1 | `gen_err1_no_fields.inp`    | Tool aborts when `field.names` is missing |
| err2 | `gen_err2_bad_type.inp`     | Tool aborts on unrecognised field type string |
| err3 | `gen_err3_bad_refname.inp`  | Tool aborts when a refinement indicator references a non-existent field |
| err4 | `gen_err4_no_criterion.inp` | Tool aborts when a refinement indicator has no recognised criterion key |

### Phase 8 — Cylindrical coordinate system

| Test | Input | What is checked |
|------|-------|-----------------|
| T14 | `gen_t14_cylinder.inp` | `cylinder_step` on all three axes and `cylinder_smooth` each produce range [0, 1] (3D serial) |
| T15 | `gen_t15_rz.inp` | `geometry.coord_sys = 1` accepted; Header records `spacedim=2` and `coord_sys=1`, matching PeleLMeX 2D RZ output (**2D serial** binary required) |

T14 uses 3D Cartesian coordinates and exercises `cylinder_step`/`cylinder_smooth` (infinite cylinder defined by radial distance from an axis line).  T15 requires the 2D binary (`generateTestPlt2d.gnu.ex`) and is skipped if that executable is absent.

## Python helpers

| Script | Purpose |
|--------|---------|
| `check_plt_value.py`   | Assert every cell of a named variable equals an expected value (within tolerance) across all AMR levels |
| `check_plt_range.py`   | Assert the observed min and max of a variable are within tolerance of expected values across all AMR levels |
| `check_plt_nlevels.py` | Assert the number of `Level_N/` directories in a plotfile equals an expected count |
| `check_plt_coord.py`   | Assert the `coord_sys` integer in the plotfile Header equals an expected value (0 = Cartesian, 1 = cylindrical/RZ) |
