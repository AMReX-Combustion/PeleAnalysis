#!/usr/bin/env bash
# run_tests.sh — build and run all generateTestPlt tests (3D serial + MPI)
#
# Usage:  ./run_tests.sh [--no-compile]
#
#   --no-compile   Skip the build step; assume executable already exists in Src/
#
# All test artefacts are written to Tests/generateTestPlt/testrun/ so the
# source tree stays clean.  Python 3 is required for value-correctness and
# level-count assertions.

set -euo pipefail

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SRC_DIR="$SCRIPT_DIR/../../Src"
RUN_DIR="$SCRIPT_DIR/testrun"
NPROC=$(nproc 2>/dev/null || sysctl -n hw.ncpu 2>/dev/null || echo 4)

# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------
COMPILE=true
for arg in "$@"; do
    case "$arg" in
        --no-compile) COMPILE=false ;;
        *) echo "Unknown argument: $arg"; exit 1 ;;
    esac
done

# ---------------------------------------------------------------------------
# Terminal colours
# ---------------------------------------------------------------------------
RED='\033[0;31m'; GREEN='\033[0;32m'; YELLOW='\033[1;33m'
BOLD='\033[1m'; NC='\033[0m'

PASS=0; FAIL=0

section() { echo -e "\n${BOLD}=== $* ===${NC}"; }
pass()    { echo -e "  ${GREEN}PASS${NC} $*"; ((PASS++)) || true; }
fail()    { echo -e "  ${RED}FAIL${NC} $*"; ((FAIL++)) || true; }
skip()    { echo -e "  ${YELLOW}SKIP${NC} $*"; }

# ---------------------------------------------------------------------------
# Assertion helpers
# ---------------------------------------------------------------------------

check_file() {   # check_file <desc> <path>
    if [[ -f "$2" ]]; then pass "$1"; else fail "$1 — missing: $2"; fi
}

check_dir() {    # check_dir <desc> <path>
    if [[ -d "$2" ]]; then pass "$1"; else fail "$1 — missing dir: $2"; fi
}

no_dir() {       # no_dir <desc> <path>  — assert directory does NOT exist
    if [[ ! -d "$2" ]]; then pass "$1"; else fail "$1 — unexpected dir: $2"; fi
}

# header_ncomp_ok <desc> <plt_dir> <expected_ncomp>
header_ncomp_ok() {
    local actual
    actual=$(sed -n '2p' "$2/Header" 2>/dev/null || echo "MISSING")
    if [[ "$actual" == "$3" ]]; then
        pass "$1  (nComp = $actual)"
    else
        fail "$1  (expected nComp=$3, got '$actual')"
    fi
}

# header_has_var <desc> <plt_dir> <var_name>
header_has_var() {
    if grep -qx "$3" "$2/Header" 2>/dev/null; then
        pass "$1"
    else
        fail "$1 — '$3' not found in $2/Header"
    fi
}

# header_no_var <desc> <plt_dir> <var_name>
header_no_var() {
    if grep -qx "$3" "$2/Header" 2>/dev/null; then
        fail "$1 — '$3' unexpectedly found in $2/Header"
    else
        pass "$1"
    fi
}

# check_plt_value <desc> <plt_dir> <var_name> <expected> [<tol>]
check_plt_value() {
    local desc="$1" plt_dir="$2" var_name="$3" expected="$4" tol="${5:-1e-10}"
    if [[ ! -d "$plt_dir" ]]; then
        fail "$desc — plotfile dir missing: $plt_dir"; return
    fi
    if ! command -v python3 &>/dev/null; then
        skip "$desc — python3 not available"; return
    fi
    local out
    if out=$(python3 "$SCRIPT_DIR/check_plt_value.py" \
                "$plt_dir" "$var_name" "$expected" "$tol" 2>&1); then
        pass "$desc  ($out)"
    else
        fail "$desc — $out"
    fi
}

# check_plt_range <desc> <plt_dir> <var_name> <min> <max> [<tol>]
check_plt_range() {
    local desc="$1" plt_dir="$2" var_name="$3" exp_min="$4" exp_max="$5" tol="${6:-1e-10}"
    if [[ ! -d "$plt_dir" ]]; then
        fail "$desc — plotfile dir missing: $plt_dir"; return
    fi
    if ! command -v python3 &>/dev/null; then
        skip "$desc — python3 not available"; return
    fi
    local out
    if out=$(python3 "$SCRIPT_DIR/check_plt_range.py" \
                "$plt_dir" "$var_name" "$exp_min" "$exp_max" "$tol" 2>&1); then
        pass "$desc  ($out)"
    else
        fail "$desc — $out"
    fi
}

# check_plt_coord <desc> <plt_dir> <expected_coord_sys>
check_plt_coord() {
    local desc="$1" plt_dir="$2" expected="$3"
    if [[ ! -d "$plt_dir" ]]; then
        fail "$desc — plotfile dir missing: $plt_dir"; return
    fi
    if ! command -v python3 &>/dev/null; then
        skip "$desc — python3 not available"; return
    fi
    local out
    if out=$(python3 "$SCRIPT_DIR/check_plt_coord.py" \
                "$plt_dir" "$expected" 2>&1); then
        pass "$desc  ($out)"
    else
        fail "$desc — $out"
    fi
}

# check_plt_nlevels <desc> <plt_dir> <expected_nlevels>
check_plt_nlevels() {
    local desc="$1" plt_dir="$2" expected="$3"
    if [[ ! -d "$plt_dir" ]]; then
        fail "$desc — plotfile dir missing: $plt_dir"; return
    fi
    if ! command -v python3 &>/dev/null; then
        skip "$desc — python3 not available"; return
    fi
    local out
    if out=$(python3 "$SCRIPT_DIR/check_plt_nlevels.py" \
                "$plt_dir" "$expected" 2>&1); then
        pass "$desc  ($out)"
    else
        fail "$desc — $out"
    fi
}

# ---------------------------------------------------------------------------
# Executable discovery
# ---------------------------------------------------------------------------
find_exe() {   # find_exe <glob>
    ls "$SRC_DIR"/$1 2>/dev/null | head -1 || true
}

# ---------------------------------------------------------------------------
# Phase 1 — Build
# ---------------------------------------------------------------------------
section "Phase 1: Build"

BUILD_OPTS="DEBUG=FALSE PRECISION=DOUBLE COMP=gnu -j$NPROC"

if $COMPILE; then
    echo "  Building generateTestPlt 3D serial..."
    if make -C "$SRC_DIR" $BUILD_OPTS EBASE=generateTestPlt DIM=3 USE_MPI=FALSE \
            > "$SCRIPT_DIR/build_genPlt3d.log" 2>&1; then
        pass "generateTestPlt 3D serial built"
    else
        fail "generateTestPlt 3D serial build failed — see build_genPlt3d.log"
    fi
    echo "  Building generateTestPlt 3D MPI..."
    if make -C "$SRC_DIR" $BUILD_OPTS EBASE=generateTestPlt DIM=3 USE_MPI=TRUE \
            > "$SCRIPT_DIR/build_genPlt3dmpi.log" 2>&1; then
        pass "generateTestPlt 3D MPI built"
    else
        fail "generateTestPlt 3D MPI build failed — see build_genPlt3dmpi.log"
    fi
    echo "  Building generateTestPlt 2D serial..."
    if make -C "$SRC_DIR" $BUILD_OPTS EBASE=generateTestPlt DIM=2 USE_MPI=FALSE \
            > "$SCRIPT_DIR/build_genPlt2d.log" 2>&1; then
        pass "generateTestPlt 2D serial built"
    else
        fail "generateTestPlt 2D serial build failed — see build_genPlt2d.log"
    fi
else
    skip "Build skipped (--no-compile)"
fi

GEN3D=$(find_exe "generateTestPlt3d.gnu.ex")
GEN3D_MPI=$(find_exe "generateTestPlt3d.gnu.MPI.ex")
GEN2D=$(find_exe "generateTestPlt2d.gnu.ex")
check_file "generateTestPlt 3D serial executable" "$GEN3D"
check_file "generateTestPlt 3D MPI executable"    "$GEN3D_MPI"
check_file "generateTestPlt 2D serial executable" "$GEN2D"

if [[ -z "$GEN3D" ]]; then
    echo -e "${RED}Cannot find serial executable — aborting test run.${NC}"
    exit 1
fi

if command -v mpirun &>/dev/null; then
    MPI_AVAILABLE=true; MPI_CMD=mpirun
elif command -v mpiexec &>/dev/null; then
    MPI_AVAILABLE=true; MPI_CMD=mpiexec
else
    MPI_AVAILABLE=false
    skip "mpirun/mpiexec not found — MPI tests will be skipped"
fi

# ---------------------------------------------------------------------------
# Phase 2 — Structural tests (all 10 field types, single-level)
# ---------------------------------------------------------------------------
section "Phase 2: Structural tests"

mkdir -p "$RUN_DIR"
cd "$RUN_DIR"

echo "  T1 — all 10 field types..."
"$GEN3D" "$SCRIPT_DIR/gen_t1_all_fields.inp" > /dev/null 2>&1
check_dir      "T1 plotfile dir"              "plt_t1_all_fields"
check_file     "T1 Header"                    "plt_t1_all_fields/Header"
header_ncomp_ok "T1 nComp = 10"              "plt_t1_all_fields" 10
header_has_var  "T1 header has const_f"      "plt_t1_all_fields" "const_f"
header_has_var  "T1 header has plane_s"      "plt_t1_all_fields" "plane_s"
header_has_var  "T1 header has plane_sm"     "plt_t1_all_fields" "plane_sm"
header_has_var  "T1 header has dplane_s"     "plt_t1_all_fields" "dplane_s"
header_has_var  "T1 header has dplane_sm"    "plt_t1_all_fields" "dplane_sm"
header_has_var  "T1 header has circle_s"     "plt_t1_all_fields" "circle_s"
header_has_var  "T1 header has circle_sm"    "plt_t1_all_fields" "circle_sm"
header_has_var  "T1 header has shell_s"      "plt_t1_all_fields" "shell_s"
header_has_var  "T1 header has shell_sm"     "plt_t1_all_fields" "shell_sm"
header_has_var  "T1 header has sine_f"       "plt_t1_all_fields" "sine_f"

# ---------------------------------------------------------------------------
# Phase 3 — Value correctness (single-level)
# ---------------------------------------------------------------------------
section "Phase 3: Value correctness"

echo "  T2 — constant field = 7.0..."
"$GEN3D" "$SCRIPT_DIR/gen_t2_const.inp" > /dev/null 2>&1
check_plt_value "T2 all cells = 7.0" "plt_t2_const" "myfield" 7.0

echo "  T3 — degenerate fields all = 5.0..."
"$GEN3D" "$SCRIPT_DIR/gen_t3_degenerate.inp" > /dev/null 2>&1
check_plt_value "T3 c1 constant = 5.0"          "plt_t3_degen" "c1" 5.0
check_plt_value "T3 c2 sphere (outside) = 5.0"  "plt_t3_degen" "c2" 5.0
check_plt_value "T3 c3 ring (outside) = 5.0"    "plt_t3_degen" "c3" 5.0
check_plt_value "T3 c4 dplane (outside) = 5.0"  "plt_t3_degen" "c4" 5.0
check_plt_value "T3 c5 sine (zero amp) = 5.0"   "plt_t3_degen" "c5" 5.0

echo "  T4 — plane_step range [0, 1] on all three axes..."
"$GEN3D" "$SCRIPT_DIR/gen_t4_planes.inp" > /dev/null 2>&1
check_plt_range "T4 px range [0,1]" "plt_t4_planes" "px" 0.0 1.0
check_plt_range "T4 py range [0,1]" "plt_t4_planes" "py" 0.0 1.0
check_plt_range "T4 pz range [0,1]" "plt_t4_planes" "pz" 0.0 1.0

# ---------------------------------------------------------------------------
# Phase 4 — AMR structural tests
# ---------------------------------------------------------------------------
section "Phase 4: AMR structural tests"

echo "  T5 — InBox → 2 levels..."
"$GEN3D" "$SCRIPT_DIR/gen_t5_amr_inbox.inp" > /dev/null 2>&1
check_plt_nlevels "T5 nlevels = 2" "plt_t5_inbox" 2

echo "  T6 — value_greater → 2 levels..."
"$GEN3D" "$SCRIPT_DIR/gen_t6_amr_vgt.inp" > /dev/null 2>&1
check_plt_nlevels "T6 nlevels = 2" "plt_t6_vgt" 2

echo "  T7 — value_less → 2 levels..."
"$GEN3D" "$SCRIPT_DIR/gen_t7_amr_vlt.inp" > /dev/null 2>&1
check_plt_nlevels "T7 nlevels = 2" "plt_t7_vlt" 2

echo "  T8 — adjacent_difference_greater → 2 levels..."
"$GEN3D" "$SCRIPT_DIR/gen_t8_amr_adj.inp" > /dev/null 2>&1
check_plt_nlevels "T8 nlevels = 2" "plt_t8_adj" 2

echo "  T9 — 3-level InBox..."
"$GEN3D" "$SCRIPT_DIR/gen_t9_amr_3level.inp" > /dev/null 2>&1
check_plt_nlevels "T9 nlevels = 3" "plt_t9_3level" 3

echo "  T10 — ref_ratio=4 → 2 levels..."
"$GEN3D" "$SCRIPT_DIR/gen_t10_amr_refrat4.inp" > /dev/null 2>&1
check_plt_nlevels "T10 nlevels = 2" "plt_t10_refrat4" 2

echo "  T11 — no tagged cells → 1 level (truncated)..."
"$GEN3D" "$SCRIPT_DIR/gen_t11_amr_notags.inp" > /dev/null 2>&1
check_plt_nlevels "T11 nlevels = 1" "plt_t11_notags" 1
no_dir            "T11 no Level_1 dir" "plt_t11_notags/Level_1"

# ---------------------------------------------------------------------------
# Phase 5 — AMR value correctness
# ---------------------------------------------------------------------------
section "Phase 5: AMR value correctness"

echo "  T12 — constant field = 3.14 on all AMR levels..."
"$GEN3D" "$SCRIPT_DIR/gen_t12_amr_const.inp" > /dev/null 2>&1
check_plt_value "T12 constant = 3.14 (all levels)" "plt_t12_amr_const" "cfield" 3.14

# ---------------------------------------------------------------------------
# Phase 6 — output_name (issue #57)
# ---------------------------------------------------------------------------
section "Phase 6: output_name (slash-in-variable-name support)"

echo "  T13 — output_name = Y(OH) overrides ParmParse prefix..."
"$GEN3D" "$SCRIPT_DIR/gen_t13_output_name.inp" > /dev/null 2>&1
header_has_var "T13 header contains 'Y(OH)'"  "plt_t13_output_name" "Y(OH)"
header_no_var  "T13 header does not contain 'myvar'" "plt_t13_output_name" "myvar"

# ---------------------------------------------------------------------------
# Phase 7 — Error handling
# ---------------------------------------------------------------------------
section "Phase 7: Error handling"

echo "  err1 — no field.names → abort"
if "$GEN3D" "$SCRIPT_DIR/gen_err1_no_fields.inp" > /dev/null 2>&1; then
    fail "err1 no fields — expected abort, tool exited 0"
else
    pass "err1 no fields — aborted correctly"
fi

echo "  err2 — invalid field type → abort"
if "$GEN3D" "$SCRIPT_DIR/gen_err2_bad_type.inp" > /dev/null 2>&1; then
    fail "err2 bad type — expected abort, tool exited 0"
else
    pass "err2 bad type — aborted correctly"
fi

echo "  err3 — refinement indicator references unknown field → abort"
if "$GEN3D" "$SCRIPT_DIR/gen_err3_bad_refname.inp" > /dev/null 2>&1; then
    fail "err3 bad refname — expected abort, tool exited 0"
else
    pass "err3 bad refname — aborted correctly"
fi

echo "  err4 — refinement indicator with no criterion → abort"
if "$GEN3D" "$SCRIPT_DIR/gen_err4_no_criterion.inp" > /dev/null 2>&1; then
    fail "err4 no criterion — expected abort, tool exited 0"
else
    pass "err4 no criterion — aborted correctly"
fi

# ---------------------------------------------------------------------------
# Phase 8 — Cylindrical coordinate system support
# ---------------------------------------------------------------------------
section "Phase 8: Cylindrical coordinate system"

echo "  T14 — cylinder_step and cylinder_smooth on all three axes..."
"$GEN3D" "$SCRIPT_DIR/gen_t14_cylinder.inp" > /dev/null 2>&1
check_plt_range "T14 cyl_z range [0,1]"   "plt_t14_cylinder" "cyl_z"    0.0 1.0
check_plt_range "T14 cyl_y range [0,1]"   "plt_t14_cylinder" "cyl_y"    0.0 1.0
check_plt_range "T14 cyl_x range [0,1]"   "plt_t14_cylinder" "cyl_x"    0.0 1.0
check_plt_range "T14 cyl_z_sm range (0,1)" "plt_t14_cylinder" "cyl_z_sm" 0.0 1.0

echo "  T15 — coord_sys=1 (cylindrical/RZ), 2D build..."
if [[ -z "$GEN2D" ]]; then
    skip "T15 RZ coord_sys=1 (2D executable not found)"
else
    "$GEN2D" "$SCRIPT_DIR/gen_t15_rz.inp" > /dev/null 2>&1
    check_plt_coord  "T15 Header coord_sys = 1" "plt_t15_rz" 1
    check_plt_value  "T15 RZ constant = 2.71828" "plt_t15_rz" "cfield" 2.71828 1e-5
fi

# ---------------------------------------------------------------------------
# Phase 9 — MPI tests
# ---------------------------------------------------------------------------
section "Phase 9: MPI tests"

if ! $MPI_AVAILABLE; then
    skip "All MPI tests (mpirun/mpiexec not found)"
elif [[ -z "$GEN3D_MPI" ]]; then
    skip "All MPI tests (generateTestPlt MPI executable not found)"
else
    # amr.max_grid_size=4 ensures ≥8 boxes on every grid size used here,
    # so all 4 ranks receive work and no MPI collective deadlocks.
    MPI_GS="amr.max_grid_size=4"

    # constant field: every cell must equal 7.0 regardless of rank count
    for NP in 2 4; do
        echo "  MPI-${NP}a — constant field, ${NP} ranks..."
        $MPI_CMD -np $NP "$GEN3D_MPI" \
            "$SCRIPT_DIR/gen_t2_const.inp" \
            plotfile_name=plt_mpi${NP}_const $MPI_GS \
            > /dev/null 2>&1 || true
        check_plt_value "MPI-${NP}a constant = 7.0" "plt_mpi${NP}_const" "myfield" 7.0
    done

    # InBox AMR: 2 levels must appear on all rank counts
    for NP in 2 4; do
        echo "  MPI-${NP}b — InBox AMR, ${NP} ranks..."
        $MPI_CMD -np $NP "$GEN3D_MPI" \
            "$SCRIPT_DIR/gen_t5_amr_inbox.inp" \
            plotfile_name=plt_mpi${NP}_inbox $MPI_GS \
            > /dev/null 2>&1 || true
        check_plt_nlevels "MPI-${NP}b nlevels = 2" "plt_mpi${NP}_inbox" 2
    done

    # InBox AMR: fewer boxes than ranks (max_grid_size default=16 → 1 box at level 0, 4 ranks)
    echo "  MPI-4b2 — InBox AMR, 4 ranks, fewer boxes than ranks..."
    $MPI_CMD -np 4 "$GEN3D_MPI" \
        "$SCRIPT_DIR/gen_t5_amr_inbox.inp" \
        plotfile_name=plt_mpi4_inbox_fewbox \
        > /dev/null 2>&1 || true
    check_plt_nlevels "MPI-4b2 nlevels = 2 (fewer boxes than ranks)" "plt_mpi4_inbox_fewbox" 2


    # value_greater AMR: 2 levels, 4 ranks
    echo "  MPI-4c — value_greater AMR, 4 ranks..."
    $MPI_CMD -np 4 "$GEN3D_MPI" \
        "$SCRIPT_DIR/gen_t6_amr_vgt.inp" \
        plotfile_name=plt_mpi4_vgt $MPI_GS \
        > /dev/null 2>&1 || true
    check_plt_nlevels "MPI-4c vgt nlevels = 2" "plt_mpi4_vgt" 2

    # value_less AMR: the boundary-cell fix must hold under MPI decomposition
    echo "  MPI-4d — value_less AMR (boundary-cell fix), 4 ranks..."
    $MPI_CMD -np 4 "$GEN3D_MPI" \
        "$SCRIPT_DIR/gen_t7_amr_vlt.inp" \
        plotfile_name=plt_mpi4_vlt $MPI_GS \
        > /dev/null 2>&1 || true
    check_plt_nlevels "MPI-4d vlt nlevels = 2" "plt_mpi4_vlt" 2

    # adjacent_difference_greater AMR: 2 levels, 4 ranks
    echo "  MPI-4e — adjacent_diff AMR, 4 ranks..."
    $MPI_CMD -np 4 "$GEN3D_MPI" \
        "$SCRIPT_DIR/gen_t8_amr_adj.inp" \
        plotfile_name=plt_mpi4_adj $MPI_GS \
        > /dev/null 2>&1 || true
    check_plt_nlevels "MPI-4e adj nlevels = 2" "plt_mpi4_adj" 2

    # 3-level AMR: all 3 levels present with 4 ranks
    echo "  MPI-4f — 3-level InBox AMR, 4 ranks..."
    $MPI_CMD -np 4 "$GEN3D_MPI" \
        "$SCRIPT_DIR/gen_t9_amr_3level.inp" \
        plotfile_name=plt_mpi4_3level $MPI_GS \
        > /dev/null 2>&1 || true
    check_plt_nlevels "MPI-4f 3-level nlevels = 3" "plt_mpi4_3level" 3

    # constant field at all AMR levels, 4 ranks
    echo "  MPI-4g — constant field value across AMR levels, 4 ranks..."
    $MPI_CMD -np 4 "$GEN3D_MPI" \
        "$SCRIPT_DIR/gen_t12_amr_const.inp" \
        plotfile_name=plt_mpi4_amr_const $MPI_GS \
        > /dev/null 2>&1 || true
    check_plt_value "MPI-4g constant = 3.14 (all levels)" "plt_mpi4_amr_const" "cfield" 3.14

    # output_name: header must contain Y(OH) on 2 ranks
    echo "  MPI-2h — output_name Y(OH), 2 ranks..."
    $MPI_CMD -np 2 "$GEN3D_MPI" \
        "$SCRIPT_DIR/gen_t13_output_name.inp" \
        plotfile_name=plt_mpi2_output_name $MPI_GS \
        > /dev/null 2>&1 || true
    header_has_var "MPI-2h header contains 'Y(OH)'" "plt_mpi2_output_name" "Y(OH)"
fi

# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------
TOTAL=$((PASS + FAIL))
echo ""
echo -e "${BOLD}Results: ${GREEN}$PASS${NC}${BOLD}/$TOTAL passed${NC}"
if (( FAIL > 0 )); then
    echo -e "${RED}$FAIL test(s) failed.${NC}"
    exit 1
else
    echo -e "${GREEN}All tests passed.${NC}"
fi
