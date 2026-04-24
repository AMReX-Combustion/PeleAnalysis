#!/usr/bin/env bash
# run_tests.sh — compile and run all jpdf tests (2D, 3D serial and MPI)
#
# Usage:  ./run_tests.sh [--no-compile]
#
#   --no-compile   Skip the build step; assume executables already exist in Src/
#
# All test artefacts are written to Tests/jpdf/testrun/ so the source tree stays
# clean. MPI tests require mpirun/mpiexec; they are skipped if neither is found.

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

section()  { echo -e "\n${BOLD}=== $* ===${NC}"; }
pass()     { echo -e "  ${GREEN}PASS${NC} $*"; ((PASS++)) || true; }
fail()     { echo -e "  ${RED}FAIL${NC} $*"; ((FAIL++)) || true; }
skip()     { echo -e "  ${YELLOW}SKIP${NC} $*"; }

# ---------------------------------------------------------------------------
# Assertion helpers
# ---------------------------------------------------------------------------

check_file() {   # check_file <desc> <path>
    if [[ -f "$2" ]]; then pass "$1"; else fail "$1 — missing: $2"; fi
}

check_dir() {    # check_dir <desc> <path>
    if [[ -d "$2" ]]; then pass "$1"; else fail "$1 — missing dir: $2"; fi
}

# pdf_sum_ok <desc> <matlab_dat_file>
# Sum all whitespace-separated values; expect result in [0.99, 1.01].
pdf_sum_ok() {
    local desc="$1" file="$2"
    if [[ ! -f "$file" ]]; then fail "$desc — missing: $file"; return; fi
    local s
    s=$(awk '{for(i=1;i<=NF;i++) s+=$i} END {printf "%.6f", s}' "$file")
    if awk -v s="$s" 'BEGIN { exit !(s >= 0.99 && s <= 1.01) }'; then
        pass "$desc  (sum = $s)"
    else
        fail "$desc  (sum = $s, expected ≈ 1.0)"
    fi
}

# gpd_sum_ok <desc> <gnuplot_gpd_file>
# Third column sums to ≈ 1.0.
gpd_sum_ok() {
    local desc="$1" file="$2"
    if [[ ! -f "$file" ]]; then fail "$desc — missing: $file"; return; fi
    local s
    s=$(awk '{s+=$3} END {printf "%.6f", s}' "$file")
    if awk -v s="$s" 'BEGIN { exit !(s >= 0.99 && s <= 1.01) }'; then
        pass "$desc  (sum = $s)"
    else
        fail "$desc  (sum = $s, expected ≈ 1.0)"
    fi
}

# gpd_close <desc> <file_serial> <file_mpi> <tol>
# Check max|col3_a − col3_b| < tol (MPI vs serial comparison).
gpd_close() {
    local desc="$1" fa="$2" fb="$3" tol="$4"
    if [[ ! -f "$fa" || ! -f "$fb" ]]; then fail "$desc — output file(s) missing"; return; fi
    local mx
    mx=$(paste "$fa" "$fb" | awk \
        '{d=$3-$6; if(d<0)d=-d; if(d>mx)mx=d} END {printf "%.2e", mx+0}')
    if awk -v mx="$mx" -v t="$tol" 'BEGIN { exit !(mx+0 < t+0) }'; then
        pass "$desc  (max diff = $mx)"
    else
        fail "$desc  (max diff = $mx, threshold = $tol)"
    fi
}

# condmean_range_ok <desc> <dat_file> <lo> <hi>
# All non-zero values in a condMean .dat matrix must lie in [lo, hi].
condmean_range_ok() {
    local desc="$1" file="$2" lo="$3" hi="$4"
    if [[ ! -f "$file" ]]; then fail "$desc — missing: $file"; return; fi
    if awk -v lo="$lo" -v hi="$hi" \
        'NF>0 {for(i=1;i<=NF;i++) if($i+0!=0 && ($i+0<lo || $i+0>hi)) exit 1}' "$file"; then
        pass "$desc"
    else
        fail "$desc — value outside [$lo, $hi]"
    fi
}

# ---------------------------------------------------------------------------
# Executable discovery — return the first match of a glob
# ---------------------------------------------------------------------------
find_exe() {
    local result
    result=$(ls "$SRC_DIR"/$1 2>/dev/null | head -1 || true)
    echo "$result"
}

# ---------------------------------------------------------------------------
# Phase 1 — Build
# ---------------------------------------------------------------------------
section "Phase 1: Build"

BUILD_OPTS="DEBUG=FALSE PRECISION=DOUBLE COMP=gnu -j$NPROC"

build_tool() {   # build_tool <label> <log> <make_args...>
    local label="$1" log="$2"; shift 2
    echo "  Building $label..."
    if make -C "$SRC_DIR" $BUILD_OPTS "$@" > "$SCRIPT_DIR/$log" 2>&1; then
        pass "$label built"
    else
        fail "$label build failed — see $log"
    fi
}

if $COMPILE; then
    build_tool "generateTestPlt 2D" build_genPlt2d.log EBASE=generateTestPlt DIM=2 USE_MPI=FALSE
    build_tool "generateTestPlt 3D" build_genPlt3d.log EBASE=generateTestPlt DIM=3 USE_MPI=FALSE
    build_tool "jpdf 2D serial"     build_jpdf2d.log   EBASE=jpdf DIM=2 USE_MPI=FALSE
    build_tool "jpdf 3D serial"     build_jpdf3d.log   EBASE=jpdf DIM=3 USE_MPI=FALSE
    build_tool "jpdf 3D MPI"        build_jpdf3dmpi.log EBASE=jpdf DIM=3 USE_MPI=TRUE
else
    skip "Build skipped (--no-compile)"
fi

GEN2D=$(find_exe "generateTestPlt2d.gnu.ex")
GEN3D=$(find_exe "generateTestPlt3d.gnu.ex")
JPDF2D=$(find_exe "jpdf2d.gnu.ex")
JPDF3D=$(find_exe "jpdf3d.gnu.ex")
JPDF3D_MPI=$(find_exe "jpdf3d.gnu.MPI.ex")

check_file "generateTestPlt 2D" "$GEN2D"
check_file "generateTestPlt 3D" "$GEN3D"
check_file "jpdf 2D serial"     "$JPDF2D"
check_file "jpdf 3D serial"     "$JPDF3D"
check_file "jpdf 3D MPI"        "$JPDF3D_MPI"

if command -v mpirun &>/dev/null; then
    MPI_AVAILABLE=true; MPI_CMD=mpirun
elif command -v mpiexec &>/dev/null; then
    MPI_AVAILABLE=true; MPI_CMD=mpiexec
else
    MPI_AVAILABLE=false
    skip "mpirun/mpiexec not found — MPI tests will be skipped"
fi

# ---------------------------------------------------------------------------
# Phase 2 — Generate test plotfiles
# ---------------------------------------------------------------------------
section "Phase 2: Generate test plotfiles"

mkdir -p "$RUN_DIR"
cd "$RUN_DIR"

echo "  Generating 3D plotfiles..."
"$GEN3D" "$SCRIPT_DIR/generate_plt_A.inp" > /dev/null 2>&1
"$GEN3D" "$SCRIPT_DIR/generate_plt_B.inp" > /dev/null 2>&1
check_dir  "plt_jpdf_A (3D)"  "plt_jpdf_A"
check_file "plt_jpdf_A/Header" "plt_jpdf_A/Header"
check_dir  "plt_jpdf_B (3D)"  "plt_jpdf_B"

echo "  Generating 2D plotfiles..."
"$GEN2D" "$SCRIPT_DIR/generate_plt_A_2d.inp" > /dev/null 2>&1
"$GEN2D" "$SCRIPT_DIR/generate_plt_B_2d.inp" > /dev/null 2>&1
check_dir  "plt_jpdf_A_2d (2D)"  "plt_jpdf_A_2d"
check_file "plt_jpdf_A_2d/Header" "plt_jpdf_A_2d/Header"
check_dir  "plt_jpdf_B_2d (2D)"  "plt_jpdf_B_2d"

# ---------------------------------------------------------------------------
# Phase 3 — 3D serial tests
# ---------------------------------------------------------------------------
section "Phase 3: 3D serial tests"

T="$SCRIPT_DIR"   # shorthand for input-file directory

echo "  Test 1 — all output formats"
"$JPDF3D" "$T/test1_all_outputs.inp" outSuffix=_3d_t1 > /dev/null 2>&1
D=plt_jpdf_A_3d_t1
check_file "T1-3D gnuplot"   "$D/Pdf_var_plane_var_sine.gpd"
check_file "T1-3D MATLAB"    "$D/Pdf_var_plane_var_sine.dat"
check_file "T1-3D Tecplot"   "$D/Pdf_var_plane_var_sine.tpd"
check_file "T1-3D FAB"       "$D/Pdf_var_plane_var_sine.fab"
check_file "T1-3D scatter"   "$D/Scatter_var_plane_var_sine.dat"
check_dir  "T1-3D plotfile"  "$D/Level_0"
pdf_sum_ok "T1-3D PDF ∑≈1"   "$D/Pdf_var_plane_var_sine.dat"

echo "  Test 2 — 2D conditional mean"
"$JPDF3D" "$T/test2_condmean.inp" outSuffix=_3d_t2 > /dev/null 2>&1
D=plt_jpdf_A_3d_t2
check_file "T2-3D condMean file" "$D/condMean_var_cond_on_var_plane_var_sine.dat"
condmean_range_ok "T2-3D condMean ≈ 0.5" \
    "$D/condMean_var_cond_on_var_plane_var_sine.dat" 0.4 0.6

echo "  Test 3 — condMean duplicate detection"
"$JPDF3D" "$T/test3_condmean_duplicate.inp" outSuffix=_3d_t3 > /dev/null 2>&1
D=plt_jpdf_A_3d_t3
check_file "T3-3D condMean plane-sine"  "$D/condMean_var_cond_on_var_plane_var_sine.dat"
check_file "T3-3D condMean plane-cond"  "$D/condMean_var_cond_on_var_plane_var_cond.dat"
check_file "T3-3D condMean sine-cond"   "$D/condMean_var_cond_on_var_sine_var_cond.dat"

echo "  Test 4 — useminmax + clamping"
"$JPDF3D" "$T/test4_useminmax.inp" outSuffix=_3d_t4 > /dev/null 2>&1
D=plt_jpdf_A_3d_t4
check_file "T4-3D plotfile Header" "$D/Header"
if grep -qE '\-[[:space:]]*5(\.0+)?[eE][-+]?01|\-0\.5' "$D/Header" 2>/dev/null; then
    pass "T4-3D Header encodes useminmax1 lower bound -0.5"
else
    fail "T4-3D Header does not contain -0.5 axis value"
fi

echo "  Test 5 — do_conditioning=1 (normalization fix)"
"$JPDF3D" "$T/test5_conditioning_1.inp" outSuffix=_3d_t5 output_gnuplot=1 > /dev/null 2>&1
D=plt_jpdf_A_3d_t5
check_file "T5-3D plotfile Header" "$D/Header"
gpd_sum_ok "T5-3D PDF ∑≈1 (norm fix)" "$D/Pdf_var_plane_var_sine.gpd"

echo "  Test 6 — do_conditioning=2"
"$JPDF3D" "$T/test6_conditioning_2.inp" outSuffix=_3d_t6 output_gnuplot=1 > /dev/null 2>&1
D=plt_jpdf_A_3d_t6
check_file "T6-3D plotfile Header" "$D/Header"
gpd_sum_ok "T6-3D PDF ∑≈1 (norm fix)" "$D/Pdf_var_plane_var_sine.gpd"

echo "  Test 7 — norm_cVal=1"
"$JPDF3D" "$T/test7_normcval.inp" outSuffix=_3d_t7 output_gnuplot=1 > /dev/null 2>&1
D=plt_jpdf_A_3d_t7
check_file "T7-3D plotfile Header" "$D/Header"
gpd_sum_ok "T7-3D PDF ∑≈1 (norm fix)" "$D/Pdf_var_plane_var_sine.gpd"

echo "  Test 8 — temporal averaging"
"$JPDF3D" "$T/test8_average.inp" outSuffix=_3d_t8 > /dev/null 2>&1
check_dir  "T8-3D per-file A"   "plt_jpdf_A_3d_t8"
check_dir  "T8-3D per-file B"   "plt_jpdf_B_3d_t8"
check_dir  "T8-3D average dir"  "JPDFAverage_3d_t8"
check_file "T8-3D average gpd"  "JPDFAverage_3d_t8/Pdf_var_plane_var_sine.gpd"
gpd_sum_ok "T8-3D average ∑≈1"  "JPDFAverage_3d_t8/Pdf_var_plane_var_sine.gpd"
gpd_close  "T8-3D average ≈ per-file" \
    "plt_jpdf_A_3d_t8/Pdf_var_plane_var_sine.gpd" \
    "JPDFAverage_3d_t8/Pdf_var_plane_var_sine.gpd" 1e-12

# ---------------------------------------------------------------------------
# Phase 4 — 2D serial tests (1, 2, 5, 8)
# ---------------------------------------------------------------------------
section "Phase 4: 2D serial tests"

echo "  Test 1 — all output formats (2D)"
"$JPDF2D" "$T/test1_all_outputs.inp" infile=plt_jpdf_A_2d outSuffix=_2d_t1 > /dev/null 2>&1
D=plt_jpdf_A_2d_2d_t1
check_file "T1-2D gnuplot"   "$D/Pdf_var_plane_var_sine.gpd"
check_file "T1-2D MATLAB"    "$D/Pdf_var_plane_var_sine.dat"
check_file "T1-2D Tecplot"   "$D/Pdf_var_plane_var_sine.tpd"
check_file "T1-2D FAB"       "$D/Pdf_var_plane_var_sine.fab"
check_file "T1-2D scatter"   "$D/Scatter_var_plane_var_sine.dat"
check_dir  "T1-2D plotfile"  "$D/Level_0"
pdf_sum_ok "T1-2D PDF ∑≈1"   "$D/Pdf_var_plane_var_sine.dat"

echo "  Test 2 — 2D conditional mean (2D)"
"$JPDF2D" "$T/test2_condmean.inp" infile=plt_jpdf_A_2d outSuffix=_2d_t2 > /dev/null 2>&1
D=plt_jpdf_A_2d_2d_t2
check_file "T2-2D condMean file" "$D/condMean_var_cond_on_var_plane_var_sine.dat"
condmean_range_ok "T2-2D condMean ≈ 0.5" \
    "$D/condMean_var_cond_on_var_plane_var_sine.dat" 0.4 0.6

echo "  Test 5 — do_conditioning=1 (2D)"
"$JPDF2D" "$T/test5_conditioning_1.inp" infile=plt_jpdf_A_2d outSuffix=_2d_t5 \
    output_gnuplot=1 > /dev/null 2>&1
D=plt_jpdf_A_2d_2d_t5
check_file "T5-2D plotfile Header" "$D/Header"
gpd_sum_ok "T5-2D PDF ∑≈1 (norm fix)" "$D/Pdf_var_plane_var_sine.gpd"

echo "  Test 8 — temporal averaging (2D)"
"$JPDF2D" "$T/test8_average_2d.inp" > /dev/null 2>&1
check_dir  "T8-2D per-file A"   "plt_jpdf_A_2d_2d_t8"
check_dir  "T8-2D per-file B"   "plt_jpdf_B_2d_2d_t8"
check_dir  "T8-2D average dir"  "JPDFAverage_2d_t8"
check_file "T8-2D average gpd"  "JPDFAverage_2d_t8/Pdf_var_plane_var_sine.gpd"
gpd_sum_ok "T8-2D average ∑≈1"  "JPDFAverage_2d_t8/Pdf_var_plane_var_sine.gpd"
gpd_close  "T8-2D average ≈ per-file" \
    "plt_jpdf_A_2d_2d_t8/Pdf_var_plane_var_sine.gpd" \
    "JPDFAverage_2d_t8/Pdf_var_plane_var_sine.gpd" 1e-12

# ---------------------------------------------------------------------------
# Phase 5 — MPI tests (3D only)
# ---------------------------------------------------------------------------
section "Phase 5: MPI tests"

if ! $MPI_AVAILABLE; then
    skip "All MPI tests (mpirun not found)"
elif [[ -z "$JPDF3D_MPI" ]]; then
    skip "All MPI tests (jpdf 3D MPI executable not found)"
else
    # Test 1: basic JPDF with 2 and 4 ranks; compare PDF to serial result
    for NP in 2 4; do
        echo "  Test 1 — basic JPDF, $NP MPI ranks"
        $MPI_CMD -np $NP "$JPDF3D_MPI" "$T/test1_all_outputs.inp" \
            outSuffix=_3d_mpi${NP}_t1 > /dev/null 2>&1
        D=plt_jpdf_A_3d_mpi${NP}_t1
        check_file "T1-MPI${NP} gnuplot"   "$D/Pdf_var_plane_var_sine.gpd"
        check_dir  "T1-MPI${NP} plotfile"  "$D/Level_0"
        gpd_sum_ok "T1-MPI${NP} PDF ∑≈1"  "$D/Pdf_var_plane_var_sine.gpd"
        gpd_close  "T1-MPI${NP} matches serial" \
            "plt_jpdf_A_3d_t1/Pdf_var_plane_var_sine.gpd" \
            "$D/Pdf_var_plane_var_sine.gpd" 1e-10
    done

    # Test 5: conditioning with 4 ranks (validates normalization fix under MPI)
    echo "  Test 5 — conditioning, 4 MPI ranks"
    $MPI_CMD -np 4 "$JPDF3D_MPI" "$T/test5_conditioning_1.inp" \
        outSuffix=_3d_mpi4_t5 output_gnuplot=1 > /dev/null 2>&1
    D=plt_jpdf_A_3d_mpi4_t5
    check_file "T5-MPI4 plotfile Header"  "$D/Header"
    gpd_sum_ok "T5-MPI4 PDF ∑≈1"         "$D/Pdf_var_plane_var_sine.gpd"
    gpd_close  "T5-MPI4 matches serial"  \
        "plt_jpdf_A_3d_t5/Pdf_var_plane_var_sine.gpd" \
        "$D/Pdf_var_plane_var_sine.gpd" 1e-10

    # Test 8: averaging with 4 ranks
    echo "  Test 8 — averaging, 4 MPI ranks"
    $MPI_CMD -np 4 "$JPDF3D_MPI" "$T/test8_average.inp" \
        outSuffix=_3d_mpi4_t8 > /dev/null 2>&1
    check_dir  "T8-MPI4 average dir" "JPDFAverage_3d_mpi4_t8"
    check_file "T8-MPI4 average gpd" "JPDFAverage_3d_mpi4_t8/Pdf_var_plane_var_sine.gpd"
    gpd_sum_ok "T8-MPI4 average ∑≈1" "JPDFAverage_3d_mpi4_t8/Pdf_var_plane_var_sine.gpd"
    gpd_close  "T8-MPI4 average ≈ serial average" \
        "JPDFAverage_3d_t8/Pdf_var_plane_var_sine.gpd" \
        "JPDFAverage_3d_mpi4_t8/Pdf_var_plane_var_sine.gpd" 1e-10
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
