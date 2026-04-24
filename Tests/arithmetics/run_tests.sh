#!/usr/bin/env bash
# run_tests.sh — compile and run all arithmetics tests (3D serial and MPI)
#
# Usage:  ./run_tests.sh [--no-compile]
#
#   --no-compile   Skip the build step; assume executables already exist in Src/
#
# All test artefacts are written to Tests/arithmetics/testrun/ so the source
# tree stays clean.  MPI tests require mpirun/mpiexec; they are skipped if
# neither is found.  Python 3 is required for value-correctness assertions.

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
    build_tool "generateTestPlt 3D" build_genPlt3d.log   EBASE=generateTestPlt DIM=3 USE_MPI=FALSE
    build_tool "arithmetics 3D serial" build_arith3d.log EBASE=arithmetics     DIM=3 USE_MPI=FALSE
    build_tool "arithmetics 3D MPI"  build_arith3dmpi.log EBASE=arithmetics    DIM=3 USE_MPI=TRUE
else
    skip "Build skipped (--no-compile)"
fi

GEN3D=$(find_exe "generateTestPlt3d.gnu.ex")
ARITH3D=$(find_exe "arithmetics3d.gnu.ex")
ARITH3D_MPI=$(find_exe "arithmetics3d.gnu.MPI.ex")

check_file "generateTestPlt 3D" "$GEN3D"
check_file "arithmetics 3D serial" "$ARITH3D"
check_file "arithmetics 3D MPI"  "$ARITH3D_MPI"

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

echo "  Generating plt_arith (constant varA=4, varB=2)..."
"$GEN3D" "$SCRIPT_DIR/generate_plt_arith.inp" > /dev/null 2>&1
check_dir  "plt_arith"        "plt_arith"
check_file "plt_arith/Header" "plt_arith/Header"

echo "  Generating plt_arith_divzero (varB = 0 on left half)..."
"$GEN3D" "$SCRIPT_DIR/generate_plt_arith_divzero.inp" > /dev/null 2>&1
check_dir  "plt_arith_divzero"        "plt_arith_divzero"
check_file "plt_arith_divzero/Header" "plt_arith_divzero/Header"

# ---------------------------------------------------------------------------
# Phase 3 — Operator structural tests
# ---------------------------------------------------------------------------
section "Phase 3: Operator structural tests"

# Each test checks: exit 0, output dir exists, nComp=3, all three variable
# names present in Header (two original + one derived).

echo "  Test 1 — add"
"$ARITH3D" "$SCRIPT_DIR/test1_add.inp" > /dev/null 2>&1
check_dir      "T1 output dir"        "plt_t1_add"
header_ncomp_ok "T1 nComp = 3"        "plt_t1_add" 3
header_has_var  "T1 Header has varA"  "plt_t1_add" "varA"
header_has_var  "T1 Header has varB"  "plt_t1_add" "varB"
header_has_var  "T1 Header has sum"   "plt_t1_add" "sum"

echo "  Test 2 — subtract"
"$ARITH3D" "$SCRIPT_DIR/test2_subtract.inp" > /dev/null 2>&1
check_dir      "T2 output dir"        "plt_t2_subtract"
header_ncomp_ok "T2 nComp = 3"        "plt_t2_subtract" 3
header_has_var  "T2 Header has diff"  "plt_t2_subtract" "diff"

echo "  Test 3 — multiply"
"$ARITH3D" "$SCRIPT_DIR/test3_multiply.inp" > /dev/null 2>&1
check_dir      "T3 output dir"        "plt_t3_multiply"
header_ncomp_ok "T3 nComp = 3"        "plt_t3_multiply" 3
header_has_var  "T3 Header has prod"  "plt_t3_multiply" "prod"

echo "  Test 4 — divide (no zeros in denominator)"
"$ARITH3D" "$SCRIPT_DIR/test4_divide.inp" > /dev/null 2>&1
check_dir      "T4 output dir"        "plt_t4_divide"
header_ncomp_ok "T4 nComp = 3"        "plt_t4_divide" 3
header_has_var  "T4 Header has quot"  "plt_t4_divide" "quot"

# ---------------------------------------------------------------------------
# Phase 4 — Divide-by-zero handling
# ---------------------------------------------------------------------------
section "Phase 4: Divide-by-zero handling"

echo "  Test 5 — divide + zeros detected (checkDivByZero=1)"
if "$ARITH3D" "$SCRIPT_DIR/test5_divzero_check.inp" > /dev/null 2>&1; then
    fail "T5 divide-by-zero check — expected abort, tool exited 0"
else
    pass "T5 divide-by-zero check — tool aborted as expected"
fi
no_dir "T5 no output written on abort" "plt_t5_divzero_check"

echo "  Test 6 — divide + zeros skipped (checkDivByZero=0)"
"$ARITH3D" "$SCRIPT_DIR/test6_divzero_skip.inp" > /dev/null 2>&1
check_dir "T6 output written despite zeros" "plt_t6_divzero_skip"

# ---------------------------------------------------------------------------
# Phase 5 — Output naming
# ---------------------------------------------------------------------------
section "Phase 5: Output naming"

echo "  Test 7 — default naming (infile_operator)"
# plt_arith + add without outfile override → plt_arith_add
"$ARITH3D" infile=plt_arith inVarAName=varA inVarBName=varB \
    outVarName=sum operator=add > /dev/null 2>&1
check_dir "T7 default output dir plt_arith_add" "plt_arith_add"

echo "  Test 7b — outfile override"
"$ARITH3D" "$SCRIPT_DIR/test7_outfile_override.inp" > /dev/null 2>&1
check_dir "T7b custom output dir exists"  "plt_custom_output"

# ---------------------------------------------------------------------------
# Phase 6 — Value correctness (round-trip)
# ---------------------------------------------------------------------------
section "Phase 6: Value correctness (round-trip)"

echo "  Test 8 — add/subtract inverse: (varA + varB) - varB = varA"
"$ARITH3D" "$SCRIPT_DIR/test8a_roundtrip_add.inp" > /dev/null 2>&1
"$ARITH3D" "$SCRIPT_DIR/test8b_roundtrip_sub.inp" > /dev/null 2>&1
check_plt_value "T8 add+subtract recovered = 4.0" "plt_t8b_sub" "recovered" 4.0

echo "  Test 9 — multiply/divide inverse: (varA * varB) / varB = varA"
"$ARITH3D" "$SCRIPT_DIR/test9a_roundtrip_mul.inp" > /dev/null 2>&1
"$ARITH3D" "$SCRIPT_DIR/test9b_roundtrip_div.inp" > /dev/null 2>&1
check_plt_value "T9 multiply+divide recovered = 4.0" "plt_t9b_div" "recovered" 4.0

# Also verify the intermediate results to catch sign / operand-order bugs
check_plt_value "T8 intermediate sum = 6.0"  "plt_t8a_add" "sum"  6.0
check_plt_value "T9 intermediate prod = 8.0" "plt_t9a_mul" "prod" 8.0

# ---------------------------------------------------------------------------
# Phase 7 — Error handling
# ---------------------------------------------------------------------------
section "Phase 7: Error handling"

echo "  Test 10 — inVarAName not in plotfile"
if "$ARITH3D" infile=plt_arith inVarAName=no_such_var inVarBName=varB \
    outVarName=out operator=add > /dev/null 2>&1; then
    fail "T10 missing inVarAName — expected abort, got exit 0"
else
    pass "T10 missing inVarAName — aborted correctly"
fi

echo "  Test 11 — inVarBName not in plotfile"
if "$ARITH3D" infile=plt_arith inVarAName=varA inVarBName=no_such_var \
    outVarName=out operator=add > /dev/null 2>&1; then
    fail "T11 missing inVarBName — expected abort, got exit 0"
else
    pass "T11 missing inVarBName — aborted correctly"
fi

echo "  Test 12 — invalid operator"
if "$ARITH3D" infile=plt_arith inVarAName=varA inVarBName=varB \
    outVarName=out operator=modulo > /dev/null 2>&1; then
    fail "T12 invalid operator — expected abort, got exit 0"
else
    pass "T12 invalid operator — aborted correctly"
fi

# ---------------------------------------------------------------------------
# Phase 8 — MPI tests
# ---------------------------------------------------------------------------
section "Phase 8: MPI tests"

if ! $MPI_AVAILABLE; then
    skip "All MPI tests (mpirun not found)"
elif [[ -z "$ARITH3D_MPI" ]]; then
    skip "All MPI tests (arithmetics 3D MPI executable not found)"
else
    # Tests 13–14: add with 2 and 4 ranks; verify derived value = 6.0
    for NP in 2 4; do
        echo "  Test — add, $NP MPI ranks"
        $MPI_CMD -np $NP "$ARITH3D_MPI" \
            infile=plt_arith inVarAName=varA inVarBName=varB \
            outVarName=sum operator=add outfile=plt_mpi${NP}_add \
            > /dev/null 2>&1
        check_dir       "T-MPI${NP} output dir"    "plt_mpi${NP}_add"
        check_plt_value "T-MPI${NP} sum = 6.0"     "plt_mpi${NP}_add" "sum" 6.0
        check_plt_value "T-MPI${NP} varA preserved" "plt_mpi${NP}_add" "varA" 4.0
    done

    # Test 15: divide + zeros with 4 ranks — all ranks must abort cleanly
    echo "  Test — divide + zero check, 4 MPI ranks"
    if $MPI_CMD -np 4 "$ARITH3D_MPI" \
            "$SCRIPT_DIR/test5_divzero_check.inp" \
            outfile=plt_mpi4_divzero_check \
            > /dev/null 2>&1; then
        fail "T-MPI4 divide-by-zero — expected abort on all ranks, got exit 0"
    else
        pass "T-MPI4 divide-by-zero — all ranks aborted"
    fi
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
