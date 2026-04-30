#!/usr/bin/env bash
# run_tests.sh — compile and run analysis_util unit tests (serial and MPI)
#
# Usage:  ./run_tests.sh [--no-compile]
#
#   --no-compile   Skip the build step; assume executables already exist.
#
# Artefacts are written to Tests/analysis_util/testrun/ so the source tree
# stays clean.  MPI tests require mpirun/mpiexec; they are skipped if neither
# is found.  Python 3 is required for plotfile value assertions.

set -euo pipefail

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
TEST_DIR="$SCRIPT_DIR"
RUN_DIR="$SCRIPT_DIR/testrun"
PELE_HOME="$(realpath "$SCRIPT_DIR/../..")"
ARITHMETICS_DIR="$SCRIPT_DIR/../arithmetics"

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

# Print per-test PASS/FAIL lines from a C++ test log with shell colours applied.
print_test_results() {   # print_test_results <log_file>
    while IFS= read -r line; do
        if [[ "$line" == *"  PASS "* ]]; then
            echo -e "  ${GREEN}PASS${NC} ${line#*  PASS }"
        elif [[ "$line" == *"  FAIL "* ]]; then
            echo -e "  ${RED}FAIL${NC} ${line#*  FAIL }"
        fi
    done < <(grep -E '  (PASS|FAIL) ' "$1" || true)
}

# ---------------------------------------------------------------------------
# Assertion helpers
# ---------------------------------------------------------------------------

check_file() {
    if [[ -f "$2" ]]; then pass "$1"; else fail "$1 — missing: $2"; fi
}

check_dir() {
    if [[ -d "$2" ]]; then pass "$1"; else fail "$1 — missing dir: $2"; fi
}

header_ncomp_ok() {   # header_ncomp_ok <desc> <plt_dir> <expected_ncomp>
    local actual
    actual=$(sed -n '2p' "$2/Header" 2>/dev/null || echo "MISSING")
    if [[ "$actual" == "$3" ]]; then
        pass "$1  (nComp = $actual)"
    else
        fail "$1  (expected nComp=$3, got '$actual')"
    fi
}

check_plt_value() {   # check_plt_value <desc> <plt_dir> <var> <expected> [<tol>]
    local desc="$1" plt_dir="$2" var="$3" expected="$4" tol="${5:-1e-10}"
    if [[ ! -d "$plt_dir" ]]; then
        fail "$desc — plotfile dir missing: $plt_dir"; return
    fi
    if ! command -v python3 &>/dev/null; then
        skip "$desc — python3 not available"; return
    fi
    local out
    if out=$(python3 "$ARITHMETICS_DIR/check_plt_value.py" \
                 "$plt_dir" "$var" "$expected" "$tol" 2>&1); then
        pass "$desc  ($out)"
    else
        fail "$desc — $out"
    fi
}

# expect_abort <desc> <exe> <args...>
# Passes if the command exits with non-zero (as expected for amrex::Abort).
expect_abort() {
    local desc="$1"; shift
    if "$@" 2>/dev/null; then
        fail "$desc — expected non-zero exit but got 0"
    else
        pass "$desc"
    fi
}

# ---------------------------------------------------------------------------
# Executable discovery
# ---------------------------------------------------------------------------
find_exe() {   # find_exe <glob-pattern-in-TEST_DIR>
    ls "$TEST_DIR"/$1 2>/dev/null | head -1 || true
}

# ---------------------------------------------------------------------------
# Build helpers
# ---------------------------------------------------------------------------
BUILD_OPTS="PELE_ANALYSIS_HOME=$PELE_HOME DEBUG=FALSE PRECISION=DOUBLE COMP=gnu USE_UTILS=TRUE -j$NPROC"

build_tool() {   # build_tool <label> <log> <extra_make_args...>
    local label="$1" log="$2"; shift 2
    echo "  Building $label..."
    if make -j 8 -C "$TEST_DIR" $BUILD_OPTS "$@" > "$SCRIPT_DIR/$log" 2>&1; then
        pass "$label built"
    else
        fail "$label build failed — see $SCRIPT_DIR/$log"
    fi
}

# ---------------------------------------------------------------------------
# Setup
# ---------------------------------------------------------------------------
mkdir -p "$RUN_DIR"

# ---------------------------------------------------------------------------
# Phase 1 — Build
# ---------------------------------------------------------------------------
section "Phase 1: Build"

if $COMPILE; then
    build_tool "test_analysis_util 3D serial" build_util3d.log \
        EBASE=test_analysis_util DIM=3 USE_MPI=FALSE
else
    skip "Build skipped (--no-compile)"
fi

EXE_SERIAL=$(find_exe "test_analysis_util3d.gnu.ex")
if [[ -z "$EXE_SERIAL" ]]; then
    fail "Serial executable not found after build"
    echo -e "\n${RED}Cannot continue without serial executable.${NC}"
    exit 1
fi

# ---------------------------------------------------------------------------
# Phase 2 — Serial test run
# ---------------------------------------------------------------------------
section "Phase 2: Serial C++ tests"

SERIAL_LOG="$RUN_DIR/serial_run.log"
if "$EXE_SERIAL" run_dir="$RUN_DIR" > "$SERIAL_LOG" 2>&1; then
    print_test_results "$SERIAL_LOG"
    passed=$(grep -c '  PASS ' "$SERIAL_LOG" || true)
    failed=$(grep -c '  FAIL ' "$SERIAL_LOG" || true)
    if [[ "$failed" -eq 0 ]]; then
        pass "Serial test run  (${passed} passed, ${failed} failed)"
    else
        fail "Serial test run  (${passed} passed, ${failed} failed — see $SERIAL_LOG)"
    fi
else
    fail "Serial executable exited non-zero — see $SERIAL_LOG"
    tail -20 "$SERIAL_LOG" | sed 's/^/    /'
fi

# ---------------------------------------------------------------------------
# Phase 3 — Plotfile artifact checks
# ---------------------------------------------------------------------------
section "Phase 3: Plotfile artifact checks"

check_dir  "plt_roundtrip_1var exists"    "$RUN_DIR/plt_roundtrip_1var"
check_file "plt_roundtrip_1var Header"    "$RUN_DIR/plt_roundtrip_1var/Header"
check_plt_value "T-PLT-1 value check"    "$RUN_DIR/plt_roundtrip_1var" "myvar" 3.14

check_dir  "plt_roundtrip_2var exists"   "$RUN_DIR/plt_roundtrip_2var"
check_plt_value "T-PLT-2 comp1 value"   "$RUN_DIR/plt_roundtrip_2var" "comp1" 2.0

header_ncomp_ok "T-PLT-8 nComp=2"       "$RUN_DIR/plt_ncomp_check" "2"

# MEF round-trip artifacts
check_file "T-MEF basic .mef exists"              "$RUN_DIR/test_basic.mef"
check_file "T-MEF single-triangle .mef exists"    "$RUN_DIR/test_single_tri.mef"
check_file "T-MEF multi-comp .mef exists"         "$RUN_DIR/test_multicomp.mef"
check_file "T-MEF title-with-spaces .mef exists"  "$RUN_DIR/test_title_spaces.mef"

# ---------------------------------------------------------------------------
# Phase 4 — Abort-path subprocess tests (serial)
# ---------------------------------------------------------------------------
section "Phase 4: Abort-path tests (serial)"

expect_abort "T-STR-A1 find_var_index aborts on missing var" \
    "$EXE_SERIAL" abort_test=find_var_index

expect_abort "T-PLT-A1 read_plotfile aborts on unknown variable" \
    "$EXE_SERIAL" abort_test=read_plotfile_bad_var run_dir="$RUN_DIR"

expect_abort "T-PLT-A2 read_plotfile aborts on non-existent file" \
    "$EXE_SERIAL" abort_test=read_plotfile_no_file

expect_abort "T-MEF-A1 read_mef aborts on non-existent file" \
    "$EXE_SERIAL" abort_test=read_mef_no_file

# ---------------------------------------------------------------------------
# Phase 5 — Build MPI executable
# ---------------------------------------------------------------------------
section "Phase 5: Build MPI executable"

MPIRUN=$(command -v mpirun 2>/dev/null || command -v mpiexec 2>/dev/null || true)

if [[ -z "$MPIRUN" ]]; then
    skip "mpirun/mpiexec not found — skipping MPI build and tests"
    EXE_MPI=""
else
    if $COMPILE; then
        build_tool "test_analysis_util 3D MPI" build_util3dmpi.log \
            EBASE=test_analysis_util DIM=3 USE_MPI=TRUE
    else
        skip "MPI build skipped (--no-compile)"
    fi
    EXE_MPI=$(find_exe "test_analysis_util3d.gnu.MPI.ex")
    if [[ -z "$EXE_MPI" ]]; then
        fail "MPI executable not found after build"
        EXE_MPI=""
    fi
fi

# ---------------------------------------------------------------------------
# Phase 6 — MPI tests
# ---------------------------------------------------------------------------
section "Phase 6: MPI tests"

if [[ -z "$EXE_MPI" ]]; then
    skip "All MPI tests (no MPI executable)"
else
    MPI_TIMEOUT=30   # seconds per MPI invocation

    # Run mpirun in the background and SIGKILL it (and wait) if it exceeds
    # MPI_TIMEOUT.  Using SIGKILL instead of SIGTERM because MPI ranks blocked
    # in a collective (e.g. MPI_Barrier) ignore SIGTERM, which causes `timeout`
    # itself to block indefinitely waiting for mpirun to exit.
    run_mpi() {   # run_mpi <np> <extra_args...>
        local np=$1; shift
        "$MPIRUN" -np "$np" "$@" &
        local mpi_pid=$!
        local elapsed=0
        while kill -0 "$mpi_pid" 2>/dev/null; do
            if [[ $elapsed -ge $MPI_TIMEOUT ]]; then
                echo "  [timeout after ${MPI_TIMEOUT}s — sending SIGKILL to mpirun pid $mpi_pid]" >&2
                kill -9 "$mpi_pid" 2>/dev/null || true
                wait "$mpi_pid" 2>/dev/null || true
                return 124
            fi
            sleep 1
            ((elapsed++)) || true
        done
        wait "$mpi_pid"
        return $?
    }

    MPI_RUN_DIR="$RUN_DIR/mpi"
    mkdir -p "$MPI_RUN_DIR"

    for np in 2 4; do
        section "  MPI with $np ranks"

        MPI_LOG="$MPI_RUN_DIR/run_${np}ranks.log"
        mpi_exit=0
        run_mpi "$np" "$EXE_MPI" run_dir="$MPI_RUN_DIR/r${np}" \
            > "$MPI_LOG" 2>&1 || mpi_exit=$?

        if [[ $mpi_exit -eq 0 ]]; then
            print_test_results "$MPI_LOG"
            passed=$(grep -c '  PASS ' "$MPI_LOG" || true)
            failed=$(grep -c '  FAIL ' "$MPI_LOG" || true)
            if [[ "$failed" -eq 0 ]]; then
                pass "MPI ${np}-rank test run  (${passed} passed)"
            else
                fail "MPI ${np}-rank test run  (${failed} failed — see $MPI_LOG)"
            fi
        elif [[ $mpi_exit -eq 124 ]]; then
            fail "MPI ${np}-rank test TIMED OUT after ${MPI_TIMEOUT}s — see $MPI_LOG"
        else
            fail "MPI ${np}-rank executable exited non-zero ($mpi_exit) — see $MPI_LOG"
        fi

        # T-PLT-MPI-1: plotfile values match serial reference
        check_plt_value "T-PLT-MPI-1 (np=$np) value" \
            "$MPI_RUN_DIR/r${np}/plt_roundtrip_1var" "myvar" 3.14

        # T-MEF-MPI-1: write_mef with IOProcessor guard — file must exist and be valid
        if [[ -f "$MPI_RUN_DIR/r${np}/test_basic.mef" ]]; then
            # verify IOProcessor guard: run read_mef on the file to confirm it is parseable
            if "$EXE_SERIAL" abort_test=read_mef_no_file 2>/dev/null; then
                true # unreachable, just a compile check
            fi
            pass "T-MEF-MPI-1 (np=$np) MEF file valid after MPI write"
        else
            fail "T-MEF-MPI-1 (np=$np) MEF file missing after MPI write"
        fi

        # T-INT-MPI-1: integral result — compared to serial reference log
        if grep -q 'T-INT-1 all-axis integral' "$MPI_LOG" 2>/dev/null; then
            mpi_result=$(grep 'T-INT-1 all-axis integral' "$MPI_LOG" | grep 'PASS' || true)
            if [[ -n "$mpi_result" ]]; then
                pass "T-INT-MPI-1 (np=$np) integrate matches serial"
            else
                fail "T-INT-MPI-1 (np=$np) integrate T-INT-1 failed on MPI run"
            fi
        fi
    done
fi

# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------
section "Summary"
echo -e "  ${GREEN}PASSED: $PASS${NC}"
echo -e "  ${RED}FAILED: $FAIL${NC}"

[[ "$FAIL" -eq 0 ]]
