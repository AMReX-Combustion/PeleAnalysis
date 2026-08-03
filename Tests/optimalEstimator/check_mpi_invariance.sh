#!/usr/bin/env bash
#
# P-invariance check for optimalEstimatorTraining.
#
# With batch_size=-1 every rank's single mini-batch is its whole local training
# set, so the union is the entire training set no matter how many ranks are
# used. The train/validation split is drawn over the global box numbering with a
# fixed seed, so that set is itself rank-count independent. Together those mean
# the whole training trajectory must be identical on 1, 2 and 4 ranks, to
# rounding.
#
# That identity holds only if each rank's loss is normalised by the mini-batch
# weight summed over ALL ranks. Under the previous scheme - average the per-rank
# mean gradients with weight 1/nProcs - it fails on a multi-level plotfile,
# because the ranks do not hold equal shares of the volume weight. Run this
# against the old build too: it should FAIL there and pass here. Note that
# convergence curves alone will not show the difference, because Adam is
# invariant to a constant rescaling of the gradient.
#
# P-invariance is claimed only for the full-batch case. With batch_size smaller
# than a rank's data the ranks partition their samples differently and the
# trajectories legitimately diverge.
#
# Usage:  ./check_mpi_invariance.sh <training-exe> <multi-level-plotfile>
# e.g.    ./check_mpi_invariance.sh ../../Src/optimalEstimatorTraining2d.gnu.MPI.ex plt_test

set -u

EXE=${1:-}
PLT=${2:-}
FEATURES=${FEATURES:-plane_smooth}
TARGETS=${TARGETS:-circle_smooth}
NEURONS=${NEURONS:-"8 8"}
RANKS=${RANKS:-"1 2 4"}
TOL=${TOL:-1e-10}

if [ -z "$EXE" ] || [ -z "$PLT" ]; then
  sed -n '2,26p' "$0"
  exit 1
fi
if [ ! -x "$EXE" ];  then echo "not executable: $EXE"; exit 1; fi
if [ ! -d "$PLT" ];  then echo "no such plotfile: $PLT"; exit 1; fi

MPIRUN=$(command -v mpirun || command -v mpiexec || true)
if [ -z "$MPIRUN" ]; then
  echo "SKIP: no mpirun/mpiexec found."
  exit 0
fi

WORK=$(mktemp -d)
trap 'rm -rf "$WORK"' EXIT

# use_double so the comparison threshold can be tight; nSteps==1 via batch_size=-1
ARGS="infile=$PLT features=$FEATURES targets=$TARGETS neurons=$NEURONS \
      nEpochs=5 batch_size=-1 use_double=1 print_every=1 \
      model_path=$WORK/m minmax_path=$WORK/mm"

for NP in $RANKS; do
  echo "--- $NP rank(s) ---"
  "$MPIRUN" -n "$NP" "$EXE" $ARGS 2>&1 \
    | grep -E "^Epoch" | sed 's/, lr:.*//' > "$WORK/out.$NP"
  if [ ! -s "$WORK/out.$NP" ]; then
    echo "FAIL: no epoch output at $NP ranks (run it by hand to see why)"
    exit 1
  fi
  cat "$WORK/out.$NP"
done

REF=$(echo "$RANKS" | awk '{print $1}')
STATUS=0
for NP in $RANKS; do
  [ "$NP" = "$REF" ] && continue
  # Compare every numeric field of every epoch line against the reference run.
  if ! paste "$WORK/out.$REF" "$WORK/out.$NP" | awk -v tol="$TOL" -v np="$NP" '
      {
        na = split($0, ab, "\t")
        n1 = split(ab[1], a, /[ ,]+/)
        n2 = split(ab[2], b, /[ ,]+/)
        for (i = 1; i <= n1; i++) {
          if (a[i] + 0 == a[i] && a[i] != "") {
            d = a[i] - b[i]; if (d < 0) d = -d
            s = a[i];        if (s < 0) s = -s
            rel = (s > 1 ? d / s : d)
            if (rel > tol) {
              printf "  epoch %d field %d: %s vs %s (rel %g)\n", NR, i, a[i], b[i], rel
              bad = 1
            }
          }
        }
      }
      END { exit bad ? 1 : 0 }'; then
    echo "FAIL: $NP ranks disagrees with $REF rank beyond $TOL"
    STATUS=1
  else
    echo "PASS: $NP ranks matches $REF rank within $TOL"
  fi
done

exit $STATUS
