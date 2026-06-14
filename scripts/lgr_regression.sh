#!/usr/bin/env bash
#
# LGR regression harness: fork (this build) vs upstream master, over the
# opm-tests/lgr deck suite.
#
# The fork rebuilds CpGrid's LGR layer and aims to stay bit-compatible with
# upstream master where the two overlap (see docs/PLAN.md). This script makes
# that measurable and regression-guards it: for every opm-tests/lgr/*.DATA it
#   1. runs the fork serial,
#   2. runs the fork in parallel (np=2),
#   3. runs the *master* binary serial as the reference,
# then reports, per deck, the tightest relative tolerance at which
#   * fork-serial  == master-serial   ("bit-compat" with upstream)
#   * fork-serial  == fork-parallel    ("ser==par" consistency)
# A baseline file records the expected status; the script exits non-zero if a
# case that used to pass regresses (or, with -u, if an xfail unexpectedly
# passes).
#
# Usage:
#   lgr_regression.sh -f FORK_FLOW -m MASTER_FLOW -c compareECL \
#                     -d opm-tests/lgr [-w WORKDIR] [-b BASELINE] [-u]
#
# Defaults assume the workspace layout in memory (builds/refined, builds/release).
set -u

WS="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
FORK="${WS}/builds/refined/opm-simulators/bin/flow"
MASTER="${WS}/builds/release/opm-simulators/bin/flow"
COMPARE="${WS}/builds/refined/opm-common/bin/compareECL"
DECKS="${WS}/opm-tests/lgr"
WORK="$(mktemp -d "${TMPDIR:-/tmp}/lgrreg.XXXXXX")"
BASELINE="$(dirname "${BASH_SOURCE[0]}")/lgr_regression.baseline"
UPDATE_XFAIL=0
NP=2
TOL_LADDER="1e-6 1e-5 1e-4 1e-3 1e-2 1e-1"

while getopts "f:m:c:d:w:b:n:u" opt; do
  case "$opt" in
    f) FORK=$OPTARG ;;
    m) MASTER=$OPTARG ;;
    c) COMPARE=$OPTARG ;;
    d) DECKS=$OPTARG ;;
    w) WORK=$OPTARG ;;
    b) BASELINE=$OPTARG ;;
    n) NP=$OPTARG ;;
    u) UPDATE_XFAIL=1 ;;
    *) echo "see header for usage"; exit 2 ;;
  esac
done

for x in "$FORK" "$MASTER" "$COMPARE"; do
  [ -x "$x" ] || { echo "ERROR: not executable: $x"; exit 2; }
done
[ -d "$DECKS" ] || { echo "ERROR: deck dir not found: $DECKS"; exit 2; }
mkdir -p "$WORK"

# tightest <ref-prefix> <new-prefix> : echo smallest matching tol, or FAIL
tightest() {
  local t
  for t in $TOL_LADDER; do
    if "$COMPARE" -t SMRY -a "$1" "$2" "$t" "$t" >/dev/null 2>&1; then
      echo "$t"; return
    fi
  done
  echo FAIL
}

declare -A want
if [ -f "$BASELINE" ]; then
  while read -r deck fmexp spexp; do
    [ -z "${deck:-}" ] && continue
    case "$deck" in \#*) continue ;; esac
    want["$deck"]="$fmexp $spexp"
  done < "$BASELINE"
fi

printf "%-42s %-7s %-7s %-13s %-11s %s\n" DECK forkS forkP "fork==master" "ser==par" status
rc=0
results=""
for f in "$DECKS"/*.DATA; do
  d=$(basename "$f" .DATA)
  fsd="$WORK/${d}_fs"; fpd="$WORK/${d}_fp"; msd="$WORK/${d}_ms"
  rm -rf "$fsd" "$fpd" "$msd"; mkdir -p "$fsd" "$fpd" "$msd"

  "$FORK"   "$f" --parsing-strictness=low --output-dir="$fsd" >"$fsd/log" 2>&1; fs=$?
  mpirun -np "$NP" "$FORK" "$f" --parsing-strictness=low --output-dir="$fpd" >"$fpd/log" 2>&1; fp=$?
  "$MASTER" "$f" --parsing-strictness=low --output-dir="$msd" >"$msd/log" 2>&1; ms=$?

  fm="-"; sp="-"
  [ $fs -eq 0 ] && [ $ms -eq 0 ] && fm=$(tightest "$msd/$d" "$fsd/$d")
  [ $fs -eq 0 ] && [ $fp -eq 0 ] && sp=$(tightest "$fsd/$d" "$fpd/$d")

  # status vs baseline: regression if a previously-passing comparison now fails
  status=ok
  exp="${want[$d]:-}"
  if [ -n "$exp" ]; then
    set -- $exp; expfm=$1; expsp=$2
    rank() { case "$1" in FAIL) echo 99;; -) echo 50;; 1e-6) echo 0;; 1e-5) echo 1;;
                          1e-4) echo 2;; 1e-3) echo 3;; 1e-2) echo 4;; 1e-1) echo 5;; *) echo 98;; esac; }
    if [ "$(rank "$fm")" -gt "$(rank "$expfm")" ] || [ "$(rank "$sp")" -gt "$(rank "$expsp")" ]; then
      status=REGRESSED; rc=1
    elif { [ "$(rank "$fm")" -lt "$(rank "$expfm")" ] || [ "$(rank "$sp")" -lt "$(rank "$expsp")" ]; }; then
      status=improved
      [ "$UPDATE_XFAIL" -eq 1 ] && rc=1
    fi
  fi

  printf "%-42s %-7s %-7s %-13s %-11s %s\n" "$d" "exit$fs" "exit$fp" "$fm" "$sp" "$status"
  results="${results}${d} ${fm} ${sp}"$'\n'
done

echo
echo "# baseline lines (deck  fork==master  ser==par):"
printf "%s" "$results" | sed 's/^/#   /'
echo
[ $rc -eq 0 ] && echo "RESULT: no regressions" || echo "RESULT: regressions detected (or improvements with -u)"
exit $rc
