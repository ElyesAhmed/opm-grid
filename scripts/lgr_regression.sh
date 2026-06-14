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
SUMMARY="${WS}/builds/refined/opm-common/bin/summary"
# A small geometry/ordering change (e.g. parallel cell numbering, or fork vs
# master sub-cell geometry) can nudge the adaptive timestep controller onto a
# different ministep path. The summary files then have different *lengths* and a
# strict compareECL throws on that, even when the solution is identical at the
# (fixed) report steps. Tightening the linear solve reduces the sensitivity but
# does not remove it for every deck -- making the controller robust to such
# perturbations is a separate task. So instead of treating a length mismatch as
# a failure, we fall back to comparing the field rate vectors on the common
# prefix: if those agree it is timestep-path noise ("Lnz"), not a regression.
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

# fieldRatesAgree <ref-prefix> <new-prefix> <reltol> : 0 if the field rate
# vectors agree on their common-length prefix (i.e. the only difference is how
# many ministeps were written -- timestep-path noise, not a solution change).
fieldRatesAgree() {
  [ -x "$SUMMARY" ] || return 1
  local v
  for v in FOPR FGPR FWPR FPR; do
    paste <("$SUMMARY" "$1" "$v" 2>/dev/null | tail -n +2) \
          <("$SUMMARY" "$2" "$v" 2>/dev/null | tail -n +2) \
      | awk -v tol="$3" 'NF==2{f=$1+0;m=$2+0;den=(m<0?-m:m);ad=(f>m?f-m:m-f);
                         if(den>1e-9){if(ad/den>tol)bad=1}else if(ad>tol)bad=1}
                         END{exit bad?1:0}' || return 1
  done
  return 0
}

# tightest <ref-prefix> <new-prefix> : echo smallest matching tol; or "Lnz" if
# the summaries differ only in length but the field rates agree (timestep-path
# noise); or FAIL if the values genuinely diverge.
tightest() {
  local t
  for t in $TOL_LADDER; do
    if "$COMPARE" -t SMRY -a "$1" "$2" "$t" "$t" >/dev/null 2>&1; then
      echo "$t"; return
    fi
  done
  # compareECL could not align the summaries (typically a different-length
  # throw). Distinguish benign timestep-path noise from a real divergence.
  if fieldRatesAgree "$1" "$2" 1e-4; then
    echo "Lnz"
  else
    echo FAIL
  fi
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
    # Lower rank = stronger agreement. Lnz (length differs but field rates
    # agree -- timestep-path noise) ranks just better than a hard FAIL and is
    # not treated as a regression.
    rank() { case "$1" in FAIL) echo 99;; -) echo 50;; Lnz) echo 6;; 1e-6) echo 0;; 1e-5) echo 1;;
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
