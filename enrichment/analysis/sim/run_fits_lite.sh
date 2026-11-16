#!/bin/zsh
# Low-footprint, resumable version of run_fits.sh.
# Usage (from the project root): zsh enrichment/analysis/sim/run_fits_lite.sh [NSIM] [JOBS]
#   JOBS = concurrent Stata processes (default 3), each capped at 1 core and niced;
#   plus one niced R process that works through scenarios 1-4 in turn.
# Skips replications already present in results/ (joint: any joint_s<k>_*.csv; R: r_s<k>.csv),
# so it can be killed and restarted at any time. Set SKIP_R=1 or SKIP_STATA=1 to run only one side.
set -u
NSIM=${1:-100}
JOBS=${2:-3}
STATA=/Applications/StataNow/StataMP.app/Contents/MacOS/StataMP
SIM=enrichment/analysis/sim
RES=$SIM/results
mkdir -p $SIM/logs $RES

done_reps() {  # $1 = file glob; prints replication numbers (column 2) already written
  local f
  for f in ${~1}(N); do tail -n +2 $f | cut -d, -f2; done | sort -n | uniq
}

# ---- Stata: missing (scenario, rep) pairs, dealt round-robin into JOBS queues ----
if [[ -z "${SKIP_STATA:-}" ]]; then
  typeset -A queue
  i=0
  for r in $(seq 1 $NSIM); do
    for s in 1 2 3 4; do
      done_reps "$RES/joint_s${s}_*.csv" | grep -qx $r && continue
      q=$(( i % JOBS + 1 )); queue[$q]+="$s:$r "; (( i++ ))
    done
  done
  echo "Stata: $i joint fits to do in $JOBS queues"
  for q in $(seq 1 $JOBS); do
    (
      mkdir -p $SIM/logs/lite_q$q && cd $SIM/logs/lite_q$q
      for sr in ${=queue[$q]:-}; do
        SCEN=${sr%%:*} REPS=${sr##*:} OUTTAG=_q$q NPROC=1 nice -n 10 $STATA -e do ../../sim_fit_joint.do
        mv -f sim_fit_joint.log sim_fit_joint_s${sr%%:*}_r${sr##*:}.log 2>/dev/null
      done
    ) &
  done
fi

# ---- R: one process, scenarios in turn, resuming after the last written rep ----
if [[ -z "${SKIP_R:-}" ]]; then
  (
    for s in 1 2 3 4; do
      last=$(done_reps "$RES/r_s$s.csv" | tail -1); last=${last:-0}
      (( last >= NSIM )) && continue
      nice -n 10 Rscript $SIM/sim_fit_r.R $s $(( last + 1 )) $NSIM >> $SIM/logs/r_s$s.log 2>&1
    done
  ) &
fi

wait
echo "fits done"
