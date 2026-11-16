#!/bin/zsh
# Scenario 5 (agency heterogeneity): generate 100 networks in 5 chunks, then fit the simple and the full model
# to each in JOBS one-core queues. Run from the project root: zsh enrichment/analysis/sim/run_s6.sh [JOBS]
set -u
JOBS=${1:-8}
NREP=${2:-60}   # fit only networks 1..NREP (resumable: skips fits already in results/)
STATA=/Applications/StataNow/StataMP.app/Contents/MacOS/StataMP
SIM=enrichment/analysis/sim
mkdir -p $SIM/logs/s6
if [[ ! -f $SIM/nets/s6_edges.csv ]]; then
  for c in 1 2 3 4 5; do
    ( mkdir -p $SIM/logs/s6/gen$c && cd $SIM/logs/s6/gen$c && cp ../../../sim_generate6.do gen.do &&
      CHUNK=$c NSIM=20 NPROC=1 nice -n 10 $STATA -e do gen.do ) &
  done
  wait
  Rscript -e 'suppressMessages(library(tidyverse)); d <- map_dfr(1:5, \(c) read_csv(sprintf("enrichment/analysis/sim/nets/s6_c%d.csv", c), show_col_types = FALSE) |> mutate(`_sim` = (c - 1) * 20 + `_sim`)); stopifnot(n_distinct(d$`_sim`) == 100); write_csv(d, "enrichment/analysis/sim/nets/s6_edges.csv")'
  echo "generated"
fi
cp $SIM/sim_fit_joint6.do $SIM/logs/s6/fit.do
done_reps() { local f; for f in $SIM/results/joint6_$1_*.csv(N); do tail -n +2 $f | cut -d, -f2; done | sort -n | uniq; }
typeset -A queue
i=0
for r in $(seq 1 $NREP); do
  for m in simple full size; do
    done_reps $m | grep -qx $r && continue
    q=$(( i % JOBS + 1 )); queue[$q]+="$m:$r "; (( i++ ))
  done
done
echo "$i fits to do in $JOBS queues"
for q in $(seq 1 $JOBS); do
  (
    mkdir -p $SIM/logs/s6/q$q && cd $SIM/logs/s6/q$q
    for mr in ${=queue[$q]}; do
      MODEL=${mr%%:*} SCEN=5 REPS=${mr##*:} OUTTAG=_q$q NPROC=1 nice -n 10 $STATA -e do ../fit.do
      mv -f fit.log fit_${mr%%:*}_r${mr##*:}.log 2>/dev/null
    done
  ) &
done
wait
echo "fits done"
