#!/bin/zsh
# Simulation study: generate networks, then fit all estimators, one parallel lane per scenario.
# Usage (from the project root): enrichment/analysis/sim/run_sim.sh [NSIM]   (default 100)
# Logs: enrichment/analysis/sim/logs/
set -u
NSIM=${1:-100}
STATA=/Applications/StataNow/StataMP.app/Contents/MacOS/StataMP
SIM=enrichment/analysis/sim
mkdir -p $SIM/logs $SIM/results

# 1. networks (interval 1e6, burn-in 5e6: lag-1 autocorrelation of the statistics about 0 at this interval)
NSIM=$NSIM SCEN="1 2 3 4" $STATA -e do $SIM/sim_generate.do
mv sim_generate.log $SIM/logs/

# 2. per scenario: Stata joint fits and R fits, scenarios in parallel (each Stata lane in its own folder,
#    because -e writes <dofile>.log to the working directory)
for s in 1 2 3 4; do
  (
    mkdir -p $SIM/logs/lane$s
    cd $SIM/logs/lane$s
    SCEN=$s REPS="$(seq -s ' ' 1 $NSIM)" $STATA -e do "$OLDPWD/$SIM/sim_fit_joint.do"
  ) &
  ( Rscript $SIM/sim_fit_r.R $s 1 $NSIM > $SIM/logs/r_s$s.log 2>&1 ) &
done
wait
echo "simulation done"
