* ============================================================
* Simulation study, step 3: fit the joint (three-layer) model to simulated networks.
* Env: SCEN = scenario (one number); REPS = replication numbers, e.g. "1 2 3" (default 1);
*      START = "cd" (contrastive divergence, default) or "true" (start at the true values; for timing checks).
* Appends one row per replication to sim/results/joint_s<k><OUTTAG>.csv (env OUTTAG, e.g. _a, for parallel lanes):
*   scen rep converged iterations interval seconds, then b_<term> and se_<term> for every coefficient.
* Run from the project root: /Applications/StataNow/StataMP.app/Contents/MacOS/StataMP -e do enrichment/analysis/sim/sim_fit_joint.do
* ============================================================
clear all
set more off
set maxvar 40000
* Env NPROC = cap on Stata/MP cores (used by run_fits_lite.sh to leave the machine usable)
local nproc : env NPROC
if "`nproc'" != "" set processors `nproc'
local root "`c(pwd)'"
local sim "`root'/enrichment/analysis/sim"
local nwergm "`root'/software/nwergm-tripartite"
adopath ++ "`nwergm'"
adopath ++ "`nwergm'/lib/plugins/macos"
quietly do "`nwergm'/unw_ergm.do"
capture mkdir "`sim'/results"

local s : env SCEN
local reps : env REPS
if "`reps'" == "" local reps 1
* Env MCTERM = MCMLE termination rule (default hummel: stop at full Hummel steps; confidence = R's multivariate test)
local term : env MCTERM
if "`term'" == "" local term hummel
local start : env START
if "`start'" == "" local start cd

global TPDIR "`sim'/world"
global TPEDGES "`sim'/world/edges.csv"
global TPCOVS pc_home ac_home
do "`sim'/tp_load.do"
tempfile attrs
save `attrs'

local terms edges bdsp1(2) bnodematch1(nat, byattr(type)) edgecov(pc_home) ///
	structural(tpcoloc(attr=key_cty) tpstarcov(x=ac_home attr=cls))
local con constraints(degrees + bd(attribs = isag iscl, maxout = bmax bmax, minout = bmin bmin))

import delimited using "`sim'/nets/s`s'_edges.csv", clear
rename (_sim _tail _head) (sim p x)
tempfile nets
save `nets'

local tag : env OUTTAG
local outf "`sim'/results/joint_s`s'`tag'.csv"
foreach r of local reps {
	use `nets', clear
	keep if sim == `r'
	mata: tp_bmat(st_data(., ("p", "x")), $TPN, $TPNP)
	capture nwdrop tp
	nwset, mat(__tpB) bipartite name(tp) nooutput
	matrix drop __tpB
	use `attrs', clear
	timer clear 1
	timer on 1
	local init ""
	if "`start'" == "cd" {
		capture noisily nwergm tp, `terms' `con' estimate(cd) cdnsteps(20) mcmleiterations(60) seed(`r') nologlik
		if !_rc {
			matrix bcd = e(b)
			forvalues k = 1/`=colsof(bcd)' {
				local init "`init' `=bcd[1,`k']'"
			}
		}
	}
	else {
		local init : env TRUEVALS
	}
	local ok 0
	if "`init'" != "" {
		capture noisily nwergm tp, `terms' `con' method(mcmle) init(`init') mcmleiterations(30) mcmletermination(`term') seed(`r') nologlik
		local ok = (_rc == 0)
	}
	timer off 1
	quietly timer list 1
	local secs = r(t1)
	preserve
	clear
	set obs 1
	gen scen = `s'
	gen rep = `r'
	gen fitted = `ok'
	gen seconds = `secs'
	if `ok' {
		gen converged = e(converged)
		gen iterations = e(mcmle_iterations)
		gen interval = e(mcmc_interval_final)
		matrix b = e(b)
		matrix V = e(V)
		local names : colnames b
		local k 0
		foreach nm of local names {
			local ++k
			local nm2 = substr("`nm'", 1, 28)
			gen b_`nm2' = b[1, `k']
			gen se_`nm2' = sqrt(V[`k', `k'])
		}
	}
	capture confirm file "`outf'"
	if _rc export delimited using "`outf'", replace
	else {
		tempfile one
		save `one'
		import delimited using "`outf'", clear
		append using `one'
		export delimited using "`outf'", replace
	}
	restore
	di as res "scenario `s' replication `r': fitted `ok', `secs' s"
}
