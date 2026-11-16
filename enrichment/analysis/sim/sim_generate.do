* ============================================================
* Simulation study, step 2: simulate networks from the joint model for each scenario.
* Env: SCEN = scenario numbers (default 1 2 3 4); NSIM = networks per scenario (default 5);
*      INTERVAL, BURNIN = MCMC settings; STATSONLY=1 = only write the statistics (calibration).
* Writes sim/nets/s<k>_edges.csv (sim, p, x) and sim/nets/s<k>_stats.csv.
* Run from the project root: /Applications/StataNow/StataMP.app/Contents/MacOS/StataMP -e do enrichment/analysis/sim/sim_generate.do
* ============================================================
clear all
set more off
set maxvar 40000
local root "`c(pwd)'"
local sim "`root'/enrichment/analysis/sim"
local nwergm "`root'/software/nwergm-tripartite"
adopath ++ "`nwergm'"
adopath ++ "`nwergm'/lib/plugins/macos"
quietly do "`nwergm'/unw_ergm.do"
capture mkdir "`sim'/nets"

global TPDIR "`sim'/world"
global TPEDGES "`sim'/world/edges.csv"
global TPCOVS pc_home ac_home
do "`sim'/tp_load.do"

local terms edges bdsp1(2) bnodematch1(nat, byattr(type)) edgecov(pc_home) ///
	structural(tpcoloc(attr=key_cty) tpstarcov(x=ac_home attr=cls))
local con constraints(degrees + bd(attribs = isag iscl, maxout = bmax bmax, minout = bmin bmin))
* a short fit only to set up the model for simulation
nwergm tp, `terms' `con' estimate(cd) cdnsteps(2) cdsamplesize(32) mcmleiterations(1) seed(1) nologlik
matrix b0 = e(b)
local names : colnames b0
di "coefficients: `names'"

* scenario parameters by term: S = b1dsp2, HA/HC = b1nodematch_nat_1/_2, L = tpcoloc_key_cty,
* TW = tpstarcov_ac_home, HOME = edgecov_pc_home. Env PARAMS overrides one scenario: "S HA HC L TW HOME".
*                    S     HA    HC    L     TW    HOME
local P0             0     0     0     0     0     0
local P1             0     0.06  0.1   0     0     1.0
local P2             0.08  0.06  0.1   0     0     1.0
local P3             0     0     0.1   0.05  0     1.0
local P4             0.08  0.05  0.1   0.02  0.5   1.0
local pover : env PARAMS
local tnames b1dsp2 b1nodematch_nat_1 b1nodematch_nat_2 tpcoloc_key_cty tpstarcov_ac_home edgecov_pc_home
capture program drop _tp_setz
program define _tp_setz
	args vals names tnames
	matrix Z = J(1, `: word count `names'', 0)
	matrix colnames Z = `names'
	forvalues k = 1/6 {
		local nm : word `k' of `tnames'
		local col = colnumb(Z, "`nm'")
		if `col' >= . {
			di as err "term `nm' not in the model"
			exit 198
		}
		matrix Z[1, `col'] = `: word `k' of `vals''
	}
end

local scen : env SCEN
if "`scen'" == "" local scen 1 2 3 4
local nsim : env NSIM
if "`nsim'" == "" local nsim 5
local iv : env INTERVAL
if "`iv'" == "" local iv 1000000
local bi : env BURNIN
if "`bi'" == "" local bi 5000000
local statsonly : env STATSONLY

foreach s of local scen {
	local vals `P`s''
	if "`pover'" != "" local vals `pover'
	_tp_setz "`vals'" "`names'" "`tnames'"
	di as txt _n "scenario `s': " _c
	matrix list Z, noheader format(%6.3f)
	timer clear 1
	timer on 1
	nwergm simulate, nsim(`nsim') coef(Z) basis(tp) output(stats) burnin(`bi') interval(`iv') seed(`=1000+`s'')
	timer off 1
	timer list 1
	matrix S = r(stats)
	matrix colnames S = `names'
	preserve
	clear
	svmat S, names(col)
	gen sim = _n
	export delimited using "`sim'/nets/s`s'_stats.csv", replace
	summarize
	restore
	if "`statsonly'" == "" {
		preserve
		nwergm simulate, nsim(`nsim') coef(Z) basis(tp) output(edgelist) clear burnin(`bi') interval(`iv') seed(`=1000+`s'')
		export delimited using "`sim'/nets/s`s'_edges.csv", replace
		restore
	}
}
