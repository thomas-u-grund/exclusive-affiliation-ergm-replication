* ============================================================
* Simulation study, scenario 5: generate networks with stacking and agency heterogeneity (sim_world5.R).
* Env CHUNK (seed and file tag), NSIM (networks in this chunk), INTERVAL, BURNIN, NPROC.
* Writes sim/nets/s5_c<CHUNK>.csv (_sim, _tail, _head); sim_combine5.R merges the chunks into s5_edges.csv.
* Run from a folder without spaces in its path.
* ============================================================
clear all
set more off
set maxvar 40000
local nproc : env NPROC
if "`nproc'" != "" set processors `nproc'
local root "`c(pwd)'"
local sim "`root'/enrichment/analysis/sim"
local nwergm "`root'/software/nwergm-tripartite"
adopath ++ "`nwergm'"
adopath ++ "`nwergm'/lib/plugins/macos"
quietly do "`nwergm'/unw_ergm.do"

global TPDIR "`sim'/world"
global TPEDGES "`sim'/world/edges.csv"
global TPCOVS pc_home ap_spec ac_home ac_mkt
do "`sim'/tp_load.do"

local terms edges bdsp1(2) bnodematch1(nat, byattr(type)) edgecov(pc_home ap_spec) ///
	structural(tpcoloc(attr=key_cty) tpstarcov(x=ac_home attr=cls) tpstarcov(x=ac_mkt attr=cls))
local con constraints(degrees + bd(attribs = isag iscl, maxout = bmax bmax, minout = bmin bmin))
nwergm tp, `terms' `con' estimate(cd) cdnsteps(2) cdsamplesize(32) mcmleiterations(1) seed(1) nologlik
matrix Z = e(b)
local names : colnames Z
matrix Z = J(1, colsof(Z), 0)
matrix colnames Z = `names'
*      S     HA    HC    L     TW   HOME  SPEC  MKT
local vals 0.08 0.03 0.1 0 0 1.0 2.0 1.5
local tn b1dsp2 b1nodematch_nat_1 b1nodematch_nat_2 tpcoloc_key_cty tpstarcov_ac_home edgecov_pc_home edgecov_ap_spec tpstarcov_ac_mkt
forvalues k = 1/8 {
	local col = colnumb(Z, "`: word `k' of `tn''")
	assert `col' < .
	matrix Z[1, `col'] = `: word `k' of `vals''
}
matrix list Z, format(%6.3f)

local chunk : env CHUNK
local nsim : env NSIM
local iv : env INTERVAL
if "`iv'" == "" local iv 1000000
local bi : env BURNIN
if "`bi'" == "" local bi 5000000
nwergm simulate, nsim(`nsim') coef(Z) basis(tp) output(edgelist) clear burnin(`bi') interval(`iv') seed(`=5000+`chunk'')
export delimited using "`sim'/nets/s5_c`chunk'.csv", replace
