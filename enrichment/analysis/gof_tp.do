* ============================================================
* GOF for the three-layer ERGM: simulate networks from a fitted stage and export them as edge lists.
* Env MODEL = stage tag of a fit saved by tp_ergm.do (output/tp/tp_stage<MODEL>.ster): 1, 2, 3, 1m, 3m.
* Env NSIM (default 100), INTERVAL (default 1e6), BURNIN (default 1e7), NPROC.
* Writes output/tp/gof/gof_<MODEL>_edges.csv (_sim, _tail, _head); summarised by gof_tp.R.
* Run from a folder without spaces in its path: StataMP -e do <relative path>/gof_tp.do
* ============================================================
clear all
set more off
set maxvar 40000
local nproc : env NPROC
if "`nproc'" != "" set processors `nproc'
local root "`c(pwd)'"
local tp "`root'/enrichment/analysis/output/tp"
local nwergm "`root'/software/nwergm-tripartite"
adopath ++ "`nwergm'"
adopath ++ "`nwergm'/lib/plugins/macos"
quietly do "`nwergm'/unw_ergm.do"
capture mkdir "`tp'/gof"

local model : env MODEL
local nsim : env NSIM
if "`nsim'" == "" local nsim 100
local iv : env INTERVAL
if "`iv'" == "" local iv 1000000
local bi : env BURNIN
if "`bi'" == "" local bi 10000000

* term sets: must match tp_ergm.do
local t1  edges bdsp1(2) bnodematch1(nat, byattr(type)) structural(tpcoloc(attr=key_cty))
local c1
local t2  edges bdsp1(2) bnodematch1(nat, byattr(type)) structural(tpcoloc(attr=key_cty) tpstarcov(x=ac_home attr=cls) tpstarcov(x=ac_val attr=cls))
local c2  ac_home ac_val
local t3  edges bdsp1(2) bnodematch1(nat, byattr(type)) edgecov(ap_corr ap_val pc_home pc_val) structural(tpcoloc(attr=key_cty) tpstarcov(x=ac_home attr=cls) tpstarcov(x=ac_val attr=cls))
local c3  ap_corr ap_val pc_home pc_val ac_home ac_val
local t1m edges bdsp1(2) bnodematch1(nat, byattr(type)) structural(tpcoloc(attr=key_cty) tpstarcov(x=ac_mkt attr=cls))
local c1m ac_mkt
local t3m edges bdsp1(2) bnodematch1(nat, byattr(type)) edgecov(ap_corr ap_val pc_home pc_val) structural(tpcoloc(attr=key_cty) tpstarcov(x=ac_home attr=cls) tpstarcov(x=ac_val attr=cls) tpstarcov(x=ac_mkt attr=cls))
local c3m ap_corr ap_val pc_home pc_val ac_home ac_val ac_mkt

local t3ms edges bdsp1(2) bnodematch1(nat, byattr(type)) edgecov(ap_corr ap_val pc_home pc_val) structural(tpcoloc(attr=key_cty) tpstarcov(x=ac_home attr=cls) tpstarcov(x=ac_val attr=cls) tpstarcov(x=ac_mkt attr=cls) tpcoloc(attr=key_clubs))
local c3ms ap_corr ap_val pc_home pc_val ac_home ac_val ac_mkt

local t3ms2 edges bnodematch1(nat, byattr(type)) edgecov(ap_corr ap_val pc_home pc_val) structural(tpcoloc(attr=key_cty) tpstarcov(x=ac_home attr=cls) tpstarcov(x=ac_val attr=cls) tpstarcov(x=ac_mkt attr=cls) tpcoloc(attr=key_clubs) tpcoloc(attr=key_clubl))
local c3ms2 ap_corr ap_val pc_home pc_val ac_home ac_val ac_mkt

global TPDIR "`tp'"
global TPEDGES "`tp'/edges.csv"
global TPCOVS `c`model''
do "`root'/enrichment/analysis/sim/tp_load.do"

local con constraints(degrees + bd(attribs = isag iscl, maxout = bmax bmax, minout = bmin bmin))
* a minimal fit only to set up the model for simulation (as in sim_generate.do)
nwergm tp, `t`model'' `con' estimate(cd) cdnsteps(2) cdsamplesize(32) mcmleiterations(1) seed(1) nologlik
matrix b0 = e(b)
local names : colnames b0
local fit : env FIT
if "`fit'" == "" local fit `model'
estimates use "`tp'/tp_stage`fit'"
matrix Z = e(b)
local fnames : colnames Z
assert "`names'" == "`fnames'"
local K = colsof(Z)

local thcsv : env THCSV
if "`thcsv'" != "" {
	preserve
	import delimited using "`thcsv'", clear
	assert _N == `K'
	forvalues k = 1/`K' {
		assert name[`k'] == "`: word `k' of `names''"
		matrix Z[1,`k'] = theta[`k']
	}
	restore
}
matrix list Z, format(%9.4f)
nwergm simulate, nsim(`nsim') coef(Z) basis(tp) output(edgelist) clear burnin(`bi') interval(`iv') seed(20261113)
export delimited using "`tp'/gof/gof_`fit'_edges.csv", replace
