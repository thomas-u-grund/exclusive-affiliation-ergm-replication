* ============================================================
* Moment-matching refinement of a three-layer ERGM fit (2026-11-13).
* Why: nwergm's MCMLE stopped at estimates whose simulated statistics do not match the observed ones (e.g. model 1:
* simulated S = 1156 and L = 16,970 against 873 and 11,800), because the MCMC interval it used was far shorter than the
* mixing time of the sampler (about 1e6 typed swaps on this network). At the MLE the expected statistics equal the observed.
* Method: Newton iterations theta <- theta + damp * Sigma^{-1} (g_obs - mean g_sim), with Sigma the covariance of the
* simulated statistics (draws at a long interval, nearly independent), then Polyak averaging of the last iterates and a
* final simulation at the average for the standard errors, standardised moment discrepancies, and a refitted .ster file.
* Env: MODEL (term set 1, 3m, 3ms2 ...), FIT (starting .ster tag), OUT (tag of the refined .ster: tp_stage<FIT><OUT>.ster),
*      ITER (default 12), NSIM (default 60), INTERVAL (default 500000), BURNIN (default 5000000), DAMP (default 0.5),
*      FINALSIM (default 120), NPROC.
* Writes output/tp/refine/trace_<FIT><OUT>.csv (iteration, theta, standardised discrepancies) and final_<FIT><OUT>.csv.
* Run from a folder without spaces: StataMP -e do <path>/refine_tp.do
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
capture mkdir "`tp'/refine"

local model : env MODEL
local fit : env FIT
local out : env OUT
if "`out'" == "" local out R
local niter : env ITER
if "`niter'" == "" local niter 12
local nsim : env NSIM
if "`nsim'" == "" local nsim 60
local iv : env INTERVAL
if "`iv'" == "" local iv 500000
local bi : env BURNIN
if "`bi'" == "" local bi 5000000
local damp : env DAMP
if "`damp'" == "" local damp 0.5
local nfin : env FINALSIM
if "`nfin'" == "" local nfin 120

local t1  edges bdsp1(2) bnodematch1(nat, byattr(type)) structural(tpcoloc(attr=key_cty))
local c1
local t3m edges bdsp1(2) bnodematch1(nat, byattr(type)) edgecov(ap_corr ap_val pc_home pc_val) structural(tpcoloc(attr=key_cty) tpstarcov(x=ac_home attr=cls) tpstarcov(x=ac_val attr=cls) tpstarcov(x=ac_mkt attr=cls))
local c3m ap_corr ap_val pc_home pc_val ac_home ac_val ac_mkt
local t3ms2 edges bnodematch1(nat, byattr(type)) edgecov(ap_corr ap_val pc_home pc_val) structural(tpcoloc(attr=key_cty) tpstarcov(x=ac_home attr=cls) tpstarcov(x=ac_val attr=cls) tpstarcov(x=ac_mkt attr=cls) tpcoloc(attr=key_clubs) tpcoloc(attr=key_clubl))
local c3ms2 ap_corr ap_val pc_home pc_val ac_home ac_val ac_mkt

global TPDIR "`tp'"
global TPEDGES "`tp'/edges.csv"
global TPCOVS `c`model''
do "`root'/enrichment/analysis/sim/tp_load.do"
local con constraints(degrees + bd(attribs = isag iscl, maxout = bmax bmax, minout = bmin bmin))
nwergm tp, `t`model'' `con' estimate(cd) cdnsteps(2) cdsamplesize(32) mcmleiterations(1) seed(1) nologlik
matrix b0 = e(b)
local names : colnames b0
mata: st_matrix("GOBS", ergm_full_stat(__nwergm_last_M, __nwergm_last_G))
estimates use "`tp'/tp_stage`fit'"
matrix TH = e(b)
local fnames : colnames TH
assert "`names'" == "`fnames'"
local K = colsof(TH)

local thcsv : env THCSV
if "`thcsv'" != "" {
	preserve
	import delimited using "`thcsv'", clear
	assert _N == `K'
	forvalues k = 1/`K' {
		assert name[`k'] == "`: word `k' of `names''"
		matrix TH[1,`k'] = theta[`k']
	}
	restore
}

* ---- Mata: Newton step from a simulated statistics matrix
mata:
void refine_step(real scalar damp, real scalar it, string scalar tracefile)
{
	real matrix S, V
	real rowvector gobs, th, free, m, d, z
	real colvector step
	real scalar K, k, nf
	string rowvector nm
	S = st_matrix("r(stats)"); gobs = st_matrix("GOBS"); th = st_matrix("TH")
	nm = tokens(st_local("names")); K = cols(th)
	free = select(1..K, nm :!= "edges")
	nf = cols(free)
	m = mean(S[., free]); V = variance(S[., free]); d = gobs[free] - m
	z = d :/ (sqrt(diagonal(V))' :/ sqrt(rows(S)))
	step = invsym(V) * d'
	if (max(abs(step)) > 0.1) step = step * (0.1 / max(abs(step)))   // cap the step: at most 0.1 per iteration (before damping)
	th[free] = th[free] + damp * step'
	st_matrix("TH", th)
	st_matrix("ZZ", z)
	st_matrix("DD", d)
	st_matrix("MM", m)
	st_matrix("VV", V)
	printf("iteration %g: max |z| = %5.1f, max |step| = %7.4f\n", it, max(abs(z)), max(abs(damp * step)))
}
end

tempname post
local cols "it"
forvalues k = 1/`K' {
	local cols "`cols' th`k'"
}
local KF = `K' - 1
forvalues k = 1/`KF' {
	local cols "`cols' z`k'"
}
postfile `post' `cols' using "`tp'/refine/trace_`fit'`out'", replace
tempname THh
forvalues it = 1/`niter' {
	nwergm simulate, nsim(`nsim') coef(TH) basis(tp) output(stats) burnin(`bi') interval(`iv') seed(`=9000+`it'')
	mata: refine_step(`damp', `it', "")
	local row "(`it')"
	forvalues k = 1/`K' {
		local row "`row' (`=TH[1,`k']')"
	}
	forvalues k = 1/`KF' {
		local row "`row' (`=ZZ[1,`k']')"
	}
	post `post' `row'
	matrix THh`it' = TH
	matrix list TH, format(%9.5f)
}
postclose `post'

* ---- Polyak average of the last iterates (burn-in half of the iterations)
local first = max(2, ceil(`niter' / 2) + 1)
matrix THA = J(1, `K', 0)
local n 0
forvalues it = `first'/`niter' {
	matrix THA = THA + THh`it'
	local ++n
}
matrix THA = THA / `n'
matrix colnames THA = `names'
matrix TH = THA
matrix list TH, format(%9.5f)

* ---- final simulation at the average: moment discrepancies and standard errors
nwergm simulate, nsim(`nfin') coef(TH) basis(tp) output(stats) burnin(`bi') interval(`iv') seed(9999)
mata: refine_step(0, 0, "")
* SE from the inverse covariance of the statistics (Fisher information at the average)
mata:
S = st_matrix("r(stats)")
nmx = tokens(st_local("names")); K = cols(S)
free = select(1..K, nmx :!= "edges")
V = variance(S[., free])
Vfull = J(K, K, 0)
Vfull[free, free] = invsym(V)
st_matrix("VCOV", Vfull)
end
matrix colnames VCOV = `names'
matrix rownames VCOV = `names'
preserve
clear
local k 0
set obs `K'
gen str40 name = ""
gen theta = .
gen se = .
gen z_moment = .
local fi 0
forvalues k = 1/`K' {
	local nm : word `k' of `names'
	replace name = "`nm'" in `k'
	replace theta = TH[1,`k'] in `k'
	replace se = sqrt(VCOV[`k',`k']) in `k'
	if "`nm'" != "edges" {
		local ++fi
		replace z_moment = ZZ[1,`fi'] in `k'
	}
}
export delimited using "`tp'/refine/final_`fit'`out'.csv", replace
list, clean
restore
