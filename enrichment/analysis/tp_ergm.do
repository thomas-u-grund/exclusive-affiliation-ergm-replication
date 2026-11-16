* ============================================================
* Three-layer ERGM (agency - player - club) with nwergm, tripartite branch
*
* Two-mode network: players (mode 1) x affiliations (agencies, clubs). Each player has one agency and one
* club; the sample space keeps every portfolio size, every squad size and the one-agency-one-club rule
* (constraints degrees + bd, typed-swap proposal). Neither layer is fixed.
*
* Run from the project root: /Applications/StataNow/StataMP.app/Contents/MacOS/StataMP -e do enrichment/analysis/tp_ergm.do
* Data: enrichment/analysis/ergm_tp_data.R -> enrichment/analysis/output/tp/
* Env QUICK=1: staged fit (env STAGE, INIT, INTERVAL; see the block below)
* ============================================================
clear all
set more off
set maxvar 40000
* Env NPROC = cap on Stata/MP cores (keeps the shared machine usable)
local nproc : env NPROC
if "`nproc'" != "" set processors `nproc'
local root "`c(pwd)'"
local tp "`root'/enrichment/analysis/output/tp"
local nwergm "`root'/software/nwergm-tripartite"
adopath ++ "`nwergm'"
adopath ++ "`nwergm'/lib/plugins/macos"
quietly do "`nwergm'/unw_ergm.do"
local quick : env QUICK

* ---------------------------------------------------------------- network
import delimited using "`tp'/nodes.csv", clear
local N = _N
quietly count if kind == 1
local NP = r(N)
mata: nodes = st_data(., ("node", "kind"))
import delimited using "`tp'/edges.csv", clear
mata:
E = st_data(., ("p", "x"))
NP = strtoreal(st_local("NP")); N = strtoreal(st_local("N"))
B = J(N - NP, NP, 0)
for (k = 1; k <= rows(E); k++) B[E[k, 2] - NP, E[k, 1]] = 1
assert(colsum(B) == J(1, NP, 2))
end
nwset, mat(B) bipartite name(tp) nooutput
mata: mata drop B

* ---------------------------------------------------------------- dyadic covariates as valued one-mode networks
mata:
real matrix tp_densify(real matrix T, real scalar n)
{
	real matrix C
	real scalar k
	C = J(n, n, 0)
	for (k = 1; k <= rows(T); k++) {
		C[T[k, 1], T[k, 2]] = T[k, 3]
		C[T[k, 2], T[k, 1]] = T[k, 3]
	}
	return(C)
}
end
local covs : env COVS
if "`covs'" == "" local covs ap_corr ap_val pc_home pc_val ac_home ac_val
foreach c of local covs {
	import delimited using "`tp'/cov_`c'.csv", clear
	mata: C = tp_densify(st_data(., ("i", "j", "v")), `N')
	nwset, mat(C) name(`c') nooutput
	mata: mata drop C
	di "covariate `c' loaded"
}

* ---------------------------------------------------------------- node attributes (rows = nodes)
import delimited using "`tp'/nodes.csv", clear
assert _N == `N' & node == _n

local terms edges bdsp1(2) bnodematch1(nat, byattr(type)) edgecov(ap_corr ap_val pc_home pc_val) ///
	structural(tpcoloc(attr=key_cty) tpstarcov(x=ac_home attr=cls) tpstarcov(x=ac_val attr=cls))
local stage : env STAGE
if "`stage'" == "1" local terms edges bdsp1(2) bnodematch1(nat, byattr(type)) structural(tpcoloc(attr=key_cty))
* Nested stages (paper Section 6.2; redefined 2026-11-13, the earlier stages 2-3 were reduced timing models):
*   1 structural core; 2 = 1 + agency-club covariates (corridor ac_home, value gap ac_val); 3 = 2 + agency-player
*   and player-club covariates = the full model above. Stage 2 needs COVS="ac_home ac_val".
if "`stage'" == "2" local terms edges bdsp1(2) bnodematch1(nat, byattr(type)) ///
	structural(tpcoloc(attr=key_cty) tpstarcov(x=ac_home attr=cls) tpstarcov(x=ac_val attr=cls))
* Market-presence variants (2026-11-13): ac_mkt = share of the agency's outside clients in the club's country
* (ergm_tp_mkt.R). 1m = stage 1 + ac_mkt (COVS="ac_mkt"); 3m = full model + ac_mkt (COVS = all seven).
if "`stage'" == "1m" local terms edges bdsp1(2) bnodematch1(nat, byattr(type)) ///
	structural(tpcoloc(attr=key_cty) tpstarcov(x=ac_mkt attr=cls))
if "`stage'" == "3m" local terms edges bdsp1(2) bnodematch1(nat, byattr(type)) edgecov(ap_corr ap_val pc_home pc_val) ///
	structural(tpcoloc(attr=key_cty) tpstarcov(x=ac_home attr=cls) tpstarcov(x=ac_val attr=cls) tpstarcov(x=ac_mkt attr=cls))
* Size-dependent stacking (2026-11-13): 3ms = 3m + extra stacking among agencies with < 30 clients, via tpcoloc with
* every club its own key and small agencies as the only centres (attr key_clubs; COVS = all seven as for 3m).
if "`stage'" == "3ms" local terms edges bdsp1(2) bnodematch1(nat, byattr(type)) edgecov(ap_corr ap_val pc_home pc_val) ///
	structural(tpcoloc(attr=key_cty) tpstarcov(x=ac_home attr=cls) tpstarcov(x=ac_val attr=cls) tpstarcov(x=ac_mkt attr=cls) tpcoloc(attr=key_clubs))
* 3ms2 = the same model in separate stacking terms (better conditioned): no common bdsp1(2); stacking among small
* agencies (key_clubs: small = -1) and among large agencies (key_clubl: large = -1), each club its own key.
if "`stage'" == "3ms2" local terms edges bnodematch1(nat, byattr(type)) edgecov(ap_corr ap_val pc_home pc_val) ///
	structural(tpcoloc(attr=key_cty) tpstarcov(x=ac_home attr=cls) tpstarcov(x=ac_val attr=cls) tpstarcov(x=ac_mkt attr=cls) tpcoloc(attr=key_clubs) tpcoloc(attr=key_clubl))
local con constraints(degrees + bd(attribs = isag iscl, maxout = bmax bmax, minout = bmin bmin))

if "`quick'" != "" {
	* Env INIT = starting values (skips CD); env INTERVAL = fixed MCMC interval (default: adaptive effective size)
	* The log-likelihood is skipped: under the degrees + bd constraint nwergm's bridge reference model is not valid.
	timer on 1
	local init : env INIT
	if "`init'" == "" {
		nwergm tp, `terms' `con' estimate(cd) cdnsteps(20) mcmleiterations(60) seed(20261107) nologlik
		matrix bcd = e(b)
		* Env INITNAMED = "name=value ..." overrides CD start values by coefficient name (e.g. structural terms
		* from a converged smaller model); covariates keep their CD values.
		local initnamed : env INITNAMED
		foreach nv of local initnamed {
			gettoken nm val : nv, parse("=")
			local val = substr("`val'", 2, .)
			local col = colnumb(bcd, "`nm'")
			if `col' >= . {
				di as err "INITNAMED: `nm' not in the model"
				exit 198
			}
			matrix bcd[1, `col'] = `val'
		}
		matrix list bcd, format(%9.4f)
		forvalues k = 1/`=colsof(bcd)' {
			local init "`init' `=bcd[1,`k']'"
		}
	}
	di "start: `init'"
	* Env TERM = MCMLE termination rule (default hummel: stop at full Hummel steps; confidence = R's multivariate test)
local term : env MCTERM
if "`term'" == "" local term hummel
local iv : env INTERVAL
	local ivopt ""
	if "`iv'" != "" local ivopt mcmcinterval(`iv') mcmcburnin(`=20*`iv'') mcmleeffectivesize(0) mcmcsamplesize(500)
	nwergm tp, `terms' `con' method(mcmle) init(`init') mcmleiterations(30) mcmletermination(`term') `ivopt' seed(20261107) verbose nologlik
	timer off 1
	timer list 1
	local tag : env TAG
	estimates save "`tp'/tp_stage`stage'`tag'", replace
	mata: st_matrix("OBS", ergm_full_stat(__nwergm_last_M, __nwergm_last_G))
	matrix colnames OBS = `: colnames e(b)'
	matrix list OBS, format(%12.3f)
	exit
}

nwergm tp, `terms' `con' method(mcmle) seed(20261107)
estimates save "`tp'/tp_m1", replace
estat mcmcdiag
