* ============================================================
* Simulation study, step 4: the conditional model - agencies modelled given the observed clubs.
* Two-mode network players (mode 1) x agencies (mode 2) with constraints(degrees): every player keeps one agency,
* every agency its portfolio size. With the clubs fixed, the joint model's agency-side statistics become:
*   S   = b1nodematch(club)     co-clients at the same club
*   H_A = b1nodematch(nat)      co-national co-clients
*   L   = b1nodematch(cty)      co-clients at clubs of the same country
*   T_W = edgecov(w_ap)         w_ap[p, a] = W[a, club(p)]
* Club-side statistics (H_C, home) are constant given the clubs and drop out.
* Env: SCEN, REPS as in sim_fit_joint.do. Appends to sim/results/cond_s<k>.csv.
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
capture mkdir "`sim'/results"

local s : env SCEN
local reps : env REPS
if "`reps'" == "" local reps 1

import delimited using "`sim'/world/nodes.csv", clear
quietly count if kind == 1
local NP = r(N)
quietly count if kind == 2
local NA = r(N)
mata: nat = st_data(1::`NP', "nat"); aghome = st_data((`NP'+1)::(`NP'+`NA'), "ag_home")
mata: clubcty = st_data((`NP'+`NA'+1)::`=_N', "club_country")

import delimited using "`sim'/nets/s`s'_edges.csv", clear
rename (_sim _tail _head) (sim p x)
tempfile nets
save `nets'

capture mata: mata drop tp_cond_setup()
mata:
// from the edge list of one simulated network: B (agencies x players), player club and country, and w_ap
void tp_cond_setup(real matrix E, real scalar NP, real scalar NA, real colvector aghome, real colvector clubcty)
{
	real matrix B, W
	real colvector club, cty
	real scalar k, p, x
	B = J(NA, NP, 0)
	club = J(NP, 1, .)
	for (k = 1; k <= rows(E); k++) {
		p = E[k, 1]; x = E[k, 2]
		if (x <= NP + NA) B[x - NP, p] = 1
		else club[p] = x - NP - NA
	}
	cty = clubcty[club]
	W = J(NP + NA, NP + NA, 0)
	for (p = 1; p <= NP; p++) {
		for (k = 1; k <= NA; k++) {
			if (aghome[k] == cty[p]) {
				W[p, NP + k] = 1
				W[NP + k, p] = 1
			}
		}
	}
	st_matrix("__B", B)
	st_matrix("__W", W)
	st_store(., "club", (club \ J(NA, 1, .)))
	st_store(., "cty", (cty \ J(NA, 1, .)))
}
end

local outf "`sim'/results/cond_s`s'.csv"
foreach r of local reps {
	use `nets', clear
	keep if sim == `r'
	mata: E = st_data(., ("p", "x"))
	clear
	set obs `=`NP'+`NA''
	gen club = .
	gen cty = .
	gen nat = .
	mata: st_store(1::`NP', "nat", nat)
	mata: tp_cond_setup(E, `NP', `NA', aghome, clubcty)
	capture nwdrop ap
	capture nwdrop w_ap
	nwset, mat(__B) bipartite name(ap) nooutput
	nwset, mat(__W) name(w_ap) nooutput
	matrix drop __B __W
	timer clear 1
	timer on 1
	capture noisily nwergm ap, edges bnodematch1(club) bnodematch1(nat) bnodematch1(cty) edgecov(w_ap) ///
		constraints(degrees) method(mcmle) mcmleiterations(30) seed(`r') nologlik
	local ok = (_rc == 0)
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
