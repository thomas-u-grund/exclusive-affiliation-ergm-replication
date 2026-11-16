* ============================================================
* Load a three-layer network in the output/tp format into nwergm.
* Arguments (globals): TPDIR = folder with nodes.csv, cov_*.csv; TPEDGES = edge file (p, x); TPCOVS = covariates.
* Leaves: network "tp", covariate networks named after the files, node attributes in memory; locals N, NP in
* globals TPN, TPNP. Needs nwergm (tripartite branch) on the adopath and unw_ergm.do loaded.
* ============================================================
import delimited using "$TPDIR/nodes.csv", clear
global TPN = _N
quietly count if kind == 1
global TPNP = r(N)
import delimited using "$TPEDGES", clear
capture mata: mata drop tp_densify()
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
void tp_bmat(real matrix E, real scalar N, real scalar NP)
{
	real matrix B
	real scalar k
	B = J(N - NP, NP, 0)
	for (k = 1; k <= rows(E); k++) B[E[k, 2] - NP, E[k, 1]] = 1
	if (colsum(B) != J(1, NP, 2)) _error("every player needs exactly one agency and one club")
	st_matrix("__tpB", B)
}
end
mata: tp_bmat(st_data(., ("p", "x")), $TPN, $TPNP)
capture nwdrop tp
nwset, mat(__tpB) bipartite name(tp) nooutput
matrix drop __tpB
foreach c of global TPCOVS {
	import delimited using "$TPDIR/cov_`c'.csv", clear
	mata: st_matrix("__tpC", tp_densify(st_data(., ("i", "j", "v")), $TPN))
	capture nwdrop `c'
	nwset, mat(__tpC) name(`c') nooutput
	matrix drop __tpC
}
import delimited using "$TPDIR/nodes.csv", clear
assert _N == $TPN & node == _n
