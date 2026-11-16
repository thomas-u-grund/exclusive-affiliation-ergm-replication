/* Dedicated sampler for exclusive two-affiliation networks (every player has one agency and one club).
 *
 * The state is the pair of assignments alpha (agency of each player) and gamma (club). All supported statistics are functions of
 * the cell counts n[a][c] (clients of agency a at club c) and of a few count tables derived from them, so a change of partner
 * updates every statistic in O(1) per term. A typed swap exchanges the agencies or the clubs of two players: four such
 * "remove player from (a,c) / add player to (a',c')" operations, accepted with probability min(1, exp(theta' dg)). The
 * proposal is a 50/50 mixture of (i) the plain typed swap and (ii) a swap with a partner drawn among the players that share
 * one of the first player's partners of the other type; both are symmetric.
 *
 * Statistics (columns of the output, in the order given by the column indices):
 *   stack      sum_{a,c} C(n_ac, 2)
 *   coloc_j    sum_a w_ja sum_k C(T_jak, 2), T_jak = clients of a at clubs with key k (key 0 = none)
 *   hom_j      sum over centres of sum_v C(H_v, 2), H_v = members with attribute v (agency or club side)
 *   acov_j     sum_p W_j[agency(p), club(p)]
 *   edgecov_j  sum_p M_j[p, agency(p)] + M_j[p, club(p)]
 */
#include <R.h>
#include <Rinternals.h>
#include <stdint.h>
#include <string.h>
#include <stdlib.h>

typedef struct { uint64_t s[4]; } rng_t;
static inline uint64_t rotl(uint64_t x, int k) { return (x << k) | (x >> (64 - k)); }
static inline uint64_t rnext(rng_t *r) {
  uint64_t *s = r->s, res = rotl(s[1] * 5, 7) * 9, t = s[1] << 17;
  s[2] ^= s[0]; s[3] ^= s[1]; s[1] ^= s[2]; s[0] ^= s[3]; s[2] ^= t; s[3] = rotl(s[3], 45); return res;
}
static inline double runif01(rng_t *r) { return (rnext(r) >> 11) * (1.0 / 9007199254740992.0); }
static inline uint32_t rbelow(rng_t *r, uint32_t n) { return (uint32_t) (((rnext(r) >> 32) * (uint64_t) n) >> 32); }
static void rseed(rng_t *r, uint64_t seed) {
  for (int i = 0; i < 4; i++) { seed += 0x9e3779b97f4a7c15ULL; uint64_t z = seed; z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL; z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL; r->s[i] = z ^ (z >> 31); }
}

typedef struct {
  int n1, A, C, n2, nstat, stack_col;
  int *alpha, *gamma, *nat;           /* per player: local agency index, local club index, attribute */
  int *agnode, *clnode;               /* mode-2 index of every agency / club */
  int *n;                             /* A x C cell counts */
  int ncol, nhom, nacov, nedge;
  int *col_col, *col_nk; double **col_w; int **col_key; int **col_T;     /* coloc terms */
  int *hom_col, *hom_side, *hom_V; double **hom_cen; int **hom_H;       /* homcentre terms */
  int *acov_col, *edge_col; double *W, *M;                              /* interleaved W[(x*n2+y)*nacov+j], M[(p*n2+x)*nedge+j] */
  int *agstart, *agList, *agPos, *clstart, *clList, *clPos;             /* member lists */
} st_t;

/* add (sgn=+1) or remove (sgn=-1) player p at (a,c); dg (may be NULL) receives the change of every statistic */
static inline void op(st_t *s, int p, int a, int c, int add, double *dg) {
  const int C = s->C;
  if (s->stack_col >= 0) {
    int *x = &s->n[a * C + c];
    if (add) { if (dg) dg[s->stack_col] += *x; (*x)++; } else { (*x)--; if (dg) dg[s->stack_col] -= *x; }
  }
  for (int j = 0; j < s->ncol; j++) {
    int k = s->col_key[j][c]; double w = s->col_w[j][a];
    if (k > 0 && w != 0) {
      int *x = &s->col_T[j][a * (s->col_nk[j] + 1) + k];
      if (add) { if (dg) dg[s->col_col[j]] += w * (*x); (*x)++; } else { (*x)--; if (dg) dg[s->col_col[j]] -= w * (*x); }
    }
  }
  for (int j = 0; j < s->nhom; j++) {
    int node = s->hom_side[j] == 1 ? a : c;
    if (s->hom_cen[j][node] != 0) {
      int *x = &s->hom_H[j][node * s->hom_V[j] + s->nat[p]];
      if (add) { if (dg) dg[s->hom_col[j]] += *x; (*x)++; } else { (*x)--; if (dg) dg[s->hom_col[j]] -= *x; }
    }
  }
  if (dg) {
    double sg = add ? 1.0 : -1.0;
    if (s->nacov) { const double *w = s->W + ((size_t) s->agnode[a] * s->n2 + s->clnode[c]) * s->nacov; for (int j = 0; j < s->nacov; j++) dg[s->acov_col[j]] += sg * w[j]; }
    if (s->nedge) {
      const double *m1 = s->M + ((size_t) p * s->n2 + s->agnode[a]) * s->nedge, *m2 = s->M + ((size_t) p * s->n2 + s->clnode[c]) * s->nedge;
      for (int j = 0; j < s->nedge; j++) dg[s->edge_col[j]] += sg * (m1[j] + m2[j]);
    }
  }
}

static inline void swap_members(int *list, int *pos, int p, int q) { int ip = pos[p], iq = pos[q]; list[ip] = q; list[iq] = p; pos[p] = iq; pos[q] = ip; }

/* one proposal; returns 1 if accepted */
static int step(st_t *s, rng_t *r, const double *theta, double *dg, double hold) {
  int p = (int) rbelow(r, s->n1), q, t = (int) rbelow(r, 2);    /* t = 0: swap agencies, 1: swap clubs */
  if (rnext(r) & 1) {
    int g = t == 0 ? s->gamma[p] : s->alpha[p];                  /* group: players sharing p's partner of the other type */
    int *list = t == 0 ? s->clList : s->agList, *pos = t == 0 ? s->clPos : s->agPos, st = t == 0 ? s->clstart[g] : s->agstart[g];
    int m = (t == 0 ? s->clstart[g + 1] : s->agstart[g + 1]) - st;
    if (m < 2) return 0;
    int rr = (int) rbelow(r, (uint32_t) (m - 1)); if (rr >= pos[p] - st) rr++;
    q = list[st + rr];
  } else { q = (int) rbelow(r, (uint32_t) (s->n1 - 1)); if (q >= p) q++; }
  int a = s->alpha[p], b = s->alpha[q], c = s->gamma[p], d = s->gamma[q];
  if (t == 0 ? a == b : c == d) return 0;
  (void) hold;
  for (int j = 0; j < s->nstat; j++) dg[j] = 0;
  if (t == 0) { op(s, p, a, c, 0, dg); op(s, q, b, d, 0, dg); op(s, p, b, c, 1, dg); op(s, q, a, d, 1, dg); }
  else        { op(s, p, a, c, 0, dg); op(s, q, b, d, 0, dg); op(s, p, a, d, 1, dg); op(s, q, b, c, 1, dg); }
  double lr = 0; for (int j = 0; j < s->nstat; j++) lr += theta[j] * dg[j];
  if (lr >= 0 || runif01(r) < exp(lr)) {
    if (t == 0) { s->alpha[p] = b; s->alpha[q] = a; swap_members(s->agList, s->agPos, p, q); }
    else        { s->gamma[p] = d; s->gamma[q] = c; swap_members(s->clList, s->clPos, p, q); }
    return 1;
  }
  if (t == 0) { op(s, p, b, c, 0, NULL); op(s, q, a, d, 0, NULL); op(s, p, a, c, 1, NULL); op(s, q, b, d, 1, NULL); }
  else        { op(s, p, a, d, 0, NULL); op(s, q, b, c, 0, NULL); op(s, p, a, c, 1, NULL); op(s, q, b, d, 1, NULL); }
  for (int j = 0; j < s->nstat; j++) dg[j] = 0;
  return 0;
}

#define GETINT(x) INTEGER(x)
#define LEN(x) ((int) XLENGTH(x))

/* arguments: spec (list, see R/fast.R), theta, nsim, burnin, interval, seed */
SEXP ex_fast_c(SEXP spec, SEXP theta_, SEXP nsim_, SEXP burnin_, SEXP interval_, SEXP seed_) {
  st_t S; memset(&S, 0, sizeof S);
  #define EL(name) VECTOR_ELT(spec, name)
  S.n1 = INTEGER(EL(0))[0]; S.A = INTEGER(EL(1))[0]; S.C = INTEGER(EL(2))[0]; S.n2 = S.A + S.C; S.nstat = INTEGER(EL(3))[0];
  int *alpha0 = INTEGER(EL(4)), *gamma0 = INTEGER(EL(5)), *nat0 = INTEGER(EL(6));
  S.agnode = INTEGER(EL(7)); S.clnode = INTEGER(EL(8));
  S.stack_col = INTEGER(EL(9))[0];
  SEXP cols = EL(10), homs = EL(11), acs = EL(12), eds = EL(13);
  S.alpha = (int *) R_alloc(S.n1, sizeof(int)); S.gamma = (int *) R_alloc(S.n1, sizeof(int)); S.nat = nat0;
  for (int p = 0; p < S.n1; p++) { S.alpha[p] = alpha0[p]; S.gamma[p] = gamma0[p]; }
  S.n = (int *) R_alloc((size_t) S.A * S.C, sizeof(int)); memset(S.n, 0, (size_t) S.A * S.C * sizeof(int));
  for (int p = 0; p < S.n1; p++) S.n[S.alpha[p] * S.C + S.gamma[p]]++;
  /* coloc */
  S.ncol = LEN(cols);
  S.col_col = (int *) R_alloc(S.ncol + 1, sizeof(int)); S.col_nk = (int *) R_alloc(S.ncol + 1, sizeof(int));
  S.col_w = (double **) R_alloc(S.ncol + 1, sizeof(double *)); S.col_key = (int **) R_alloc(S.ncol + 1, sizeof(int *)); S.col_T = (int **) R_alloc(S.ncol + 1, sizeof(int *));
  for (int j = 0; j < S.ncol; j++) {
    SEXP e = VECTOR_ELT(cols, j); S.col_col[j] = INTEGER(VECTOR_ELT(e, 0))[0]; S.col_w[j] = REAL(VECTOR_ELT(e, 1)); S.col_key[j] = INTEGER(VECTOR_ELT(e, 2)); S.col_nk[j] = INTEGER(VECTOR_ELT(e, 3))[0];
    size_t sz = (size_t) S.A * (S.col_nk[j] + 1); S.col_T[j] = (int *) R_alloc(sz, sizeof(int)); memset(S.col_T[j], 0, sz * sizeof(int));
    for (int p = 0; p < S.n1; p++) { int k = S.col_key[j][S.gamma[p]]; if (k > 0) S.col_T[j][S.alpha[p] * (S.col_nk[j] + 1) + k]++; }
  }
  /* homcentre */
  S.nhom = LEN(homs);
  S.hom_col = (int *) R_alloc(S.nhom + 1, sizeof(int)); S.hom_side = (int *) R_alloc(S.nhom + 1, sizeof(int)); S.hom_V = (int *) R_alloc(S.nhom + 1, sizeof(int));
  S.hom_cen = (double **) R_alloc(S.nhom + 1, sizeof(double *)); S.hom_H = (int **) R_alloc(S.nhom + 1, sizeof(int *));
  for (int j = 0; j < S.nhom; j++) {
    SEXP e = VECTOR_ELT(homs, j); S.hom_col[j] = INTEGER(VECTOR_ELT(e, 0))[0]; S.hom_side[j] = INTEGER(VECTOR_ELT(e, 1))[0]; S.hom_cen[j] = REAL(VECTOR_ELT(e, 2)); S.hom_V[j] = INTEGER(VECTOR_ELT(e, 3))[0];
    int nn = S.hom_side[j] == 1 ? S.A : S.C; size_t sz = (size_t) nn * S.hom_V[j]; S.hom_H[j] = (int *) R_alloc(sz, sizeof(int)); memset(S.hom_H[j], 0, sz * sizeof(int));
    for (int p = 0; p < S.n1; p++) S.hom_H[j][(S.hom_side[j] == 1 ? S.alpha[p] : S.gamma[p]) * S.hom_V[j] + S.nat[p]]++;
  }
  /* acov / edgecov (interleaved) */
  S.nacov = LEN(acs); S.nedge = LEN(eds);
  if (S.nacov) { S.acov_col = INTEGER(VECTOR_ELT(acs, 0)); S.W = REAL(VECTOR_ELT(acs, 1)); S.nacov = LEN(VECTOR_ELT(acs, 0)); }
  if (S.nedge) { S.edge_col = INTEGER(VECTOR_ELT(eds, 0)); S.M = REAL(VECTOR_ELT(eds, 1)); S.nedge = LEN(VECTOR_ELT(eds, 0)); }
  /* member lists */
  S.agstart = (int *) R_alloc(S.A + 1, sizeof(int)); S.clstart = (int *) R_alloc(S.C + 1, sizeof(int));
  S.agList = (int *) R_alloc(S.n1, sizeof(int)); S.clList = (int *) R_alloc(S.n1, sizeof(int)); S.agPos = (int *) R_alloc(S.n1, sizeof(int)); S.clPos = (int *) R_alloc(S.n1, sizeof(int));
  { int *fa = (int *) R_alloc(S.A, sizeof(int)), *fc = (int *) R_alloc(S.C, sizeof(int)); memset(fa, 0, S.A * sizeof(int)); memset(fc, 0, S.C * sizeof(int));
    for (int p = 0; p < S.n1; p++) { fa[S.alpha[p]]++; fc[S.gamma[p]]++; }
    S.agstart[0] = 0; for (int a = 0; a < S.A; a++) S.agstart[a + 1] = S.agstart[a] + fa[a];
    S.clstart[0] = 0; for (int c = 0; c < S.C; c++) S.clstart[c + 1] = S.clstart[c] + fc[c];
    memset(fa, 0, S.A * sizeof(int)); memset(fc, 0, S.C * sizeof(int));
    for (int p = 0; p < S.n1; p++) { int i = S.agstart[S.alpha[p]] + fa[S.alpha[p]]++; S.agList[i] = p; S.agPos[p] = i; int k = S.clstart[S.gamma[p]] + fc[S.gamma[p]]++; S.clList[k] = p; S.clPos[p] = k; } }
  /* initial statistics from the tables */
  double *g = (double *) R_alloc(S.nstat, sizeof(double)), *dg = (double *) R_alloc(S.nstat, sizeof(double)); memset(g, 0, S.nstat * sizeof(double));
  if (S.stack_col >= 0) for (long i = 0; i < (long) S.A * S.C; i++) g[S.stack_col] += 0.5 * S.n[i] * (S.n[i] - 1.0);
  for (int j = 0; j < S.ncol; j++) for (int a = 0; a < S.A; a++) for (int k = 1; k <= S.col_nk[j]; k++) { double x = S.col_T[j][a * (S.col_nk[j] + 1) + k]; g[S.col_col[j]] += S.col_w[j][a] * 0.5 * x * (x - 1); }
  for (int j = 0; j < S.nhom; j++) { int nn = S.hom_side[j] == 1 ? S.A : S.C; for (int u = 0; u < nn; u++) for (int v = 0; v < S.hom_V[j]; v++) { double x = S.hom_H[j][u * S.hom_V[j] + v]; if (S.hom_cen[j][u] != 0) g[S.hom_col[j]] += 0.5 * x * (x - 1); } }
  for (int p = 0; p < S.n1; p++) {
    for (int j = 0; j < S.nacov; j++) g[S.acov_col[j]] += S.W[((size_t) S.agnode[S.alpha[p]] * S.n2 + S.clnode[S.gamma[p]]) * S.nacov + j];
    for (int j = 0; j < S.nedge; j++) g[S.edge_col[j]] += S.M[((size_t) p * S.n2 + S.agnode[S.alpha[p]]) * S.nedge + j] + S.M[((size_t) p * S.n2 + S.clnode[S.gamma[p]]) * S.nedge + j];
  }
  int nsim = INTEGER(nsim_)[0]; double burnin = REAL(burnin_)[0], interval = REAL(interval_)[0]; const double *theta = REAL(theta_);
  rng_t rng; rseed(&rng, (uint64_t) REAL(seed_)[0]);
  SEXP out = PROTECT(allocMatrix(REALSXP, nsim, S.nstat)), fa_ = PROTECT(allocVector(INTSXP, S.n1)), fg_ = PROTECT(allocVector(INTSXP, S.n1));
  double *o = REAL(out); double acc = 0, tot = 0;
  for (double it = 0; it < burnin; it++) { if (step(&S, &rng, theta, dg, 0)) { acc++; for (int j = 0; j < S.nstat; j++) g[j] += dg[j]; } tot++; if (((long) it & 0xFFFFF) == 0) R_CheckUserInterrupt(); }
  for (int d = 0; d < nsim; d++) {
    if (d > 0) for (double it = 0; it < interval; it++) { if (step(&S, &rng, theta, dg, 0)) { acc++; for (int j = 0; j < S.nstat; j++) g[j] += dg[j]; } tot++; if (((long) it & 0xFFFFF) == 0) R_CheckUserInterrupt(); }
    for (int j = 0; j < S.nstat; j++) o[d + (size_t) j * nsim] = g[j];
  }
  for (int p = 0; p < S.n1; p++) { INTEGER(fa_)[p] = S.alpha[p]; INTEGER(fg_)[p] = S.gamma[p]; }
  SEXP res = PROTECT(allocVector(VECSXP, 4)); SET_VECTOR_ELT(res, 0, out); SET_VECTOR_ELT(res, 1, fa_); SET_VECTOR_ELT(res, 2, fg_); SET_VECTOR_ELT(res, 3, ScalarReal(acc / (tot > 0 ? tot : 1)));
  UNPROTECT(4); return res;
}
