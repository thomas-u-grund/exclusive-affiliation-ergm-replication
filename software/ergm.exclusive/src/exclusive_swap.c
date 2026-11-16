/* Typed-swap Metropolis-Hastings proposal for exclusive multi-affiliation networks.
 *
 * Bipartite representation: mode 1 = actors (players), mode 2 = affiliations of K types (agencies, clubs, ...).
 * Every actor has exactly one tie to an affiliation of each type. The proposal draws an actor p uniformly, a second
 * actor q uniformly among the others and a type t uniformly, and exchanges the type-t partners of p and q (four dyad
 * toggles). It is symmetric, preserves every affiliation size and the one-of-each-type rule, and leaves the
 * state unchanged (no toggles) when p and q have the same partner of type t.
 *
 * MH inputs: types of the mode-2 nodes (1..K), in node order.
 */
#include "ergm_MHproposal.h"
#include "ergm_MHstorage.h"
#include "ergm_edgetree.h"
#include "ergm_changestat.h"
#include <R.h>

MH_I_FN(Mi_ExclusiveSwap) {
  MHp->ntoggles = 4;
  /* number of types, computed once (the proposal runs millions of times) */
  MH_ALLOC_STORAGE(1, int, Kp);
  *Kp = 0;
  for (unsigned int i = 0; i < MHp->ninputs; i++) if ((int) MHp->inputs[i] > *Kp) *Kp = (int) MHp->inputs[i];
}

/* partner of actor p of type t: loop over p's out-edges (p is a mode-1 node, its ties are stored as out-edges) */
static Vertex partner_of_type(Vertex p, int t, Network *nwp, Vertex bip, double *types) {
  Edge e; Vertex v;
  STEP_THROUGH_OUTEDGES(p, e, v) {
    if ((int) types[v - bip - 1] == t) return v;
  }
  return 0;
}

MH_P_FN(MH_ExclusiveSwap) {
  Vertex bip = nwp->bipartite;
  Vertex p, q, ap, aq;
  int t;
  MH_GET_STORAGE(int, Kp);
  int K = *Kp;
  p = 1 + (Vertex) R_unif_index((double) bip);
  t = 1 + (int) R_unif_index((double) K);
  if (K >= 2 && unif_rand() < 0.5) {
    /* grouped swap: q is drawn uniformly among the other actors that share p's partner of another type t2 (t2 uniform among
     * the types other than t). The reverse move has the same probability (same t, t2, group), so the kernel is symmetric. */
    int t2 = 1 + (int) R_unif_index((double) (K - 1));
    if (t2 >= t) t2++;
    Vertex v = partner_of_type(p, t2, nwp, bip, MHp->inputs);
    unsigned int m = (v == 0) ? 0 : nwp->indegree[v];
    q = 0;
    if (m >= 2) {
      unsigned int r = (unsigned int) R_unif_index((double) (m - 1)), k = 0;
      Edge e; Vertex w;
      STEP_THROUGH_INEDGES(v, e, w) {
        if (w == p) continue;
        if (k++ == r) { q = w; break; }
      }
    }
    if (q == 0) { /* p is alone in its group: stay */
      ap = partner_of_type(p, t, nwp, bip, MHp->inputs);
      MHp->ntoggles = 2; Mtail[0] = p; Mhead[0] = ap; Mtail[1] = p; Mhead[1] = ap; MHp->logratio = 0; return;
    }
  } else {
    q = 1 + (Vertex) R_unif_index((double) (bip - 1));
    if (q >= p) q++;
  }
  ap = partner_of_type(p, t, nwp, bip, MHp->inputs);
  aq = partner_of_type(q, t, nwp, bip, MHp->inputs);
  if (ap == 0 || aq == 0) { /* state violates the constraint */
    MHp->toggletail[0] = MH_FAILED; MHp->togglehead[0] = MH_CONSTRAINT; return;
  }
  if (ap == aq) { /* null move: same partner, nothing to exchange */
    /* a valid no-op: toggle the same tie twice (net change and net change statistics are exactly zero) */
    MHp->ntoggles = 2;
    Mtail[0] = p; Mhead[0] = ap; Mtail[1] = p; Mhead[1] = ap;
    MHp->logratio = 0; return;
  }
  MHp->ntoggles = 4;
  Mtail[0] = p; Mhead[0] = ap;   /* p leaves ap   */
  Mtail[1] = q; Mhead[1] = aq;   /* q leaves aq   */
  Mtail[2] = p; Mhead[2] = aq;   /* p joins aq    */
  Mtail[3] = q; Mhead[3] = ap;   /* q joins ap    */
  MHp->logratio = 0;
}

MH_F_FN(Mf_ExclusiveSwap) {
}
