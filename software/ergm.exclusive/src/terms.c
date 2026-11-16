/* Cross-layer change statistics for exclusive multi-affiliation networks (bipartite representation: mode 1 = actors,
 * mode 2 = affiliations). A toggled dyad has tail p (an actor) and head x (an affiliation); the partners of an actor are
 * its out-edges, the partners of an affiliation its in-edges.
 *
 * coloc:      sum over "centre" affiliations a of sum over unordered pairs {p, q} of a's partners of the number of
 *             pairs (y, z), y a partner of p and z a partner of q, with key[y] == key[z] > 0.
 *             (co-location by club country with agencies as centres; stacking among the centres chosen when every
 *             club is its own key.)
 * homcentre:  sum over centre affiliations a of the number of pairs {p, q} of a's partners with attr[p] == attr[q].
 * acov:       sum over actors p of sum over unordered pairs {x, y} of p's affiliation partners of W[x, y]
 *             (affiliation-affiliation covariate summed over the actors that link them).
 */
#include "ergm_changestat.h"

/* inputs: [0] = number of actors; [1 .. n2] = centre weight (0 = not a centre; 1 = plain count); [n2+1 .. 2 n2] = key;
 * [2 n2+1] = number of keys G; [2 n2+2 .. 2 n2+2+G] = offsets of the key groups into the member list (key k occupies
 * members[off[k-1] .. off[k]-1]); then the member list (vertex ids of the nodes of each key).
 *
 * For a pair of a centre's actors the count needs, for a given key, the number of (q, z) with q in the centre, z a partner of q
 * and key[z] = key. Two ways to get it: scan the centre's actors and their partners, or scan the nodes of that key and test
 * each of their actors for membership in the centre. The cheaper one is used (cost compared per call), which matters because
 * agencies with hundreds of clients are toggled most often. */
static double coloc_count_d(Vertex x, double k, Vertex p, double *mem, Edge from, Edge to, double *key, Vertex bip, Network *nwp) {
  /* number of pairs (q, z): q in in(x), q != p, z in out(q), key[z] == k, by the cheaper of two scans */
  double cnt = 0;
  Edge e, g;
  Vertex q, z, w;
  double costA = (double) nwp->indegree[x] * 2.0, costB = 0;
  if (costA > 16) {
    Edge m;
    for (m = from; m < to; m++) costB += (double) nwp->indegree[(Vertex) mem[m]];
  } else costB = costA + 1;
  if (costA <= costB) {
    STEP_THROUGH_INEDGES(x, e, q) {
      if (q == p) continue;
      STEP_THROUGH_OUTEDGES(q, g, z) if (key[z - bip - 1] == k) cnt += 1;
    }
  } else {
    Edge m;
    for (m = from; m < to; m++) {
      w = (Vertex) mem[m];
      STEP_THROUGH_INEDGES(w, e, q) if (q != p && GetEdge(q, x, nwp)) cnt += 1;
    }
  }
  return cnt;
}

C_CHANGESTAT_FN(c_coloc) {
  Vertex bip = (Vertex) INPUT_PARAM[0];
  Vertex n2 = nwp->nnodes - bip;
  double *centre = INPUT_PARAM + 1, *key = INPUT_PARAM + 1 + n2;
  Edge G = (Edge) INPUT_PARAM[1 + 2 * n2];
  double *off = INPUT_PARAM + 2 + 2 * n2;
  double *memd = off + G + 1;
  Vertex p = tail, x = head, y, a;
  Edge f;
  double change = 0;
  if (centre[x - bip - 1] != 0) {
    double wx = centre[x - bip - 1];
    STEP_THROUGH_OUTEDGES(p, f, y) {
      double ky = key[y - bip - 1];
      if (ky <= 0) continue;
      change += wx * coloc_count_d(x, ky, p, memd, (Edge) off[(int) ky - 1], (Edge) off[(int) ky], key, bip, nwp);
    }
  } else if (key[x - bip - 1] > 0) {
    double kx = key[x - bip - 1];
    STEP_THROUGH_OUTEDGES(p, f, a) {
      double wa = centre[a - bip - 1];
      if (wa == 0) continue;
      change += wa * coloc_count_d(a, kx, p, memd, (Edge) off[(int) kx - 1], (Edge) off[(int) kx], key, bip, nwp);
    }
  }
  CHANGE_STAT[0] += edgestate ? -change : change;
}

/* inputs: [0] = number of actors; [1 .. n2] = centre indicator; [n2+1 .. n2+n1] = actor attribute */
C_CHANGESTAT_FN(c_homcentre) {
  Vertex bip = (Vertex) INPUT_PARAM[0];
  Vertex n2 = nwp->nnodes - bip;
  double *centre = INPUT_PARAM + 1, *attr = INPUT_PARAM + 1 + n2;
  Vertex p = tail, x = head, q;
  Edge e;
  double change = 0;
  if (centre[x - bip - 1] != 0) {
    STEP_THROUGH_INEDGES(x, e, q) if (q != p && attr[q - 1] == attr[p - 1]) change += 1;
  }
  CHANGE_STAT[0] += edgestate ? -change : change;
}

/* inputs: [0] = number of actors; [1 .. n2*n2] = W, row x (0-based), column y */
C_CHANGESTAT_FN(c_acov) {
  Vertex bip = (Vertex) INPUT_PARAM[0];
  Vertex n2 = nwp->nnodes - bip;
  double *W = INPUT_PARAM + 1;
  Vertex p = tail, x = head, y;
  Edge f;
  double change = 0;
  STEP_THROUGH_OUTEDGES(p, f, y) if (y != x) change += W[(size_t)(x - bip - 1) * n2 + (y - bip - 1)];
  CHANGE_STAT[0] += edgestate ? -change : change;
}

/* stack: number of pairs of actors that share their partner of every affiliation type (for two types: the same agency and
 * the same club; equal to b1dsp(2) on the constrained space, but computed from the smaller of the partner sets).
 * inputs: [0] = number of actors; [1] = number of types K. Actors hold at most one partner per type in every state the
 * typed-swap proposal passes through, so the toggled actor's other partners are the K-1 partners of the other types. */
C_CHANGESTAT_FN(c_stack) {
  int K = (int) INPUT_PARAM[1], no = 0, i, j;
  Vertex p = tail, x = head, y, q, best = x, others[16];
  Edge f, e;
  double change = 0;
  STEP_THROUGH_OUTEDGES(p, f, y) if (y != x && no < 16) others[no++] = y;
  if (no == K - 1) {
    unsigned int bd = nwp->indegree[x];
    for (i = 0; i < no; i++) if (nwp->indegree[others[i]] < bd) { best = others[i]; bd = nwp->indegree[best]; }
    STEP_THROUGH_INEDGES(best, e, q) {
      int ok = 1;
      if (q == p) continue;
      if (best != x && !GetEdge(q, x, nwp)) continue;
      for (j = 0; j < no && ok; j++) if (others[j] != best && !GetEdge(q, others[j], nwp)) ok = 0;
      if (ok) change += 1;
    }
  }
  CHANGE_STAT[0] += edgestate ? -change : change;
}


/* dstack: number of (centre affiliation, other affiliation) cells holding at least two actors, i.e. stacked cells counted once
 * however large they are (a saturating counterpart of the pair count in coloc). inputs: [0] = number of actors; [1 .. n2] =
 * centre indicator (agencies). Every non-centre affiliation is a "club". */
C_CHANGESTAT_FN(c_dstack) {
  Vertex bip = (Vertex) INPUT_PARAM[0];
  double *centre = INPUT_PARAM + 1;
  Vertex p = tail, x = head, y, q, small, other;
  Edge f, e;
  double change = 0;
  int xc = centre[x - bip - 1] != 0;
  STEP_THROUGH_OUTEDGES(p, f, y) {
    if (y == x || (centre[y - bip - 1] != 0) == xc) continue;   /* the other-type partner of p */
    small = nwp->indegree[x] <= nwp->indegree[y] ? x : y; other = small == x ? y : x;
    unsigned int m = 0;
    STEP_THROUGH_INEDGES(small, e, q) if (q != p && GetEdge(q, other, nwp)) { if (++m > 1) break; }
    if (m == 1) change += 1;
  }
  CHANGE_STAT[0] += edgestate ? -change : change;
}
