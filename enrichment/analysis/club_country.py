"""Club -> country map from the cached career data (raw/perf/*.json.gz), standard library only.

Domestic league matches (competitionTypeId 1-6 = tier, not national-team games) link a competition to both clubs. Seed: the clubs of
The snapshot leagues (player_attributes.csv: snapshot_club_id, market_id such as germany_2). Propagation: a
competition's country is the majority country of its clubs with a known country; every club playing league games in
a competition with a known country gets that country (majority over its competitions). Repeated until stable.
Output: enrichment/output/club_country.csv (club_id, country, n_league_games, source).
Run from the project root: python3 -I enrichment/analysis/club_country.py [workers]
"""
import csv, glob, gzip, json, os, sys
from collections import Counter, defaultdict
from multiprocessing import Pool

ROOT = os.getcwd()
RAW = os.path.join(ROOT, "enrichment", "raw", "perf")
OUT = os.path.join(ROOT, "enrichment", "output")


def parse(path):
    pairs = Counter()
    try:
        d = json.loads(gzip.decompress(open(path, "rb").read()))
    except Exception:
        return pairs
    for g in (d.get("data") or {}).get("performance") or []:
        gi = g.get("gameInformation") or {}
        if gi.get("competitionTypeId") not in (1, 2, 3, 4, 5, 6) or gi.get("isNationalGame"):
            continue
        comp = gi.get("competitionId")
        ci = g.get("clubsInformation") or {}
        for side in ("club", "opponent"):
            cid = (ci.get(side) or {}).get("clubId")
            if comp and cid:
                pairs[(comp, str(cid))] += 1
    return pairs


def main():
    workers = int(sys.argv[1]) if len(sys.argv) > 1 else 8
    files = sorted(glob.glob(os.path.join(RAW, "*.json.gz")))
    total = Counter()
    with Pool(workers) as pool:
        for c in pool.imap_unordered(parse, files, chunksize=20):
            total.update(c)
    seed = {}
    with open(os.path.join(OUT, "player_attributes.csv"), newline="") as f:
        for r in csv.DictReader(f):
            if r["snapshot_club_id"] and r["market_id"]:
                seed[r["snapshot_club_id"]] = r["market_id"].rsplit("_", 1)[0]
    club_comps = defaultdict(Counter); comp_clubs = defaultdict(set)
    for (comp, club), n in total.items():
        club_comps[club][comp] += n; comp_clubs[comp].add(club)
    country = dict(seed); source = {c: "seed" for c in seed}
    for _ in range(20):
        comp_cty = {}
        for comp, clubs in comp_clubs.items():
            votes = Counter(country[c] for c in clubs if c in country)
            if votes:
                comp_cty[comp] = votes.most_common(1)[0][0]
        changed = 0
        for club, comps in club_comps.items():
            if club in seed:
                continue
            votes = Counter()
            for comp, n in comps.items():
                if comp in comp_cty:
                    votes[comp_cty[comp]] += n
            if votes:
                new = votes.most_common(1)[0][0]
                if country.get(club) != new:
                    country[club] = new; source[club] = "propagated"; changed += 1
        if not changed:
            break
    with open(os.path.join(OUT, "club_country.csv"), "w", newline="") as f:
        w = csv.writer(f); w.writerow(["club_id", "country", "n_league_games", "source"])
        for club in sorted(set(club_comps) | set(seed), key=lambda x: int(x) if x.isdigit() else 0):
            w.writerow([club, country.get(club, ""), sum(club_comps[club].values()), source.get(club, "")])
    known = sum(1 for c in club_comps if c in country)
    print(f"{len(files)} files, {len(comp_clubs)} competitions, {len(club_comps)} clubs, {known} with country;",
          Counter(country.values()).most_common(12))
    for comp, cty in sorted(comp_cty.items()):
        if cty in ("england", "germany", "spain", "italy", "france", "portugal", "poland"):
            print(cty, comp, len(comp_clubs[comp]))


if __name__ == "__main__":
    main()
