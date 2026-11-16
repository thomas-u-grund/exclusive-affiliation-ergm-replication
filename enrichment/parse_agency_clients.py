"""
Parse the global client lists of all 1,093 agencies (raw/agency/{id}.html.gz = page 1,
raw/agencypage/{id}_p{n}.html.gz = pages 2..) into output/agency_clients.csv.

One row per agency x client (as of the scrape date, 2026-09-23):
source_agent_id, player_id, nationality1, nationality2, age, club_id, club_name, mv_eur, page
"""
import csv, glob, gzip, html, os, re

ROOT = os.path.dirname(os.path.abspath(__file__))
RAW = os.path.join(ROOT, "raw")
OUT = os.path.join(ROOT, "output", "agency_clients.csv")

ROW = re.compile(r'<tr class="(?:odd|even)">(.*?)</tr>\s*(?=<tr class="(?:odd|even)">|</tbody>)', re.S)
PLAYER = re.compile(r'href="/[^"]*/profil/spieler/(\d+)"')
FLAG = re.compile(r'title="([^"]*)" alt="[^"]*" class="flaggenrahmen"')
CLUB = re.compile(r'<a title="([^"]*)" href="/[^"]*/startseite/verein/(\d+)"')
CELLS = re.compile(r'<td class="zentriert(?: hauptlink)?">([^<]*)</td>')


def mv(s):
    s = s.strip().replace("€", "")
    m = re.match(r"([\d.]+)(m|k)", s)
    if not m:
        return ""
    return str(int(round(float(m.group(1)) * (1e6 if m.group(2) == "m" else 1e3))))


def parse(path):
    text = gzip.open(path, "rt", encoding="utf-8", errors="replace").read()
    i = text.find('class="items"')
    if i < 0:
        return []
    body = text[i:]
    j = body.find('<div class="pager"')
    body = body[:j] if j > 0 else body
    rows = []
    for r in ROW.findall(body):
        p = PLAYER.search(r)
        if not p:
            continue
        flags = [html.unescape(f) for f in FLAG.findall(r)]
        club = CLUB.search(r)
        cells = CELLS.findall(r)  # age, contract, market value (club cell contains a link, so not matched)
        rows.append({
            "player_id": p.group(1),
            "nationality1": flags[0] if flags else "",
            "nationality2": flags[1] if len(flags) > 1 else "",
            "age": cells[0].strip() if cells else "",
            "club_id": club.group(2) if club else "",
            "club_name": html.unescape(club.group(1)) if club else "",
            "mv_eur": mv(cells[-1]) if cells else "",
        })
    return rows


def main():
    files = [(os.path.basename(f).split(".")[0], 1, f) for f in glob.glob(os.path.join(RAW, "agency", "*.html.gz"))]
    for f in glob.glob(os.path.join(RAW, "agencypage", "*.html.gz")):
        a, n = os.path.basename(f).split(".")[0].split("_p")
        files.append((a, int(n), f))
    out, seen = [], set()
    for a, n, f in sorted(files, key=lambda x: (int(x[0]), x[1])):
        for r in parse(f):
            if (a, r["player_id"]) in seen:
                continue
            seen.add((a, r["player_id"]))
            out.append({"source_agent_id": a, **r, "page": n})
    with open(OUT, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(out[0].keys()))
        w.writeheader()
        w.writerows(out)
    print(f"{len(out)} agency-client rows from {len(files)} pages, {len({r['source_agent_id'] for r in out})} agencies")


if __name__ == "__main__":
    main()
