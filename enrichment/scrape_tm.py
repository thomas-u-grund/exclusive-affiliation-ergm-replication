"""
Collect player- and agency-level enrichment data from Transfermarkt for the
players in Milan Dittrich's BA snapshot (snapshot_2026-08-24).

Sources (one cached raw file per request; reruns skip what is cached):
  player    tmapi-alpha.transfermarkt.technology/player/{id}                 -> DOB, nationalities, contract, joined, agency, MV current/previous/peak
  mvhist    tmapi-alpha.transfermarkt.technology/player/{id}/market-value-history -> full market value history (value at snapshot date)
  perf      www.transfermarkt.com/ceapi/performance-game/{id}                -> per-match career record (club + national team)
  agency    www.transfermarkt.com/-/beraterfirma/berater/{id}                -> agency HQ, global client count, total MV
  agencypage www.transfermarkt.com/-/beraterfirma/berater/{agency}/page/{n}  -> client list pages 2.. (page 1 is in agency)

Usage:  python3 scrape_tm.py [player|mvhist|perf|agency|agencypage|all]
"""
import csv, gzip, json, os, random, sys, time, threading, urllib.request, urllib.error
from concurrent.futures import ThreadPoolExecutor

ROOT = os.path.dirname(os.path.abspath(__file__))
REPRO = os.path.join(ROOT, "..", "work", "Bachelorarbeit_R_Reproduktionsmaterial")
PLAYERS_CSV = os.path.join(REPRO, "data/clean/snapshot_2026-08-24/players_clean.csv")
AGENTS_CSV = os.path.join(REPRO, "data/metadata/snapshot_2026-08-24/cleaning/agent_reference.csv")
RAW = os.path.join(ROOT, "raw")
LOG = os.path.join(ROOT, "raw", "scrape_log.csv")

UA = "Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/128.0 Safari/537.36"
THREADS = 2
PAUSE = (0.6, 1.2)  # seconds per request per thread

ENDPOINTS = {
    "player": ("https://tmapi-alpha.transfermarkt.technology/player/{id}", "json"),
    "mvhist": ("https://tmapi-alpha.transfermarkt.technology/player/{id}/market-value-history", "json"),
    "perf":   ("https://www.transfermarkt.com/ceapi/performance-game/{id}", "json"),
    "agency": ("https://www.transfermarkt.com/-/beraterfirma/berater/{id}", "html"),
    "agencypage": ("https://www.transfermarkt.com/-/beraterfirma/berater/{id}", "html"),  # id = "<agency>/page/<n>"
}

log_lock = threading.Lock()


def ids_for(kind):
    if kind == "agencypage":
        # client lists show 25 players per page; needs output/agency_metadata.csv (build_tables.py) for the counts
        with open(os.path.join(ROOT, "output", "agency_metadata.csv"), newline="", encoding="utf-8") as f:
            return [f"{r['source_agent_id']}/page/{n}" for r in csv.DictReader(f)
                    for n in range(2, -(-int(r["n_players_global"] or 0) // 25) + 1)]
    if kind == "agency":
        with open(AGENTS_CSV, newline="", encoding="utf-8") as f:
            return sorted({r["source_agent_id"] for r in csv.DictReader(f) if r["source_agent_id"]}, key=int)
    with open(PLAYERS_CSV, newline="", encoding="utf-8") as f:
        return sorted({r["player_id"] for r in csv.DictReader(f)}, key=int)


def cache_path(kind, id_):
    ext = ENDPOINTS[kind][1]
    return os.path.join(RAW, kind, f"{id_.replace('/page/', '_p')}.{ext}.gz")


def fetch(url, tries=5):
    for attempt in range(1, tries + 1):
        req = urllib.request.Request(url, headers={"User-Agent": UA, "Accept-Encoding": "gzip",
                                                   "Accept": "application/json,text/html;q=0.9"})
        try:
            with urllib.request.urlopen(req, timeout=60) as r:
                body = r.read()
                if r.headers.get("Content-Encoding") == "gzip":
                    body = gzip.decompress(body)
                return r.status, body, attempt
        except urllib.error.HTTPError as e:
            if e.code == 404:
                return 404, e.read(), attempt
            wait = 30 * attempt if e.code in (403, 429) else 5 * attempt
        except Exception:
            wait = 5 * attempt
        time.sleep(wait)
    return None, b"", tries


def worker(kind, id_):
    path = cache_path(kind, id_)
    if os.path.exists(path):
        return "cached"
    url = ENDPOINTS[kind][0].format(id=id_)
    status, body, attempts = fetch(url)
    ok = status == 200 and (ENDPOINTS[kind][1] != "json" or body.lstrip().startswith(b"{"))
    if ok:
        tmp = path + ".tmp"
        with gzip.open(tmp, "wb") as f:
            f.write(body)
        os.replace(tmp, path)
    with log_lock:
        new = not os.path.exists(LOG)
        with open(LOG, "a", newline="") as f:
            w = csv.writer(f)
            if new:
                w.writerow(["kind", "id", "url", "status", "ok", "attempts", "bytes", "retrieved_at_utc"])
            w.writerow([kind, id_, url, status, ok, attempts, len(body),
                        time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())])
    time.sleep(random.uniform(*PAUSE))
    return "ok" if ok else f"fail:{status}"


def run(kind):
    os.makedirs(os.path.join(RAW, kind), exist_ok=True)
    ids = ids_for(kind)
    todo = [i for i in ids if not os.path.exists(cache_path(kind, i))]
    print(f"[{kind}] {len(ids)} ids, {len(todo)} to fetch", flush=True)
    done = fails = 0
    with ThreadPoolExecutor(THREADS) as ex:
        for res in ex.map(lambda i: worker(kind, i), todo):
            done += 1
            fails += res.startswith("fail")
            if done % 200 == 0:
                print(f"[{kind}] {done}/{len(todo)} fetched, {fails} failures", flush=True)
    print(f"[{kind}] finished: {done} attempted, {fails} failures", flush=True)


if __name__ == "__main__":
    kinds = sys.argv[1:] or ["all"]
    if kinds == ["all"]:
        kinds = ["player", "agency", "mvhist", "perf"]
    for k in kinds:
        run(k)
