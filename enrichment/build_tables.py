"""
Build analysis tables from the raw Transfermarkt cache written by scrape_tm.py.

Outputs (enrichment/output/):
  player_attributes.csv    one row per player_id (4,800): demographics, contract, market value (at snapshot + peak), status changes since snapshot
  player_performance.csv   one row per player_id: 2025/26 season, 2026/27 up to the snapshot date, 2023-2026 window, and career totals up to the snapshot date
  agency_metadata.csv      one row per Transfermarkt agency id: HQ location, global client count, first-tier clients, total portfolio value
  national_team_clubs.csv  lookup used to separate senior from youth national-team matches
  build_report.txt         coverage and consistency checks
"""
import csv, gzip, html, json, os, re, statistics, time, urllib.request
from collections import Counter, defaultdict
from datetime import date

ROOT = os.path.dirname(os.path.abspath(__file__))
RAW = os.path.join(ROOT, "raw")
OUT = os.path.join(ROOT, "output")
REPRO = os.path.join(ROOT, "..", "work", "Bachelorarbeit_R_Reproduktionsmaterial")
PLAYERS_CSV = os.path.join(REPRO, "data/clean/snapshot_2026-08-24/players_clean.csv")
AGENTS_CSV = os.path.join(REPRO, "data/metadata/snapshot_2026-08-24/cleaning/agent_reference.csv")
SNAPSHOT = date(2026, 8, 24)
UA = "Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) AppleWebKit/537.36 Chrome/128.0 Safari/537.36"

# Transfermarkt competitionTypeId (inferred from competition ids observed in the data)
LEAGUE_TYPES = {1, 2, 3, 4, 5}          # domestic league tiers (1 = first tier of any country)
DOMESTIC_CUP_TYPES = {8, 9, 13, 14}     # national cup, super cups, league cup
INTL_CLUB_TYPES = {10}                  # UCL/UEL/UECL incl. qualifying, Club World Cup
PLAYED = {"played"}


def load_json(kind, id_):
    p = os.path.join(RAW, kind, f"{id_}.json.gz")
    if not os.path.exists(p):
        return None
    with gzip.open(p, "rt", encoding="utf-8") as f:
        d = json.load(f)
    return d.get("data") if d.get("success") else None


def load_html(kind, id_):
    p = os.path.join(RAW, kind, f"{id_}.html.gz")
    if not os.path.exists(p):
        return None
    with gzip.open(p, "rt", encoding="utf-8", errors="ignore") as f:
        return f.read()


def d_iso(s):
    try:
        return date.fromisoformat(s[:10]) if s else None
    except ValueError:
        return None


def age_on(dob, ref):
    if not dob:
        return None
    return ref.year - dob.year - ((ref.month, ref.day) < (dob.month, dob.day))


# ---------------------------------------------------------------- snapshot
with open(PLAYERS_CSV, newline="", encoding="utf-8") as f:
    snap = list(csv.DictReader(f))
snap_by_id = {r["player_id"]: r for r in snap}
player_ids = sorted(snap_by_id, key=int)

# ---------------------------------------------------------------- nationality id -> name
# Transfermarkt's API returns country ids only. Learn names from snapshot players with a
# single citizenship (their text label is unambiguous), then apply to all ids.
nat_votes = defaultdict(Counter)
players_json = {pid: load_json("player", pid) for pid in player_ids}
for pid, pj in players_json.items():
    if not pj:
        continue
    nats = pj["nationalityDetails"]["nationalities"]
    if not nats.get("secondNationalityId"):
        nat_votes[nats["nationalityId"]][snap_by_id[pid]["citizenship"].strip()] += 1
for pid, pj in players_json.items():  # second pass: dual nationals, strip the known first name
    if not pj:
        continue
    nats = pj["nationalityDetails"]["nationalities"]
    n1, n2 = nats["nationalityId"], nats.get("secondNationalityId")
    txt = snap_by_id[pid]["citizenship"].strip()
    if n2 and n1 in nat_votes and n2 not in nat_votes:
        first = nat_votes[n1].most_common(1)[0][0]
        if txt.startswith(first + " "):
            nat_votes[n2][txt[len(first) + 1:]] += 1
    elif n2 and n2 in nat_votes and n1 not in nat_votes:
        second = nat_votes[n2].most_common(1)[0][0]
        if txt.endswith(" " + second):
            nat_votes[n1][txt[:-len(second) - 1]] += 1
country_name = {k: v.most_common(1)[0][0] for k, v in nat_votes.items()}

# ---------------------------------------------------------------- player attributes
market_country = {"england_1": "England", "france_1": "France", "germany_1": "Germany", "germany_2": "Germany",
                  "germany_3": "Germany", "italy_1": "Italy", "poland_1": "Poland", "portugal_1": "Portugal",
                  "spain_1": "Spain"}

attr_rows = []
for pid in player_ids:
    s = snap_by_id[pid]
    snap_agent = "" if s["source_agent_id"] in ("", "NA") else s["source_agent_id"]
    pj = players_json[pid]
    mvh = load_json("mvhist", pid)
    row = {"player_id": pid, "market_id": s["market_id"], "snapshot_club_id": s["club_id"],
           "snapshot_agent_id": snap_agent, "snapshot_agent_status": s["agent_status"], "api_player_found": pj is not None}
    if pj:
        a = pj["attributes"]
        nats = pj["nationalityDetails"]["nationalities"]
        dob = d_iso(pj["lifeDates"]["dateOfBirth"]) if not pj["lifeDates"]["isDateOfBirthUnknown"] else None
        cur = [c for c in pj.get("clubAssignments", []) if c.get("type") == "current"]
        lr = a.get("lastContractRenewal") or {}
        lr_date = (f"{lr['year']:04d}-{lr['month']:02d}-{lr['day']:02d}"
                   if lr.get("year") and lr.get("month") and lr.get("day") else "")
        n1, n2 = nats.get("nationalityId"), nats.get("secondNationalityId") or None
        row.update({
            "date_of_birth": dob.isoformat() if dob else "",
            "age_at_snapshot": age_on(dob, SNAPSHOT),
            "age_at_snapshot_years": round((SNAPSHOT - dob).days / 365.25, 2) if dob else "",
            "nationality1_id": n1, "nationality1": country_name.get(n1, ""),
            "nationality2_id": n2 or "", "nationality2": country_name.get(n2, "") if n2 else "",
            "dual_national": bool(n2),
            "country_of_birth_id": pj["birthPlaceDetails"].get("countryOfBirthId") or "",
            "country_of_birth": country_name.get(pj["birthPlaceDetails"].get("countryOfBirthId"), ""),
            "height_m": a.get("height") or "",
            "foot": (a.get("preferredFoot") or {}).get("name", ""),
            "position_group": a.get("positionGroup", ""),
            "position_id": (a.get("position") or {}).get("id", ""),
            "position": (a.get("position") or {}).get("name", ""),
            "position_short": (a.get("position") or {}).get("shortName", ""),
            "contract_until": a.get("contractUntil") or "",
            "contract_option_id": a.get("contractOptionId") or "",
            "last_contract_renewal": lr_date,
            "current_club_id_at_scrape": cur[0]["clubId"] if cur else "",
            "joined_current_club": cur[0].get("start", "") if cur else "",
            "debut_current_club": cur[0].get("debut", "") if cur else "",
            "agency_id_at_scrape": a.get("consultantAgencyId") or "",
            "agency_verified_at_scrape": ((a.get("consultantAgency") or {}).get("verificationStatus") == "verified"),
            "mv_at_scrape_eur": ((pj.get("marketValueDetails") or {}).get("current") or {}).get("value", ""),
            "mv_at_scrape_determined": ((pj.get("marketValueDetails") or {}).get("current") or {}).get("determined", ""),
        })
        row["club_changed_since_snapshot"] = bool(cur) and str(cur[0]["clubId"]) != s["club_id"]
        row["agency_changed_since_snapshot"] = bool(snap_agent) and str(a.get("consultantAgencyId") or "") != snap_agent
        # contract/joined refer to the scrape date; blank them for players who moved after the snapshot
        if row["club_changed_since_snapshot"]:
            row["contract_until_valid_for_snapshot_club"] = False
        else:
            row["contract_until_valid_for_snapshot_club"] = True
    # market value history, evaluated at the snapshot date
    hist = sorted(((d_iso(h["marketValue"]["determined"]), h["marketValue"]["value"])
                   for h in (mvh or {}).get("history", []) if h.get("marketValue")), key=lambda x: x[0] or date.min)
    pre = [(d, v) for d, v in hist if d and d <= SNAPSHOT]
    row.update({
        "mv_snapshot_original_eur": s["player_valuation"] if s["player_valuation"] != "NA" else "",
        "mv_at_snapshot_eur": pre[-1][1] if pre else "",
        "mv_at_snapshot_determined": pre[-1][0].isoformat() if pre else "",
        "mv_peak_to_snapshot_eur": max(v for _, v in pre) if pre else "",
        "mv_peak_to_snapshot_date": max(pre, key=lambda x: x[1])[0].isoformat() if pre else "",
        "mv_12m_before_snapshot_eur": next((v for d, v in reversed(pre) if d <= date(2025, 8, 24)), ""),
        "mv_first_recorded_date": pre[0][0].isoformat() if pre else "",
        "mv_n_valuations_to_snapshot": len(pre),
    })
    attr_rows.append(row)

# ---------------------------------------------------------------- national-team club lookup (senior vs youth)
nt_clubs = set()
perf_cache = {}
for pid in player_ids:
    pf = load_json("perf", pid)
    perf_cache[pid] = pf["performance"] if pf else None
    for g in perf_cache[pid] or []:
        if g["gameInformation"]["isNationalGame"]:
            nt_clubs.add(g["clubsInformation"]["club"]["clubId"])

nt_path = os.path.join(OUT, "national_team_clubs.csv")
nt_info = {}
if os.path.exists(nt_path):
    with open(nt_path, newline="", encoding="utf-8") as f:
        nt_info = {r["club_id"]: r for r in csv.DictReader(f) if r.get("name")}  # retry unresolved ids
for cid in sorted(nt_clubs - set(nt_info)):
    try:
        req = urllib.request.Request(f"https://tmapi-alpha.transfermarkt.technology/club/{cid}",
                                     headers={"User-Agent": UA, "Accept": "application/json"})
        with urllib.request.urlopen(req, timeout=60) as r:
            d = json.load(r)["data"]
        name, main = d["name"], str(d["baseDetails"].get("mainClubId") or "")
    except Exception:
        name, main = "", ""
    # youth/olympic sides point to the senior national team via mainClubId
    youth = (bool(main) and main != str(cid)) or bool(re.search(r"\bU-?\s?\d{2}\b|Olymp|Youth|Juniors", name, re.I))
    nt_info[cid] = {"club_id": cid, "name": name, "main_club_id": main, "youth": youth}
    time.sleep(0.5)
os.makedirs(OUT, exist_ok=True)
with open(nt_path, "w", newline="", encoding="utf-8") as f:
    w = csv.DictWriter(f, ["club_id", "name", "main_club_id", "youth"])
    w.writeheader()
    w.writerows(sorted(nt_info.values(), key=lambda r: int(r["club_id"])))
senior_nt = {cid for cid, r in nt_info.items() if str(r["youth"]) in ("False", "false", "0", "") and r["name"]}

# ---------------------------------------------------------------- performance aggregation
WINDOWS = {
    "s2526": lambda g, d: g["seasonId"] == 2025,
    "s2627_to_snapshot": lambda g, d: g["seasonId"] == 2026 and d is not None and d <= SNAPSHOT,
    "s2326": lambda g, d: 2023 <= g["seasonId"] <= 2025,
    "career_to_snapshot": lambda g, d: d is not None and d <= SNAPSHOT,
}
SCOPES = {
    "club": lambda g, gi, cid: not gi["isNationalGame"],
    "league": lambda g, gi, cid: not gi["isNationalGame"] and gi["competitionTypeId"] in LEAGUE_TYPES,
    "toptier": lambda g, gi, cid: not gi["isNationalGame"] and gi["competitionTypeId"] == 1,
    "intlclub": lambda g, gi, cid: not gi["isNationalGame"] and gi["competitionTypeId"] in INTL_CLUB_TYPES,
    "nt_senior": lambda g, gi, cid: gi["isNationalGame"] and cid in senior_nt,
    "nt_youth": lambda g, gi, cid: gi["isNationalGame"] and cid not in senior_nt,
}
KEEP = {  # (window, scope) combinations written out
    ("s2526", "club"), ("s2526", "league"), ("s2526", "toptier"), ("s2526", "intlclub"), ("s2526", "nt_senior"),
    ("s2627_to_snapshot", "club"),
    ("s2326", "club"), ("s2326", "toptier"), ("s2326", "intlclub"),
    ("career_to_snapshot", "club"), ("career_to_snapshot", "toptier"), ("career_to_snapshot", "intlclub"),
    ("career_to_snapshot", "nt_senior"), ("career_to_snapshot", "nt_youth"),
}


def agg(games):
    played = [x for x in games if x["statistics"]["generalStatistics"]["participationState"] in PLAYED]
    st = lambda x, a, b: (x["statistics"].get(a) or {}).get(b) or 0
    minutes = sum(st(x, "playingTimeStatistics", "playedMinutes") for x in played)
    possible = sum((x["gameInformation"].get("gameDuration") or 90) for x in games)
    return {
        "games_listed": len(games),
        "apps": len(played),
        "starts": sum(bool(st(x, "playingTimeStatistics", "isStarting")) for x in played),
        "minutes": minutes,
        "goals": sum(st(x, "goalStatistics", "goalsScoredTotalOfficial") for x in played),
        "assists": sum(st(x, "goalStatistics", "assistsOfficial") for x in played),
        "yellow": sum(st(x, "cardStatistics", "yellowCardNet") for x in played),
        "minutes_share": round(minutes / possible, 4) if possible else "",
        "injured_games": sum(x["statistics"]["generalStatistics"]["participationState"] == "injured" for x in games),
    }


perf_rows = []
for pid in player_ids:
    games = perf_cache[pid]
    row = {"player_id": pid, "api_perf_found": games is not None}
    if games is not None:
        tagged = []
        for g in games:
            gi = g["gameInformation"]
            d = d_iso(((gi.get("date") or {}).get("dateTimeUTC")) or "")
            tagged.append((g, gi, d, g["clubsInformation"]["club"]["clubId"]))
        for (wname, scope) in sorted(KEEP):
            sel = [g for g, gi, d, cid in tagged if WINDOWS[wname](gi, d) and SCOPES[scope](g, gi, cid)]
            for k, v in agg(sel).items():
                if scope.startswith("nt") and k in ("games_listed", "minutes_share", "injured_games"):
                    continue  # NT game lists only contain call-ups, share is not meaningful
                row[f"{wname}_{scope}_{k}"] = v
        # main league and club in 2025/26 (by minutes, then games listed)
        lg = Counter(); cl = Counter()
        for g, gi, d, cid in tagged:
            if gi["seasonId"] == 2025 and not gi["isNationalGame"] and gi["competitionTypeId"] in LEAGUE_TYPES:
                w = 1 + ((g["statistics"].get("playingTimeStatistics") or {}).get("playedMinutes") or 0)
                lg[(gi["competitionId"], gi["competitionTypeId"])] += w
                cl[cid] += w
        if lg:
            (cid_, tier_), _ = lg.most_common(1)[0]
            row["s2526_main_league_id"], row["s2526_main_league_tier"] = cid_, tier_
            row["s2526_main_club_id"] = cl.most_common(1)[0][0]
        first = [d for g, gi, d, cid in tagged if d and not gi["isNationalGame"]
                 and g["statistics"]["generalStatistics"]["participationState"] == "played"]
        row["first_senior_club_app_date"] = min(first).isoformat() if first else ""
        row["n_clubs_career_to_snapshot"] = len({cid for g, gi, d, cid in tagged
                                                 if d and d <= SNAPSHOT and not gi["isNationalGame"]})
    perf_rows.append(row)

# ---------------------------------------------------------------- agencies


def txt(s):
    s = re.sub(r"<script.*?</script>|<style.*?</style>", " ", s, flags=re.S)
    s = re.sub(r"<[^>]+>", " ", s)
    return re.sub(r"\s+", " ", html.unescape(s)).strip()


def money(s):
    m = re.match(r"€\s*([\d.,]+)\s*(bn|m|k|Th\.)?", s or "")
    if not m:
        return ""
    v = float(m.group(1).replace(",", ""))
    return int(round(v * {"bn": 1e9, "m": 1e6, "k": 1e3, "Th.": 1e3, None: 1}[m.group(2)]))


with open(AGENTS_CSV, newline="", encoding="utf-8") as f:
    agent_ref = {r["source_agent_id"]: r for r in csv.DictReader(f) if r["source_agent_id"]}
ag_rows = []
for aid in sorted(agent_ref, key=int):
    h = load_html("agency", aid)
    row = {"source_agent_id": aid, "agent_name": agent_ref[aid]["agent_name"], "page_found": h is not None}
    if h:
        t = txt(h)
        grab = lambda pat: (v if (m := re.search(pat, t)) and (v := m.group(1).strip().strip(",").strip()) not in ("-", "") else "")
        row.update({
            "verified_since": grab(r"\(verified since (\d{2}/\d{2}/\d{4})\)"),
            "premium": " Premium Address " in f" {t} " or bool(re.search(r"\) Premium ", t)),
            "hq_street": grab(r"Street: (.*?) ZIP, location:"),
            "hq_zip_city": grab(r"ZIP, location: (.*?) Country:"),
            "hq_country": grab(r"Country: (.*?) Contact"),
            "website": grab(r"Website: (\S+)"),
            "n_players_global": grab(r"Key Facts Players: (\d+)"),
            "n_players_first_tier_global": grab(r"Players 1st tier: (\d+)"),
            "total_mv_global_eur": money(grab(r"Total market value: (€\s*[\d.,]+\s*(?:bn|m|k|Th\.)?)")),
        })
    ag_rows.append(row)

# ---------------------------------------------------------------- write


def write(name, rows):
    cols = []
    for r in rows:
        for k in r:
            if k not in cols:
                cols.append(k)
    with open(os.path.join(OUT, name), "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, cols)
        w.writeheader()
        w.writerows(rows)


write("player_attributes.csv", attr_rows)
write("player_performance.csv", perf_rows)
write("agency_metadata.csv", ag_rows)

# ---------------------------------------------------------------- report
rep = []
n = len(attr_rows)
cov = lambda rows, k: sum(1 for r in rows if r.get(k) not in (None, "", "NA")) / len(rows)
rep.append(f"players: {n}")
for k in ["api_player_found", "date_of_birth", "nationality1", "nationality2", "country_of_birth", "height_m", "position",
          "contract_until", "joined_current_club", "last_contract_renewal", "mv_at_snapshot_eur",
          "mv_snapshot_original_eur", "mv_peak_to_snapshot_eur"]:
    val = sum(1 for r in attr_rows if r.get(k) not in (None, "", "NA", False)) / n
    rep.append(f"  {k:32s} {val:6.1%}")
rep.append(f"  nationality ids without a learned name: {sorted({r.get('nationality1_id') for r in attr_rows if r.get('nationality1_id') and not r.get('nationality1')} | {r.get('nationality2_id') for r in attr_rows if r.get('nationality2_id') and not r.get('nationality2')})}")
both = [r for r in attr_rows if r["mv_at_snapshot_eur"] != "" and r["mv_snapshot_original_eur"] != ""]
agree = sum(float(r["mv_at_snapshot_eur"]) == float(r["mv_snapshot_original_eur"]) for r in both)
rep.append(f"MV at snapshot (from history) equals MV in original snapshot: {agree}/{len(both)}")
rep.append(f"club changed since snapshot: {sum(bool(r.get('club_changed_since_snapshot')) for r in attr_rows)}")
rep.append(f"agency changed since snapshot (among players with snapshot agency): {sum(bool(r.get('agency_changed_since_snapshot')) for r in attr_rows)}")
rep.append(f"performance found: {cov(perf_rows, 's2526_club_apps'):.1%}; with >=1 club app in 2025/26: "
           f"{sum(1 for r in perf_rows if (r.get('s2526_club_apps') or 0) > 0) / n:.1%}")
rep.append(f"agencies: {len(ag_rows)}; page found {cov(ag_rows, 'hq_country') if ag_rows else 0:.1%} with HQ country; "
           f"{cov(ag_rows, 'n_players_global') if ag_rows else 0:.1%} with global client count")
rep.append(f"HQ countries (top 15): {Counter(r.get('hq_country') for r in ag_rows).most_common(15)}")
rep.append(f"senior national-team club ids: {len(senior_nt)}, youth: {len(nt_info) - len(senior_nt)}")
open(os.path.join(OUT, "build_report.txt"), "w").write("\n".join(rep) + "\n")
print("\n".join(rep))
