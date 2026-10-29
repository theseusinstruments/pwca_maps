"""Map regions and their groups for a release's region catalogue R<n>.reg (Jira THES-263; Hermes - 5's
HERMES5_MAP_REGIONS_DESIGN.md and HERMES5_MAP_REGION_GROUPS.md, approved at zoom 6).

A region is a state, province or territory (OSM admin_level=4 with an ISO3166-2 code: US-, CA-; Mexico as one
region, the owner 2026-10-27), or Greenland (admin_level=2). Its cells are the zoom-6 Web Mercator cells its boundary
covers or reaches.

THE GROUP RULE (the owner, 2026-10-27, through PC-33: "Yes, thats what I want."): most regions stand alone; a region
too small to own a map cell joins the neighbour owning its largest cell, and that group selects and deselects
together. "Too small" = it covers less than 10% of every cell it touches (SMALL_COVER). Each cell's owner is the
region covering most of it. On release 2026102201 this gives Mid-Atlantic (VA + DC, DE, NJ) and New England (ME + RI,
CT); everything else stands alone. Shared border cells stay on a card while any selected region needs them (Hermes).

    python tools/region_groups.py --pbf sources/north-america-latest.osm.pbf --manifest R2026102201.man \\
        --cache <dir> --out regions/2026102201.json [--reg-body R2026102201.reg.body] [--processes 5]

THE STORE REGIONS (Branding - 5's BRANDING5_WEB_SHOP_ITEM.md section 2, "Store regions"; the owner, 2026-10-29: "Okay,
this is perfect."): the shop sells named groups of the release's codes, regions/<release>.store.json, written by hand
and reviewed with the release: {release, regions: [{code, name, country, codes}]}. --check-store checks one against the
release's table: every code in exactly one store region, no code the table lacks, every card group whole in one, the
store codes in the database's key rule and the release the same. No extract needed:

    python tools/region_groups.py --check-store regions/2026102201.store.json --table regions/2026102201.json

THE PLACES INDEX (Branding - 5's BRANDING5_WEB_SHOP_ITEM.md section 2, "The regions card (the ?)" and "What the search
knows"): regions/<release>.places.json, which the shop's "Where do you ride?" search reads (pwca_api
app/core/map_regions.py search_places; the API adds the states, provinces and store regions itself):

    python tools/region_groups.py --places-out regions/2026102201.places.json --table regions/2026102201.json \\
        --cache <dir> --release-dir C:/pwca_maps_releases/R2026102201 --water-names out/na-full/water/names.tsv

{release, places: [{name, kind, states, rank, lat, lon}]}, one place a line; states are the table's codes; rank 0 is
best; lat and lon (degrees, 4 decimals) the point a map centres on (the group ride map's search reads the same file).
  - Cities, towns and villages: the release's places layer (layer 5 pieces, the level cut at zoom 12; hamlets left
    out). A point: its one code. Rank: city 0, town 1, village 2.
  - Marinas and ramps: the named marinas, fuel docks and boat ramps of the release's symbols layer (layer 6 pieces,
    the level for zoom 14, where names are written). A point. Kind "ramp" (rank 1: riders search for their launch)
    for a slipway, and for a marina or fuel dock whose name says launch, ramp or access site (RAMP_WORDS: OSM maps
    many public water access sites as marinas); else "marina" (rank 2). A name of GENERIC_WORDS only ("Boat
    Launch", "Public Access", "Marina") is left out.
  - Rivers: names.tsv's label points of named rivers (one every 1.5 km) and named river areas; the points of one name
    within RIVER_LINK of each other are one stretch, and stretches of one name sharing a code are one river (a lake on
    it breaks its points), listed in every code its points fall in. Rank by length. Its point: the label point
    nearest the mean of its points.
  - Lakes and reservoirs of MIN_LAKE_M2 or more: names.tsv's label point is its point and gives a code; a lake of
    OUTLINE_ALWAYS_M2 or more, or one with a state line within NEAR_FACTOR equal-area radii of its label point, is
    traced on the release's water layer (layer 1 pieces, its polygons rastered RES x RES a tile and followed across
    tiles from the label point, at the coarsest level that draws it 4 pixels in radius), and listed in every code its
    outline and water lie in (Lake Havasu: US-AZ and US-CA). Rank by area.
  - Bays and straits: names.tsv's named bay areas. A bay is not drawn (the sea under it is), so the water layer has
    no outline of it: its codes are those of its label point and of 8 points round it at its equal-area radius.
    Rank by area. A bay mapped as a point only (Tampa Bay, Georgian Bay) is the places layer's BAY_POINT_CLASS less
    the names names.tsv has (that class also holds mid-sized lakes): its point's code, rank 2.
  - Not yet: OHV areas (kind "ohv"). No layer of the release or the build carries them; they need a pass over the
    extract. The self-check looks for Hungry Valley once the index holds any.
  - Left out: hamlets, unnamed things, canals, streams, drains, harbour basins, street addresses and ZIP codes. One
    entry per name, kind and codes (the API shows no more), with the best rank and the point of its best (largest)
    instance.
The boundaries come from --cache (with --pbf only when the cache is not whole: that is the heavy part). The step ends
with a self-check (the locator on two squares; Lake Havasu City in US-AZ, Lake Havasu in US-AZ and US-CA, ...), which
--check-places runs on a written index:

    python tools/region_groups.py --check-places regions/2026102201.places.json

--reg-body writes the catalogue's unsigned body (Hermes' crates/hermes-core/src/map/regions.rs reads it); Web/API
signs it with the map release key over SHA-256("PWC-MAP-REGIONS-V1" || body) and appends the 64-byte signature.
The three passes over the extract (relations, their ways, their nodes) are cached in --cache.
"""
import argparse
import collections
import hashlib
import json
import math
import os
import pickle
import re
import struct
import sys
import time
import unicodedata

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import osmpbf  # noqa: E402
from build_places import STAND_INS  # noqa: E402

CELL_ZOOM = 6
N = 1 << CELL_ZOOM
SMALL_COVER = 0.10
SAMPLES = 24
COUNTRIES = [(1, "UNITED STATES", "US"), (2, "CANADA", "CA"), (3, "MEXICO", "MX"), (4, "GREENLAND", "GL")]
REG_MAGIC = b"PREG"


# ---- the boundaries (three passes over the extract) -------------------------------------------------------------------

def wanted(tags):
    iso = tags.get("ISO3166-2", "")
    if tags.get("boundary") == "administrative" and tags.get("admin_level") == "4" and iso[:3] in ("US-", "CA-", "MX-"):
        return True
    return tags.get("boundary") == "administrative" and tags.get("admin_level") == "2" and tags.get("ISO3166-1") == "GL"


def no_ways(tags):
    return None


CACHE_FILES = ("admin1_rel.pkl", "admin1_ways.pkl", "admin1_nodes.npz")


def extract(pbf, cache, processes):
    """(relations, {way id: node refs}, (node ids, lat, lon)), each pass cached in `cache`. Without `pbf` the cache
    must be whole (reading the extract is the heavy part: about 70 minutes on 5 processes)."""
    if not pbf and not all(os.path.exists(os.path.join(cache, n)) for n in CACHE_FILES):
        raise SystemExit(f"{cache} does not hold {', '.join(CACHE_FILES)}: --pbf is needed to read the extract")
    os.makedirs(cache, exist_ok=True)
    rel_path = os.path.join(cache, "admin1_rel.pkl")
    if not os.path.exists(rel_path):
        pickle.dump(osmpbf.map_blocks(pbf, osmpbf._relations_in, (wanted,), processes), open(rel_path, "wb"))
    rels = pickle.load(open(rel_path, "rb"))
    ids_path = os.path.join(cache, "admin1_wayids.npy")
    if pbf:
        np.save(ids_path, np.unique(np.array([m for _, _, ms in rels for t, m, r in ms if t == 1], np.int64)))
    ways_path = os.path.join(cache, "admin1_ways.pkl")
    if not os.path.exists(ways_path):
        pickle.dump(osmpbf.map_blocks(pbf, osmpbf._ways_in, (no_ways, ids_path), processes), open(ways_path, "wb"))
    ways = {wid: refs for wid, _, refs in pickle.load(open(ways_path, "rb"))}
    nodes_path = os.path.join(cache, "admin1_nodes.npz")
    if not os.path.exists(nodes_path):
        needed = np.unique(np.concatenate(list(ways.values())).astype(np.int64))
        need_path = os.path.join(cache, "admin1_nodeids.npy")
        np.save(need_path, needed)
        lat = np.full(needed.size, np.nan)
        lon = np.full(needed.size, np.nan)
        for pos, la, lo, _ in osmpbf.map_blocks(pbf, osmpbf._nodes_in, (need_path, None), processes):
            lat[pos] = la
            lon[pos] = lo
        np.savez(nodes_path, ids=needed, lat=lat, lon=lon)
    nd = np.load(nodes_path)
    return rels, ways, (nd["ids"], nd["lat"], nd["lon"])


def tile_x(lon):
    return (lon + 180.0) / 360.0 * N


def tile_y(lat):
    s = np.sin(np.radians(np.clip(lat, -85.0511, 85.0511)))
    return (0.5 - np.log((1 + s) / (1 - s)) / (4 * np.pi)) * N


def boundaries(rels, ways, nodes):
    """{code: (name, [polylines in tile units])}; the Mexican states as one region MX."""
    ids, lat, lon = nodes
    out = {}
    for _, tags, members in rels:
        code = tags.get("ISO3166-2") or tags.get("ISO3166-1")
        if not code or (code in out and not code.startswith("MX-")):
            continue
        name = tags.get("name:en") or tags.get("name") or code
        if code.startswith("MX-"):
            code, name = "MX", "Mexico"
        segs = []
        for t, m, _ in members:
            if t != 1 or m not in ways:
                continue
            pos = np.searchsorted(ids, ways[m])
            pos = pos[pos < ids.size]
            x, y = tile_x(lon[pos]), tile_y(lat[pos])
            good = ~np.isnan(x) & ~np.isnan(y)
            q = np.round(np.stack([x[good], y[good]], 1) * 2000) / 2000
            if len(q) < 2:
                continue
            keep = np.ones(len(q), bool)
            keep[1:] = np.any(q[1:] != q[:-1], 1)
            segs.append(q[keep])
        if segs:
            out.setdefault(code, (name, []))[1].extend(segs)
    return out


def coverage(segs, cells):
    """{cell: share of the cell inside the boundary}, even-odd over every edge (shared borders of merged states
    appear twice and cancel); a cell the boundary only reaches gets a tiny share (the region needs it drawn)."""
    edges = np.concatenate([np.stack([s[:-1], s[1:]], 1) for s in segs if len(s) > 1])
    x0, y0, x1, y1 = edges[:, 0, 0], edges[:, 0, 1], edges[:, 1, 0], edges[:, 1, 1]
    pts = np.concatenate(segs)
    touched = {(int(a), int(b)) for a, b in np.floor(pts).astype(int)}
    g = (np.arange(SAMPLES) + 0.5) / SAMPLES
    out = {}
    for cx in range(int(np.floor(pts[:, 0].min())), int(np.floor(pts[:, 0].max())) + 1):
        for cy in range(int(np.floor(pts[:, 1].min())), int(np.floor(pts[:, 1].max())) + 1):
            if (cx, cy) not in cells:
                continue
            sel = ~((np.minimum(y0, y1) > cy + 1) | (np.maximum(y0, y1) < cy))
            if not sel.any():
                continue
            ex0, ey0, ex1, ey1 = x0[sel], y0[sel], x1[sel], y1[sel]
            px = np.tile(cx + g, SAMPLES)
            py = np.repeat(cy + g, SAMPLES)
            inside = np.zeros(px.size, bool)
            for i in range(0, px.size, 64):
                X, Y = px[i:i + 64, None], py[i:i + 64, None]
                with np.errstate(divide="ignore", invalid="ignore"):
                    xi = ex0 + (Y - ey0) * (ex1 - ex0) / (ey1 - ey0)
                inside[i:i + 64] = (np.sum(((ey0 > Y) != (ey1 > Y)) & (X < xi), 1) % 2) == 1
            f = float(inside.mean())
            if f > 0 or (cx, cy) in touched:
                out[(cx, cy)] = max(f, 1e-4)
    return out


# ---- the rule --------------------------------------------------------------------------------------------------------

def groups(cov):
    """{code: group code}: each cell's owner covers most of it; a region covering < SMALL_COVER of every cell it
    touches joins the owner of its largest cell; everyone else leads its own group (no chains: a small region's host
    is never itself clubbed)."""
    by_cell = collections.defaultdict(dict)
    for code, d in cov.items():
        for cell, f in d.items():
            by_cell[cell][code] = f
    owner = {cell: max(d, key=lambda c: (d[c], c)) for cell, d in by_cell.items()}
    small = {code for code, d in cov.items() if d and max(d.values()) < SMALL_COVER}
    out = {}
    for code, d in cov.items():
        out[code] = code
        if code in small:
            best = max(d, key=lambda cell: (d[cell], cell))
            host = owner[best]
            if host != code and host not in small:
                out[code] = host
    return out


# ---- the table and the catalogue's body -------------------------------------------------------------------------------

def table(release, cov, group):
    """The computed table: countries, then regions (stable ids in code order) with their group and cells."""
    codes = sorted(cov, key=lambda c: (next(i for i, _, p in COUNTRIES if c.startswith(p)), c))
    ids = {c: i + 1 for i, c in enumerate(codes)}
    regions = []
    for c in codes:
        country = next(i for i, _, p in COUNTRIES if c.startswith(p))
        regions.append({"id": ids[c], "code": c, "name": c and NAMES.get(c, c), "country": country, "group": ids[group[c]],
                        "best_cover": round(max(cov[c].values()), 4), "cells": sorted([list(k) for k in cov[c]])})
    return {"release": release, "cell_zoom": CELL_ZOOM, "small_cover": SMALL_COVER,
            "countries": [{"id": i, "name": n} for i, n, _ in COUNTRIES], "regions": regions}


NAMES = {}


def reg_body(tbl, manifest_sha256):
    """R<n>.reg without its signature (Hermes' map::regions::Catalogue reads it)."""
    regions = tbl["regions"]
    countries = [c for c in tbl["countries"] if any(r["country"] == c["id"] for r in regions)]
    b = bytearray(REG_MAGIC + bytes([1, 0]) + struct.pack("<HI", 48, tbl["release"]) + manifest_sha256
                  + struct.pack("<HH", len(countries), len(regions)))
    order = []
    for c in countries:
        mine = [r for r in regions if r["country"] == c["id"]]
        b += struct.pack("<H", c["id"]) + c["name"].encode("ascii")[:32].ljust(32, b"\0") + struct.pack("<HH", len(order), len(mine))
        order += mine
    for r in order:
        b += struct.pack("<HH", r["id"], r["country"]) + r["code"].encode("ascii")[:8].ljust(8, b"\0")
        b += r["name"].encode("ascii", "replace")[:32].ljust(32, b"\0") + struct.pack("<HH", r["group"], len(r["cells"]))
        for x, y in r["cells"]:
            b += struct.pack("<HH", x, y)
    return bytes(b)


def manifest_cells(path):
    data = open(path, "rb").read()
    count = struct.unpack_from("<I", data, 20)[0]
    cells = set()
    for i in range(count):
        x, y = struct.unpack_from("<HH", data, 48 + 56 * i + 42)
        if x != 0xFFFF:
            cells.add((x, y))
    return struct.unpack_from("<I", data, 8)[0], cells, data


STORE_CODE = re.compile(r"^[a-z][a-z0-9_]{1,31}$")  # the shop database's key rule (migration 0193)


def store_problems(tbl, store):
    """What is wrong with a store regions file against the release's region table; [] when nothing is."""
    problems = []
    if store.get("release") != tbl.get("release"):
        problems.append(f"release {store.get('release')} is not the table's {tbl.get('release')}")
    known = {r["code"]: r for r in tbl["regions"]}
    seen = {}
    names = set()
    for s in store.get("regions", []):
        code = s.get("code", "")
        if not STORE_CODE.match(code):
            problems.append(f"store code {code!r} is not lowercase letters, digits and underscores")
        if not s.get("name") or s["name"] in names:
            problems.append(f"{code}: no name, or a name used twice")
        names.add(s.get("name"))
        if not s.get("country"):
            problems.append(f"{code}: no country")
        if not s.get("codes"):
            problems.append(f"{code}: no codes")
        for c in s.get("codes", []):
            if c not in known:
                problems.append(f"{code}: {c} is not a region of the release")
            elif c in seen:
                problems.append(f"{c} is in {seen[c]} and {code}")
            seen[c] = code
    for c in known:
        if c not in seen:
            problems.append(f"{c} is in no store region")
    groups_of = collections.defaultdict(set)
    for c, r in known.items():
        groups_of[r["group"]].add(c)
    for members in groups_of.values():
        homes = {seen.get(c) for c in members}
        if len(members) > 1 and len(homes) > 1:
            problems.append(f"card group {' + '.join(sorted(members))} is split across {sorted(h or '-' for h in homes)}")
    return problems


def check_store(store_path, table_path):
    tbl = json.load(open(table_path, encoding="utf-8"))
    store = json.load(open(store_path, encoding="utf-8"))
    problems = store_problems(tbl, store)
    for p in problems:
        print("PROBLEM:", p)
    if not problems:
        print(f"{len(store['regions'])} store regions hold the release's {len(tbl['regions'])} codes, each once,"
              " every card group whole: OK")
    return not problems


# ---- the places index (the shop's "Where do you ride?" search) -------------------------------------------------------

LOCATE_BINS = 64  # the point locator's rows per zoom-6 cell
EARTH_CIRC = 40_075_016.686
# build_places.py's classes in the release's places layer: city (and the cities of 250,000 and 1,000,000 people),
# town, village. Hamlets (4) are left out of the search.
PLACE_KINDS = {1: "city", 6: "city", 7: "city", 2: "town", 3: "village"}
PLACE_RANK = {"city": 0, "town": 1, "village": 2}
# A bay or strait mapped as a point (no size): build_places.py's class WATER0 + BAY_NODE_TIER, which it shares with
# water areas of 25 ha to 5 km2; a name of that class that names.tsv does not have is such a bay.
BAY_POINT_CLASS = 20
BAY_POINT_RANK = 2
# build_extras.py's symbol classes: 1 fuel dock, 2 marina, 3 boat ramp (leisure=slipway or a seamark slipway), named
# ones only. A slipway is kind "ramp"; a marina or fuel dock is "marina", unless its name says it is a launch (OSM
# maps many public water access sites as leisure=marina).
SYMBOL_KINDS = {1: "marina", 2: "marina", 3: "ramp"}
SYMBOL_RAMP = 3
RAMP_WORDS = re.compile(r"\b(boat ramp|ramp|launch|launching|slipway|slip way|water access|access site|boat access|"
                        r"public access|fishing access|put-in|put in|take-out)\b", re.I)
KIND_RANK = {"marina": 2, "ramp": 1}
# A marina or ramp whose name has only these words is left out (it names no place).
GENERIC_WORDS = set("""a and the of at to for in out on boat boats launch launching ramp ramps marina marinas public
access water site sites slipway slip slips way fishing fish dock docks landing fuel gas kayak canoe carry put take
area small craft hand trailer parking lot entrance new old main north south east west upper lower
private community city county state public""".split())
# names.tsv (build_pack.py, the water-full build): A = a named water area (class 2 lake, 3 reservoir, 4 river area,
# 5 canal, harbour or marina basin: left out), B = a named bay or strait (its area), L = a label point of a named
# waterway every 1.5 km (16 river; canals, streams and drains are left out).
LAKE_CLASSES = (2, 3)
RIVER_AREA_CLASS = 4
RIVER_LINE_CLASS = 16
MIN_LAKE_M2 = 100_000.0  # 10 ha (25 acres): smaller named water is a pond (the lead, 2026-10-29: from 4 ha)
# rank: (from size, rank), largest first; below the last, rank 3. Lakes and bays by area (m2), rivers by length (m).
LAKE_RANKS = [(50e6, 0), (5e6, 1), (0.5e6, 2)]
BAY_RANKS = [(500e6, 0), (20e6, 1), (1e6, 2)]
RIVER_RANKS = [(300e3, 0), (50e3, 1), (5e3, 2)]
OUTLINE_ALWAYS_M2 = 5e6  # a lake this large is always traced (a long reservoir reaches far from its label point)
NEAR_FACTOR = 4.0  # a smaller lake is traced when a state line is within this many of its equal-area radii
RIVER_LINK = 0.016  # zoom-6 units (10 km at the equator, 5 km at 60 N): label points of one name this near are one river
RES = 512  # raster side of a water tile, for tracing a lake
TRACE_MAX_TILES = 400
# A known place and the codes it must have (the self-check of --places-out).
# Known places: (name, kinds, the codes it must have, (lat, lon) it must be near or None, how near in degrees). The
# self-check of --places-out; an "ohv" entry is checked only once the index holds OHV areas.
PLACES_EXPECT = [
    ("Lake Havasu City", ("city", "town"), {"US-AZ"}, (34.48, -114.33), 0.05),
    ("Lake Havasu", ("lake",), {"US-AZ", "US-CA"}, (34.41, -114.29), 0.2),
    ("Lake Okeechobee", ("lake",), {"US-FL"}, (26.94, -80.80), 0.3),
    ("Lake of the Ozarks", ("lake",), {"US-MO"}, None, 0),
    ("Lake of the Woods", ("lake",), {"CA-MB", "CA-ON", "US-MN"}, None, 0),
    ("Lake Simcoe", ("lake",), {"CA-ON"}, None, 0),
    ("Lake Havasu Marina", ("marina",), {"US-AZ"}, (34.461, -114.351), 0.01),
    ("Lees Ferry Boat Launch", ("ramp",), {"US-AZ"}, (36.866, -111.586), 0.01),
    ("Hungry Valley", ("ohv",), {"US-CA"}, (34.75, -118.85), 0.1),  # an ohv name is matched by its start
]


class Locator:
    """Which region each of many points lies in: even-odd over the region's boundary edges with a ray to the east.
    Points and boundaries in zoom-6 cell units; -1 for a point in no region."""

    def __init__(self, bounds, codes):
        self.codes = list(codes)
        parts = []
        for i, c in enumerate(self.codes):
            for s in bounds[c][1]:
                a, b = s[:-1], s[1:]
                keep = a[:, 1] != b[:, 1]
                parts.append(np.column_stack([a[keep], b[keep], np.full(int(keep.sum()), float(i))]))
        e = np.concatenate(parts)
        b0 = np.floor(np.minimum(e[:, 1], e[:, 3]) * LOCATE_BINS).astype(np.int64)
        b1 = np.floor(np.maximum(e[:, 1], e[:, 3]) * LOCATE_BINS).astype(np.int64)
        n = b1 - b0 + 1
        idx = np.repeat(np.arange(len(e)), n)
        bins = np.repeat(b0, n) + (np.arange(int(n.sum())) - np.repeat(np.cumsum(n) - n, n))
        order = np.lexsort((e[idx, 4], bins))
        self.bins, self.e = bins[order], e[idx[order]]

    def locate(self, x, y):
        x, y = np.asarray(x, float).ravel(), np.asarray(y, float).ravel()
        out = np.full(x.size, -1, np.int64)
        if not x.size:
            return out
        pb = np.floor(y * LOCATE_BINS).astype(np.int64)
        order = np.argsort(pb, kind="stable")
        ub, starts = np.unique(pb[order], return_index=True)
        ends = np.append(starts[1:], order.size)
        lo, hi = np.searchsorted(self.bins, ub, "left"), np.searchsorted(self.bins, ub, "right")
        for s, t, a, b in zip(starts, ends, lo, hi):
            if a == b:
                continue
            e = self.e[a:b]
            codes = e[:, 4].astype(np.int64)
            first = np.concatenate([[0], np.flatnonzero(np.diff(codes)) + 1])
            ex0, ey0, ex1, ey1 = e[:, 0], e[:, 1], e[:, 2], e[:, 3]
            pts = order[s:t]
            step = max(1, 4_000_000 // len(e))
            for i in range(0, len(pts), step):
                p = pts[i:i + step]
                X, Y = x[p, None], y[p, None]
                xi = ex0 + (Y - ey0) * (ex1 - ex0) / (ey1 - ey0)
                cross = (((ey0 > Y) != (ey1 > Y)) & (X < xi)).view(np.uint8)
                odd = (np.add.reduceat(cross, first, axis=1, dtype=np.uint8) & 1).astype(bool)
                hit = odd.any(1)
                out[p[hit]] = codes[first[odd.argmax(1)[hit]]]
        return out


def locator_problems():
    """The locator on two squares side by side (no data needed)."""
    sq = lambda x0: (x0, [np.array([[x0, 1.0], [x0 + 1, 1.0], [x0 + 1, 2.0], [x0, 2.0], [x0, 1.0]])])
    loc = Locator({"A": sq(1.0), "B": sq(2.0)}, ["A", "B"])
    got = loc.locate([1.5, 2.5, 3.5, 1.999, 2.001], [1.5, 1.5, 1.5, 1.2, 1.9]).tolist()
    return [] if got == [0, 1, -1, 0, 1] else [f"the locator on two squares gave {got}"]


def to_cell_units(lat_deg, lon_deg):
    return tile_x(np.asarray(lon_deg, float)), tile_y(np.asarray(lat_deg, float))


def cell_lat_lon(x, y):
    """(lat, lon) in degrees, 4 decimals (about 11 m: enough to centre a map), of a point in zoom-6 units."""
    lon = x / N * 360.0 - 180.0
    lat = math.degrees(math.atan(math.sinh(math.pi * (1 - 2 * y / N))))
    return round(lat, 4), round(lon, 4)


def generic_symbol_name(name):
    """True for a marina's or ramp's name made only of generic words ("Boat Launch", "Public Access", "Marina"): it
    names no place a rider could search for."""
    words = re.findall(r"[a-z]+", name.lower())
    return bool(words) and all(w in GENERIC_WORDS for w in words)


def symbol_kind(cls, name):
    """marina or ramp: a slipway is a ramp; so is a marina or fuel dock whose name says it is a launch."""
    if cls == SYMBOL_RAMP or RAMP_WORDS.search(name):
        return "ramp"
    return "marina"


def rank_of(size, table):
    return next((r for m, r in table if size >= m), 3)


def release_pieces(release_dir):
    """(release, {layer id: {(cell x, cell y): piece path}}) from the release's unsigned manifest."""
    import manifest
    mans = [n for n in os.listdir(release_dir) if n.endswith(".man.unsigned") or n.endswith(".man")]
    if not mans:
        raise SystemExit(f"{release_dir}: no R<release>.man.unsigned")
    data = open(os.path.join(release_dir, sorted(mans)[-1]), "rb").read()
    body = manifest.parse_body(data[:manifest.HEADER + manifest.RECORD * struct.unpack_from("<I", data, 20)[0]])
    out = collections.defaultdict(dict)
    for r in body["records"]:
        out[r["layer"]][(r["cell_x"], r["cell_y"])] = os.path.join(release_dir, "pieces", r["sha256"].hex() + ".pmt")
    return body["release"], out


def piece_points(paths, tile_zoom, zoom_min, classes):
    """[(class, x, y in zoom-6 units, name)] of the named points of one level of a point layer's pieces."""
    import pmt
    out = []
    for (cx, _cy), path in sorted(paths.items()):
        if cx == 0xFFFF:
            continue
        f = pmt.open_pmt_mapped(path)
        level = next((i for i, lv in enumerate(f.levels) if lv.tile_zoom == tile_zoom and lv.zoom_min == zoom_min), None)
        if level is None:
            raise SystemExit(f"{path}: no level cut at zoom {tile_zoom} for zoom {zoom_min}")
        bits, div = f.levels[level].coord_bits, float(1 << (tile_zoom - CELL_ZOOM))
        for entry in f.entries(level):
            tx, ty = pmt.unmorton(entry[0])
            for _g, cls, _a, (px, py, name) in f.tile(level, entry):
                if cls in classes and name:
                    out.append((cls, (tx + px / (1 << bits)) / div, (ty + py / (1 << bits)) / div,
                                name.decode("utf-8")))
        f.data.close()
    return out


def rasterize(feats, bits):
    """A RES x RES raster of a water tile: each pixel the number (1-based) of the lake or reservoir polygon over its
    centre, 0 for none (even-odd within a polygon)."""
    scale = RES / float(1 << bits)
    cols = []
    for k, (_g, _c, _a, parts) in enumerate(feats, 1):
        for part in parts:
            a = np.asarray(part, float) * scale
            b = np.roll(a, -1, 0)
            keep = a[:, 1] != b[:, 1]
            cols.append(np.column_stack([a[keep], b[keep], np.full(int(keep.sum()), float(k))]))
    lab = np.zeros((RES, RES), np.int32)
    if not cols:
        return lab
    e = np.concatenate(cols)
    ymin, ymax = np.minimum(e[:, 1], e[:, 3]), np.maximum(e[:, 1], e[:, 3])
    r0 = np.clip(np.ceil(ymin - 0.5), 0, RES).astype(np.int64)
    r1 = np.clip(np.ceil(ymax - 0.5), 0, RES).astype(np.int64)
    n = r1 - r0
    idx = np.repeat(np.arange(len(e)), n)
    rows = np.repeat(r0, n) + (np.arange(int(n.sum())) - np.repeat(np.cumsum(n) - n, n))
    if not rows.size:
        return lab
    ee = e[idx]
    xc = ee[:, 0] + (rows + 0.5 - ee[:, 1]) * (ee[:, 2] - ee[:, 0]) / (ee[:, 3] - ee[:, 1])
    order = np.lexsort((xc, rows, ee[:, 4]))
    rows, xc, k = rows[order], xc[order], ee[order, 4].astype(np.int64)
    c0 = np.clip(np.ceil(xc[0::2] - 0.5), 0, RES).astype(np.int64)
    c1 = np.clip(np.ceil(xc[1::2] - 0.5), 0, RES).astype(np.int64)
    diff = np.zeros((RES, RES + 1), np.int64)
    np.add.at(diff, (rows[0::2], c0), k[0::2])
    np.add.at(diff, (rows[0::2], c1), -k[0::2])
    lab = np.cumsum(diff[:, :RES], 1)
    lab[(lab < 0) | (lab > len(feats))] = 0  # two polygons over one pixel: neither
    return lab.astype(np.int32)


class WaterTiles:
    """The release's water layer (layer 1) as rasters of its lakes and reservoirs, tile by tile, for tracing a lake's
    outline across tiles."""

    def __init__(self, paths):
        import pmt
        self.pmt, self.paths, self.files, self.cache = pmt, paths, {}, collections.OrderedDict()
        f = self.open(next(c for c in sorted(paths) if c[0] != 0xFFFF))
        from build_pack import LEVELS_FULL
        self.levels = [(lv.tile_zoom, lv.coord_bits) for lv in f.levels]
        if [z for z, _ in self.levels] != [lv[0] for lv in LEVELS_FULL]:
            raise SystemExit("the release's water levels are not build_pack.py's LEVELS_FULL")
        self.min_area = [lv[6] for lv in LEVELS_FULL]

    def open(self, cell):
        if cell not in self.files:
            self.files[cell] = self.pmt.open_pmt_mapped(self.paths[cell]) if cell in self.paths else None
        return self.files[cell]

    def tile(self, level, tx, ty):
        """(raster, the lake and reservoir polygons) of one tile; an absent tile has none."""
        key = (level, tx, ty)
        if key in self.cache:
            self.cache.move_to_end(key)
            return self.cache[key]
        tz, bits = self.levels[level]
        f = self.open((tx >> (tz - CELL_ZOOM), ty >> (tz - CELL_ZOOM)))
        entry = f.find(level, tx, ty) if f is not None and 0 <= tx < (1 << tz) and 0 <= ty < (1 << tz) else None
        feats = [ft for ft in (f.tile(level, entry) if entry else []) if ft[0] == 1 and ft[1] in LAKE_CLASSES]
        out = (rasterize(feats, bits), feats)
        self.cache[key] = out
        if len(self.cache) > 300:
            self.cache.popitem(last=False)
        return out

    def level_for(self, area, lat):
        """The coarsest level that keeps a lake of this area and draws it at least 4 pixels in radius."""
        r = math.sqrt(area / math.pi)
        for i, (tz, _bits) in enumerate(self.levels):
            px = EARTH_CIRC * math.cos(math.radians(lat)) / ((1 << tz) * RES)
            if area >= self.min_area[i] and r >= 4 * px:
                return i
        return len(self.levels) - 1

    def trace(self, level, x6, y6):
        """Points on and in the outline of the lake under (x6, y6) (zoom-6 units), followed across tile edges; []
        when no lake is there at this level."""
        tz, bits = self.levels[level]
        wx, wy = x6 * (1 << (tz - CELL_ZOOM)), y6 * (1 << (tz - CELL_ZOOM))
        tx, ty = int(wx), int(wy)
        lab, _ = self.tile(level, tx, ty)
        px, py = min(RES - 1, int((wx - tx) * RES)), min(RES - 1, int((wy - ty) * RES))
        win = lab[max(0, py - 3):py + 4, max(0, px - 3):px + 4]
        k = int(lab[py, px]) or (int(np.bincount(win[win > 0]).argmax()) if (win > 0).any() else 0)
        if not k:
            return []
        seen, todo = {(tx, ty, k)}, [(tx, ty, k)]
        sides = [(1, 0, (slice(None), -1), (slice(None), 0)), (-1, 0, (slice(None), 0), (slice(None), -1)),
                 (0, 1, (-1, slice(None)), (0, slice(None))), (0, -1, (0, slice(None)), (-1, slice(None)))]
        xs, ys = [], []
        div, side = float(1 << (tz - CELL_ZOOM)), float(1 << bits)
        while todo:
            tx, ty, k = todo.pop()
            lab, feats = self.tile(level, tx, ty)
            mine = lab == k
            for part in feats[k - 1][3]:
                v = np.asarray(part, float)
                v = v[(v[:, 0] >= 0) & (v[:, 0] <= side) & (v[:, 1] >= 0) & (v[:, 1] <= side)]
                xs.append((tx + v[:, 0] / side) / div)
                ys.append((ty + v[:, 1] / side) / div)
            gy, gx = np.nonzero(mine[::8, ::8])
            xs.append((tx + (gx * 8 + 0.5) / RES) / div)
            ys.append((ty + (gy * 8 + 0.5) / RES) / div)
            if len(seen) >= TRACE_MAX_TILES:
                continue
            for dx, dy, here, there in sides:
                edge = mine[here]
                if not edge.any():
                    continue
                nlab, _ = self.tile(level, tx + dx, ty + dy)
                for v in np.unique(nlab[there][edge]):
                    node = (tx + dx, ty + dy, int(v))
                    if v and node not in seen:
                        seen.add(node)
                        todo.append(node)
        return [np.concatenate(xs), np.concatenate(ys)] if xs else []


def read_water_names(path):
    """(the rows of names.tsv we use: (type, class, size, lat, lon, name), the names of all its water areas and bays
    as the places layer writes them)."""
    rows, areas = [], set()
    with open(path, encoding="utf-8") as f:
        for line in f:
            typ, cls, size, la, lo, name = line.rstrip("\n").split("\t", 5)
            cls = int(cls)
            if typ != "L":
                areas.add(unicodedata.normalize("NFC", " ".join(name.split())).translate(STAND_INS))
            if (typ == "A" and cls in LAKE_CLASSES + (RIVER_AREA_CLASS,)) or typ == "B" or \
                    (typ == "L" and cls == RIVER_LINE_CLASS):
                if name.strip():
                    rows.append((typ, cls, float(size), int(la) / 1e7, int(lo) / 1e7, name))
    return rows, areas


def rivers(points, loc_codes):
    """[(name, codes, total metres, x, y)]: a river's label points joined where they lie within RIVER_LINK of each
    other (one river per name and stretch, so two Mill Rivers in states apart stay two), and stretches sharing a code
    joined. (x, y): its middle label point, the one nearest the mean of its points (zoom-6 units)."""
    by_name = collections.defaultdict(list)
    for (name, x, y, metres), code in zip(points, loc_codes):
        by_name[name].append((x, y, metres, code))
    out = []
    for name, pts in by_name.items():
        cell_of = {}
        for i, (x, y, _m, _c) in enumerate(pts):
            cell_of.setdefault((int(x / RIVER_LINK), int(y / RIVER_LINK)), []).append(i)
        parent = {c: c for c in cell_of}

        def find(c):
            while parent[c] != c:
                parent[c] = parent[parent[c]]
                c = parent[c]
            return c
        for (cx, cy) in cell_of:
            for dx in (-1, 0, 1):
                for dy in (-1, 0, 1):
                    n = (cx + dx, cy + dy)
                    if n in parent:
                        a, b = find((cx, cy)), find(n)
                        if a != b:
                            parent[a] = b
        groups = collections.defaultdict(lambda: [set(), 0.0, []])
        for c, members in cell_of.items():
            g = groups[find(c)]
            for i in members:
                if pts[i][3]:
                    g[0].add(pts[i][3])
                g[1] += pts[i][2]
                g[2].append(i)
        # Stretches of one name that share a code are one river (a river broken by a lake on it: the Snake in ID and
        # in ID, OR and WA); the search could not tell them apart anyway.
        merged = []
        for codes, metres, members in sorted(groups.values(), key=lambda g: -g[1]):
            if not codes:
                continue
            for m in merged:
                if m[0] & codes:
                    m[0] |= codes
                    m[1] += metres
                    m[2] += members
                    break
            else:
                merged.append([set(codes), metres, list(members)])
        joined = True
        while joined:
            joined = False
            for i in range(len(merged)):
                for j in range(i + 1, len(merged)):
                    if merged[i][0] & merged[j][0]:
                        merged[i][0] |= merged[j][0]
                        merged[i][1] += merged[j][1]
                        merged[i][2] += merged[j][2]
                        del merged[j]
                        joined = True
                        break
                if joined:
                    break
        for codes, metres, members in merged:
            xy = np.array([pts[i][:2] for i in members])
            mid = xy[int(np.argmin(((xy - xy.mean(0)) ** 2).sum(1)))]
            out.append((name, codes, metres, float(mid[0]), float(mid[1])))
    return out


def build_places(table_path, bounds, release_dir, water_names, log=print):
    """The places index {release, places: [{name, kind, states, rank}]} (see the module docstring)."""
    tbl = json.load(open(table_path, encoding="utf-8"))
    codes = [r["code"] for r in tbl["regions"]]
    missing = [c for c in codes if c not in bounds]
    if missing:
        raise SystemExit(f"the boundaries lack {missing}")
    loc = Locator(bounds, codes)
    release, pieces = release_pieces(release_dir)
    if release != tbl["release"]:
        raise SystemExit(f"{release_dir} is release {release}, the table {tbl['release']}")
    t0 = time.time()
    entries = {}  # (name, kind, codes) -> (rank, -size, x, y): the best instance's rank and point

    def put(name, kind, states, rank, x, y, size=0.0):
        key = (name, kind, tuple(sorted(states)))
        val = (rank, -size, float(x), float(y))
        if key not in entries or val[:2] < entries[key][:2]:
            entries[key] = val

    rows, area_names = read_water_names(water_names)
    log(f"names.tsv: {len(rows)} rows of lakes, river areas, bays and river label points ({time.time() - t0:.0f} s)")
    # Towns and cities: the places layer's level cut at zoom 12 for zooms 12-13 (cities, towns and villages). With
    # them its bays mapped as a point: the class it shares with mid-sized water areas, less the names names.tsv has.
    pts = piece_points(pieces[5], 12, 12, set(PLACE_KINDS) | {BAY_POINT_CLASS})
    at = loc.locate([p[1] for p in pts], [p[2] for p in pts])
    n_bays = 0
    for (cls, x, y, name), c in zip(pts, at):
        if c < 0:
            continue
        if cls == BAY_POINT_CLASS:
            if name not in area_names:
                n_bays += 1
                put(name, "bay", [codes[c]], BAY_POINT_RANK, x, y)
        else:
            put(name, PLACE_KINDS[cls], [codes[c]], PLACE_RANK[PLACE_KINDS[cls]], x, y)
    log(f"places layer: {sum(1 for p in pts if p[0] != BAY_POINT_CLASS)} cities, towns and villages, {n_bays} bays "
        f"mapped as a point ({time.time() - t0:.0f} s)")
    # Marinas, fuel docks and boat ramps with a name: the symbols layer's level for zoom 14 (names from 14).
    pts = piece_points(pieces[6], 12, 14, set(SYMBOL_KINDS))
    at = loc.locate([p[1] for p in pts], [p[2] for p in pts])
    n_generic = 0
    for (cls, x, y, name), c in zip(pts, at):
        if c < 0:
            continue
        if generic_symbol_name(name):
            n_generic += 1
            continue
        kind = symbol_kind(cls, name)
        put(name, kind, [codes[c]], KIND_RANK[kind], x, y)
    log(f"symbols layer: {len(pts)} named marinas, fuel docks and boat ramps, {int((at >= 0).sum())} in a region, "
        f"{n_generic} of them with a generic name only ({time.time() - t0:.0f} s)")
    # Water: names.tsv.
    # A river area's size is its area, not a length: it places the river, the label points measure it.
    river_pts = [(name, *to_cell_units(la, lo), size if typ == "L" else 0.0) for typ, cls, size, la, lo, name in rows
                 if typ == "L" or (typ == "A" and cls == RIVER_AREA_CLASS)]
    at = loc.locate([p[1] for p in river_pts], [p[2] for p in river_pts])
    river_list = rivers([(n, float(x), float(y), m) for n, x, y, m in river_pts],
                        [codes[c] if c >= 0 else None for c in at])
    for name, states, metres, x, y in river_list:
        put(name, "river", states, rank_of(metres, RIVER_RANKS), x, y, metres)
    log(f"rivers: {len(river_pts)} label points, {len(river_list)} rivers in a region ({time.time() - t0:.0f} s)")
    # Bays: a disc of the bay's own area round its label point (a bay is not drawn, so the water tiles have no
    # outline of it).
    ang = np.linspace(0, 2 * np.pi, 8, endpoint=False)
    bays = [r for r in rows if r[0] == "B"]
    bx, by = [], []
    for _t, _c, size, la, lo, _n in bays:
        rad = math.sqrt(size / math.pi) / (EARTH_CIRC * math.cos(math.radians(la)) / (1 << CELL_ZOOM))
        x, y = to_cell_units(la, lo)
        bx.append(np.concatenate([[x], x + rad * np.cos(ang)]))
        by.append(np.concatenate([[y], y + rad * np.sin(ang)]))
    at = loc.locate(np.concatenate(bx), np.concatenate(by)).reshape(-1, 9) if bays else np.zeros((0, 9), int)
    n_bays = 0
    for (_t, _c, size, la, lo, name), row in zip(bays, at):
        states = {codes[c] for c in row if c >= 0}
        if states:
            n_bays += 1
            put(name, "bay", states, rank_of(size, BAY_RANKS), *to_cell_units(la, lo), size)
    log(f"bays: {len(bays)}, {n_bays} in a region ({time.time() - t0:.0f} s)")
    # Lakes and reservoirs: the label point's region; a lake near a state line (or large) traced on the water tiles.
    lakes = [r for r in rows if r[0] == "A" and r[1] in LAKE_CLASSES and r[2] >= MIN_LAKE_M2]
    lx, ly, rad = [], [], []
    for _t, _c, size, la, lo, _n in lakes:
        x, y = to_cell_units(la, lo)
        lx.append(x)
        ly.append(y)
        rad.append((NEAR_FACTOR * math.sqrt(size / math.pi) + 300) / (EARTH_CIRC * math.cos(math.radians(la)) / (1 << CELL_ZOOM)))
    lx, ly, rad = np.array(lx), np.array(ly), np.array(rad)
    sx = np.concatenate([lx[:, None], lx[:, None] + rad[:, None] * np.cos(ang)], 1)
    sy = np.concatenate([ly[:, None], ly[:, None] + rad[:, None] * np.sin(ang)], 1)
    at = loc.locate(sx, sy).reshape(-1, 9)
    big = np.array([r[2] >= OUTLINE_ALWAYS_M2 for r in lakes], bool)
    near = (at != at[:, :1]).any(1)
    trace = np.flatnonzero(big | near)
    log(f"lakes: {len(lakes)} of {MIN_LAKE_M2 / 1e4:.0f} ha or more; {int(big.sum())} large, {int(near.sum())} near a "
        f"state line: {trace.size} to trace ({time.time() - t0:.0f} s)")
    water = WaterTiles(pieces[1])
    traced = {}
    levels = [water.level_for(lakes[i][2], lakes[i][3]) for i in trace]
    for n, (i, lev) in enumerate(sorted(zip(trace, levels), key=lambda t: (t[1], int(lx[t[0]] * 64), int(ly[t[0]] * 64)))):
        pts = water.trace(lev, lx[i], ly[i])
        if pts:
            got = loc.locate(pts[0], pts[1])
            traced[i] = {codes[c] for c in np.unique(got[got >= 0])}
        if n % 2000 == 1999:
            log(f"  traced {n + 1} of {trace.size} ({time.time() - t0:.0f} s)")
    n_lakes = 0
    for i, ((_t, _c, size, _la, _lo, name), row) in enumerate(zip(lakes, at)):
        states = set(traced.get(i, ()))
        if row[0] >= 0:
            states.add(codes[row[0]])
        if states:
            n_lakes += 1
            put(name, "lake", states, rank_of(size, LAKE_RANKS), lx[i], ly[i], size)
    log(f"lakes: {len(traced)} traced, {n_lakes} in a region ({time.time() - t0:.0f} s)")
    places = []
    for (n, k, s), (r, _size, x, y) in sorted(entries.items(), key=lambda e: (e[1][0], e[0][0], e[0][1], e[0][2])):
        lat, lon = cell_lat_lon(x, y)
        places.append({"name": n, "kind": k, "states": list(s), "rank": r, "lat": lat, "lon": lon})
    return {"release": release, "places": places}


def places_problems(index):
    """The known places that are not where they must be, and any place without a point in North America."""
    problems = locator_problems()
    by_name = collections.defaultdict(list)
    for p in index["places"]:
        by_name[p["name"]].append(p)
    no_point = [p["name"] for p in index["places"] if not (isinstance(p.get("lat"), float) and
                                                          isinstance(p.get("lon"), float) and 5 <= p["lat"] <= 85
                                                          and (-180 <= p["lon"] <= -10 or p["lon"] >= 170))]
    # (the western Aleutians, US-AK, lie past 180 degrees)
    if no_point:
        problems.append(f"{len(no_point)} places without a point in North America, e.g. {no_point[:5]}")
    has_ohv = any(p["kind"] == "ohv" for p in index["places"])
    for name, kinds, states, near, within in PLACES_EXPECT:
        if "ohv" in kinds and not has_ohv:
            continue
        if "ohv" in kinds:
            named = [p for p in index["places"] if p["kind"] == "ohv" and p["name"].startswith(name)]
        else:
            named = by_name.get(name, [])
        found = [p for p in named if p["kind"] in kinds and set(p["states"]) == states]
        if not found:
            problems.append(f"{name} ({'/'.join(kinds)}) is not in {' + '.join(sorted(states))}: "
                            f"{[(p['kind'], p['states']) for p in named] or 'not found'}")
        elif near and not any(abs(p["lat"] - near[0]) <= within and abs(p["lon"] - near[1]) <= within
                              for p in found):
            problems.append(f"{name} is at {[(p['lat'], p['lon']) for p in found]}, not within {within} of {near}")
    return problems


def write_places(path, index):
    """One place a line: diffable, and small."""
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    with open(path, "w", encoding="utf-8", newline="\n") as f:
        f.write('{"release": %d,\n "places": [\n' % index["release"])
        f.write(",\n".join(json.dumps(p, ensure_ascii=False, separators=(",", ":")) for p in index["places"]))
        f.write("\n]}\n")


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--check-store", help="a regions/<release>.store.json to check against --table, then stop")
    ap.add_argument("--table", help="the release's region table, for --check-store")
    ap.add_argument("--pbf")
    ap.add_argument("--manifest")
    ap.add_argument("--cache")
    ap.add_argument("--out")
    ap.add_argument("--reg-body")
    ap.add_argument("--processes", type=int, default=5)
    ap.add_argument("--places-out", help="write the places index regions/<release>.places.json, then stop")
    ap.add_argument("--release-dir", help="the release's folder (R<release>.man.unsigned and pieces/), for --places-out")
    ap.add_argument("--water-names", help="the water-full build's names.tsv, for --places-out")
    ap.add_argument("--check-places", help="a places index to self-check, then stop")
    a = ap.parse_args(argv)
    if a.check_store:
        if not a.table:
            ap.error("--check-store needs --table")
        raise SystemExit(0 if check_store(a.check_store, a.table) else 1)
    if a.check_places:
        problems = places_problems(json.load(open(a.check_places, encoding="utf-8")))
        for p in problems:
            print("PROBLEM:", p)
        print("the places index's self-check:", "FAILED" if problems else "OK")
        raise SystemExit(1 if problems else 0)
    if a.places_out:
        if not (a.table and a.cache and a.release_dir and a.water_names):
            ap.error("--places-out needs --table, --cache, --release-dir and --water-names")
        t0 = time.time()
        log = lambda m: print(f"[{time.time() - t0:6.0f} s] {m}", flush=True)
        bounds = boundaries(*extract(a.pbf, a.cache, a.processes))
        log(f"boundaries: {len(bounds)} regions")
        index = build_places(a.table, bounds, a.release_dir, a.water_names, log)
        write_places(a.places_out, index)
        kinds = collections.Counter(p["kind"] for p in index["places"])
        log(f"{len(index['places'])} places ({', '.join(f'{k} {n}' for k, n in sorted(kinds.items()))}), "
            f"{os.path.getsize(a.places_out)} B -> {a.places_out}")
        problems = places_problems(index)
        for p in problems:
            print("PROBLEM:", p)
        print("the places index's self-check:", "FAILED" if problems else "OK")
        raise SystemExit(1 if problems else 0)
    if not (a.pbf and a.manifest and a.cache and a.out):
        ap.error("--pbf, --manifest, --cache and --out are needed to build a table")
    release, cells, man = manifest_cells(a.manifest)
    rels, ways, nodes = extract(a.pbf, a.cache, a.processes)
    bounds = boundaries(rels, ways, nodes)
    NAMES.update({c: n for c, (n, _) in bounds.items()})
    cov = {c: coverage(segs, cells) for c, (_, segs) in bounds.items()}
    cov = {c: d for c, d in cov.items() if d}
    tbl = table(release, cov, groups(cov))
    os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)
    with open(a.out, "w", encoding="utf-8", newline="\n") as f:
        json.dump(tbl, f, indent=1)
        f.write("\n")
    if a.reg_body:
        open(a.reg_body, "wb").write(reg_body(tbl, hashlib.sha256(man).digest()))
    led = collections.defaultdict(list)
    for r in tbl["regions"]:
        led[r["group"]].append(r["code"])
    for members in led.values():
        if len(members) > 1:
            print("group:", " + ".join(members))
    print(len(tbl["regions"]), "regions,", sum(1 for m in led.values() if len(m) > 1), "groups ->", a.out)


if __name__ == "__main__":
    main()
