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

--reg-body writes the catalogue's unsigned body (Hermes' crates/hermes-core/src/map/regions.rs reads it); Web/API
signs it with the map release key over SHA-256("PWC-MAP-REGIONS-V1" || body) and appends the 64-byte signature.
The three passes over the extract (relations, their ways, their nodes) are cached in --cache.
"""
import argparse
import collections
import hashlib
import json
import os
import pickle
import struct
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import osmpbf  # noqa: E402

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


def extract(pbf, cache, processes):
    """(relations, {way id: node refs}, (node ids, lat, lon)), each pass cached in `cache`."""
    os.makedirs(cache, exist_ok=True)
    rel_path = os.path.join(cache, "admin1_rel.pkl")
    if not os.path.exists(rel_path):
        pickle.dump(osmpbf.map_blocks(pbf, osmpbf._relations_in, (wanted,), processes), open(rel_path, "wb"))
    rels = pickle.load(open(rel_path, "rb"))
    ids_path = os.path.join(cache, "admin1_wayids.npy")
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


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--pbf", required=True)
    ap.add_argument("--manifest", required=True)
    ap.add_argument("--cache", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--reg-body")
    ap.add_argument("--processes", type=int, default=5)
    a = ap.parse_args(argv)
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
