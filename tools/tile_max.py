"""The largest tiles of PMT files, exactly, without decoding every tile (THES-260 / the display's PSRAM budget).

Usage: python tools/tile_max.py <file.pmt> [<file.pmt> ...]

The display fills a tile's polygons through an edge table (one edge per polygon point, the closing edges included) and
decodes a tile into fixed scratch (pwca_panoptes components/panoptes_map_view/map_view.c: map_draw_init). This prints,
per file and level, the most points and the most polygon points in one tile, and the largest blob.

Exact, and fast on a 2 GB file: a point is at least two bytes of a tile's blob (two varints), so a tile cannot hold
more points than half its blob. The tiles are taken largest blob first and decoded until the next blob's bound is
below the most points found: no tile left can hold more.
"""

import sys

import pmt


def level_max(f, li):
    entries = sorted(f.entries(li), key=lambda e: e[2], reverse=True)
    lv = f.levels[li]
    best_pts = best_poly = decoded = 0
    for e in entries:
        if e[2] // 2 < best_pts:
            break  # no tile left can hold more points (polygon points are among them)
        feats = pmt.decode_tile(f.tile_blob(li, e), lv.coord_bits, lv.buffer, f.attr_count, f.kind)
        pts = poly = 0
        for geom, _cls, _attr, parts in feats:
            if geom == pmt.POINT:
                pts += 1
                continue
            if f.kind == pmt.KIND_STREETS:  # a named line: (name, points)
                pts += len(parts[1])
                continue
            n = sum(len(p) for p in parts)
            pts += n
            if geom == pmt.POLYGON:
                poly += n
        best_pts = max(best_pts, pts)
        best_poly = max(best_poly, poly)
        decoded += 1
    return len(entries), (entries[0][2] if entries else 0), best_pts, best_poly, decoded


def main(argv):
    if len(argv) < 2:
        print(__doc__.splitlines()[2])
        return 2
    worst_pts = worst_poly = 0
    for path in argv[1:]:
        f = pmt.open_pmt_mapped(path)
        print(f"{path} (kind {f.kind})")
        for li, lv in enumerate(f.levels):
            n, blob, pts, poly, dec = level_max(f, li)
            print(f"  level {li}: tz {lv.tile_zoom} z{lv.zoom_min}-{lv.zoom_max}: {n} tiles; largest blob {blob} B; "
                  f"most points {pts}; most polygon points {poly} ({dec} tiles decoded)")
            worst_pts = max(worst_pts, pts)
            worst_poly = max(worst_poly, poly)
    print(f"ALL: most points in a tile {worst_pts}; most polygon points in a tile {worst_poly}")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
