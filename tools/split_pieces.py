"""Cuts a built pack into the pieces of a map release (THES-260; pwca_panoptes
docs/design/PANOPTES_MAP_UPDATES_CONTRACT.md sections 1 and 2) and writes the release's unsigned manifest.

Usage:
    python tools/split_pieces.py OUT_DIR --release N --pack-name NAME [--min-fw X.Y.Z] [--source-date YYYYMMDD]
                                 --layer water=a.pmt,b.pmt --layer roads=r.pmt [--layer ...]

A piece is one layer x one zoom-6 Web Mercator cell: a PMT file holding exactly the tiles, of every level, whose
zoom-6 ancestor is the cell (with the coverage grid's zoom-8 cells inside it). Levels whose tile zoom is under 6 do
not fit a cell: they go, uncut, into one more piece of the layer, cell (0xFFFF, 0xFFFF). Several files of one layer
(north-america-01..03) are one layer: their tiles are merged per level, and their levels, strings and attributes must
agree. Tiles are copied as they are (their encoded bytes), never decoded again.

DETERMINISTIC: the same tiles give the same piece bytes in every release, so an unchanged cell keeps its SHA-256 and
is never sent to a display again. So a piece holds no build time, no build string and no source date (the manifest's
record carries the source date); zones keep theirs (the display shows it beside the zones).

Writes OUT_DIR/pieces/<sha256 hex>.pmt and OUT_DIR/R<release, 8 digits>.man.unsigned (the header and the records:
sign it with sign_manifest.py). OUT_DIR must not exist or must be empty.
"""

import argparse
import sys
import time
from pathlib import Path

import manifest
import pmt

CELL_ZOOM = manifest.CELL_ZOOM


class SplitError(Exception):
    pass


def level_params(lv):
    return (lv.tile_zoom, lv.zoom_min, lv.zoom_max, lv.coord_bits, lv.buffer, lv.tolerance_dm, lv.full_class,
            lv.grid is not None)


def plain_strings(f, kind):
    """The file's strings with the build string, and (but for zones) the source date, emptied: their indices stay,
    so every tile reads as before, and the bytes no longer change from one build to the next."""
    strings = [s.decode("utf-8") for s in f.strings]
    blank = ["build"] + ([] if kind == pmt.KIND_ZONES else ["source_date"])
    for name in blank:
        sid = f.ids.get(name, pmt.NO_STRING)
        if sid != pmt.NO_STRING:
            strings[sid] = ""
    return strings


def layer_pieces(name, paths):
    """{(cell_x, cell_y): pmt bytes} for one layer (UNCUT for its low levels), and the layer's PMT kind."""
    layer_id, want_kind = manifest.LAYERS[name]
    files = [pmt.open_pmt_mapped(p) for p in paths]
    first = files[0]
    if any(f.kind != want_kind for f in files):
        raise SplitError(f"{name}: a file is not of layer kind {want_kind}")
    params = [level_params(lv) for lv in first.levels]
    strings = plain_strings(first, first.kind)
    for f, p in zip(files[1:], paths[1:]):
        if [level_params(lv) for lv in f.levels] != params:
            raise SplitError(f"{name}: {p} has other levels than {paths[0]}")
        if plain_strings(f, f.kind) != strings or f.attrs != first.attrs or f.ids != first.ids:
            raise SplitError(f"{name}: {p} has other strings or attributes than {paths[0]}")
    ids = {k: (None if v == pmt.NO_STRING else v) for k, v in first.ids.items()}
    attrs = list(first.attrs)
    keep_date = first.kind == pmt.KIND_ZONES

    # Every tile of every file, by level and cell, as index ranges: the index is in Morton order, and the Morton key's
    # high bits are the zoom-6 ancestor's, so a cell's tiles are one run of the index (memory stays small for a 2 GB
    # file with 7 million tiles).
    by_cell = {}  # cell -> [ [(file, first index, end index), ...] per level ]
    nlev = len(params)
    for f in files:
        for i, lv in enumerate(f.levels):
            tz = lv.tile_zoom
            run_cell, run_start = None, 0
            for j in range(lv.entry_count + 1):
                if j < lv.entry_count:
                    x, y = pmt.unmorton(f._entry(lv, j)[0])
                    cell = ((x >> (tz - CELL_ZOOM), y >> (tz - CELL_ZOOM)) if tz >= CELL_ZOOM
                            else (manifest.UNCUT, manifest.UNCUT))
                else:
                    cell = None
                if cell != run_cell:
                    if run_cell is not None:
                        by_cell.setdefault(run_cell, [[] for _ in range(nlev)])[i].append((f, run_start, j))
                    run_cell, run_start = cell, j
    # Cells that only the coverage grid names (all water, no tile): they still need their piece.
    for i, lv in enumerate(first.levels):
        if lv.grid is None:
            continue
        for f in files:
            g = f.levels[i].grid
            for idx in range(pmt.GRID_SIDE * pmt.GRID_SIDE):
                if (g[idx // 4] >> (2 * (idx % 4))) & 3 not in (0, 3):
                    cx, cy = idx % pmt.GRID_SIDE, idx // pmt.GRID_SIDE
                    by_cell.setdefault((cx >> 2, cy >> 2), [[] for _ in range(nlev)])

    pieces = {}
    for cell in sorted(by_cell):
        uncut = cell[0] == manifest.UNCUT
        slots = by_cell[cell]
        levels = []
        for i, (tz, zmin, zmax, bits, buf, tol, full, has_grid) in enumerate(params):
            if (tz < CELL_ZOOM) != uncut:
                continue
            tiles = {}
            for f, start, end in slots[i]:
                lv = f.levels[i]
                for j in range(start, end):
                    e = f._entry(lv, j)
                    xy = pmt.unmorton(e[0])
                    if xy in tiles:
                        raise SplitError(f"{name}: tile {xy[0]},{xy[1]} of level {i} is in two files")
                    tiles[xy] = f.tile_blob(i, e)
            grid = None
            if has_grid:
                grid = {}
                for cx in range(cell[0] * 4, cell[0] * 4 + 4):
                    for cy in range(cell[1] * 4, cell[1] * 4 + 4):
                        for f in files:
                            v = f.coverage(i, cx << (tz - 8), cy << (tz - 8))
                            if v and grid.get((cx, cy), v) != v:
                                raise SplitError(f"{name}: zoom-8 cell {cx},{cy} has two coverage values")
                            if v:
                                grid[(cx, cy)] = v
            levels.append({"tile_zoom": tz, "zoom_min": zmin, "zoom_max": zmax, "coord_bits": bits, "buffer": buf,
                           "tolerance_dm": tol, "full_class": full, "grid": grid, "tiles": tiles})
        if not levels or not any(lv["tiles"] or lv["grid"] for lv in levels):
            continue
        data = pmt.build_pmt(layer_kind=first.kind, levels=levels, strings=strings, ids=ids, attrs=attrs,
                             build_time=0, data_time=first.data_time if keep_date else 0)
        pieces[cell] = data
    return layer_id, first.kind, pieces


def main(argv):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("out_dir")
    ap.add_argument("--release", type=int, required=True)
    ap.add_argument("--pack-name", required=True)
    ap.add_argument("--min-fw", default="")
    ap.add_argument("--source-date", type=int, default=int(time.strftime("%Y%m%d", time.gmtime())))
    ap.add_argument("--layer", action="append", required=True, help="NAME=file.pmt[,file.pmt...]")
    a = ap.parse_args(argv[1:])
    out = Path(a.out_dir)
    if out.exists() and any(out.iterdir()):
        raise SplitError(f"{out} is not empty")
    (out / "pieces").mkdir(parents=True, exist_ok=True)
    records, seen = [], set()
    for spec in a.layer:
        name, _, files = spec.partition("=")
        if name not in manifest.LAYERS or not files or name in seen:
            raise SplitError(f"--layer {spec}: NAME=files, NAME one of {', '.join(manifest.LAYERS)}, once each")
        seen.add(name)
        layer_id, kind, pieces = layer_pieces(name, files.split(","))
        for (cx, cy), data in pieces.items():
            path = out / "pieces" / "tmp.pmt"
            path.write_bytes(data)
            sha = manifest.sha256_file(path)
            final = out / "pieces" / (sha.hex() + ".pmt")
            if final.exists():
                path.unlink()  # the same bytes in two layers or cells: stored once
            else:
                path.rename(final)
            records.append({"sha256": sha, "size": len(data), "layer": layer_id, "pmt_minor": pmt.kind_minor(kind),
                            "cell_x": cx, "cell_y": cy, "source_date": a.source_date})
        print(f"{name}: {len(pieces)} pieces, {sum(len(d) for d in pieces.values())} bytes")
    body = manifest.encode_body(release=a.release, min_fw=manifest.pack_version(a.min_fw), built=int(time.time()),
                                pack_name=a.pack_name, records=records)
    target = out / f"R{a.release:08d}.man.unsigned"
    target.write_bytes(body)
    print(f"{target}: {len(records)} pieces, {sum(r['size'] for r in records)} bytes; sign it with sign_manifest.py")
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main(sys.argv))
    except (SplitError, manifest.ManifestError, pmt.PmtError) as e:
        print(f"split_pieces: {e}", file=sys.stderr)
        sys.exit(1)
