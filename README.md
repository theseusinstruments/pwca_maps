# pwca_maps: offline map packs

Build recipe for the offline North American water and shoreline map data used on our helm display: coastlines,
lakes, reservoirs and rivers. Riders see where they are without a network connection.

The map packs are a **Derivative Database of OpenStreetMap** and are licensed under the
[Open Database License (ODbL) 1.0](https://opendatacommons.org/licenses/odbl/1-0/). This repository is the
published method of making them (ODbL section 4.6). With it and the source extracts named in each pack's
`SOURCES.json`, anyone can rebuild a pack.

## Sources

| Layer | Source | Licence |
|---|---|---|
| Detail (z10-16): shorelines, lakes, reservoirs, rivers | [OpenStreetMap](https://www.openstreetmap.org/copyright) extracts from [Geofabrik](https://download.geofabrik.de/north-america.html); sea polygons from [osmdata.openstreetmap.de](https://osmdata.openstreetmap.de/data/water-polygons.html) | ODbL 1.0 |
| Overview (z0-9): land, coast, large lakes and rivers | [Natural Earth](https://www.naturalearthdata.com/) 1:10m, with the North America supplements | Public domain |
| Speed-restricted / no-wake areas (planned) | Government datasets, e.g. Florida FWC State Boating Safety Zones | Per source |

Attribution shown on the device and in its documentation:

> Map data © OpenStreetMap contributors, available under the Open Database License (ODbL):
> openstreetmap.org/copyright. Overview map made with Natural Earth. Not for navigation.

Data from different sources is kept in separate layers (a Collective Database), so non-OSM layers keep their own
terms.

## Pack format

The file format is frozen as **PMT 1.0**: [FORMAT.md](FORMAT.md). `tools/pmt.py` is its reference encoder and
decoder, and `tests/vectors/` pins it (`python -m unittest discover -s tests`). The design below is how the pack
builder uses it.

- Web Mercator z12 tiles. Each tile holds three geometry levels:

  | Level | Zoom | Simplification tolerance | Units per tile side |
  |---|---|---|---|
  | L1 | z10-11 | about 40 m | 1024 |
  | L2 | z12-13 | about 10 m | 4096 |
  | L3 | z14-16 | about 4 m | 65536 |

- Coordinates are tile-local, stored as zigzag deltas encoded as varints.
- Closed water bodies under 1 ha are dropped.
- Packs are FAT32-friendly: no file over 2 GiB. A large extract (all of North America) is written as several
  files, whole regions each, every file with the coverage grid of its own regions.
- Each pack root carries `ATTRIBUTION.txt`, `LICENSE.txt` (ODbL 1.0 notice), `SOURCES.json` (source URLs,
  timestamps, SHA-256, tool version) and a SHA-256 manifest.
- Packs are never encrypted or locked (ODbL section 4.7).

## Tools

- `tools/measure_water.py <extract.osm.pbf> [--json out.json]`: reads an OSM extract and reports how large its water
  layers become in the pack encoding, per level. It needs Python 3 and numpy only; it includes a minimal PBF reader.
- `tools/build_pack.py <extract.osm.pbf> <out_dir>`: builds the detail water layer of an extract (water areas as
  ways and multipolygons, river and canal lines; intermittent water dropped; areas under 1 ha dropped), simplified per
  level, cut into buffered tiles, written as one PMT file with `SOURCES.json`, `ATTRIBUTION.txt`, `LICENSE.txt` and
  `MANIFEST.sha256`. For a coastal extract, `--sea water-polygons-split-4326.zip` adds the sea from the osmdata
  water polygons; `--processes N` runs the passes in N worker processes.
- `tools/build_zones.py <out_dir> --boating B.geojson --manatee M.geojson --date YYYY-MM`: the speed-restricted
  and no-wake zones layer from the Florida FWC State Boating Safety Zones (FAC 68D-24) and State Manatee
  Protection Zones (FAC 68C-22), mapped cautiously (seasons by whole months, partial times flagged; anchoring and
  propeller rules left out); the coverage grid holds only the cells with zones.
- `tools/fetch_fwc_zones.py <out_dir>`: fetches those two FWC layers as GeoJSON (WGS 84) in pages by OBJECTID (the
  server cuts large responses short), checks every feature arrived, and records URLs, time and SHA-256.
- `tools/split_pieces.py OUT_DIR --release N --pack-name NAME --layer water=a.pmt,b.pmt ...`: cuts a built pack
  into the pieces of a map release (one layer x one zoom-6 cell, tiles copied as they are, deterministic so an
  unchanged cell keeps its SHA-256) and writes the release's unsigned manifest; `tools/sign_manifest.py` signs it
  (ECDSA P-256, low-S). The format and the rules: pwca_panoptes `docs/design/PANOPTES_MAP_UPDATES_CONTRACT.md`
  (THES-260, map updates over Wi-Fi and by Hermes).
- `tools/pmt.py dump|check FILE`: prints a PMT file's canonical dump, or checks every tile.
- `tools/pmt_crop.py`: cuts a region out of a PMT file by tile ranges (small test packs).
- `tools/build_overview.py <natural_earth_dir> <out.pmt>`: the overview (base layer, z2-9) from the six Natural
  Earth 1:10m zips named in the script.
- `tools/make_test_vectors.py`, `tools/make_synthetic_pack.py`: the test vectors, and invented geometry for renderer
  benches (labelled as such in the file).
- Library modules: `osmpbf.py` (a minimal PBF reader), `packgeom.py` (projection, simplification, tile cutting),
  `featurestore.py` (the builder's on-disk feature store, memory-mapped by the worker processes).

Tests: `python -m unittest discover -s tests`.

Source extracts go in `sources/` and outputs in `out/`. Both are git-ignored.

## Licence

- Map data produced by these tools: ODbL 1.0 (see above).
- The tool code, format specification and test vectors in this repository: Apache License 2.0 ([LICENSE](LICENSE)).
  The map data the tools produce is not covered by it: it stays under the ODbL (and the terms of each non-OSM
  layer's source).
