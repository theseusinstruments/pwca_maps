"""Tests for the map release pieces and their manifest (tools/split_pieces.py, tools/manifest.py,
tools/sign_manifest.py; THES-260, pwca_panoptes docs/design/PANOPTES_MAP_UPDATES_CONTRACT.md)."""

import struct
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "tools"))

import manifest  # noqa: E402
import pmt  # noqa: E402
import sign_manifest  # noqa: E402
import split_pieces  # noqa: E402


def square(n, cls=1):
    return [(pmt.POLYGON, cls, None, [[(10, 10), (n, 10), (n, n), (10, n)]])]


def water_file(path, tiles12, tiles4=(), grid=None, build="abc1234", source_date="2026-09"):
    """A water PMT: level 0 tile zoom 4 (no grid), level 1 tile zoom 8 (grid), level 2 tile zoom 12."""
    t8 = {}
    for (x, y) in tiles12:
        t8[(x >> 4, y >> 4)] = square(200)
    levels = [
        {"tile_zoom": 4, "zoom_min": 4, "zoom_max": 7, "coord_bits": 12, "buffer": 64, "tolerance_dm": 4000,
         "tiles": {xy: square(300) for xy in tiles4}},
        {"tile_zoom": 8, "zoom_min": 8, "zoom_max": 11, "coord_bits": 12, "buffer": 64, "tolerance_dm": 400,
         "grid": grid if grid is not None else {}, "tiles": t8},
        {"tile_zoom": 12, "zoom_min": 12, "zoom_max": 16, "coord_bits": 12, "buffer": 64, "tolerance_dm": 40,
         "tiles": {xy: square(100 + (xy[0] % 7)) for xy in tiles12}},
    ]
    strings = ["osm-water", "(c) OpenStreetMap contributors", "ODbL 1.0", source_date, build]
    ids = {"layer": 0, "credit": 1, "license": 2, "source": None, "source_date": 3, "build": 4}
    data = pmt.build_pmt(layer_kind=pmt.KIND_WATER, levels=levels, strings=strings, ids=ids, build_time=1234,
                         data_time=5678)
    Path(path).write_bytes(data)
    return path


def cell_of(x, y, tz):
    return x >> (tz - 6), y >> (tz - 6)


# Tiles at zoom 12 in three zoom-6 cells: (10, 20) with two tiles, (11, 20) with one, (40, 30) with one.
T12 = [(10 * 64 + 3, 20 * 64 + 5), (10 * 64 + 60, 20 * 64 + 1), (11 * 64, 20 * 64 + 63), (40 * 64 + 9, 30 * 64 + 9)]


def split(out, *layer_args, release=7):
    args = ["split_pieces.py", str(out), "--release", str(release), "--pack-name", "NA-TEST", "--min-fw", "0.4.34",
            "--source-date", "20261101"]
    for a in layer_args:
        args += ["--layer", a]
    assert split_pieces.main(args) == 0
    return (out / f"R{release:08d}.man.unsigned").read_bytes()


class Split(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.d = Path(self.tmp.name)

    def tearDown(self):
        self.tmp.cleanup()

    def pieces(self, out, body):
        info = manifest.parse_body(body)
        got = {}
        for r in info["records"]:
            p = out / "pieces" / (r["sha256"].hex() + ".pmt")
            self.assertEqual(manifest.sha256_file(p), r["sha256"])
            self.assertEqual(p.stat().st_size, r["size"])
            got[(r["layer"], r["cell_x"], r["cell_y"])] = pmt.open_pmt(p)
        return info, got

    def test_cells_partition_the_tiles(self):
        src = water_file(self.d / "w.pmt", T12, tiles4=[(2, 5), (3, 5)], grid={(10 * 4 + 1, 20 * 4 + 2): 2})
        body = split(self.d / "out", f"water={src}")
        info, got = self.pieces(self.d / "out", body)
        self.assertEqual(info["release"], 7)
        self.assertEqual(info["min_fw"], 0x00040022)
        self.assertEqual(info["pack_name"], "NA-TEST")
        self.assertEqual(set(got), {(1, 10, 20), (1, 11, 20), (1, 40, 30), (1, 0xFFFF, 0xFFFF)})
        original = pmt.open_pmt(src)
        for (layer, cx, cy), f in got.items():
            self.assertEqual(f.kind, pmt.KIND_WATER)
            self.assertEqual((f.build_time, f.data_time), (0, 0), "no times in a piece")
            self.assertEqual(f.string(f.ids["build"]), "", "no build string in a piece")
            self.assertEqual(f.string(f.ids["source_date"]), "", "no source date in a piece but zones")
            self.assertEqual(f.string(f.ids["layer"]), "osm-water")
            if cx == 0xFFFF:
                self.assertEqual([lv.tile_zoom for lv in f.levels], [4])
                self.assertEqual(sorted(pmt.unmorton(e[0]) for e in f.entries(0)), [(2, 5), (3, 5)])
                continue
            self.assertEqual([lv.tile_zoom for lv in f.levels], [8, 12])
            for i, lv in enumerate(f.levels):
                for e in f.entries(i):
                    x, y = pmt.unmorton(e[0])
                    self.assertEqual(cell_of(x, y, lv.tile_zoom), (cx, cy))
                    # The tile's bytes are the source's.
                    self.assertEqual(f.tile_blob(i, e), original.tile_blob(i + 1, original.find(i + 1, x, y)))
        # Every source tile is in exactly one piece.
        for i, tz in ((1, 8), (2, 12)):
            for e in original.entries(i):
                x, y = pmt.unmorton(e[0])
                cell = cell_of(x, y, tz)
                f = got[(1,) + cell]
                self.assertIsNotNone(f.find(i - 1, x, y))
        # The coverage grid: the zoom-8 cell inside (10, 20) keeps its value; nothing outside the cell has one.
        f = got[(1, 10, 20)]
        self.assertEqual(f.coverage(0, (41 << 0), (82 << 0)), 2)
        self.assertEqual(f.coverage(0, 0, 0), 0)
        self.assertEqual(got[(1, 40, 30)].coverage(0, 41, 82), 0)

    def test_a_cell_only_the_grid_names_gets_a_piece(self):
        src = water_file(self.d / "w.pmt", T12, grid={(50 * 4, 25 * 4): 1})
        body = split(self.d / "out", f"water={src}")
        info, got = self.pieces(self.d / "out", body)
        self.assertIn((1, 50, 25), got)
        self.assertEqual(got[(1, 50, 25)].coverage(0, 200, 100), 1)

    def test_deterministic_across_builds(self):
        """A new build (other build time, commit and source date) with the same tiles: the same piece bytes; one
        changed tile changes only its own cell's piece."""
        a = water_file(self.d / "a.pmt", T12)
        b = water_file(self.d / "b.pmt", T12, build="def5678", source_date="2026-10")
        changed = [(x, y) for (x, y) in T12 if (x, y) != T12[3]] + [(T12[3][0] + 1, T12[3][1])]
        c = water_file(self.d / "c.pmt", changed)
        ra = {(r["cell_x"], r["cell_y"]): r["sha256"] for r in manifest.parse_body(split(self.d / "oa", f"water={a}"))["records"]}
        rb = {(r["cell_x"], r["cell_y"]): r["sha256"] for r in manifest.parse_body(split(self.d / "ob", f"water={b}"))["records"]}
        rc = {(r["cell_x"], r["cell_y"]): r["sha256"] for r in manifest.parse_body(split(self.d / "oc", f"water={c}"))["records"]}
        self.assertEqual(ra, rb)
        self.assertEqual({k for k in ra if ra[k] != rc[k]}, {(40, 30)})

    def test_files_of_one_layer_merge(self):
        a = water_file(self.d / "a.pmt", T12[:2])
        b = water_file(self.d / "b.pmt", T12[2:])
        whole = water_file(self.d / "w.pmt", T12)
        merged = manifest.parse_body(split(self.d / "o1", f"water={a},{b}"))
        single = manifest.parse_body(split(self.d / "o2", f"water={whole}"))
        self.assertEqual([r["sha256"] for r in merged["records"]], [r["sha256"] for r in single["records"]])
        dup = water_file(self.d / "d.pmt", T12[:1])
        with self.assertRaises(split_pieces.SplitError):
            split(self.d / "o3", f"water={a},{dup}")

    def test_layers_and_kinds(self):
        src = water_file(self.d / "w.pmt", T12)
        with self.assertRaises(split_pieces.SplitError):
            split(self.d / "o1", f"roads={src}")  # a water file named roads
        with self.assertRaises(split_pieces.SplitError):
            split(self.d / "o2", f"rivers={src}")


class Manifest(unittest.TestCase):
    def setUp(self):
        from cryptography.hazmat.primitives.asymmetric import ec
        from cryptography.hazmat.primitives import serialization

        self.key = ec.generate_private_key(ec.SECP256R1())
        self.pub = self.key.public_key().public_bytes(serialization.Encoding.X962,
                                                      serialization.PublicFormat.UncompressedPoint)
        self.recs = [{"sha256": bytes([i]) * 32, "size": 1000 + i, "layer": 1 + i % 3, "pmt_minor": 0,
                      "cell_x": i * 3 % 64, "cell_y": i * 5 % 64, "source_date": 20261101} for i in range(20)]

    def body(self, **kw):
        args = {"release": 12, "min_fw": manifest.pack_version("0.4.34"), "built": 1790000000, "pack_name": "NA",
                "records": self.recs}
        args.update(kw)
        return manifest.encode_body(**args)

    def test_layout(self):
        b = self.body()
        self.assertEqual(len(b), 48 + 56 * 20)
        self.assertEqual(b[:4], b"PMAN")
        self.assertEqual(struct.unpack_from("<BBHIIII", b, 4), (1, 6, 48, 12, 0x00040022, 1790000000, 20))
        keys = [(r["layer"], r["cell_y"], r["cell_x"]) for r in manifest.parse_body(b)["records"]]
        self.assertEqual(keys, sorted(keys))
        self.assertEqual(manifest.encode_body(release=12, min_fw=0x00040022, built=1790000000, pack_name="NA",
                                              records=list(reversed(self.recs))), b, "the order given does not matter")

    def test_sign_and_verify(self):
        b = self.body()
        for _ in range(16):  # low-S every time
            sig = manifest.low_s_signature(self.key, b)
            self.assertLessEqual(int.from_bytes(sig[32:], "big"), manifest.P256_N // 2)
            self.assertTrue(manifest.verify_signature(self.pub, b, sig))
        info = manifest.parse(b + sig, [self.pub])
        self.assertEqual(info["release"], 12)
        # High-S (the same signature's twin) is refused, as the display refuses it.
        s = int.from_bytes(sig[32:], "big")
        high = sig[:32] + (manifest.P256_N - s).to_bytes(32, "big")
        self.assertFalse(manifest.verify_signature(self.pub, b, high))
        # Any changed byte of the body, or another key: refused.
        for pos in (0, 9, 30, 48, 48 + 40, len(b) - 1):
            bad = bytearray(b + sig)
            bad[pos] ^= 1
            with self.assertRaises(manifest.ManifestError):
                manifest.parse(bytes(bad), [self.pub])
        from cryptography.hazmat.primitives.asymmetric import ec
        from cryptography.hazmat.primitives import serialization
        other = ec.generate_private_key(ec.SECP256R1()).public_key().public_bytes(
            serialization.Encoding.X962, serialization.PublicFormat.UncompressedPoint)
        with self.assertRaises(manifest.ManifestError):
            manifest.parse(b + sig, [other])

    def test_refusals(self):
        b = self.body()
        with self.assertRaises(manifest.ManifestError):
            manifest.parse_body(b + b"\0")  # length
        for pos, val in ((4, 2), (5, 7), (6, 47), (8, 0), (48 + 46, 1), (48 + 52, 1), (30, 0x01)):
            bad = bytearray(b)
            if pos == 8:
                bad[8:12] = b"\0\0\0\0"
            else:
                bad[pos] = val
            with self.assertRaises(manifest.ManifestError, msg=f"byte {pos}"):
                manifest.parse_body(bytes(bad))
        swapped = bytearray(b)  # two records out of order
        swapped[48:48 + 56], swapped[48 + 56:48 + 112] = b[48 + 56:48 + 112], b[48:48 + 56]
        with self.assertRaises(manifest.ManifestError):
            manifest.parse_body(bytes(swapped))
        with self.assertRaises(manifest.ManifestError):
            self.body(records=self.recs + [dict(self.recs[0], sha256=b"\x99" * 32)])  # one layer and cell twice
        with self.assertRaises(manifest.ManifestError):
            self.body(pack_name="NA\n")
        with self.assertRaises(manifest.ManifestError):
            self.body(records=[])

    def test_sign_tool(self):
        from cryptography.hazmat.primitives import serialization

        with tempfile.TemporaryDirectory() as t:
            d = Path(t)
            src = water_file(d / "w.pmt", T12)
            body = split(d / "out", f"water={src}")
            (d / "key.pem").write_bytes(self.key.private_bytes(serialization.Encoding.PEM,
                                                               serialization.PrivateFormat.PKCS8,
                                                               serialization.NoEncryption()))
            unsigned = d / "out" / "R00000007.man.unsigned"
            out = d / "R00000007.man"
            self.assertEqual(sign_manifest.main(["x", str(unsigned), str(d / "key.pem"), str(out), "--pieces",
                                                 str(d / "out" / "pieces")]), 0)
            data = out.read_bytes()
            self.assertEqual(data[:-64], body)
            manifest.parse(data, [self.pub])
            # A piece that is not the one named: refused before signing.
            first = sorted((d / "out" / "pieces").iterdir())[0]
            first.write_bytes(first.read_bytes() + b"x")
            with self.assertRaises(manifest.ManifestError):
                sign_manifest.main(["x", str(unsigned), str(d / "key.pem"), str(d / "R2.man"), "--pieces",
                                    str(d / "out" / "pieces")])


if __name__ == "__main__":
    unittest.main()
