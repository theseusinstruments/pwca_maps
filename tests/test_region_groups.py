"""tools/region_groups.py: the group rule (THES-263) and the catalogue's body as Hermes reads it."""
import hashlib
import os
import struct
import sys
import unittest

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "tools"))
import region_groups as rg  # noqa: E402


class GroupRule(unittest.TestCase):
    def test_a_region_under_the_small_cover_joins_the_owner_of_its_largest_cell(self):
        cov = {
            "US-VA": {(18, 24): 0.25, (17, 24): 0.20},
            "US-PA": {(18, 24): 0.19, (18, 23): 0.15},   # not small: stands alone, though it owns no cell here
            "US-DC": {(18, 24): 0.002},
            "US-NJ": {(18, 24): 0.089, (19, 23): 0.01},
            "US-RI": {(19, 23): 0.017},
            "US-ME": {(19, 23): 0.255},
        }
        g = rg.groups(cov)
        self.assertEqual(g["US-DC"], "US-VA")
        self.assertEqual(g["US-NJ"], "US-VA")
        self.assertEqual(g["US-RI"], "US-ME")
        for code in ("US-VA", "US-PA", "US-ME"):
            self.assertEqual(g[code], code)

    def test_ten_percent_is_the_line_and_there_are_no_chains(self):
        cov = {"US-A": {(1, 1): 0.10, (2, 1): 0.5}, "US-B": {(1, 1): 0.5}, "US-C": {(3, 3): 0.05}, "US-D": {(3, 3): 0.06}}
        g = rg.groups(cov)
        self.assertEqual(g["US-A"], "US-A", "10% is not under 10%")
        # C and D are both small: neither leads the other (no chain into a clubbed host).
        self.assertEqual((g["US-C"], g["US-D"]), ("US-C", "US-D"))

    def test_the_body_is_the_layout_hermes_reads(self):
        cov = {"US-AA": {(10, 22): 0.6}, "US-DD": {(10, 22): 0.01}, "CA-CC": {(11, 22): 0.9}}
        rg.NAMES.update({"US-AA": "Alpha", "US-DD": "Delta", "CA-CC": "Charlie"})
        tbl = rg.table(2026102201, cov, rg.groups(cov))
        sha = hashlib.sha256(b"manifest").digest()
        body = rg.reg_body(tbl, sha)
        self.assertEqual(body[:4], b"PREG")
        self.assertEqual(struct.unpack_from("<HI", body, 6), (48, 2026102201))
        self.assertEqual(body[12:44], sha)
        self.assertEqual(struct.unpack_from("<HH", body, 44), (2, 3))
        # countries 38 bytes each; then regions 48 + 4 per cell
        self.assertEqual(len(body), 48 + 2 * 38 + 3 * (48 + 4))
        by_code = {r["code"]: r for r in tbl["regions"]}
        self.assertEqual(by_code["US-DD"]["group"], by_code["US-AA"]["id"])
        self.assertEqual(by_code["US-AA"]["group"], by_code["US-AA"]["id"])


if __name__ == "__main__":
    unittest.main()
