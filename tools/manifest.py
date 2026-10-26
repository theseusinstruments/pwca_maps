"""The map release manifest (THES-260): pwca_panoptes docs/design/PANOPTES_MAP_UPDATES_CONTRACT.md section 2.

A manifest lists every piece of a release (one layer x one zoom-6 cell, a PMT file) by SHA-256 and size, and is
signed once: ECDSA P-256 over SHA-256 of the header and the records, r || s, low-S. Little-endian.

    header 48 B | piece record 56 B x count | signature 64 B

This module encodes, signs and checks it; split_pieces.py writes the unsigned body, sign_manifest.py signs it.
"""

import hashlib
import struct

MAGIC = b"PMAN"
FORMAT = 1
CELL_ZOOM = 6
HEADER = 48
RECORD = 56
SIGNATURE = 64
MAX_PIECES = 4096
MAX_PIECE_SIZE = (1 << 31) - 1
UNCUT = 0xFFFF  # cell x and y of a piece that is not cut by cell
PACK_NAME = 24

# The contract's layer ids (section 1) and the PMT layer kind each one's pieces state.
LAYERS = {
    "water": (1, 2),
    "roads": (2, 4),
    "minor": (3, 4),
    "streets": (4, 8),
    "places": (5, 5),
    "symbols": (6, 7),
    "zones": (7, 3),
    "base": (8, 1),
    "land": (9, 6),
    "buildings": (10, 9),
}

# P-256's group order, for the low-S rule.
P256_N = 0xFFFFFFFF00000000FFFFFFFFFFFFFFFFBCE6FAADA7179E84F3B9CAC2FC632551


class ManifestError(Exception):
    pass


def pack_version(text):
    """'0.4.33' -> major << 24 | minor << 16 | patch; '' or None -> 0 (any)."""
    if not text:
        return 0
    parts = text.split(".")
    if len(parts) != 3 or not all(p.isdigit() for p in parts):
        raise ManifestError(f"not a version: {text!r}")
    major, minor, patch = (int(p) for p in parts)
    if major > 255 or minor > 255 or patch > 65535:
        raise ManifestError(f"version out of range: {text!r}")
    return major << 24 | minor << 16 | patch


def record_key(r):
    return (r["layer"], r["cell_y"], r["cell_x"])


def encode_body(*, release, min_fw, built, pack_name, records):
    """The header and the records (what is signed). records: dicts with sha256 (32 bytes), size, layer, pmt_minor,
    cell_x, cell_y, source_date (yyyymmdd). They are sorted here; a duplicate layer and cell is refused."""
    if not 0 < release <= 0xFFFFFFFF:
        raise ManifestError("release must be 1..2^32-1")
    if not 1 <= len(records) <= MAX_PIECES:
        raise ManifestError(f"1..{MAX_PIECES} pieces")
    name = pack_name.encode("ascii")
    if len(name) > PACK_NAME or any(c < 0x20 or c > 0x7E for c in name):
        raise ManifestError("pack name: at most 24 printable ASCII characters")
    recs = sorted(records, key=record_key)
    for a, b in zip(recs, recs[1:]):
        if record_key(a) == record_key(b):
            raise ManifestError(f"two pieces for layer {a['layer']} cell {a['cell_x']},{a['cell_y']}")
    out = bytearray(struct.pack("<4sBBHIIII", MAGIC, FORMAT, CELL_ZOOM, HEADER, release, min_fw, built, len(recs)))
    out += name.ljust(PACK_NAME, b"\0")
    for r in recs:
        if len(r["sha256"]) != 32 or not 1 <= r["size"] <= MAX_PIECE_SIZE:
            raise ManifestError("bad piece record")
        out += struct.pack("<32sQBBHHHII", r["sha256"], r["size"], r["layer"], r["pmt_minor"], r["cell_x"],
                           r["cell_y"], 0, r["source_date"], 0)
    return bytes(out)


def low_s_signature(private_key, body):
    """r || s (64 bytes) of an ECDSA P-256 / SHA-256 signature, s in the lower half (the display takes low-S only)."""
    from cryptography.hazmat.primitives import hashes
    from cryptography.hazmat.primitives.asymmetric import ec
    from cryptography.hazmat.primitives.asymmetric.utils import decode_dss_signature

    r, s = decode_dss_signature(private_key.sign(body, ec.ECDSA(hashes.SHA256())))
    if s > P256_N // 2:
        s = P256_N - s
    return r.to_bytes(32, "big") + s.to_bytes(32, "big")


def verify_signature(public_sec1, body, sig):
    """True when sig (r || s) is a low-S ECDSA P-256 / SHA-256 signature of body by the uncompressed SEC1 key."""
    from cryptography.exceptions import InvalidSignature
    from cryptography.hazmat.primitives import hashes
    from cryptography.hazmat.primitives.asymmetric import ec
    from cryptography.hazmat.primitives.asymmetric.utils import encode_dss_signature

    if len(sig) != SIGNATURE:
        return False
    r, s = int.from_bytes(sig[:32], "big"), int.from_bytes(sig[32:], "big")
    if not (0 < r < P256_N and 0 < s <= P256_N // 2):
        return False
    key = ec.EllipticCurvePublicKey.from_encoded_point(ec.SECP256R1(), public_sec1)
    try:
        key.verify(encode_dss_signature(r, s), body, ec.ECDSA(hashes.SHA256()))
        return True
    except InvalidSignature:
        return False


def parse(data, public_keys=()):
    """Checks a manifest as a display does (the contract's section 2.4, without the release and firmware
    comparisons, which need the unit) and returns its fields. public_keys: SEC1 keys any of which may have signed it;
    empty: the signature is not checked (the builder's own read-back of an unsigned body is parse_body)."""
    if len(data) < HEADER + SIGNATURE:
        raise ManifestError("too short")
    count = struct.unpack_from("<I", data, 20)[0]
    if not 1 <= count <= MAX_PIECES or len(data) != HEADER + RECORD * count + SIGNATURE:
        raise ManifestError("length is not 48 + 56 x count + 64")
    body, sig = data[:-SIGNATURE], data[-SIGNATURE:]
    if public_keys and not any(verify_signature(k, body, sig) for k in public_keys):
        raise ManifestError("signature by no key that counts")
    return parse_body(body)


def parse_body(body):
    magic, fmt, zoom, hlen, release, min_fw, built, count = struct.unpack_from("<4sBBHIIII", body, 0)
    if magic != MAGIC or fmt != FORMAT or zoom != CELL_ZOOM or hlen != HEADER:
        raise ManifestError("bad header")
    if not 1 <= count <= MAX_PIECES or len(body) != HEADER + RECORD * count:
        raise ManifestError("bad count or length")
    if release == 0:
        raise ManifestError("release 0")
    name = body[24:48]
    text = name.rstrip(b"\0")
    if any(c < 0x20 or c > 0x7E for c in text) or any(name[len(text):]):
        raise ManifestError("bad pack name")
    records = []
    for i in range(count):
        sha, size, layer, minor, cx, cy, res1, date, res2 = struct.unpack_from("<32sQBBHHHII", body,
                                                                               HEADER + RECORD * i)
        if not 1 <= size <= MAX_PIECE_SIZE or res1 or res2:
            raise ManifestError(f"record {i}: bad size or reserved field")
        records.append({"sha256": sha, "size": size, "layer": layer, "pmt_minor": minor, "cell_x": cx,
                        "cell_y": cy, "source_date": date})
    for a, b in zip(records, records[1:]):
        if not record_key(a) < record_key(b):
            raise ManifestError("records not strictly sorted")
    return {"release": release, "min_fw": min_fw, "built": built, "pack_name": text.decode("ascii"),
            "records": records}


def sha256_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 22), b""):
            h.update(chunk)
    return h.digest()
