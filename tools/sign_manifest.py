"""Signs a map release's manifest (THES-260; pwca_panoptes docs/design/PANOPTES_MAP_UPDATES_CONTRACT.md section 2.3).

Usage: python tools/sign_manifest.py R00000012.man.unsigned KEY.pem R00000012.man [--pieces DIR]

KEY.pem is an ECDSA P-256 private key (PEM, unencrypted, or set MAP_KEY_PASSWORD). The release key itself lives on
the server (PWC-MAP-RELEASE, the dev one PWC-MAP-RELEASE-DEV; the owner's choice 2); this tool is for the dev key and
the tests. Before signing, the body is read back as a display reads it, and with --pieces every piece it names is
checked against its SHA-256 and size. The signature is low-S (the display takes low-S only), and is checked once
more with the key's public half before the file is written. The public key (uncompressed SEC1, hex) is printed: it
is what goes into the display's trust table.
"""

import os
import sys
from pathlib import Path

import manifest


def main(argv):
    args = [a for a in argv[1:] if not a.startswith("--")]
    pieces = None
    if "--pieces" in argv:
        pieces = Path(argv[argv.index("--pieces") + 1])
        args = [a for a in args if Path(a) != pieces]
    if len(args) != 3:
        print(__doc__, file=sys.stderr)
        return 2
    body_path, key_path, out_path = (Path(a) for a in args)
    from cryptography.hazmat.primitives import serialization
    from cryptography.hazmat.primitives.asymmetric import ec

    body = body_path.read_bytes()
    info = manifest.parse_body(body)
    if pieces is not None:
        for r in info["records"]:
            p = pieces / (r["sha256"].hex() + ".pmt")
            if not p.is_file() or p.stat().st_size != r["size"] or manifest.sha256_file(p) != r["sha256"]:
                raise manifest.ManifestError(f"{p}: missing, or not the piece the manifest names")
    password = os.environ.get("MAP_KEY_PASSWORD")
    key = serialization.load_pem_private_key(key_path.read_bytes(), password.encode() if password else None)
    if not isinstance(key, ec.EllipticCurvePrivateKey) or key.curve.name != "secp256r1":
        raise manifest.ManifestError("the key is not ECDSA P-256")
    pub = key.public_key().public_bytes(serialization.Encoding.X962, serialization.PublicFormat.UncompressedPoint)
    sig = manifest.low_s_signature(key, body)
    if not manifest.verify_signature(pub, body, sig):
        raise manifest.ManifestError("the signature does not check with the key's own public half")
    out_path.write_bytes(body + sig)
    manifest.parse(body + sig, [pub])
    print(f"{out_path}: release {info['release']}, {len(info['records'])} pieces, signed by {pub.hex()}")
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main(sys.argv))
    except manifest.ManifestError as e:
        print(f"sign_manifest: {e}", file=sys.stderr)
        sys.exit(1)
