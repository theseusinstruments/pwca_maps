"""Map tile files (PMT) version 1.3: the reference encoder, decoder and canonical dump. The format is FORMAT.md
(with section 7.7, the land cover of version 1.2) and, for the layer kinds 7 and 8 of version 1.3,
FORMAT-1.3-SYMBOLS-STREET-NAMES.md.

Usage:
  python tools/pmt.py dump FILE      print the canonical dump (FORMAT.md section 9)
  python tools/pmt.py check FILE     open FILE and decode every tile; exit status 1 on the first error

Standard library only. The decoder applies every reader rule in FORMAT.md section 8, so it also serves as the
executable form of those rules.
"""

import mmap
import struct
import sys
import zlib

MAGIC = b"PMTF"
VERSION_MAJOR = 1
VERSION_MINOR = 4  # the newest minor this tool writes and reads; a base, water or zones file is still written as 1.0
FIXED_HEADER = 64
LEVEL_RECORD = 48
INDEX_ENTRY = 16
PAGE_ENTRIES = 256
ALIGN = 4096
MAX_HEADER = 16 << 20
MAX_FILE = (1 << 31) - 1
MAX_TILE = 8 << 20
MAX_LEVELS = 8
MAX_FEATURES = 65535
MAX_PARTS = 65535
MAX_POINTS = 1 << 20
MAX_TILE_ZOOM = 16
MAX_DISPLAY_ZOOM = 22
NO_STRING = 0xFFFF
GRID_SIDE = 256
GRID_BYTES = GRID_SIDE * GRID_SIDE // 4

POLYGON, LINE, POINT = 1, 2, 3
KIND_BASE, KIND_WATER, KIND_ZONES = 1, 2, 3
KIND_ROADS, KIND_PLACES = 4, 5  # FORMAT 1.1 (sections 7.4 and 7.5)
KIND_LAND = 6  # FORMAT 1.2 (section 7.7): land cover
# The geometries a layer kind may hold (FORMAT.md section 6). A kind that is not here is not decoded at all.
# FORMAT 1.3 (FORMAT-1.3-SYMBOLS-STREET-NAMES.md): symbols (points with a class and an optional name) and street names
# (lines that carry their name). Kind 6 (minor version 2) is the land cover, FORMAT.md section 7.7.
KIND_SYMBOLS, KIND_STREETS = 7, 8
# FORMAT 1.4 (FORMAT-1.4-BUILDINGS.md, THES-83): buildings by the water, polygons only.
KIND_BUILDINGS = 9
KIND_GEOMETRIES = {KIND_BASE: (POLYGON, LINE), KIND_WATER: (POLYGON, LINE), KIND_ZONES: (POLYGON, LINE),
                   KIND_ROADS: (LINE,), KIND_PLACES: (POINT,), KIND_LAND: (POLYGON,), KIND_SYMBOLS: (POINT,),
                   KIND_STREETS: (LINE,), KIND_BUILDINGS: (POLYGON,)}
NAME_MAX = 127  # bytes of a name (a point's, a named line's)
LEVEL_FLAG_GRID = 1
ID_FIELDS = ("layer", "credit", "license", "source", "source_date", "build")

# Zone attribute tags (FORMAT.md section 7.3).
TAG_NAME, TAG_SPEED_DKMH, TAG_MONTHS, TAG_DAYS, TAG_FLAGS, TAG_SOURCE_REF = 1, 2, 3, 4, 5, 6


class PmtError(ValueError):
    """The input breaks a rule of FORMAT.md."""


def kind_minor(layer_kind):
    """The minor version a file of this layer kind states: 1 for the kinds 1.1 added, 2 for the land cover (1.2),
    3 for the kinds 1.3 added, 4 for the buildings (1.4), else 0."""
    if layer_kind == KIND_BUILDINGS:
        return 4
    if layer_kind in (KIND_SYMBOLS, KIND_STREETS):
        return 3
    if layer_kind == KIND_LAND:
        return 2
    return 1 if layer_kind in (KIND_ROADS, KIND_PLACES) else 0


def name_ok(raw):
    """Whether a point's name may be drawn: 1 to NAME_MAX bytes of strict UTF-8 with no control character (C0, DEL
    or C1)."""
    if not 1 <= len(raw) <= NAME_MAX:
        return False
    try:
        text = bytes(raw).decode("utf-8", "strict")
    except UnicodeDecodeError:
        return False
    return not any(ord(c) < 0x20 or 0x7F <= ord(c) <= 0x9F for c in text)


def crc32(data):
    return zlib.crc32(data) & 0xFFFFFFFF


def align(n, a=ALIGN):
    return (n + a - 1) // a * a


def morton(x, y):
    key = 0
    for i in range(MAX_TILE_ZOOM):
        key |= ((x >> i) & 1) << (2 * i) | ((y >> i) & 1) << (2 * i + 1)
    return key


def unmorton(key):
    x = y = 0
    for i in range(MAX_TILE_ZOOM):
        x |= ((key >> (2 * i)) & 1) << i
        y |= ((key >> (2 * i + 1)) & 1) << i
    return x, y


def zigzag(v):
    if not -(1 << 31) <= v < (1 << 31):
        raise PmtError(f"delta {v} does not fit 32 bits")
    return ((v << 1) ^ (v >> 31)) & 0xFFFFFFFF


def unzigzag(u):
    return (u >> 1) ^ -(u & 1)


def put_varint(out, v):
    if not 0 <= v < (1 << 32):
        raise PmtError(f"varint value {v} out of range")
    while v >= 0x80:
        out.append((v & 0x7F) | 0x80)
        v >>= 7
    out.append(v)


def get_varint(buf, pos, end):
    value = 0
    for i in range(5):
        if pos >= end:
            raise PmtError("varint runs past the end")
        c = buf[pos]
        pos += 1
        value |= (c & 0x7F) << (7 * i)
        if c < 0x80:
            if value >= 1 << 32:
                raise PmtError("varint over 32 bits")
            return value, pos
    raise PmtError("varint longer than 5 bytes")


# ---- tiles ---------------------------------------------------------------------------------------------------------

def encode_tile(features, coord_bits, buffer, layer_kind=KIND_WATER):
    """features: a list of (geometry, class, attr or None, parts); parts: a list of lists of (x, y) tile units. A
    point (layer kinds 5 and 7) is (POINT, class, None, (x, y, name)) with the name as str or bytes; in kind 7 the name
    may be empty. A named line (layer kind 8) is (LINE, class, None, (name, [(x, y), ...]))."""
    if len(features) > MAX_FEATURES:
        raise PmtError("too many features")
    allowed = KIND_GEOMETRIES.get(layer_kind, ())
    lo, hi = -buffer, (1 << coord_bits) + buffer
    out = bytearray()
    px = py = points = 0
    for geom, cls, attr, parts in features:
        if geom not in allowed or not 0 <= cls < (1 << 29):
            raise PmtError(f"bad feature geometry {geom} or class {cls} for layer kind {layer_kind}")
        if geom == POINT:
            x, y, name = parts
            name = name.encode("utf-8") if isinstance(name, str) else bytes(name)
            if attr is not None or not (name_ok(name) or (layer_kind == KIND_SYMBOLS and not name)):
                raise PmtError(f"bad point feature (an attribute, or the name {name!r})")
            if not (lo <= x <= hi and lo <= y <= hi):
                raise PmtError(f"point ({x}, {y}) outside the buffered tile")
            points += 1
            if points > MAX_POINTS:
                raise PmtError("too many points in the tile")
            put_varint(out, geom | (cls << 3))
            put_varint(out, zigzag(x - px))
            put_varint(out, zigzag(y - py))
            out.append(len(name))
            out += name
            px, py = x, y
            continue
        if layer_kind == KIND_STREETS:
            name, part = parts
            name = name.encode("utf-8") if isinstance(name, str) else bytes(name)
            if attr is not None or not name_ok(name) or len(part) < 2:
                raise PmtError(f"bad named line (an attribute, under 2 points, or the name {name!r})")
            points += len(part)
            if points > MAX_POINTS:
                raise PmtError("too many points in the tile")
            put_varint(out, geom | (cls << 3))
            out.append(len(name))
            out += name
            put_varint(out, len(part))
            for x, y in part:
                if not (lo <= x <= hi and lo <= y <= hi):
                    raise PmtError(f"point ({x}, {y}) outside the buffered tile")
                put_varint(out, zigzag(x - px))
                put_varint(out, zigzag(y - py))
                px, py = x, y
            continue
        if not 1 <= len(parts) <= MAX_PARTS:
            raise PmtError("bad part count")
        put_varint(out, geom | (4 if attr is not None else 0) | (cls << 3))
        if attr is not None:
            put_varint(out, attr)
        put_varint(out, len(parts))
        for part in parts:
            if len(part) < (3 if geom == POLYGON else 2):
                raise PmtError("too few points in a part")
            points += len(part)
            if points > MAX_POINTS:
                raise PmtError("too many points in the tile")
            put_varint(out, len(part))
            for x, y in part:
                if not (lo <= x <= hi and lo <= y <= hi):
                    raise PmtError(f"point ({x}, {y}) outside the buffered tile")
                put_varint(out, zigzag(x - px))
                put_varint(out, zigzag(y - py))
                px, py = x, y
    if len(out) > MAX_TILE:
        raise PmtError("tile over 8 MiB")
    return bytes(out)


def decode_tile(blob, coord_bits, buffer, attr_count, layer_kind=KIND_WATER):
    """The inverse of encode_tile, with every check of FORMAT.md section 6. A point comes back as (POINT, class,
    None, (x, y, name bytes)); a point whose name fails name_ok comes back with the name None (a reader skips that
    label and keeps the tile); a symbol (kind 7) without a name comes back with the name b"". A named line (kind 8)
    comes back as (LINE, class, None, (name bytes or None, [(x, y), ...]))."""
    allowed = KIND_GEOMETRIES.get(layer_kind)
    if not allowed:
        raise PmtError(f"layer kind {layer_kind} is not one this reader decodes")
    lo, hi = -buffer, (1 << coord_bits) + buffer
    features = []
    pos, end = 0, len(blob)
    px = py = points = 0
    while pos < end:
        if len(features) == MAX_FEATURES:
            raise PmtError("too many features")
        head, pos = get_varint(blob, pos, end)
        geom, has_attr, cls = head & 3, (head >> 2) & 1, head >> 3
        if geom not in allowed:
            raise PmtError(f"geometry type {geom} in layer kind {layer_kind}")
        if geom == POINT:
            if has_attr:
                raise PmtError("a point with an attribute")
            points += 1
            if points > MAX_POINTS:
                raise PmtError("too many points in the tile")
            u, pos = get_varint(blob, pos, end)
            v, pos = get_varint(blob, pos, end)
            px += unzigzag(u)
            py += unzigzag(v)
            if not (lo <= px <= hi and lo <= py <= hi):
                raise PmtError(f"point ({px}, {py}) outside the buffered tile")
            if pos >= end:
                raise PmtError("a point without its name")
            n = blob[pos]
            pos += 1
            if not (0 if layer_kind == KIND_SYMBOLS else 1) <= n <= NAME_MAX or n > end - pos:
                raise PmtError("a name of length 0, over the limit or past the tile's end")
            name = bytes(blob[pos:pos + n])
            pos += n
            features.append((POINT, cls, None, (px, py, name if not n or name_ok(name) else None)))
            continue
        if layer_kind == KIND_STREETS:
            if has_attr:
                raise PmtError("a named line with an attribute")
            if pos >= end:
                raise PmtError("a named line without its name")
            n = blob[pos]
            pos += 1
            if not 1 <= n <= NAME_MAX or n > end - pos:
                raise PmtError("a name of length 0, over the limit or past the tile's end")
            name = bytes(blob[pos:pos + n])
            pos += n
            n, pos = get_varint(blob, pos, end)
            if n < 2:
                raise PmtError("too few points in a named line")
            points += n
            if points > MAX_POINTS:
                raise PmtError("too many points in the tile")
            part = []
            for _ in range(n):
                u, pos = get_varint(blob, pos, end)
                v, pos = get_varint(blob, pos, end)
                px += unzigzag(u)
                py += unzigzag(v)
                if not (lo <= px <= hi and lo <= py <= hi):
                    raise PmtError(f"point ({px}, {py}) outside the buffered tile")
                part.append((px, py))
            features.append((LINE, cls, None, (name if name_ok(name) else None, part)))
            continue
        attr = None
        if has_attr:
            attr, pos = get_varint(blob, pos, end)
            if attr >= attr_count:
                raise PmtError(f"attribute index {attr} out of range")
        nparts, pos = get_varint(blob, pos, end)
        if not 1 <= nparts <= MAX_PARTS:
            raise PmtError("bad part count")
        parts = []
        for _ in range(nparts):
            n, pos = get_varint(blob, pos, end)
            if n < (3 if geom == POLYGON else 2):
                raise PmtError("too few points in a part")
            points += n
            if points > MAX_POINTS:
                raise PmtError("too many points in the tile")
            part = []
            for _ in range(n):
                u, pos = get_varint(blob, pos, end)
                v, pos = get_varint(blob, pos, end)
                px += unzigzag(u)
                py += unzigzag(v)
                if not (lo <= px <= hi and lo <= py <= hi):
                    raise PmtError(f"point ({px}, {py}) outside the buffered tile")
                part.append((px, py))
            parts.append(part)
        features.append((geom, cls, attr, parts))
    return features


def encode_attr(items):
    """items: a list of (tag, value); value is str (UTF-8), bytes, or int (varint)."""
    out = bytearray()
    for tag, value in items:
        if isinstance(value, int):
            body = bytearray()
            put_varint(body, value)
        elif isinstance(value, str):
            body = value.encode("utf-8")
        else:
            body = bytes(value)
        put_varint(out, tag)
        put_varint(out, len(body))
        out += body
    return bytes(out)


def decode_attr(record):
    items, pos, end = [], 0, len(record)
    while pos < end:
        tag, pos = get_varint(record, pos, end)
        n, pos = get_varint(record, pos, end)
        if n > end - pos:
            raise PmtError("attribute item runs past the record")
        items.append((tag, bytes(record[pos:pos + n])))
        pos += n
    return items


# ---- writer --------------------------------------------------------------------------------------------------------

def grid_bytes(cells):
    """cells: {(cx, cy): value} at zoom 8; missing cells are 0."""
    out = bytearray(GRID_BYTES)
    for (cx, cy), v in cells.items():
        if not (0 <= cx < GRID_SIDE and 0 <= cy < GRID_SIDE and 0 <= v <= 3):
            raise PmtError("bad grid cell")
        i = cy * GRID_SIDE + cx
        out[i // 4] |= v << (2 * (i % 4))
    return bytes(out)


def table_bytes(records):
    """A table of `u32 off[n + 1]` relative to the table start, then the records back to back."""
    head = 4 * (len(records) + 1)
    offs, pos = [], head
    for r in records:
        offs.append(pos)
        pos += len(r)
    offs.append(pos)
    return struct.pack(f"<{len(offs)}I", *offs) + b"".join(records)


def build_pmt(*, layer_kind, levels, strings, ids, attrs=(), build_time=0, data_time=0):
    """Returns the bytes of a PMT file.

    levels: a list of dicts with tile_zoom, zoom_min, zoom_max, coord_bits, buffer, tolerance_dm, tiles
    ({(x, y): features}, or the tile's bytes from encode_tile; an empty list makes an empty tile) and optionally
    full_class and grid ({(cx, cy): value}).
    strings: a list of str. ids: {field: string index or None} for the fields in ID_FIELDS.
    attrs: a list of attribute records (bytes, e.g. from encode_attr).
    """
    if not 1 <= len(levels) <= MAX_LEVELS:
        raise PmtError("1 to 8 levels")
    if len(strings) >= NO_STRING:
        raise PmtError("too many strings")
    # Tile blobs, per level, in key order; identical blobs are stored once.
    encoded = []
    for lv in levels:
        tz, bits, buf = lv["tile_zoom"], lv["coord_bits"], lv["buffer"]
        if not (0 <= tz <= MAX_TILE_ZOOM and 8 <= bits <= 16 and 0 <= buf < (1 << (bits - 3))):
            raise PmtError("bad level parameters")
        if not tz <= lv["zoom_min"] <= lv["zoom_max"] <= MAX_DISPLAY_ZOOM:
            raise PmtError("bad level zoom range")
        if lv.get("grid") is not None and tz < 8:
            raise PmtError("a coverage grid needs tile_zoom >= 8")
        tiles = []
        for (x, y), feats in lv["tiles"].items():
            if not (0 <= x < (1 << tz) and 0 <= y < (1 << tz)):
                raise PmtError(f"tile ({x}, {y}) outside zoom {tz}")
            if isinstance(feats, (bytes, bytearray)):  # already encoded (encode_tile with this level's parameters)
                tiles.append((morton(x, y), bytes(feats)))
                continue
            for _, _, attr, _ in feats:
                if attr is not None and attr >= len(attrs):
                    raise PmtError("attribute index out of range")
            tiles.append((morton(x, y), encode_tile(feats, bits, buf, layer_kind)))
        tiles.sort()
        encoded.append(tiles)

    # Header block: fixed header, level table, grids, page tables, attribute table, string table.
    pos = FIXED_HEADER + LEVEL_RECORD * len(levels)
    grid_off, page_off, parts = [], [], []
    for lv in levels:
        if lv.get("grid") is not None:
            grid_off.append(pos)
            parts.append(grid_bytes(lv["grid"]))
            pos += GRID_BYTES
        else:
            grid_off.append(0)
    for tiles in encoded:
        keys = [k for k, _ in tiles[::PAGE_ENTRIES]]
        page_off.append(pos)
        parts.append(struct.pack(f"<{len(keys)}I", *keys))
        pos += 4 * len(keys)
    attr_table = table_bytes(list(attrs)) if attrs else b""
    attr_off = pos if attrs else 0
    parts.append(attr_table)
    pos += len(attr_table)
    str_table = table_bytes([s.encode("utf-8") + b"\0" for s in strings])
    str_off = pos
    parts.append(str_table)
    pos += len(str_table)
    header_size = align(pos)
    if header_size > MAX_HEADER:
        raise PmtError("header block over 16 MiB")

    # Indexes (each on a 4096 boundary), then the data regions.
    pos = header_size
    index_off = []
    for tiles in encoded:
        pos = align(pos)
        index_off.append(pos)
        pos += INDEX_ENTRY * len(tiles)
    data, index, data_off, data_size = [], [], [], []
    for tiles in encoded:
        data_off.append(pos)
        region, seen, entries = bytearray(), {}, []
        for key, blob in tiles:
            if not blob:
                entries.append((key, pos, 0, 0))
                continue
            if blob not in seen:
                seen[blob] = pos + len(region)
                region += blob
            entries.append((key, seen[blob], len(blob), crc32(blob)))
        data.append(bytes(region))
        index.append(b"".join(struct.pack("<4I", *e) for e in entries))
        data_size.append(len(region))
        pos += len(region)
    file_size = pos
    if file_size > MAX_FILE:
        raise PmtError("file over 2 GiB")

    def sid(name):
        v = ids.get(name)
        if v is None:
            return NO_STRING
        if not 0 <= v < len(strings):
            raise PmtError(f"string id for {name} out of range")
        return v

    head = bytearray(struct.pack(
        "<4sHHIIIBBHIIHHHHHHHHIIII", MAGIC, VERSION_MAJOR, kind_minor(layer_kind), header_size, 0, file_size,
        layer_kind, len(levels), 0, build_time, data_time, len(strings), *(sid(n) for n in ID_FIELDS), 0, str_off,
        len(attrs), attr_off, 0))
    assert len(head) == FIXED_HEADER
    for i, (lv, tiles) in enumerate(zip(levels, encoded)):
        keys = [k for k, _ in tiles]
        head += struct.pack(
            "<BBBBHBBIIIIIIIIHHI", lv["tile_zoom"], lv["zoom_min"], lv["zoom_max"], lv["coord_bits"], lv["buffer"],
            LEVEL_FLAG_GRID if grid_off[i] else 0, lv.get("full_class", 0), keys[0] if keys else 0,
            keys[-1] if keys else 0, len(tiles), page_off[i], index_off[i], data_off[i], data_size[i], grid_off[i],
            lv["tolerance_dm"], 0, 0)
    for p in parts:
        head += p
    head += bytes(header_size - len(head))
    struct.pack_into("<I", head, 12, crc32(head))

    out = bytearray(head)
    for i in range(len(levels)):
        out += bytes(index_off[i] - len(out))
        out += index[i]
    for i in range(len(levels)):
        assert len(out) == data_off[i]
        out += data[i]
    assert len(out) == file_size
    return bytes(out)


# ---- reader --------------------------------------------------------------------------------------------------------

class Level:
    pass


class PmtFile:
    """An opened PMT file, validated per FORMAT.md section 8. Tiles are checked when read."""

    def __init__(self, data):
        # A memory map stays one (open_pmt_mapped: a 2 GB file is not copied); anything else becomes bytes.
        self.data = data = data if isinstance(data, mmap.mmap) else bytes(data)
        if len(data) < FIXED_HEADER:
            raise PmtError("shorter than the fixed header")
        (magic, self.major, self.minor, hsize, hcrc, fsize, self.kind, nlev, flags, self.build_time,
         self.data_time, nstr, *rest) = struct.unpack_from("<4sHHIIIBBHIIHHHHHHHHIIII", data, 0)
        ids, (_res, str_off, self.attr_count, attr_off, _res2) = rest[:6], rest[6:]
        if magic != MAGIC or self.major != VERSION_MAJOR:
            raise PmtError("not a PMT version 1 file")
        if self.minor < kind_minor(self.kind):
            raise PmtError(f"layer kind {self.kind} in a version 1.{self.minor} file (it needs "
                           f"1.{kind_minor(self.kind)})")
        if hsize % ALIGN or not FIXED_HEADER <= hsize <= MAX_HEADER or hsize > len(data):
            raise PmtError("bad header size")
        head = bytearray(data[:hsize])
        struct.pack_into("<I", head, 12, 0)
        if crc32(head) != hcrc:
            raise PmtError("header CRC mismatch")
        if fsize != len(data):
            raise PmtError("file size mismatch (truncated or padded copy)")
        if not 1 <= nlev <= MAX_LEVELS or FIXED_HEADER + LEVEL_RECORD * nlev > hsize:
            raise PmtError("bad level count")
        self.header_size = hsize

        def inside_header(off, size):
            return FIXED_HEADER + LEVEL_RECORD * nlev <= off and off + size <= hsize

        self.strings = self._table(data, str_off, nstr, inside_header, "string")
        for s in self.strings:
            if not s.endswith(b"\0"):
                raise PmtError("string without its 0 byte")
        self.strings = [s[:-1] for s in self.strings]
        self.ids = dict(zip(ID_FIELDS, ids))
        for name, v in self.ids.items():
            if v != NO_STRING and v >= nstr:
                raise PmtError(f"string id for {name} out of range")
        if self.attr_count:
            self.attrs = self._table(data, attr_off, self.attr_count, inside_header, "attribute")
        else:
            self.attrs = []

        self.levels = []
        for i in range(nlev):
            lv = Level()
            (lv.tile_zoom, lv.zoom_min, lv.zoom_max, lv.coord_bits, lv.buffer, lv.flags, lv.full_class, lv.key_min,
             lv.key_max, lv.entry_count, lv.page_table_offset, lv.index_offset, lv.data_offset, lv.data_size,
             lv.grid_offset, lv.tolerance_dm, _r1, _r2) = struct.unpack_from(
                "<BBBBHBBIIIIIIIIHHI", data, FIXED_HEADER + LEVEL_RECORD * i)
            if not (lv.tile_zoom <= MAX_TILE_ZOOM and 8 <= lv.coord_bits <= 16
                    and lv.buffer < (1 << (lv.coord_bits - 3))
                    and lv.tile_zoom <= lv.zoom_min <= lv.zoom_max <= MAX_DISPLAY_ZOOM and lv.flags & ~1 == 0):
                raise PmtError(f"level {i}: bad parameters")
            lv.page_count = (lv.entry_count + PAGE_ENTRIES - 1) // PAGE_ENTRIES
            if not inside_header(lv.page_table_offset, 4 * lv.page_count):
                raise PmtError(f"level {i}: page table outside the header block")
            lv.pages = struct.unpack_from(f"<{lv.page_count}I", data, lv.page_table_offset)
            if any(a >= b for a, b in zip(lv.pages, lv.pages[1:])):
                raise PmtError(f"level {i}: page table not ascending")
            if lv.flags & LEVEL_FLAG_GRID:
                if lv.tile_zoom < 8 or not inside_header(lv.grid_offset, GRID_BYTES):
                    raise PmtError(f"level {i}: bad coverage grid")
                lv.grid = data[lv.grid_offset:lv.grid_offset + GRID_BYTES]
            else:
                lv.grid = None
            if (lv.index_offset % ALIGN or lv.index_offset < hsize
                    or lv.index_offset + INDEX_ENTRY * lv.entry_count > fsize):
                raise PmtError(f"level {i}: index outside the file")
            if lv.data_offset < hsize or lv.data_offset + lv.data_size > fsize:
                raise PmtError(f"level {i}: data region outside the file")
            self.levels.append(lv)

    @staticmethod
    def _table(data, off, n, inside_header, what):
        if not inside_header(off, 4 * (n + 1)):
            raise PmtError(f"{what} table outside the header block")
        offs = struct.unpack_from(f"<{n + 1}I", data, off)
        if offs[0] < 4 * (n + 1) or any(a > b for a, b in zip(offs, offs[1:])) or not inside_header(off, offs[-1]):
            raise PmtError(f"bad {what} table")
        return [data[off + a:off + b] for a, b in zip(offs, offs[1:])]

    def string(self, sid):
        return None if sid == NO_STRING else self.strings[sid].decode("utf-8", "replace")

    def coverage(self, level, x, y):
        """The grid value of the zoom-8 cell holding tile (x, y) of the level; 1 without a grid."""
        lv = self.levels[level]
        if lv.grid is None:
            return 1
        s = lv.tile_zoom - 8
        i = (y >> s) * GRID_SIDE + (x >> s)
        v = (lv.grid[i // 4] >> (2 * (i % 4))) & 3
        return 0 if v == 3 else v

    def entries(self, level):
        lv = self.levels[level]
        for i in range(lv.entry_count):
            yield self._entry(lv, i)

    def _entry(self, lv, i):
        key, off, length, crc = struct.unpack_from("<4I", self.data, lv.index_offset + INDEX_ENTRY * i)
        return key, off, length, crc

    def find(self, level, x, y):
        """The index entry (key, offset, length, crc) of tile (x, y), found as a device would; None if absent."""
        lv = self.levels[level]
        key = morton(x, y)
        if lv.entry_count == 0 or not lv.key_min <= key <= lv.key_max:
            return None
        lo, hi = 0, lv.page_count - 1  # last page whose first key <= key
        while lo < hi:
            mid = (lo + hi + 1) // 2
            if lv.pages[mid] <= key:
                lo = mid
            else:
                hi = mid - 1
        first = lo * PAGE_ENTRIES
        lo, hi = first, min(first + PAGE_ENTRIES, lv.entry_count) - 1
        while lo <= hi:
            mid = (lo + hi) // 2
            e = self._entry(lv, mid)
            if e[0] == key:
                return e
            if e[0] < key:
                lo = mid + 1
            else:
                hi = mid - 1
        return None

    def tile_blob(self, level, entry):
        lv = self.levels[level]
        _key, off, length, crc = entry
        if length > MAX_TILE or off < lv.data_offset or off + length > lv.data_offset + lv.data_size:
            raise PmtError("tile outside its data region")
        blob = self.data[off:off + length]
        if crc32(blob) != crc:
            raise PmtError("tile CRC mismatch")
        return blob

    def tile(self, level, entry):
        lv = self.levels[level]
        return decode_tile(self.tile_blob(level, entry), lv.coord_bits, lv.buffer, self.attr_count, self.kind)


def open_pmt(path):
    with open(path, "rb") as f:
        return PmtFile(f.read())


def open_pmt_mapped(path):
    """open_pmt for files too large to read into memory (split_pieces.py): the file is memory-mapped, read-only."""
    with open(path, "rb") as f:
        return PmtFile(mmap.mmap(f.fileno(), 0, access=mmap.ACCESS_READ))


# ---- canonical dump ------------------------------------------------------------------------------------------------

def quote(raw):
    out = bytearray(b'"')
    for c in raw:
        if c == 0x5C:
            out += b"\\\\"
        elif c == 0x22:
            out += b'\\"'
        elif c < 0x20:
            out += b"\\x%02X" % c
        else:
            out.append(c)
    return bytes(out + b'"')


def dump(pf):
    """The canonical dump as bytes (FORMAT.md section 9). Strings are printed byte for byte."""
    lines = []

    def w(text):
        lines.append(text.encode("ascii"))

    w(f"PMTF {pf.major}.{pf.minor} kind={pf.kind} levels={len(pf.levels)} strings={len(pf.strings)} "
      f"attrs={pf.attr_count} build_time={pf.build_time} data_time={pf.data_time}")
    w("ids " + " ".join(f"{n}={pf.ids[n]}" for n in ID_FIELDS))
    for i, s in enumerate(pf.strings):
        lines.append(b"string %d " % i + quote(s))
    for i, a in enumerate(pf.attrs):
        w(f"attr {i} {a.hex()}")
    for i, lv in enumerate(pf.levels):
        w(f"level {i} tz={lv.tile_zoom} zoom={lv.zoom_min}-{lv.zoom_max} bits={lv.coord_bits} buffer={lv.buffer} "
          f"flags={lv.flags} full_class={lv.full_class} keys={lv.key_min}-{lv.key_max} entries={lv.entry_count} "
          f"tol_dm={lv.tolerance_dm}")
        if lv.grid is not None:
            counts = [0, 0, 0, 0]
            for b in lv.grid:
                for k in range(4):
                    counts[(b >> (2 * k)) & 3] += 1
            w(f"grid {i} " + " ".join(f"v{k}={counts[k]}" for k in range(4)))
    for i, lv in enumerate(pf.levels):
        for entry in pf.entries(i):
            key, _off, length, crc = entry
            x, y = unmorton(key)
            w(f"tile {i} {x} {y} len={length} crc={crc:08x}")
            for geom, cls, attr, parts in pf.tile(i, entry):
                if geom == POINT:
                    px, py, name = parts
                    lines.append(f"feature point class={cls} {px},{py} name=".encode("ascii")
                                 + (quote(name) if name is not None else b"-"))
                    continue
                if pf.kind == KIND_STREETS:
                    name, part = parts
                    lines.append(f"feature nline class={cls} name=".encode("ascii")
                                 + (quote(name) if name is not None else b"-"))
                    w(f"part {len(part)} " + " ".join(f"{x},{y}" for x, y in part))
                    continue
                w(f"feature {'polygon' if geom == POLYGON else 'line'} class={cls} "
                  f"attr={'-' if attr is None else attr} parts={len(parts)}")
                for part in parts:
                    w(f"part {len(part)} " + " ".join(f"{x},{y}" for x, y in part))
    return b"\n".join(lines) + b"\n"


def main(argv):
    if len(argv) != 3 or argv[1] not in ("dump", "check"):
        print(__doc__, file=sys.stderr)
        return 2
    try:
        pf = open_pmt(argv[2])
        if argv[1] == "dump":
            sys.stdout.buffer.write(dump(pf))
        else:
            n = 0
            for i in range(len(pf.levels)):
                for entry in pf.entries(i):
                    pf.tile(i, entry)
                    n += 1
            print(f"ok: {len(pf.levels)} levels, {n} tiles")
    except PmtError as e:
        print(f"error: {e}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
