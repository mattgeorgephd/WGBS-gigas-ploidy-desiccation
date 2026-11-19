import sys, zlib, struct, bz2, lzma
def head(path):
    raw = open(path, 'rb').read()
    magic = raw[:6]
    data = None
    if raw[:2] == b'\x1f\x8b':
        d = zlib.decompressobj(31); data = d.decompress(raw); comp = 'gzip'
    elif raw[:3] == b'BZh':
        d = bz2.BZ2Decompressor(); data = d.decompress(raw); comp = 'bzip2'
    elif raw[:6] == b'\xfd7zXZ\x00':
        d = lzma.LZMADecompressor(); data = d.decompress(raw); comp = 'xz'
    else:
        data = raw; comp = 'none'
    fmt = data[:5]
    kind = data[5:7]
    ver, writer, minr = struct.unpack('>iii', data[7:19])
    def v(x): return '%d.%d.%d' % (x >> 16, (x >> 8) & 255, x & 255)
    enc = ''
    if ver >= 3:
        n = struct.unpack('>i', data[19:23])[0]
        enc = data[23:23+n].decode('ascii', 'replace')
    return dict(compression=comp, magic=fmt, kind=kind, serialization_version=ver, written_by_R=v(writer), min_reader=v(minr), native_encoding=enc)
for p in sys.argv[1:]:
    try:
        print(p.split('/')[-1], head(p))
    except Exception as e:
        print(p.split('/')[-1], 'ERROR', repr(e))
