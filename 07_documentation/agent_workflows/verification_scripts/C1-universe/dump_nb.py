import json, sys
path = sys.argv[1]
out = sys.argv[2]
nb = json.load(open(path))
with open(out, 'w') as fh:
    for i, c in enumerate(nb['cells']):
        fh.write("=" * 100 + "\n")
        fh.write("CELL %d  type=%s  exec_count=%s\n" % (i, c['cell_type'], c.get('execution_count')))
        fh.write("-" * 100 + "\n")
        src = c['source']
        if isinstance(src, list):
            src = "".join(src)
        fh.write(src + "\n")
        for o in c.get('outputs', []):
            fh.write("-" * 40 + " OUTPUT type=%s name=%s\n" % (o.get('output_type'), o.get('name')))
            if 'text' in o:
                t = o['text']
                if isinstance(t, list):
                    t = "".join(t)
                fh.write(t + "\n")
            elif 'data' in o:
                for k, v in o['data'].items():
                    if isinstance(v, list):
                        v = "".join(v)
                    fh.write("[%s]\n%s\n" % (k, v if k.startswith('text') else '<%d chars>' % len(v)))
            elif 'ename' in o:
                fh.write("ERROR %s: %s\n" % (o.get('ename'), o.get('evalue')))
                fh.write("\n".join(o.get('traceback', [])) + "\n")
