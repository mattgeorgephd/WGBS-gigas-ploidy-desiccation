"""One-time rewrite of the two manuscript maps for the settings layout of 2026-11-20.

07_documentation/agent_workflows/scope_inputs/settings_build_spec.md, section 6.13.
Run once from the repository root:  python3 -I 07_documentation/toggles_migration/rewrite_maps.py

manuscript_values_map.csv
  - column expression_corrected is renamed expression_optional;
  - in its expressions, a table of the analysis named tab("<step>", "<name>_corrected.<ext>")
    becomes tab("<step>", paste0("<name>_", A, ".<ext>")) and the literal "corrected"
    (always in mode == "corrected") becomes A, the analysis id that step 13 sets
    for each optional analysis; the two rows that read the paper's step 09
    tables (D4.3-008, D4.7-005) read them with paper_tab(). With A = "corrected"
    every expression gives the old value.
    Tables of 06_diagnostics/ are not touched.
  - every other column is unchanged.
manuscript_items_map.csv
  - items of the corrected mode and of step 10 become templates (analysis "optional",
    {A} in source and destination, under 04_output/not_adopted/{A}/ and
    05_manuscript/not_adopted/{A}/), or rows of the comparison folder
    (analysis "comparisons"); the bias figure of step 10 had a row per mode and
    gets one template.

Both files keep CRLF line endings. Written with the csv module (QUOTE_MINIMAL,
lineterminator CRLF); before writing, the unchanged rows are written to memory
and must give the original bytes (round-trip check).
"""
import csv, io, re, sys

def read(path):
    raw = open(path, 'rb').read()
    rows = list(csv.reader(io.StringIO(raw.decode('utf-8'), newline='')))
    return raw, rows

def write_bytes(rows):
    buf = io.StringIO(newline='')
    w = csv.writer(buf, quoting=csv.QUOTE_MINIMAL, lineterminator='\r\n')
    w.writerows(rows)
    return buf.getvalue().encode('utf-8')

def roundtrip(path, raw, rows):
    if write_bytes(rows) != raw:
        sys.exit('round-trip check failed for ' + path + ': the csv writer does not give the original bytes')

# ---- values map ----------------------------------------------------------------
vpath = '03_input_files/manuscript/manuscript_values_map.csv'
raw, rows = read(vpath)
roundtrip(vpath, raw, rows)
header = rows[0]
if 'expression_corrected' not in header:
    sys.exit('the values map has no column expression_corrected (already rewritten?)')
k = header.index('expression_corrected')
header[k] = 'expression_optional'
tab_re = re.compile(r'tab\("([^"]+)",\s*"([^"]*)_corrected(\.[A-Za-z.]+)"\)')
n_files = n_lit = n_rows = 0
for r in rows[1:]:
    e = r[k]
    if not e:
        continue
    n_rows += 1
    e, nf = tab_re.subn(lambda m: 'tab("%s", paste0("%s_", A, "%s"))' % (m.group(1), m.group(2), m.group(3)), e)
    # file names built with sprintf: sprintf("genes_with_dml_%s_corrected.csv", s)
    e, ns = re.subn(r'sprintf\("([^"]*)_corrected(\.[A-Za-z.]+)",\s*([A-Za-z_.]+)\)', r'sprintf("\1_%s\2", \3, A)', e)
    nf += ns
    # the column of the analysis in legacy_top10_at_gene_level.csv of step 10 (q_wallenius_eligible_<id>_dml)
    if 'q_wallenius_eligible_corrected_dml' in e:
        e = e.replace('[, all(!is.na(q_wallenius_eligible_corrected_dml) & q_wallenius_eligible_corrected_dml < config$go_fdr)]',
                      '[, {q <- get(if (A == "legacy") "q_wallenius_eligible" else paste0("q_wallenius_eligible_", A, "_dml")); all(!is.na(q) & q < config$go_fdr)}]')
        e = e.replace('[ontology == "BP", all(!is.na(q_wallenius_eligible_corrected_dml) & q_wallenius_eligible_corrected_dml < config$go_fdr)]',
                      '[ontology == "BP", {q <- get(if (A == "legacy") "q_wallenius_eligible" else paste0("q_wallenius_eligible_", A, "_dml")); all(!is.na(q) & q < config$go_fdr)}]')
        nf += 1
    # tables of the paper's step 09 (the draft's procedure on the paper's DML), read as such in every analysis
    if r[0] in ('D4.3-008', 'D4.7-005'):
        e = e.replace('tab("09_go_enrichment_legacy"', 'paper_tab("09_go_enrichment_legacy"')
    nl = e.count('"corrected"')
    e = e.replace('"corrected"', 'A')
    n_files += nf; n_lit += nl
    if '_corrected' in re.sub(r'diagnostic\(\s*"[^"]*"\)', '', e) or '"corrected"' in e:
        sys.exit('left over in ' + r[0] + ': ' + e)
    r[k] = e
open(vpath, 'wb').write(write_bytes(rows))
raw2, rows2 = read(vpath)
roundtrip(vpath, raw2, rows2)
old = list(csv.reader(io.StringIO(raw.decode('utf-8'), newline='')))
for i, (a, b) in enumerate(zip(old, rows2)):
    for j in range(len(a)):
        if j != k and a[j] != b[j] and i > 0:
            sys.exit('a column other than the expression changed in row %d' % i)
print('values map: %d expressions rewritten; %d table names, %d literals; %d rows, CRLF %d' %
      (n_rows, n_files, n_lit, len(rows2) - 1, raw2.count(b'\r\n')))

# ---- items map -----------------------------------------------------------------
ipath = '03_input_files/manuscript/manuscript_items_map.csv'
raw, rows = read(ipath)
roundtrip(ipath, raw, rows)
h = rows[0]
ia, isrc, idst = h.index('analysis'), h.index('source'), h.index('destination')
out = [h]; seen = set()
s10 = '04_output/10_go_enrichment_gene_level/'
for r in rows[1:]:
    r = list(r)
    src, dst, an = r[isrc], r[idst], r[ia]
    if src.startswith(s10 + 'figures/figure_go_bias_'):
        ext = src.rsplit('.', 1)[1]
        r[ia] = 'optional'
        r[isrc] = '04_output/not_adopted/{A}/10_go_enrichment_gene_level/figures/figure_go_bias_{A}.' + ext
        r[idst] = 'not_adopted/{A}/figures/GO_bias_by_CpG_count.' + ext
    elif src.startswith(s10):
        r[ia] = 'comparisons'
        r[isrc] = src.replace('04_output/', '04_output/not_adopted/comparisons/', 1)
        name = dst.split('/')[-1]
        sub = 'figures' if dst.startswith('figures/') else 'supplement'
        r[idst] = 'not_adopted/comparisons/%s/%s' % (sub, name)
    elif an == 'corrected':
        r[ia] = 'optional'
        step_part = src[len('04_output/'):]
        r[isrc] = '04_output/not_adopted/{A}/' + step_part.replace('_corrected.', '_{A}.')
        r[idst] = 'not_adopted/{A}/' + dst.replace('corrected/', '', 1)
    if r[idst] in seen:
        continue                                   # the bias figure: one template for both modes
    seen.add(r[idst]); out.append(r)
open(ipath, 'wb').write(write_bytes(out))
raw2, rows2 = read(ipath)
roundtrip(ipath, raw2, rows2)
from collections import Counter
print('items map: %d rows (was %d); by analysis: %s; CRLF %d' %
      (len(rows2) - 1, len(rows) - 1, dict(Counter(r[ia] for r in rows2[1:])), raw2.count(b'\r\n')))
