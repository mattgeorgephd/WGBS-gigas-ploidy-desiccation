"""Write the column expression_paper_recipe of manuscript_values_map.csv.

Decision of the author, 2026-11-20 (07_documentation/reproducibility_plan.md, section 7,
table "Decided"): for every optional analysis, step 13 evaluates the paper's own
computation (expression_legacy) on that analysis's tables, so that a difference
from the paper's value comes only from the option. The corrected computation
(expression_optional) is kept and reported in a column of its own.

expression_paper_recipe is expression_legacy written as a template: A stands for
the analysis id, as in expression_optional. With A = "legacy" every template gives
the paper's value (tests/check_paper_recipe_legacy.R). The column is derived, never
edited by hand: after a change to expression_legacy, run from the repository root
    python3 -I 07_documentation/toggles_migration/paper_recipe_column.py
and review the diff; with --check it only reports whether the stored column is
what the rules give (exit status 1 if not).

Rules, applied in this order:
  1. Rows treated by hand (HAND below, with the legacy expression each was written
     for: the script stops when that expression has changed, so that the template
     is reviewed again).
  2. Terms named in Results 3.4 (rows that read in_legacy_p_list of the paper's
     named_terms_results_3_4.csv): the GO IDs named in the text are the paper's
     (paper_tab); whether each is in the p < 0.01 list of the draft's test
     (below_p_cutoff, step 09) is read from the analysis's own Fisher tables.
     Step 11 builds in_legacy_p_list the same way for the paper.
  3. The per-CpG table of step 05: tab("05_differential_methylation",
     "united_cpg_legacy.csv.gz")[rows] reads, for an optional analysis, its rows
     at the DML, united_at_dml_<A>.csv, which is tracked where the per-CpG tables
     are not (as expression_optional does); for A = "legacy" it reads the paper's
     table, which has no such file.
  4. A table of the analysis named "<name>_legacy.<ext>" becomes
     paste0("<name>_", A, ".<ext>"), and the filter mode == "legacy" becomes
     mode == A. Tables that do not depend on the analysis (steps 01 to 03 and 12,
     library_coverage.csv, gene_go_map.csv.gz, the files of step 09 named
     legacy_<...>, which an optional analysis writes under the same name) and
     columns named legacy_<...> are not touched.
  Rows that read a table that only the paper run makes (the "as drawn" clusters and
  the REVIGO exports of step 11) are left as they are: in an optional run the table
  is absent from the analysis folder and step 13 reports "not run".

Written with the csv module (QUOTE_MINIMAL, CRLF), like rewrite_maps.py. Before
writing, the file as read must give back its own bytes (round-trip check), and
every column other than expression_paper_recipe must be unchanged afterwards.
"""
import csv, io, os, re, sys

PATH = os.path.relpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', '03_input_files', 'manuscript',
                                    'manuscript_values_map.csv'))
COL = 'expression_paper_recipe'

def read(path):
    raw = open(path, 'rb').read()
    return raw, list(csv.reader(io.StringIO(raw.decode('utf-8'), newline='')))

def write_bytes(rows):
    buf = io.StringIO(newline='')
    csv.writer(buf, quoting=csv.QUOTE_MINIMAL, lineterminator='\r\n').writerows(rows)
    return buf.getvalue().encode('utf-8')

# ---- 1. Rows treated by hand ------------------------------------------------------------------
# Table 1 group means (06_diagnostics/table1_group_means.R) for the Table 1 of the analysis:
# the read-pooled percent methylation of the four groups at each locus of table_1_<A>.csv,
# from the counts that step 05 tested, and the differences the diagnostic writes (rounded to
# two decimals, as there). Replaces diagnostic("table1_group_means.csv"), which is made for
# the paper's Table 1 only.
T1_MEANS = ('local({t <- tab("09_go_enrichment_legacy", paste0("table_1_", A, ".csv")); '
            'g <- tab("02_alignment_qc", "bismark_by_library.csv")[order(legacy_sample_index)]; '
            'u <- (if (A == "legacy") tab("05_differential_methylation", "united_cpg_legacy.csv.gz") else '
            'tab("05_differential_methylation", paste0("united_at_dml_", A, ".csv")))[match(paste(t$chromosome, t$start), paste(chr, start))]; '
            'stopifnot(!anyNA(u$chr)); '
            'm <- function(pl, tr) {k <- which(g$ploidy == pl & g$treatment == tr); cs <- rowSums(u[, paste0("numCs", k), with = FALSE]); '
            '100 * cs / (cs + rowSums(u[, paste0("numTs", k), with = FALSE]))}; '
            'data.table(dml_set = t$dml_set, gene_id = t$gene_id, start = t$start, table1_value = t$methylation_difference, '
            'triploid_minus_diploid_control_cohort = round(m("triploid", "control") - m("diploid", "control"), 2), '
            'triploid_minus_diploid_heat_cohort = round(m("triploid", "heat") - m("diploid", "heat"), 2), '
            'heat_minus_control_diploids = round(m("diploid", "heat") - m("diploid", "control"), 2), '
            'heat_minus_control_triploids = round(m("triploid", "heat") - m("triploid", "control"), 2))})')
T1_ROWS = ['R3.5-029', 'R3.5-031', 'R3.5-033', 'R3.5-041', 'R3.5-043', 'R3.5-046', 'R3.5-047', 'T1-003',
           'D4.3-016', 'D4.3-020', 'D4.9-008', 'D4.9-010', 'D4.9-012']

def fisher(set_expr):
    """The analysis's Fisher tables of step 09 (draft's test), the three ontologies bound."""
    return ('rbindlist(lapply(c("BP", "CC", "MF"), function(o) tab("09_go_enrichment_legacy", '
            'sprintf("legacy_fisher_%%s_%%s.csv", %s, o))))' % set_expr)

HAND = {
    # term_lists_compared.csv of step 11 counts the terms with p below the cutoff in the paper's Fisher tables.
    'D4.3-006': ('tab("11_semantic_similarity", "term_lists_compared.csv")[set == "shared", sum(terms_legacy_p_below_cutoff) > 0]',
                 'sum(%s$below_p_cutoff) > 0' % fisher('"shared"')),
    'D4.8-003': ('tab("11_semantic_similarity", "term_lists_compared.csv")[, sum(terms_legacy_p_below_cutoff) > 0]',
                 'sum(vapply(c("ploidy", "heat", "shared"), function(s) sum(%s$below_p_cutoff), numeric(1))) > 0' % fisher('s')),
    # The two AK9 loci: NA (not comparable) when they are not both among the analysis's loci in both lists,
    # as in expression_optional; the paper's expression stops there instead.
    'D4.9-014': ('{g <- tab("02_alignment_qc", "bismark_by_library.csv")[order(legacy_sample_index)]; s <- tab("05_differential_methylation", "dml_shared_legacy.csv")[chr == "NC_047563.1" & start %in% c(5569229, 5569232)]; u <- tab("05_differential_methylation", "united_cpg_legacy.csv.gz")[s$united_row]; stopifnot(u$chr == s$chr, u$start == s$start); pct <- function(p, t) {k <- which(g$ploidy == p & g$treatment == t); cs <- rowSums(u[, paste0("numCs", k), with = FALSE]); 100 * cs / (cs + rowSums(u[, paste0("numTs", k), with = FALSE]))}; stopifnot(nrow(s) == 2); hd <- pct("diploid", "heat") - pct("diploid", "control"); ht <- pct("triploid", "heat") - pct("triploid", "control"); all(abs(hd - ht) > config$dml_difference)}',
                 '{g <- tab("02_alignment_qc", "bismark_by_library.csv")[order(legacy_sample_index)]; s <- tab("05_differential_methylation", paste0("dml_shared_", A, ".csv"))[chr == "NC_047563.1" & start %in% c(5569229, 5569232)]; if (nrow(s) != 2) NA else {u <- (if (A == "legacy") tab("05_differential_methylation", "united_cpg_legacy.csv.gz")[, united_row := .I] else tab("05_differential_methylation", paste0("united_at_dml_", A, ".csv")))[match(s$united_row, united_row)]; stopifnot(u$chr == s$chr, u$start == s$start); pct <- function(p, t) {k <- which(g$ploidy == p & g$treatment == t); cs <- rowSums(u[, paste0("numCs", k), with = FALSE]); 100 * cs / (cs + rowSums(u[, paste0("numTs", k), with = FALSE]))}; hd <- pct("diploid", "heat") - pct("diploid", "control"); ht <- pct("triploid", "heat") - pct("triploid", "control"); all(abs(hd - ht) > config$dml_difference)}}'),
}
# Rows left as they are because they read a table that only the paper run makes (step 11, "as drawn"
# clusters and REVIGO exports); in an optional run step 13 reports them "not run".
PAPER_ONLY = ['R3.4-007', 'R3.4-008', 'R3.4-009', 'R3.4-013', 'R3.4-024', 'R3.4-035']
NAMED_ROWS = ['R3.4-010', 'R3.4-011', 'R3.4-012', 'R3.4-014', 'R3.4-015', 'R3.4-016', 'R3.4-018', 'R3.4-019',
              'R3.4-020', 'R3.4-021', 'R3.4-022', 'R3.4-023', 'R3.4-025', 'R3.4-026', 'R3.4-027', 'R3.4-030',
              'R3.4-032', 'R3.4-033', 'R3.4-034', 'R3.4-036', 'R3.4-037', 'R3.4-038', 'R3.4-043', 'R3.4-044']
UNITED_ROWS = ['R3.2-029', 'R3.2-030', 'R3.2-031', 'D4.2-004', 'D4.3-005', 'D3.5x-006']

# ---- 2. Terms named in Results 3.4 ---------------------------------------------------------
named_re = re.compile(r'tab\("11_semantic_similarity", "named_terms_results_3_4\.csv"\)\[(set == "(ploidy|heat|shared)"[^\]]*), '
                      r'(in_legacy_p_list|any\(in_legacy_p_list\)|all\(in_legacy_p_list\))\]')
def named_rule(e):
    m = named_re.fullmatch(e)
    if not m:
        return None
    filt, s, j = m.group(1), m.group(2), m.group(3)
    agg = {'in_legacy_p_list': '%s', 'any(in_legacy_p_list)': 'any(%s)', 'all(in_legacy_p_list)': 'all(%s)'}[j]
    return ('local({n <- paper_tab("11_semantic_similarity", "named_terms_results_3_4.csv")[%s, go_id]; '
            'p <- %s[below_p_cutoff == TRUE, go_id]; %s})' % (filt, fisher('"%s"' % s), agg % 'n %in% p'))

# ---- 3. The per-CpG table of step 05 -------------------------------------------------------
UNITED = 'tab("05_differential_methylation", "united_cpg_legacy.csv.gz")['
UNITED_NEW = ('(if (A == "legacy") tab("05_differential_methylation", "united_cpg_legacy.csv.gz")[, united_row := .I] else '
              'tab("05_differential_methylation", paste0("united_at_dml_", A, ".csv")))[match(%s, united_row)]')
PLACEHOLDER = '\x00UNITED\x00'
def united_rule(e):
    out = ''
    while UNITED in e:
        i = e.index(UNITED); j = i + len(UNITED); depth = 1; k = j
        while depth:
            depth += {'[': 1, ']': -1}.get(e[k], 0); k += 1
        out += e[:i] + UNITED_NEW % e[j:k - 1]
        e = e[k:]
    return out + e

# ---- 4. Tables of the analysis and the mode filter -----------------------------------------
suffix_re = re.compile(r'"([A-Za-z0-9_]+)_legacy(\.[A-Za-z.]+)"')
def mechanical(e):
    e = e.replace('"united_cpg_legacy.csv.gz"', PLACEHOLDER)
    e = suffix_re.sub(lambda m: 'paste0("%s_", A, "%s")' % (m.group(1), m.group(2)), e)
    e = e.replace('mode == "legacy"', 'mode == A')
    return e.replace(PLACEHOLDER, '"united_cpg_legacy.csv.gz"')

def template(rid, e):
    if rid in HAND:
        want, t = HAND[rid]
        if e != want:
            sys.exit('%s: expression_legacy has changed since its template was written by hand; review HAND in %s' % (rid, sys.argv[0]))
        return t, 'hand'
    if rid in T1_ROWS:
        if e.count('diagnostic("table1_group_means.csv")') != 1:
            sys.exit('%s: expected one diagnostic("table1_group_means.csv")' % rid)
        return mechanical(e.replace('diagnostic("table1_group_means.csv")', T1_MEANS)), 'table1'
    if rid in PAPER_ONLY:
        return e, 'paper_only'
    n = named_rule(e)
    if n is not None:
        return n, 'named'
    if UNITED in e:
        return mechanical(united_rule(e)), 'united'
    return mechanical(e), 'mechanical'

def main():
    check = '--check' in sys.argv[1:]
    raw, rows = read(PATH)
    if write_bytes(rows) != raw:
        sys.exit('round-trip check failed: the csv writer does not give the bytes of ' + PATH)
    header = rows[0]
    L = header.index('expression_legacy')
    if COL in header:
        k = header.index(COL)
    else:
        if check:
            sys.exit('the map has no column ' + COL)
        k = L + 1
        for r in rows:
            r.insert(k, COL if r is header else '')
    kinds = {}
    new = [list(r) for r in rows]
    for r in new[1:]:
        e = r[L]
        if not e:
            r[k] = ''
            continue
        t, kind = template(r[0], e)
        kinds.setdefault(kind, []).append(r[0])
        # Nothing of the paper's analysis may be left, except in the branch for A == "legacy" and the rows read as they are.
        rest = t.replace('(if (A == "legacy") ', '').replace('"united_cpg_legacy.csv.gz"', '')
        if kind != 'paper_only' and (suffix_re.search(rest) or re.search(r'mode\s*==\s*"legacy"', rest)):
            sys.exit('left over in %s: %s' % (r[0], t))
        r[k] = t
    # The rows of each rule are the rows they were written for.
    for kind, want in (('table1', T1_ROWS), ('paper_only', PAPER_ONLY), ('named', NAMED_ROWS), ('united', UNITED_ROWS),
                       ('hand', sorted(HAND))):
        if sorted(kinds.get(kind, [])) != sorted(want):
            sys.exit('rule %s applies to %s, expected %s' % (kind, sorted(kinds.get(kind, [])), sorted(want)))
    out = write_bytes(new)
    for i, (a, b) in enumerate(zip(rows, new)):
        if [x for j, x in enumerate(a) if j != k] != [x for j, x in enumerate(b) if j != k]:
            sys.exit('row %d: a column other than %s would change' % (i, COL))
    counts = ', '.join('%s %d' % (x, len(kinds.get(x, []))) for x in ('mechanical', 'united', 'named', 'table1', 'hand', 'paper_only'))
    if check:
        if out != raw:
            sys.exit('%s: the stored %s is not what the rules give; run this script without --check' % (PATH, COL))
        print('%s is up to date (%d templates: %s)' % (COL, sum(map(len, kinds.values())), counts))
        return
    open(PATH, 'wb').write(out)
    raw2, rows2 = read(PATH)
    if write_bytes(rows2) != raw2 or raw2 != out:
        sys.exit('written file does not read back')
    print('%s written: %d templates (%s); %d rows, CRLF %d' % (COL, sum(map(len, kinds.values())), counts, len(rows2) - 1, raw2.count(b'\r\n')))

if __name__ == '__main__':
    main()
