#!/usr/bin/env python3
"""Write one text file of findings per group of the code pass, from the code reviews (R1 to R6),
their verifications, the check of library_methylation_level.R (S3, when it has reported) and the
list of extra items kept by the session."""
import json, os, re
REPO = '/home/user/WGBS-gigas-ploidy-desiccation/'
S = '/tmp/claude-0/-home-user-WGBS-gigas-ploidy-desiccation/ed665b6a-3457-5aca-a521-7838dbc31cfa/scratchpad'
OUT = S + '/wf/codepass'
os.makedirs(OUT, exist_ok=True)
d = json.load(open(REPO + '07_documentation/agent_workflows/results/audit-verification.json'))
by = {a['agent']: a['result'] for a in d['results']}
GROUP_FILES = {
 'A': ['04_cpg_coverage', '05_differential_methylation', 'methylkit_helpers', 'sample_sheet', 'dml_sensitivity', 'table1_group_means'],
 'B': ['06_dml_figures', '07_genome_features', 'genomic_overlaps', 'genome_features/README', 'feature_odds_matched', 'plot_helpers'],
 'C': ['08_gene_annotation', '09_go_enrichment_legacy', '02_R_functions/annotation', 'annotation/README', 'legacy_reader_mechanism', 'universe_file_cut'],
 'D': ['10_go_enrichment_gene_level', 'go_tests', 'go_random_cpg_control'],
 'E': ['run_pipeline', 'run_config', 'logging.R', 'packages.R', 'paths.R', 'legacy_go_db', '01_fetch_inputs', '11_semantic_similarity', '12_qpcr', '13_manuscript_values', '.gitattributes', 'write_table', 'load_packages'],
 'F': ['library_methylation_level'],
}
texts = {g: [] for g in GROUP_FILES}
def clean(s): return str(s).replace(REPO, '')
n_total = 0
for r in ['R1-coverage-dml', 'R2-figures-features', 'R3-annotation-legacy-go', 'R4-gene-level-go', 'R5-clean-checkout', 'R6-documentation']:
    rev = by['review:' + r]['findings']; ver = by['verify:' + r]['verdicts']
    assert len(rev) == len(ver), (r, len(rev), len(ver))
    for i, (f, v) in enumerate(zip(rev, ver), 1):
        fid = f"{r.split('-')[0]}.{i}"
        blob = ' '.join(clean(x) for x in [f.get('file'), f.get('summary'), f.get('suggested_fix'), v.get('file'), v.get('summary'), v.get('fix')])
        groups = [g for g, keys in GROUP_FILES.items() if any(k in blob for k in keys)]
        entry = [f"{'=' * 100}", f"FINDING {fid}   (reviewer severity {f.get('severity')}; second agent: {v.get('verdict')}, severity {v.get('severity')})",
                 f"FILE (reviewer): {clean(f.get('file'))}:{f.get('line')}    FILE (second agent, current tree): {clean(v.get('file'))}:{v.get('line')}",
                 f"REVIEWER SUMMARY: {clean(f.get('summary'))}", f"REVIEWER FAILURE SCENARIO: {clean(f.get('failure_scenario'))}",
                 f"REVIEWER EVIDENCE: {clean(f.get('evidence'))}", f"REVIEWER SUGGESTED FIX: {clean(f.get('suggested_fix'))}",
                 f"SECOND AGENT SUMMARY: {clean(v.get('summary'))}", f"SECOND AGENT REASON: {clean(v.get('reason'))}", f"SECOND AGENT TESTED FIX: {clean(v.get('fix'))}", '']
        for g in groups: texts[g] += entry
        n_total += 1
        if not groups: texts.setdefault('none', []).extend([fid + ' | ' + clean(f.get('file')) + ' | ' + clean(v.get('verdict')) + ' | ' + clean(v.get('summary'))[:160]])
# extra items kept by the session
extra = open(OUT + '/extra_items.txt').read()
paras = [p.strip() for p in re.split(r'\n\s*\n', extra) if p.strip()]
for p in paras:
    for g in GROUP_FILES:
        if re.search(r'\[group %s[,\]]' % g, p) or p.startswith('[all groups]') or p.startswith('OPTIONAL OUTPUT ADDITIONS'):
            pass
for g in GROUP_FILES:
    mine = []
    for p in paras:
        lines = p.split('\n')
        keep = [l for l in lines if re.search(r'\[group %s[,\]]' % g, l) or l.startswith('[all groups]')]
        if p.startswith('[group %s' % g): mine.append(p)
        elif keep and not p.startswith('[group '): mine.append('\n'.join(([lines[0]] if not lines[0].startswith('[') else []) + keep))
    if mine:
        texts[g] += ['=' * 100, 'EXTRA ITEMS (from the checks of the issue log and from the traceability table)', ''] + ['\n' + m for m in mine]
# S3: code review of library_methylation_level.R
f3 = REPO + '07_documentation/agent_workflows/results/check-sensitivity-text.json'
s3 = [a['result'] for a in json.load(open(f3))['results'] if a['agent'] == 'check:S3-code-review-library-level']
if s3:
    r = s3[0]
    texts['F'] += ['=' * 100, 'CODE REVIEW OF 06_diagnostics/library_methylation_level.R BY AN INDEPENDENT AGENT (check S3)', '']
    for i, x in enumerate(r.get('code_findings') or [], 1):
        texts['F'] += [f"S3.CF{i} [{x['severity']}] {clean(x['file'])}:{x['line']} changes a committed number: {x['changes_a_committed_number']}", f"PROBLEM: {clean(x['problem'])}", f"FAILURE SCENARIO: {clean(x['failure_scenario'])}", f"PROPOSED FIX: {clean(x['proposed_fix'])}", '']
    for i, x in enumerate(r['discrepancies'], 1):
        texts['F'] += [f"S3.D{i} [{x['severity']}] {clean(x['where'])}", f"QUOTE: {x['quote']}", f"PROBLEM: {clean(x['problem'])}", f"EVIDENCE: {clean(x['evidence'])}", f"PROPOSED: {x['proposed_text']}", '']
    for i, x in enumerate(r['additions'], 1):
        texts['F'] += [f"S3.A{i} [{x['severity']}] {clean(x['item'])}", f"WHY: {clean(x['why_it_matters'])}", f"WHERE: {clean(x['where_to_add'])}", f"PROPOSED: {x['proposed_text']}", '']
    texts['F'] += ['S3 NOTES: ' + clean(r['notes']), '']
else:
    texts['F'] += ['=' * 100, 'The independent code review of the script (check S3) had not reported when this file was written.', '']
HEAD = ("Findings for group {g} of the code pass. Each finding gives the reviewer's text and the verdict of a second agent who re-examined it,\n"
        "with the fix that agent tested. The reviews are some hours old: check the current file before editing. Documentation statements in\n"
        "issues.md, README.md, the plan and other README files were already corrected by the parent session (commits b8e3ddf, dd648ca);\n"
        "where a finding has a documentation part and a code part, the code part is yours.\n\n")
for g in GROUP_FILES:
    open(f'{OUT}/group_{g}.txt', 'w').write(HEAD.format(g=g) + '\n'.join(texts[g]) + '\n')
    print(g, len(texts[g]), 'lines,', sum(1 for l in texts[g] if l.startswith('FINDING ')), 'findings,', os.path.getsize(f'{OUT}/group_{g}.txt'), 'bytes')
print('findings in the reviews:', n_total, '| assigned to no group (documentation only):', len(texts.get('none', [])))
for l in texts.get('none', []): print('   ', l)
