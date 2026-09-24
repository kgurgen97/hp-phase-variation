#!/usr/bin/env python3
"""G5 main-process merge: reference panel + catalogue (workstream B) + literature (workstream A) -> metadata/g5/*.
KNOWN_PV assignment: literature row (gene/locus tag AND repeat motif class AND stated context) matched to catalogue members;
applied to every member of the same locus_id (ortholog group). Disease-blind."""
import csv, gzip, re, os, collections, shutil, json
W = 'metadata/workstreams/'
rd = lambda f: list(csv.DictReader(open(f), delimiter='\t'))
cat = rd(W + 'g5_reference_catalogue/g5b_repeat_catalogue.tsv')
panel = rd(W + 'g5_reference_catalogue/g5b_reference_panel.tsv')
orth = rd(W + 'g5_reference_catalogue/g5b_locus_orthology.tsv')
lit = rd(W + 'g5_literature_and_known_pv/g5a_known_pv_loci.tsv')
# old_locus_tag -> RefSeq locus_tag per reference (from the GFFs)
old2new = {}
pth = {r['short_code']: r for r in panel}
for short, p in pth.items():
    acc = p['assembly_accession_version']; d = f'data/reference/g5_panel/{acc}/'
    for fn in os.listdir(d):
        if fn.endswith('.gff.gz'):
            for l in gzip.open(d + fn, 'rt'):
                if l.startswith('#'): continue
                c = l.rstrip('\n').split('\t')
                if len(c) < 9 or c[2] != 'gene': continue
                a = dict(x.split('=', 1) for x in c[8].split(';') if '=' in x)
                if 'old_locus_tag' in a and 'locus_tag' in a:
                    for ol in a['old_locus_tag'].replace('%2C', ',').split(','): old2new[(short, ol)] = a['locus_tag']
# literature -> match rules. Each key = literature gene_name (exact string in g5a file) ; entries: refs/tags, gene names, motif canonicals, contexts
MATCH = [
 ('futA/futB', dict(tags=[('26695', 'HP0379'), ('26695', 'HP0651')], genes=[], motifs={'C'}, contexts={'CDS'})),
 ('sabA', dict(tags=[('26695', 'HP0725')], genes=['sabA'], motifs={'AG', 'A'}, contexts={'CDS', 'PROMOTER_OR_INTERGENIC_UPSTREAM'})),
 ('sabB', dict(tags=[], genes=['sabB'], motifs={'AG'}, contexts={'CDS'})),
 ('hopZ', dict(tags=[], genes=['hopZ'], motifs={'AG'}, contexts={'CDS'})),
 ('oipA (hopH)', dict(tags=[('26695', 'HP0638')], genes=['oipA', 'hopH'], motifs={'AG'}, contexts={'CDS'})),
 ('babA', dict(tags=[], genes=['babA'], motifs={'AG'}, contexts={'CDS'})),
 ('modH', dict(tags=[('26695', 'HP1522')], genes=[], motifs={'C'}, contexts={'CDS'})),
 ('res/mod', dict(tags=[('J99', 'jhp1297')], genes=[], motifs={'C'}, contexts={'CDS'})),
]
def lit_row(prefix):
    for r in lit:
        if r['gene_name'].startswith(prefix): return r
    raise KeyError(prefix)
known_members = collections.defaultdict(list)  # locus_id -> [(lit gene, basis)]
unmatched = []
for key, m in MATCH:
    hit = False
    for r in cat:
        by_tag = any(old2new.get((s, t)) == r['locus_tag'] and s == r['reference_short'] for s, t in m['tags'])
        by_gene = r['gene_name'] in m['genes'] and r['reference_short'] in ('26695', 'J99')  # paralog gene names are inconsistent across strains (sabA family): use the two curated references only
        if not (by_tag or by_gene): continue
        if r['motif_canonical'] in m['motifs'] and r['context'] in m['contexts'] and r['interruption_flag'] in ('NONE', 'NA', ''):
            known_members[r['locus_id']].append((key, ('locus_tag_via_old_locus_tag' if by_tag else 'gene_name') + f":{r['reference_short']}:{r['locus_tag']}"))
            hit = True
    if not hit: unmatched.append(key)
# apply evidence class
for r in cat:
    if r['locus_id'] in known_members:
        k = known_members[r['locus_id']]
        r['evidence_class'] = 'KNOWN_PV'
        r['evidence_source'] = 'literature:' + ';'.join(sorted({x[0] for x in k})) + ' (' + ';'.join(sorted({x[1] for x in k}))[:200] + ')'
# functional interpretation flag per FS rules (justification of ON/OFF frame logic; state still requires a callable allele)
for r in cat:
    if r['caller_eligible'].startswith('NO'): r['functional_interpretation_justified'] = 'NO'
    elif r['context'] == 'CDS' and int(r['motif_length']) % 3 != 0: r['functional_interpretation_justified'] = 'YES'
    elif r['context'] == 'CDS': r['functional_interpretation_justified'] = 'NO'  # FS03: ALLELE_STATE
    else: r['functional_interpretation_justified'] = 'NO'  # FS05/FS06: REGULATORY_STATE / ALLELE_STATE
cols = list(cat[0])
os.makedirs('metadata/g5', exist_ok=True)
def wr(fn, rows, c):
    with open(fn, 'w', newline='') as f:
        w = csv.DictWriter(f, c, delimiter='\t', lineterminator='\n', restval='NA'); w.writeheader(); w.writerows(rows)
wr('metadata/g5/g5_repeat_catalogue.tsv', cat, cols)
wr('metadata/g5/g5_reference_panel.tsv', panel, list(panel[0]))
wr('metadata/g5/g5_locus_orthology.tsv', orth, list(orth[0]))
wr('metadata/g5/g5_repeat_flanks.tsv', [dict(locus_id=r['locus_id'], member_id=r['member_id'], reference_short=r['reference_short'], left_flank_40=r['left_flank_40'], right_flank_40=r['right_flank_40'], left_flank_sha256=r['left_flank_sha256'], right_flank_sha256=r['right_flank_sha256']) for r in cat], ['locus_id', 'member_id', 'reference_short', 'left_flank_40', 'right_flank_40', 'left_flank_sha256', 'right_flank_sha256'])
# known loci table: all literature rows + catalogue match
kn = []
for r in lit:
    key = next((k for k, _ in MATCH if r['gene_name'].startswith(k.split(' ')[0].split('/')[0]) or (k == 'res/mod' and r['gene_name'].startswith('res/mod'))), None)
    ids = sorted({lid for lid, v in known_members.items() if any(x[0] == key for x in v)}) if key else []
    d = dict(r); d['catalogue_locus_ids_matched'] = ';'.join(ids) or 'NONE'
    d['catalogue_match_basis'] = ('gene/locus tag plus motif class plus context (RD11)' if ids else ('literature has no motif or no locus tag: not matched, no KNOWN_PV flag assigned' if r['evidence_strength'] != 'NOT_FOUND' else 'NOT_FOUND in literature search'))
    kn.append(d)
wr('metadata/g5/g5_known_pv_loci.tsv', kn, list(kn[0]))
ev = collections.Counter(r['evidence_class'] for r in cat)
elig = sum(1 for r in cat if r['caller_eligible'].startswith('YES'))
loci = len({r['locus_id'] for r in cat})
print('members', len(cat), 'loci', loci, dict(ev), 'eligible', elig, 'unmatched literature keys', unmatched, 'known loci ids', len(known_members))
json.dump(dict(members=len(cat), loci=loci, evidence_members=dict(ev), eligible=elig, known_locus_ids=sorted(known_members), unmatched=unmatched), open('metadata/g5/g5_merge_summary.json', 'w'), indent=1)
