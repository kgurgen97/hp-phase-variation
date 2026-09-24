#!/usr/bin/env python3
"""G6 genome-wide disease-blind catalogue boundary audit per frozen rules AR01 to AR09 (metadata/g6/g6_boundary_audit_rules.tsv, sha256 in .sha256).
Inputs: catalogue members, orthology table. Nothing is altered."""
import csv, collections, itertools, subprocess
assert subprocess.run(['shasum', '-a', '256', '-c', 'metadata/g6/g6_boundary_audit_rules.sha256'], capture_output=True).returncode == 0, 'audit rules changed'
rd = lambda f: list(csv.DictReader(open(f), delimiter='\t'))
cat = rd('metadata/g5/g5_repeat_catalogue.tsv'); orth = {r['member_id']: r for r in rd('metadata/g5/g5_locus_orthology.tsv')}
COMP = str.maketrans('ACGT', 'TGCA'); rc = lambda s: s.translate(COMP)[::-1]
loci = collections.defaultdict(list)
for r in cat: loci[r['locus_id']].append(r)
def ident_end(a, b):
    """identity of a and b aligned at their END (tract-adjacent side), ungapped, shifts -3..3; returns best fraction over >=20 overlapping bases"""
    best = 0.0
    for sh in range(-3, 4):
        x = a[:len(a) - sh] if sh > 0 else a; y = b if sh > 0 else b[:len(b) + sh] if sh < 0 else b
        if sh > 0: y = b
        n = min(len(x), len(y))
        if n < 20: continue
        xs, ys = x[len(x) - n:], y[len(y) - n:]
        best = max(best, sum(1 for p, q in zip(xs, ys) if p == q) / n)
    return best
def anchors(m):
    return m['left_flank_40'][-25:], m['right_flank_40'][:25]
def pair(mi, mj):
    Li, Ri = anchors(mi); Lj, Rj = anchors(mj)
    same = (ident_end(Li, Lj), ident_end(rc(Ri), rc(Rj)))            # left ends at tract; right compared as reverse complements (end at tract)
    opp = (ident_end(Li, rc(Rj)), ident_end(rc(Ri), Lj))
    best = same if min(same) >= min(opp) else opp
    return best
OVR = {'HPL00546': ('BOUNDARY_AMBIGUOUS', 'independent alignment: 1 bp insertion 3 bp upstream of the tract and a 2 bp deletion at the tract start in nearly all reads (HPAG1 member context), so the tract length depends on where the flank/tract boundary is drawn; anchor-based short and HiFi calls 7 bp, alignment-based 6 bp (post-hoc diagnostic)'),
       'HPL00438': ('BOUNDARY_AMBIGUOUS', 'singleton; independent alignment shows a complex indel pattern inside and adjacent to the tract (insertions at +3, +6, +7 and a deletion at +13 relative to the tract start) in most reads; anchor-based 7 bp versus alignment-based 13 bp (post-hoc diagnostic)'),
       'HPL00690': ('ORTHOLOGY_AMBIGUOUS', 'members J99 and India7 are grouped by product/neighbor only; anchors are non-homologous; in the same reads the two members give different tract contexts (J99 T8 versus India7 14 bp non-homopolymer segment); see G6 independent-check report')}
rows = []
for lid, ms in sorted(loci.items()):
    n = len(ms); motifs = {m['motif_canonical'] for m in ms}; ks = {int(m['motif_length']) for m in ms}
    k = min(ks); tl = [int(m['tract_length_bp']) for m in ms]
    motif_ok = len(motifs) == 1 and (k == 1 or all((a - b) % k == 0 for a, b in itertools.combinations(tl, 2)))
    ctx_ok = len({m['context'] == 'CDS' for m in ms}) == 1
    refs = [m['reference_short'] for m in ms]; dup = len(refs) != len(set(refs))
    oflags = [orth[m['member_id']]['ambiguity_flag'] for m in ms]; obasis = [orth[m['member_id']]['basis'] for m in ms]
    orth_amb = any(f != 'NONE' for f in oflags) or any(('FLANK' not in b and b != 'SINGLETON') for b in obasis) or dup
    short_anchor = any(len(anchors(m)[0]) < 25 or len(anchors(m)[1]) < 25 or 'N' in ''.join(anchors(m)) for m in ms)
    ctx_cons = 'NA_SINGLE_MEMBER'; graph_note = ''
    if n >= 2 and not short_anchor:
        par = list(range(n))
        def f(x):
            while par[x] != x: par[x] = par[par[x]]; x = par[x]
            return x
        one_sided = 0; homol = 0; mins = []
        for i, j in itertools.combinations(range(n), 2):
            a, b = pair(ms[i], ms[j]); mins.append(min(a, b))
            if min(a, b) >= 0.72: par[f(i)] = f(j); homol += 1
            elif max(a, b) >= 0.72: one_sided += 1
        comp = len({f(x) for x in range(n)})
        ctx_cons = 'YES' if comp == 1 else 'NO'
        graph_note = f'components {comp}; homologous pairs {homol}/{n*(n-1)//2}; one-sided pairs {one_sided}; min pair identity {min(mins):.2f}'
    boundary_ok = motif_ok and ctx_ok and ctx_cons in ('YES', 'NA_SINGLE_MEMBER')
    if lid in OVR: status, note = OVR[lid]
    elif orth_amb: status = 'ORTHOLOGY_AMBIGUOUS'; note = 'orthology ambiguity flag, non-flank orthology basis, or duplicate reference genome'
    elif n >= 2 and not (motif_ok and ctx_ok and ctx_cons == 'YES'): status = 'BOUNDARY_AMBIGUOUS'; note = f'motif_ok {motif_ok}; context_ok {ctx_ok}; {graph_note}'
    elif n == 1 or short_anchor: status = 'INSUFFICIENT_TO_RESOLVE'; note = 'single member: no cross-member evidence' if n == 1 else 'anchor shorter than 25 bp or contains N'
    else: status = 'BOUNDARY_STABLE'; note = graph_note
    rows.append(dict(locus_id=lid, member_count=n, motif=';'.join(sorted(motifs)), boundary_consistent_across_members=('YES' if (motif_ok and ctx_ok) else 'NO') if n >= 2 else 'NA_SINGLE_MEMBER', tract_context_consistent=ctx_cons,
                     paralog_orthology_ambiguity='YES' if orth_amb else 'NO', audit_status=status,
                     primary_technical_eligibility='ELIGIBLE_BOUNDARY_STABLE' if status == 'BOUNDARY_STABLE' else 'INELIGIBLE_' + status, evidence_class=ms[0]['evidence_class'] if len({m['evidence_class'] for m in ms}) == 1 else 'MIXED', notes=note))
with open('metadata/g6/g6_catalogue_boundary_audit.tsv', 'w', newline='') as fh:
    w = csv.DictWriter(fh, list(rows[0]), delimiter='\t', lineterminator='\n'); w.writeheader(); w.writerows(rows)
print(len(rows), collections.Counter(r['audit_status'] for r in rows))
print('multi-member', collections.Counter(r['audit_status'] for r in rows if r['member_count'] >= 2))
