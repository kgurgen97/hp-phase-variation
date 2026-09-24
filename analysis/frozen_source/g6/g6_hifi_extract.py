#!/usr/bin/env python3
"""G6 HiFi (CCS) locus evidence per frozen rules H05. Fast exact-seed version: k=16 seeds sampled every 4 bases, anchors (25 bp) verified by Hamming
distance (at most 2 mismatches; HiFi reads are near-error-free in flanks), both orientations, arbitration between members of different loci as in
the Stage A extractor. usage: g6_hifi_extract.py BIOSAMPLE FASTQ.gz OUT.tsv.gz"""
import sys, gzip, collections
sys.path.insert(0, 'scripts/g5')
from g5caller.common import load_catalogue, revcomp, iter_fastq
A, K, STRIDE, MAXMM, MAXTRACT = 25, 16, 4, 2, 150
def build(members):
    aid = {}; alist = []           # unique anchor string -> id
    slots = collections.defaultdict(list)   # anchor id -> [(mi, side, orient)]
    for mi, m in enumerate(members):
        la, ra = m.left40[-A:], m.right40[:A]
        if len(la) < A or len(ra) < A or 'N' in la or 'N' in ra: continue
        for orient, side, seq in (('+', 'L', la), ('+', 'R', ra), ('-', 'R', revcomp(ra)), ('-', 'L', revcomp(la))):
            if seq not in aid: aid[seq] = len(alist); alist.append(seq)
            slots[aid[seq]].append((mi, side, orient))
    idx = collections.defaultdict(list)
    for a, seq in enumerate(alist):
        for j in range(A - K + 1): idx[seq[j:j + K]].append((a, j))
    return alist, slots, idx
def mm(a, b):
    n = 0
    for x, y in zip(a, b):
        if x != y:
            n += 1
            if n > MAXMM: return n
    return n
def process(sample, fq, out_path, catalogue='metadata/g5/g5_repeat_catalogue.tsv'):
    members = load_catalogue(catalogue); alist, slots, idx = build(members)
    out = gzip.open(out_path, 'wt'); out.write('\t'.join(['sample', 'locus_id', 'member_id', 'decoy', 'read_id', 'strand', 'tract_len', 'edits', 'read_len']) + '\n')
    nreads = nspan = namb = 0; get = idx.get
    for rid, s, q in iter_fastq(fq):
        nreads += 1; L = len(s); hits = {}
        for i in range(0, L - K + 1, STRIDE):
            h = get(s[i:i + K])
            if h:
                for a, j in h:
                    st = i - j
                    if st < 0 or st + A > L or (a, st) in hits: continue
                    d = mm(s[st:st + A], alist[a])
                    if d <= MAXMM: hits[(a, st)] = d
        if not hits: continue
        found = collections.defaultdict(dict)    # (mi, orient) -> {side: (mm, start, end)}
        for (a, st), d in hits.items():
            for mi, side, orient in slots[a]:
                cur = found[(mi, orient)].get(side)
                if cur is None or d < cur[0]: found[(mi, orient)][side] = (d, st, st + A)
        cands = []
        for (mi, orient), dd in found.items():
            if 'L' in dd and 'R' in dd:
                tl = dd['R'][1] - dd['L'][2] if orient == '+' else dd['L'][1] - dd['R'][2]
                if 0 <= tl <= MAXTRACT: cands.append((mi, orient, tl, dd['L'][0] + dd['R'][0], dd))
        drop = set()
        for a in range(len(cands)):
            for b in range(a + 1, len(cands)):
                ma, mb = cands[a], cands[b]
                if members[ma[0]].locus_id == members[mb[0]].locus_id or ma[1] != mb[1]: continue
                pa = min(ma[4]['L'][1], ma[4]['R'][1]); pb = min(mb[4]['L'][1], mb[4]['R'][1])
                if abs(pa - pb) <= A:
                    if ma[3] + 2 <= mb[3]: drop.add(b)
                    elif mb[3] + 2 <= ma[3]: drop.add(a)
                    else: drop.update((a, b)); namb += 1
        for k, (mi, orient, tl, ed, dd) in enumerate(cands):
            if k in drop: continue
            m = members[mi]
            out.write('\t'.join([sample, m.locus_id, m.member_id, '1' if not m.caller_eligible else '0', rid.split()[0], orient, str(tl), str(ed), str(L)]) + '\n'); nspan += 1
    out.close(); print(sample, 'reads', nreads, 'spanning rows', nspan, 'ambiguous paralog events', namb, flush=True)
if __name__ == '__main__':
    process(sys.argv[1], sys.argv[2], sys.argv[3])
