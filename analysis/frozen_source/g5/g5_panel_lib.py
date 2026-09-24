"""Helpers to read the reference panel annotation and scan sub-threshold tracts (used for known-PV reconciliation; disease-blind)."""
import csv, gzip, os
COMP = str.maketrans('ACGT', 'TGCA')
def rc(s): return s.translate(COMP)[::-1]
def load_panel():
    P = {}
    for r in csv.DictReader(open('metadata/g5/g5_reference_panel.tsv'), delimiter='\t'):
        acc = r['assembly_accession_version']; d = f'data/reference/g5_panel/{acc}/'
        fa = [d + f for f in os.listdir(d) if f.endswith('.fna.gz')][0]; gff = [d + f for f in os.listdir(d) if f.endswith('.gff.gz')][0]
        seqs = {}; name = None
        for l in gzip.open(fa, 'rt'):
            if l.startswith('>'): name = l[1:].split()[0]; seqs[name] = []
            else: seqs[name].append(l.strip().upper())
        seqs = {k: ''.join(v) for k, v in seqs.items()}
        genes = {}
        for l in gzip.open(gff, 'rt'):
            if l.startswith('#'): continue
            c = l.rstrip('\n').split('\t')
            if len(c) < 9 or c[2] != 'gene': continue
            a = dict(x.split('=', 1) for x in c[8].split(';') if '=' in x)
            g = dict(contig=c[0], start=int(c[3]), end=int(c[4]), strand=c[6], locus_tag=a.get('locus_tag'), gene=a.get('gene', ''),
                     old=[x for x in a.get('old_locus_tag', '').replace('%2C', ',').split(',') if x], biotype=a.get('gene_biotype', ''))
            genes[g['locus_tag']] = g
        P[r['short_code']] = dict(seqs=seqs, genes=genes, acc=acc)
    return P
def tracts(seq, lo, hi, min_mono=5, min_copies=3):
    """Sub-threshold pure tracts in seq[lo-1:hi] (1-based inclusive); mono >= min_mono bp, di/tri >= min_copies copies."""
    out = []; s = seq[lo - 1:hi]; n = len(s)
    for k in (1, 2, 3):
        i = 0
        while i < n - k:
            if s[i] == s[i + k]:
                j = i
                while j < n - k and s[j] == s[j + k]: j += 1
                L = j - i + k; m = s[i:i + k]
                if (k == 1 and L >= min_mono) or (k > 1 and L // k >= min_copies and len(set(m)) > 1 and m != m[:1] * k):
                    out.append(dict(start=lo + i, end=lo + i + L - 1, k=k, motif=m, length=L))
                i = j + 1
            else: i += 1
    return out
def canon(m):
    mm = [m[i:] + m[:i] for i in range(len(m))] + [rc(m)[i:] + rc(m)[:i] for i in range(len(m))]
    return min(mm)
