"""Synthetic read simulator with truth tables (Illumina-like substitution errors, SE/PE, mixtures, stutter)."""
import gzip
import math
import random

from .common import revcomp, motif_class, write_tsv

TRUTH_FIELDS = ["sample_id", "locus_id", "member_id", "layout", "depth", "mix_label", "stutter", "seed",
                "true_alleles_copies", "true_alleles_bp", "major_copies", "major_bp", "ref_copies",
                "n_fragments", "n_fragments_by_allele"]

LAYOUTS = {"SE150": dict(read_len=150, paired=False, flen_mean=350, flen_sd=60),
           "PE250": dict(read_len=250, paired=True, flen_mean=420, flen_sd=70),
           "PE150": dict(read_len=150, paired=True, flen_mean=300, flen_sd=50)}


class StutterSpec:
    """Per-fragment stutter: dominant tract shifts by -1/+1 unit with probability a+b*copies (down) and
    up_ratio times that (up); shifts of 2 with k2_ratio times. Truth model of the simulator only."""

    def __init__(self, a=None, b=None, up_ratio=0.5, k2_ratio=0.25):
        self.a = a or {"MONO": 0.005, "DI": 0.01, "POLY": 0.003}
        self.b = b or {"MONO": 0.0025, "DI": 0.004, "POLY": 0.001}
        self.up_ratio, self.k2_ratio = up_ratio, k2_ratio

    def draw_shift(self, rng, mclass, copies):
        r = max(0.0, min(0.5, self.a[mclass] + self.b[mclass] * copies))
        u = rng.random()
        for k, pd, pu in ((1, r, r * self.up_ratio), (2, r * self.k2_ratio, r * self.up_ratio * self.k2_ratio)):
            if u < pd:
                return -k
            u -= pd
            if u < pu:
                return k
            u -= pu
        return 0


def make_tract(m, copies_delta):
    """Tract for reference copies + delta (units of the motif; interrupted loci: extend/trim at the end)."""
    ref_tract = m.ref_tract
    n = len(ref_tract) + copies_delta * m.motif_len
    if n <= 0:
        raise ValueError("allele would have no tract")
    if copies_delta >= 0:
        return (ref_tract + m.motif_plus * copies_delta) if copies_delta else ref_tract
    return ref_tract[:n]


def attach_ref_tracts(members, refs):
    for m in members:
        m.ref_tract = refs[m.contig][m.start - 1:m.end]


def _mutate(rng, seq, err_scale, mate2=False):
    L = len(seq)
    out, qs = [], []
    for i, b in enumerate(seq):
        q = int(round(37 - 14 * (i / L) ** 3 - (2 if mate2 else 0) + rng.gauss(0, 2)))
        q = max(2, min(41, q))
        if rng.random() < err_scale * 10 ** (-q / 10.0):
            b = rng.choice([x for x in "ACGT" if x != b])
        out.append(b)
        qs.append(chr(q + 33))
    return "".join(out), "".join(qs)


def simulate_sample(members, refs, alleles, out_prefix, layout="SE150", depth=30, seed=1, stutter=None,
                    orphan_fraction=0.0, err_scale=1.0, sample_id="SIM", mix_label="", window=700):
    """alleles: {member_id: [(delta_units, fraction), ...]} (default reference allele). Writes FASTQ.gz:
    SE -> <prefix>.se.fastq.gz; PE -> <prefix>_R1/_R2.fastq.gz plus <prefix>.orphan.fastq.gz when
    orphan_fraction > 0 (a mate of that fraction of pairs is dropped and kept as an unpaired read).
    Returns (dict of file paths, truth rows)."""
    lay = LAYOUTS[layout]
    L = lay["read_len"]
    rng = random.Random(seed)
    if members and members[0].ref_tract is None:
        attach_ref_tracts(members, refs)
    files = {}
    if lay["paired"]:
        f1 = gzip.open(out_prefix + "_R1.fastq.gz", "wt")
        f2 = gzip.open(out_prefix + "_R2.fastq.gz", "wt")
        fo = gzip.open(out_prefix + ".orphan.fastq.gz", "wt") if orphan_fraction > 0 else None
        files = {"pe": [(out_prefix + "_R1.fastq.gz", out_prefix + "_R2.fastq.gz")],
                 "orphan": [out_prefix + ".orphan.fastq.gz"] if fo else []}
    else:
        f1 = gzip.open(out_prefix + ".se.fastq.gz", "wt")
        f2 = fo = None
        files = {"se": [out_prefix + ".se.fastq.gz"]}
    truth, nread = [], 0
    for m in members:
        spec = alleles.get(m.member_id, [(0, 1.0)])
        ref = refs[m.contig]
        lo = max(0, m.start - 1 - window)
        hi = m.end + window
        mc = motif_class(m.motif_len)
        haps = {}
        # neighbouring catalogue loci inside the window carry their major simulated allele (no phase model)
        neigh = []
        for o in members:
            if o is not m and o.contig == m.contig and o.start - 1 >= lo and o.end <= hi:
                osp = alleles.get(o.member_id, [(0, 1.0)])
                neigh.append((o, max(osp, key=lambda x: x[1])[0]))

        def hap(delta):
            if delta not in haps:
                edits = sorted([(m.start - 1, m.end, make_tract(m, delta))] +
                               [(o.start - 1, o.end, make_tract(o, od)) for o, od in neigh])
                parts, pos = [], lo
                for a, b, t in edits:
                    parts.append(ref[pos:a])
                    parts.append(t)
                    pos = b
                parts.append(ref[pos:hi])
                haps[delta] = "".join(parts)
            return haps[delta]

        left_ctx = ref[lo:m.start - 1]
        right_ctx = ref[m.end:hi]
        hl = hi - lo + 0
        mates = 2 if lay["paired"] else 1
        n_frag = int(round(depth * (hl - lay["flen_mean"]) / (L * mates)))
        counts = {}
        cum, acc = [], 0.0
        for d, f in spec:
            acc += f
            cum.append((acc, d))
        for _ in range(n_frag):
            u = rng.random() * acc
            d = next(dd for c, dd in cum if u <= c + 1e-12)
            counts[d] = counts.get(d, 0) + 1
            eff = d
            if stutter:
                sh = stutter.draw_shift(rng, mc, m.ref_copies + d)
                if sh and m.ref_copies + d + sh >= 1:
                    eff = d + sh
            h = hap(eff)
            flen = max(L, int(rng.gauss(lay["flen_mean"], lay["flen_sd"])))
            flen = min(flen, len(h))
            st = rng.randint(0, len(h) - flen)
            frag = h[st:st + flen]
            if rng.random() < 0.5:
                frag = revcomp(frag)
            nread += 1
            name = "@sim%d" % nread
            if not lay["paired"]:
                s, q = _mutate(rng, frag[:L], err_scale)
                f1.write("%s\n%s\n+\n%s\n" % (name, s, q))
            else:
                s1, q1 = _mutate(rng, frag[:L], err_scale)
                s2, q2 = _mutate(rng, revcomp(frag[-L:]), err_scale, True)
                if fo and rng.random() < orphan_fraction:
                    keep_first = rng.random() < 0.5
                    s, q = (s1, q1) if keep_first else (s2, q2)
                    fo.write("%s/%d\n%s\n+\n%s\n" % (name, 1 if keep_first else 2, s, q))
                else:
                    f1.write("%s/1\n%s\n+\n%s\n" % (name, s1, q1))
                    f2.write("%s/2\n%s\n+\n%s\n" % (name, s2, q2))
        major = max(spec, key=lambda x: x[1])[0]
        truth.append({
            "sample_id": sample_id, "locus_id": m.locus_id, "member_id": m.member_id, "layout": layout,
            "depth": depth, "mix_label": mix_label, "stutter": "YES" if stutter else "NO", "seed": seed,
            "true_alleles_copies": ";".join("%g:%.3f" % (m.ref_copies + d, f) for d, f in spec),
            "true_alleles_bp": ";".join("%d:%.3f" % (m.tract_len + d * m.motif_len, f) for d, f in spec),
            "major_copies": "%g" % (m.ref_copies + major), "major_bp": m.tract_len + major * m.motif_len,
            "ref_copies": "%g" % m.ref_copies, "n_fragments": n_frag,
            "n_fragments_by_allele": ";".join("%d:%d" % kv for kv in sorted(counts.items()))})
    for f in (f1, f2, fo):
        if f:
            f.close()
    return files, truth
