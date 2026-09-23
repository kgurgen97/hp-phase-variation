"""Stage 1: anchor matching and spanning-read evidence extraction (no aligner)."""
import gzip
import json
from collections import defaultdict

from .common import revcomp, iter_fastq, load_catalogue, Params

EVIDENCE_FIELDS = ["kind", "sample_id", "locus_id", "member_id", "fragment_id", "source", "strand",
                   "read_len", "tract", "anchor_q", "l_mm", "r_mm", "variants", "n"]


def _low_complexity(kmer):
    if len(set(kmer)) == 1:
        return True
    return kmer[:2] * (len(kmer) // 2) == kmer[: 2 * (len(kmer) // 2)] and len(kmer) % 2 == 0


class AnchorIndex:
    """k-mer seed index over left/right anchors of caller-eligible members."""

    def __init__(self, members, params):
        self.p = params
        self.A = params.anchor_len
        self.k = params.seed_k
        self.stride = params.seed_stride
        self.members = []
        self.decoy = []  # non-eligible members are indexed only to arbitrate paralog cross-talk; never reported
        self.left, self.right = [], []
        self.seeds = defaultdict(list)
        for m in members:
            if len(m.left40) < self.A or len(m.right40) < self.A:
                continue
            la = m.left40[-self.A:]
            ra = m.right40[: self.A]
            if "N" in la or "N" in ra:
                continue
            mi = len(self.members)
            self.members.append(m)
            self.decoy.append(not m.caller_eligible)
            self.left.append(la)
            self.right.append(ra)
            for side, a in ((0, la), (1, ra)):
                for j in range(self.A - self.k + 1):
                    km = a[j:j + self.k]
                    if _low_complexity(km):
                        continue
                    self.seeds[km].append((mi, side, j))

    def find(self, seq, qual):
        """Return {member_idx: {'L': (start, mm, mmlist), 'R': ...}} for one orientation."""
        A, k, stride = self.A, self.k, self.stride
        n = len(seq)
        seen = set()
        best = {}
        seeds = self.seeds
        for i in range(0, n - k + 1, stride):
            hits = seeds.get(seq[i:i + k])
            if not hits:
                continue
            for mi, side, j in hits:
                st = i - j
                if st < 0 or st + A > n or (mi, side, st) in seen:
                    continue
                seen.add((mi, side, st))
                ref = self.right[mi] if side else self.left[mi]
                sub = seq[st:st + A]
                mmpos = [x for x in range(A) if sub[x] != ref[x]]
                if len(mmpos) > self.p.max_anchor_mismatch:
                    continue
                key = (mi, side)
                cur = best.get(key)
                if cur is None or len(mmpos) < len(cur[1]):
                    best[key] = (st, mmpos)
        out = {}
        for (mi, side), (st, mmpos) in best.items():
            ref = self.right[mi] if side else self.left[mi]
            sub = seq[st:st + A]
            var = ["%s%d:%s>%s" % ("R" if side else "L", x, ref[x], sub[x]) for x in mmpos]
            out.setdefault(mi, {})["R" if side else "L"] = (st, len(mmpos), var)
        return out


def _qmean(qual, a, b):
    return sum(ord(c) for c in qual[a:b]) / max(1, b - a) - 33


def eval_read(index, seq, qual):
    """Evaluate one read on both strands. Returns {member_idx: record} where record is
    ('SPAN', strand, tract, qmean, lmm, rmm, variants) or ('PART', set_of_sides)."""
    A = index.A
    if len(seq) < index.p.min_read_len:
        return {}
    res = {}
    for strand in ("+", "-"):
        if strand == "+":
            s, q = seq, qual
        else:
            s, q = revcomp(seq), qual[::-1]
        found = index.find(s, q)
        for mi, d in found.items():
            rec = None
            if "L" in d and "R" in d:
                ls, lmm, lv = d["L"]
                rs, rmm, rv = d["R"]
                tl = rs - (ls + A)
                if 0 <= tl <= index.p.max_tract_scan_len:
                    qm = (_qmean(q, ls, ls + A) + _qmean(q, rs, rs + A)) / 2.0
                    rec = ("SPAN", strand, s[ls + A:rs], qm, lmm, rmm, lv + rv, ls, rs)
                else:
                    rec = ("PART", {"L", "R"})
            else:
                rec = ("PART", set(d.keys()))
            old = res.get(mi)
            if old is None or (rec[0] == "SPAN" and (old[0] != "SPAN" or rec[4] + rec[5] < old[4] + old[5])):
                res[mi] = rec
    res = _arbitrate_paralogs(index, res)
    return {mi: r for mi, r in res.items() if not index.decoy[mi]}


def _arbitrate_paralogs(index, res):
    """Cross-locus anchor arbitration (G5 main-process addition). Two members of DIFFERENT loci whose left and right anchors
    hit the same read region (paralogous flanks tolerated by the mismatch allowance) compete for the read: the member with at least
    2 fewer anchor mismatches keeps it; otherwise the read is dropped for both (ambiguous). Members of the same locus (orthologs)
    are never in competition. Reads spanning two truly distinct nearby loci are unaffected because their anchor positions differ."""
    A = index.A
    spans = [(mi, r) for mi, r in res.items() if r[0] == "SPAN"]
    if len(spans) < 2:
        return res
    drop = set()
    for i in range(len(spans)):
        for j in range(i + 1, len(spans)):
            (a, ra), (b, rb) = spans[i], spans[j]
            if ra[1] != rb[1] or index.members[a].locus_id == index.members[b].locus_id:
                continue
            if abs(ra[7] - rb[7]) <= A and abs(ra[8] - rb[8]) <= A:
                ca, cb = ra[4] + ra[5], rb[4] + rb[5]
                if ca + 2 <= cb:
                    drop.add(b)
                elif cb + 2 <= ca:
                    drop.add(a)
                else:
                    drop.update((a, b))
                index.ambiguous_paralog = getattr(index, "ambiguous_paralog", 0) + 1
    for mi in drop:
        res[mi] = ("PART", {"L", "R"})
    return res


def _flip(strand):
    return "-" if strand == "+" else "+"


def extract_sample(catalogue, sample_id, out_path, se_files=(), pe_pairs=(), orphan_files=(),
                   params=None, max_fragments=None):
    """Stream all files of ONE sample; write evidence TSV(.gz) plus <out>.meta.json. Returns stats dict.

    se_files: single-end FASTQ (source 'single'); pe_pairs: [(r1, r2)] (sources mate1/mate2, merged per
    fragment); orphan_files: unpaired reads left by trimming/pair-repair (source 'unpaired'; they add
    evidence to the SAME sample, each read its own fragment)."""
    params = params or Params()
    members = catalogue if isinstance(catalogue, list) else load_catalogue(catalogue)
    index = AnchorIndex(members, params)
    opener = gzip.open if str(out_path).endswith(".gz") else open
    partial = defaultdict(int)
    stats = defaultdict(int)
    fid = [0]
    out = opener(out_path, "wt")
    out.write("\t".join(EVIDENCE_FIELDS) + "\n")

    def emit(mi, source, strand, rl, tract, qm, lmm, rmm, var, nmates=1):
        m = index.members[mi]
        out.write("\t".join([
            "SPAN", sample_id, m.locus_id, m.member_id, str(fid[0]), source, strand, str(rl), tract,
            "%.2f" % qm, str(lmm), str(rmm), ",".join(var) if var else ".", str(nmates)]) + "\n")
        stats["span_rows"] += 1

    def handle_single(seq, qual, source):
        fid[0] += 1
        stats["fragments"] += 1
        stats["reads"] += 1
        for mi, rec in eval_read(index, seq, qual).items():
            if rec[0] == "SPAN":
                emit(mi, source, rec[1], len(seq), rec[2], rec[3], rec[4], rec[5], rec[6])
            else:
                for side in rec[1]:
                    partial[(mi, side)] += 1

    def handle_pair(s1, q1, s2, q2):
        fid[0] += 1
        stats["fragments"] += 1
        stats["reads"] += 2
        r1, r2 = eval_read(index, s1, q1), eval_read(index, s2, q2)
        for mi in set(r1) | set(r2):
            a, b = r1.get(mi), r2.get(mi)
            sa = a is not None and a[0] == "SPAN"
            sb = b is not None and b[0] == "SPAN"
            if sa and sb:
                if a[2] == b[2]:
                    emit(mi, "mate1+mate2", a[1], len(s1), a[2], (a[3] + b[3]) / 2, a[4] + b[4],
                         a[5] + b[5], sorted(set(a[6]) | set(b[6])), 2)
                    stats["merged_overlap_fragments"] += 1
                else:
                    stats["mate_conflict"] += 1
                    if params.pe_conflict_policy == "BEST":
                        w, src, strand = (a, "mate1", a[1]) if a[3] >= b[3] else (b, "mate2", _flip(b[1]))
                        emit(mi, src, strand, len(s1), w[2], w[3], w[4], w[5], w[6])
                    else:
                        partial[(mi, "CONFLICT")] += 1
            elif sa:
                emit(mi, "mate1", a[1], len(s1), a[2], a[3], a[4], a[5], a[6])
            elif sb:
                emit(mi, "mate2", _flip(b[1]), len(s2), b[2], b[3], b[4], b[5], b[6])
            else:
                sides = (a[1] if a else set()) | (b[1] if b else set())
                for side in sides:
                    partial[(mi, side)] += 1

    def limit():
        return max_fragments is not None and stats["fragments"] >= max_fragments

    for f in se_files:
        for _, s, q in iter_fastq(f):
            handle_single(s, q, "single")
            if limit():
                break
    for r1, r2 in pe_pairs:
        it2 = iter_fastq(r2)
        for n1, s1, q1 in iter_fastq(r1):
            try:
                n2, s2, q2 = next(it2)
            except StopIteration:
                raise ValueError("mate2 file shorter than mate1: %s / %s" % (r1, r2))
            handle_pair(s1, q1, s2, q2)
            if limit():
                break
        else:
            if next(it2, None) is not None:
                raise ValueError("mate2 file longer than mate1: %s / %s" % (r1, r2))
    for f in orphan_files:
        for _, s, q in iter_fastq(f):
            handle_single(s, q, "unpaired")
            if limit():
                break
    for (mi, side), n in sorted(partial.items()):
        m = index.members[mi]
        kind = "CONFLICT" if side == "CONFLICT" else "PARTIAL_" + side
        out.write("\t".join([kind, sample_id, m.locus_id, m.member_id, "", "", "", "", "", "", "", "", "", str(n)]) + "\n")
    out.close()
    stats = dict(stats)
    stats["indexed_members"] = len(index.members)
    stats["ambiguous_paralog_events"] = getattr(index, "ambiguous_paralog", 0)
    with open(str(out_path) + ".meta.json", "w") as fh:
        json.dump({"sample_id": sample_id, "stats": stats, "params": dict(params),
                   "se_files": list(se_files), "pe_pairs": [list(x) for x in pe_pairs],
                   "orphan_files": list(orphan_files)}, fh, indent=1)
    return stats
