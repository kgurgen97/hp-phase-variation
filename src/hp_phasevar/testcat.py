"""Tiny TEST catalogue built from a reference FASTA. CDS coordinates here are pseudo-annotations
(longest ATG..stop ORF, >=100 aa, that fully contains the tract) used ONLY for tests/smoke benchmarks;
they are not real gene annotations."""
import re

from .common import (read_fasta, revcomp, canonical_motif, sha256, write_tsv)
from .functional import CODON

CATALOGUE_FIELDS = [
    "locus_id", "member_id", "reference_id", "reference_short", "contig", "start", "end", "strand",
    "motif_observed_plus", "motif_canonical", "motif_length", "tract_length_bp", "copy_number", "purity",
    "interruption_flag", "context", "locus_tag", "gene_name", "product", "cds_id", "cds_start", "cds_end",
    "tract_offset_in_cds", "left_flank_40", "right_flank_40", "left_flank_sha256", "right_flank_sha256",
    "uniquely_anchorable", "caller_eligible", "ortholog_group", "ortholog_basis", "evidence_class",
    "evidence_source", "functional_interpretation_justified"]

_MONO = re.compile(r"([ACGT])\1{7,}")


def _orf_around(seq, ts, te, min_aa=100):
    """Longest ORF (ATG..stop, plus strand of seq) containing [ts,te); returns (start0, end0_exclusive) or None."""
    best = None
    for f in range(3):
        i = f
        start = None
        while i + 3 <= len(seq):
            cod = seq[i:i + 3]
            if start is None and cod == "ATG":
                start = i
            if CODON.get(cod) == "*":
                if start is not None and start < ts and te <= i + 3 and (i + 3 - start) // 3 - 1 >= min_aa:
                    if best is None or (i + 3 - start) > (best[1] - best[0]):
                        best = (start, i + 3)
                start = None
            i += 3
    return best


def find_tracts(genome):
    out = []
    for m in _MONO.finditer(genome):
        out.append((m.start(), m.end(), m.group(1)))
    for pat in (re.compile(r"((?:[ACGT]{2}))\1{4,}"), re.compile(r"((?:[ACGT]{3}))\1{3,}")):
        for m in pat.finditer(genome):
            mot = m.group(1)
            if len(set(mot)) == 1:
                continue
            out.append((m.start(), m.end(), mot))
    return sorted(out)


def build_tiny_catalogue(fasta, out_tsv, n_mono_cds=4, seed_offset=0, contig=None, ref_id="GCF_000008525.1"):
    """Select a small varied set of loci: mono tracts in CDS on both strands, a dinucleotide CDS locus,
    a trinucleotide CDS locus, a regulatory-labelled tract, a REPEAT_ONLY tract, a non-eligible one."""
    seqs = read_fasta(fasta)
    contig = contig or next(iter(seqs))
    g = seqs[contig]
    grc = revcomp(g)
    cand = find_tracts(g)
    picked, used_kinds = [], {"mono+": 0, "mono-": 0, "di": 0, "tri": 0, "reg": 0, "ro": 0, "inel": 0}
    for (s, e, mot) in cand:
        if s < 3000 or e > len(g) - 3000:
            continue
        tl = e - s
        k = len(mot)
        # 25-mer uniqueness of both flanks on both strands
        l25, r25 = g[s - 25:s], g[e:e + 25]
        if any(g.count(x) + grc.count(x) != 1 for x in (l25, r25)):
            continue
        if any(b in g[s - 40:s] + g[e:e + 40] for b in ("N",)):
            continue
        # avoid the flank itself extending the repeat (maximal tract)
        if g[s - k:s] == mot or g[e:e + k] == mot:
            continue
        win_lo, win_hi = s - 3000, e + 3000
        win = g[win_lo:win_hi]
        orf_p = _orf_around(win, s - win_lo, e - win_lo)
        wrc = revcomp(win)
        ts_rc = len(win) - (e - win_lo)
        orf_m = _orf_around(wrc, ts_rc, ts_rc + tl)
        info = None
        if orf_p:
            info = ("+", win_lo + orf_p[0] + 1, win_lo + orf_p[1])
        elif orf_m:
            info = ("-", win_lo + (len(win) - orf_m[1]) + 1, win_lo + (len(win) - orf_m[0]))
        kind = None
        if k == 1 and info and info[0] == "+" and used_kinds["mono+"] < (n_mono_cds + 1) // 2:
            kind = "mono+"
        elif k == 1 and info and info[0] == "-" and used_kinds["mono-"] < n_mono_cds // 2:
            kind = "mono-"
        elif k == 2 and info and used_kinds["di"] < 1:
            kind = "di"
        elif k == 3 and info and used_kinds["tri"] < 1:
            kind = "tri"
        elif k == 1 and not info and used_kinds["reg"] < 1:
            kind = "reg"
        elif k == 1 and not info and used_kinds["ro"] < 1 and used_kinds["reg"] == 1:
            kind = "ro"
        elif k == 1 and info and used_kinds["inel"] < 1 and used_kinds["mono+"] >= 1 and used_kinds["mono-"] >= 1:
            kind = "inel"
        if kind is None:
            continue
        used_kinds[kind] += 1
        picked.append((s, e, mot, info, kind))
    rows = []
    for i, (s, e, mot, info, kind) in enumerate(picked, 1):
        lid = "HPL%04d" % i
        left, right = g[s - 40:s], g[e:e + 40]
        row = dict(locus_id=lid, member_id=lid + ".M1", reference_id=ref_id, reference_short="26695",
                   contig=contig, start=s + 1, end=e, motif_observed_plus=mot[:len(mot)],
                   motif_canonical=canonical_motif(mot), motif_length=len(mot), tract_length_bp=e - s,
                   copy_number=(e - s) // len(mot), purity="1.0", interruption_flag="NO",
                   left_flank_40=left, right_flank_40=right, left_flank_sha256=sha256(left),
                   right_flank_sha256=sha256(right), uniquely_anchorable="YES", caller_eligible="YES",
                   ortholog_group=lid, ortholog_basis="TEST", evidence_class="STRONG_PV_CANDIDATE",
                   evidence_source="TEST_ONLY", functional_interpretation_justified="YES",
                   strand="NA", context="INTERGENIC_OTHER", cds_start="NA", cds_end="NA",
                   tract_offset_in_cds="NA", cds_id="NA", locus_tag="TEST_%s" % lid, gene_name="NA", product="TEST pseudo-ORF")
        if info:
            strand, cs, ce = info
            row.update(strand=strand, context="CDS", cds_start=cs, cds_end=ce, cds_id="TESTCDS_" + lid)
            row["tract_offset_in_cds"] = (s + 1 - cs + 1) if strand == "+" else (ce - e + 1)
        if kind == "reg":
            row.update(context="PROMOTER_OR_INTERGENIC_UPSTREAM", strand="+")
        if kind == "ro":
            row.update(evidence_class="REPEAT_ONLY", functional_interpretation_justified="NO")
        if kind == "inel":
            row.update(caller_eligible="NO reason=TEST_NOT_ELIGIBLE")
        if kind == "reg" and info is None:
            pass
        rows.append(row)
    write_tsv(out_tsv, rows, CATALOGUE_FIELDS)
    return rows
