"""Functional state (FS01-FS09) from a called tract allele by substituting it into the annotated CDS."""
from .common import revcomp

_B = "TCAG"
_AA = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
CODON = {a + b + c: _AA[i] for i, (a, b, c) in enumerate((x, y, z) for x in _B for y in _B for z in _B)}
REGULATORY = ("PROMOTER_OR_INTERGENIC_UPSTREAM", "INTERGENIC_OTHER")


def translate_orf_len(seq):
    """Number of amino acids from the first codon up to (excluding) the first in-frame stop.
    Returns (aa_count, stopped)."""
    n = 0
    for i in range(0, len(seq) - 2, 3):
        aa = CODON.get(seq[i:i + 3], "X")
        if aa == "*":
            return n, True
        n += 1
    return n, False


def _apply_variants(seq_part, part_start_1based, variants, anchor_len, m, side):
    """variants: list of 'L12:A>G' strings; only those of the given side; positions on plus strand."""
    if not variants:
        return seq_part
    s = list(seq_part)
    for v in variants:
        if v[0] != side:
            continue
        off, ch = v[1:].split(":")
        off = int(off)
        new = ch.split(">")[1]
        pos = (m.start - anchor_len + off) if side == "L" else (m.end + 1 + off)
        i = pos - part_start_1based
        if 0 <= i < len(s):
            s[i] = new
    return "".join(s)


def substituted_cds(m, ref, tract, params, variants=None, anchor_len=25):
    """Gene-orientation sequence: substituted CDS plus reference bases appended past the CDS end."""
    cs, ce = m.cds_start, m.cds_end
    left = ref[cs - 1:m.start - 1]
    right = ref[m.end:ce]
    left = _apply_variants(left, cs, variants, anchor_len, m, "L")
    right = _apply_variants(right, m.end + 1, variants, anchor_len, m, "R")
    ext = params.orf_downstream_extend
    if m.strand == "-":
        plus = ref[max(0, cs - 1 - ext):cs - 1] + left + tract + right
        return revcomp(plus)
    return left + tract + right + ref[ce:ce + ext]


def annotated_aa(m):
    return (m.cds_end - m.cds_start + 1) // 3 - 1


def allele_state(m, ref, tract, params, variants=None, known_relations=None, anchor_len=25):
    """Return (state, basis) for one allele tract (plus-strand string). Implements FS02-FS06, FS08."""
    if not m.caller_eligible:
        return "NOT_APPLICABLE", "FS08;not caller eligible (%s)" % m.eligible_reason
    if m.context == "RNA_GENE" or m.evidence_class == "REPEAT_ONLY" or m.justified == "NO":
        return "ALLELE_STATE", "FS06;evidence_class=%s;context=%s;interpretation_justified=%s" % (
            m.evidence_class, m.context, m.justified)
    if m.context in REGULATORY:
        copies = len(tract) / m.motif_len
        if known_relations:
            for (cmin, cmax, st) in known_relations.get(m.locus_id, []):
                if cmin <= copies <= cmax:
                    return st, "FS05;source-supported relation copies %g in [%g,%g]" % (copies, cmin, cmax)
        return "REGULATORY_STATE", "FS05;regulatory context %s;no source-supported count-to-expression relation;copies=%g" % (
            m.context, copies)
    if m.context != "CDS" or m.cds_start is None or m.cds_end is None or ref is None:
        return "AMBIGUOUS", "FS02;CDS annotation or reference sequence unavailable"
    delta = len(tract) - m.tract_len
    ann = annotated_aa(m)
    thr = params.cds_on_min_orf_fraction

    def evaluate(vars_):
        seq = substituted_cds(m, ref, tract, params, vars_, anchor_len)
        aa, stopped = translate_orf_len(seq)
        ratio = aa / ann if ann > 0 else 0.0
        if m.motif_len % 3 == 0:
            st = "ALLELE_STATE" if ratio >= thr else "AMBIGUOUS"
            rule = "FS03"
        elif ratio >= thr:
            st, rule = "ON", "FS02"
        elif delta % 3 != 0:
            st, rule = "OFF", "FS02"
        else:
            st, rule = "AMBIGUOUS", "FS02"
        return st, rule, aa, ratio

    st, rule, aa, ratio = evaluate(None)
    basis = "%s;orf_aa=%d;annotated_aa=%d;orf_fraction=%.3f;tract_delta_bp=%d;frameshift=%s" % (
        rule, aa, ann, ratio, delta, "YES" if delta % 3 else "NO")
    if rule == "FS03" and st == "AMBIGUOUS":
        basis += ";stop introduced, not verified in surrounding ORF"
    if rule == "FS02" and st == "AMBIGUOUS":
        basis += ";stop without frameshift consequence"
    if variants:
        st2, _, aa2, ratio2 = evaluate(variants)
        if st2 != st:
            return "AMBIGUOUS", "FS04;disruptive sample variant(s) in anchoring reads %s change state %s->%s (orf_fraction %.3f->%.3f);%s" % (
                ",".join(variants), st, st2, ratio, ratio2, basis)
        basis += ";sample variants in anchors %s do not change state" % ",".join(variants)
    basis += ";FS04:sample CDS outside spanning reads NOT assessed (reference CDS used)"
    return st, basis
