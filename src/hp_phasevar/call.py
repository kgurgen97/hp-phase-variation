"""Stage 2: apply caller rules to cached evidence; per sample x locus report."""
import csv
import gzip
from collections import Counter, defaultdict

from .common import (Params, load_catalogue, motif_class, read_fasta, write_tsv)
from .functional import allele_state, REGULATORY
from .noise import NoiseModel, binom_sf, wilson

REPORT_FIELDS = [
    "sample_id", "locus_id", "member_id", "callable", "uncallable_reason", "spanning_reads", "partial_left",
    "partial_right", "mate_conflicts", "filtered_low_quality", "filtered_discordant_tract", "offunit_reads",
    "ref_copies", "allele_counts", "offunit_counts", "dominant_bp", "dominant_copies", "dominant_fraction",
    "dominant_ci_low", "dominant_ci_high", "secondary_bp", "secondary_copies", "secondary_fraction",
    "strand_plus", "strand_minus", "dominant_strand_plus", "dominant_strand_minus", "anchor_q_mean",
    "repeat_length_min", "repeat_length_max", "dominant_purity", "mixed_flag", "mixed_detail",
    "strand_flag", "confidence_class", "uncertainty_notes", "functional_state", "functional_state_basis",
    "secondary_state"]


def read_evidence(path):
    """Return {member_id: {'span': [rows], 'PARTIAL_L': n, 'PARTIAL_R': n, 'CONFLICT': n}} and sample ids."""
    ev = defaultdict(lambda: {"span": [], "PARTIAL_L": 0, "PARTIAL_R": 0, "CONFLICT": 0})
    samples = set()
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            samples.add(r["sample_id"])
            e = ev[r["member_id"]]
            if r["kind"] == "SPAN":
                e["span"].append(r)
            else:
                e[r["kind"]] += int(r["n"])
    return ev, samples


def _fmt(x):
    return ("%g" % x) if isinstance(x, float) else str(x)


def summarize_alleles(m, rows, params):
    """Filter reads and build allele distribution. Returns dict."""
    ref_len, ml = m.tract_len, m.motif_len
    lowq = disc = 0
    aligned = defaultdict(list)   # length -> list of rows
    off = Counter()
    for r in rows:
        if float(r["anchor_q"]) < params.min_anchor_mean_q or int(r["l_mm"]) + int(r["r_mm"]) > params.max_total_anchor_mismatch:
            lowq += 1
            continue
        t = r["tract"]
        L = len(t)
        if (L - ref_len) % ml != 0:
            off[L] += 1
            continue
        if not m.interrupted and L > 0:
            exp = (m.motif_plus * (L // ml + 1))[:L]
            mm = sum(1 for a, b in zip(t, exp) if a != b)
            if mm > params.max_tract_mismatch:
                disc += 1
                continue
        aligned[L].append(r)
    return {"aligned": aligned, "off": off, "lowq": lowq, "disc": disc}


def _allele_consensus(rows, m):
    c = Counter(r["tract"] for r in rows)
    return c.most_common(1)[0][0], c


def call_member(m, ev, sample_id, params, noise, ref_seqs=None, known_relations=None):
    rep = {f: "" for f in REPORT_FIELDS}
    rep.update(sample_id=sample_id, locus_id=m.locus_id, member_id=m.member_id, ref_copies=_fmt(m.ref_copies))
    ev = ev or {"span": [], "PARTIAL_L": 0, "PARTIAL_R": 0, "CONFLICT": 0}
    rep.update(partial_left=ev["PARTIAL_L"], partial_right=ev["PARTIAL_R"], mate_conflicts=ev["CONFLICT"])

    def uncallable(reason, state="UNCALLABLE", basis=None):
        rep.update({"callable": "NO", "uncallable_reason": reason, "functional_state": state,
                    "functional_state_basis": basis or ("FS01;" + reason if state == "UNCALLABLE" else "FS08;" + reason),
                    "mixed_flag": "NA", "confidence_class": "NA"})
        return rep

    if not m.caller_eligible:
        return uncallable("NOT_CALLER_ELIGIBLE:" + m.eligible_reason, "NOT_APPLICABLE")
    if m.tract_len > params.max_tract_len_call:
        return uncallable("TRACT_TOO_LONG_FOR_CALLING", "NOT_APPLICABLE")
    S = summarize_alleles(m, ev["span"], params)
    aligned, off = S["aligned"], S["off"]
    n_al = sum(len(v) for v in aligned.values())
    n_off = sum(off.values())
    n_span = n_al + n_off
    rep.update(spanning_reads=n_span, filtered_low_quality=S["lowq"], filtered_discordant_tract=S["disc"],
               offunit_reads=n_off)
    rep["allele_counts"] = ";".join("%dbp/%gc:%d" % (L, L / m.motif_len, len(aligned[L])) for L in sorted(aligned))
    rep["offunit_counts"] = ";".join("%dbp:%d" % (L, off[L]) for L in sorted(off))
    if n_span == 0:
        return uncallable("NO_SPANNING_READS")
    allrows = [r for v in aligned.values() for r in v]
    if n_span < params.min_spanning_reads:
        return uncallable("INSUFFICIENT_SPANNING_READS(%d<%d)" % (n_span, params.min_spanning_reads))
    if n_off / n_span > params.max_offunit_fraction:
        return uncallable("EXCESS_NON_UNIT_LENGTH_READS(%.2f)" % (n_off / n_span))
    if n_al < params.min_spanning_reads:
        return uncallable("INSUFFICIENT_UNIT_ALIGNED_READS(%d<%d)" % (n_al, params.min_spanning_reads))
    order = sorted(aligned, key=lambda L: (-len(aligned[L]), abs(L - m.tract_len), L))
    dom = order[0]
    nd = len(aligned[dom])
    if nd / n_al < params.min_dominant_fraction:
        return uncallable("NO_DOMINANT_ALLELE(%.2f)" % (nd / n_al))
    top2 = nd + (len(aligned[order[1]]) if len(order) > 1 else 0)
    if top2 / n_al < params.min_top2_fraction:
        return uncallable("DIFFUSE_ALLELE_DISTRIBUTION(top2=%.2f)" % (top2 / n_al))

    plus = sum(1 for r in allrows if r["strand"] == "+")
    minus = n_al - plus
    dplus = sum(1 for r in aligned[dom] if r["strand"] == "+")
    dminus = nd - dplus
    lo, hi = wilson(nd, n_al)
    mc = motif_class(m.motif_len)
    dom_copies = dom / m.motif_len
    dom_tract, dom_seqs = _allele_consensus(aligned[dom], m)
    pure_expect = (m.motif_plus * (dom // m.motif_len + 1))[:dom]
    rep.update(
        dominant_bp=dom, dominant_copies=_fmt(dom_copies), dominant_fraction="%.4f" % (nd / n_al),
        dominant_ci_low="%.4f" % lo, dominant_ci_high="%.4f" % hi, strand_plus=plus, strand_minus=minus,
        dominant_strand_plus=dplus, dominant_strand_minus=dminus,
        anchor_q_mean="%.2f" % (sum(float(r["anchor_q"]) for r in allrows) / n_al),
        repeat_length_min=min(aligned), repeat_length_max=max(aligned),
        dominant_purity="PURE" if dom_tract == pure_expect else "INTERRUPTED_OR_IMPURE")

    # ---- mixture test against the noise model
    mixed_hits, cand_hits, detail = [], [], []
    for L in order[1:]:
        c = len(aligned[L])
        k = (L - dom) // m.motif_len
        p = min(0.99, noise.rate(mc, k, dom_copies) * params.noise_safety_factor)
        pval = binom_sf(c, n_al, p)
        frac = c / n_al
        sp = sum(1 for r in aligned[L] if r["strand"] == "+") / c
        strand_ok = c >= 2 and min(sp, 1 - sp) >= params.mixed_min_minor_strand_frac
        detail.append("%dbp:n=%d,frac=%.3f,expected_noise=%.4f,p=%.2g,plus_frac=%.2f" % (L, c, frac, p, pval, sp))
        full = (frac >= params.mixed_min_minor_fraction and c >= params.mixed_min_minor_reads
                and n_al >= params.mixed_min_total_reads and pval < params.mixed_alpha and strand_ok)
        if full:
            mixed_hits.append(L)
        elif (frac >= params.candidate_min_fraction and c >= params.candidate_min_reads) or \
                (pval < params.mixed_candidate_alpha and c >= 2):
            cand_hits.append(L)
    sec = None
    if mixed_hits:
        rep["mixed_flag"] = "MIXED"
        sec = mixed_hits[0]
    elif cand_hits:
        rep["mixed_flag"] = "MIXTURE_CANDIDATE"
        sec = cand_hits[0]
    else:
        rep["mixed_flag"] = "NONE"
        sec = order[1] if len(order) > 1 else None
    rep["mixed_detail"] = "|".join(detail) if detail else "."
    if sec is not None:
        rep.update(secondary_bp=sec, secondary_copies=_fmt(sec / m.motif_len),
                   secondary_fraction="%.4f" % (len(aligned[sec]) / n_al))

    # ---- strand bias
    notes = []
    sflag = "OK"
    if nd >= params.strand_bias_min_reads and min(dplus, dminus) / nd < params.strand_bias_min_frac:
        sflag = "STRAND_BIASED"
        notes.append("dominant allele strand bias (%d+/%d-)" % (dplus, dminus))
        if params.strand_bias_action == "UNCALLABLE":
            return uncallable("DOMINANT_ALLELE_STRAND_BIAS")
        if params.strand_bias_action == "IGNORE":
            sflag = "STRAND_BIASED_IGNORED"
    rep["strand_flag"] = sflag

    # ---- confidence
    level = 2 if n_span >= params.conf_high_min_reads else (1 if n_span >= params.conf_medium_min_reads else 0)
    if sflag == "STRAND_BIASED":
        level -= 1
        notes.append("downgraded: strand bias")
    if rep["mixed_flag"] == "NONE" and lo < params.conf_min_dominant_ci_low:
        level -= 1
        notes.append("downgraded: dominant fraction CI low %.2f" % lo)
    if rep["mixed_flag"] == "MIXTURE_CANDIDATE":
        notes.append("minor allele consistent with or not distinguishable from noise: MIXTURE_CANDIDATE only")
    rep["confidence_class"] = ["LOW", "MEDIUM", "HIGH"][max(0, level)]
    rep["uncertainty_notes"] = ";".join(notes) if notes else "."
    rep.update({"callable": "YES", "uncallable_reason": "."})

    # ---- functional state (dominant, secondary; anchor variants seen in >=fraction of dominant reads)
    ref = (ref_seqs or {}).get(m.contig)
    var_counts = Counter()
    for r in aligned[dom]:
        if r["variants"] != ".":
            for v in set(r["variants"].split(",")):
                var_counts[v] += 1
    variants = [v for v, c in var_counts.items()
                if c >= params.variant_min_reads and c / nd >= params.variant_min_fraction]
    dstate, dbasis = allele_state(m, ref, dom_tract, params, variants, known_relations)
    state, basis = dstate, dbasis
    if sec is not None and rep["mixed_flag"] in ("MIXED", "MIXTURE_CANDIDATE"):
        st_tract, _ = _allele_consensus(aligned[sec], m)
        sstate, sbasis = allele_state(m, ref, st_tract, params, variants, known_relations)
        rep["secondary_state"] = sstate
        if rep["mixed_flag"] == "MIXED":
            if dstate in ("ON", "OFF") and sstate in ("ON", "OFF") and dstate != sstate:
                state = "MIXED"
                basis = "FS07;dominant %s [%s] vs secondary %s [%s]" % (dstate, dbasis, sstate, sbasis)
            elif m.context == "CDS" and dstate in ("ON", "OFF", "ALLELE_STATE", "AMBIGUOUS"):
                basis = "FS07;mixed flag but states not both defined and different (secondary %s); dominant state reported;%s" % (sstate, dbasis)
            elif dstate == "REGULATORY_STATE":
                state = "ALLELE_STATE"
                basis = "FS07;mixed non-coding regulatory locus reported as ALLELE_STATE;" + dbasis
        else:
            basis = "FS07;MIXTURE_CANDIDATE only (not MIXED); dominant state reported;" + dbasis
    if state in ("ON", "OFF", "MIXED") and rep["confidence_class"] == "LOW":
        basis += ";confidence LOW"
    rep["functional_state"] = state
    rep["functional_state_basis"] = basis
    return rep


def call_sample(catalogue, evidence_path, out_path, params=None, noise=None, ref_fasta=None,
                known_relations=None, sample_id=None):
    params = params or Params()
    noise = noise or NoiseModel(params)
    members = catalogue if isinstance(catalogue, list) else load_catalogue(catalogue)
    ev, samples = read_evidence(evidence_path)
    if sample_id is None:
        if len(samples) > 1:
            raise ValueError("evidence contains several sample_ids: %s" % sorted(samples))
        sample_id = next(iter(samples), "NA")
    refs = read_fasta(ref_fasta) if isinstance(ref_fasta, str) else ref_fasta
    rows = [call_member(m, ev.get(m.member_id), sample_id, params, noise, refs, known_relations) for m in members]
    if out_path:
        write_tsv(out_path, rows, REPORT_FIELDS)
    return rows


def pure_observations(catalogue, evidence_path, params, min_dominant_fraction=0.0, truth_pure=None):
    """Observations for noise estimation: one per member with enough unit-aligned reads. If truth_pure
    (set of member_ids known pure, e.g. synthetic) is given only those are used; else members whose
    dominant fraction >= min_dominant_fraction (biases stutter downward; document when used)."""
    members = catalogue if isinstance(catalogue, list) else load_catalogue(catalogue)
    ev, _ = read_evidence(evidence_path)
    obs = []
    for m in members:
        e = ev.get(m.member_id)
        if not e or not m.caller_eligible:
            continue
        if truth_pure is not None and m.member_id not in truth_pure:
            continue
        S = summarize_alleles(m, e["span"], params)
        al = S["aligned"]
        n = sum(len(v) for v in al.values())
        noff = sum(S["off"].values())
        if n < params.min_spanning_reads:
            continue
        dom = max(al, key=lambda L: len(al[L]))
        if truth_pure is None and len(al[dom]) / n < min_dominant_fraction:
            continue
        cnt = lambda k: len(al.get(dom + k * m.motif_len, []))
        near = sum(cnt(k) for k in (-2, -1, 1, 2))
        obs.append({"member_id": m.member_id, "mclass": motif_class(m.motif_len), "copies": dom / m.motif_len,
                    "n": n, "m1": cnt(-1), "p1": cnt(1), "m2": cnt(-2), "p2": cnt(2),
                    "far": n - len(al[dom]) - near, "offunit": noff, "n_all": n + noff})
    return obs
