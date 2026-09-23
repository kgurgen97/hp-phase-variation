"""Benchmark driver: simulate -> extract -> call -> compare with truth. Use on a catalogue of your choice
(smoke: tiny test catalogue; the main process runs the real catalogue after merge)."""
import os
import time
from collections import defaultdict

from .call import call_sample, REPORT_FIELDS
from .common import Params, load_catalogue, read_fasta, write_tsv
from .extract import extract_sample
from .functional import allele_state
from .noise import NoiseModel
from .simulate import simulate_sample, attach_ref_tracts, StutterSpec, TRUTH_FIELDS, make_tract

MIXES = {"100:0": 0.0, "90:10": 0.10, "70:30": 0.30, "50:50": 0.50}
MINOR_DELTAS = [1, -1, 2, -3, 3, -2]      # rotated across loci (units of motif)
MAJOR_DELTAS = [0, -1, 1, -2, 2, 3]       # copy-number states around frame transitions


def build_alleles(members, mix_frac, rot=0, only_ids=None):
    spec = {}
    for i, m in enumerate(members):
        if not m.caller_eligible or (only_ids is not None and m.member_id not in only_ids):
            continue
        dom = MAJOR_DELTAS[(i + rot) % len(MAJOR_DELTAS)]
        while m.ref_copies + dom < 3:
            dom += 1
        if mix_frac <= 0:
            spec[m.member_id] = [(dom, 1.0)]
        else:
            mn = MINOR_DELTAS[(i + rot) % len(MINOR_DELTAS)]
            while m.ref_copies + dom + mn < 3 or mn == 0:
                mn += 1
            spec[m.member_id] = [(dom, 1.0 - mix_frac), (dom + mn, mix_frac)]
    return spec


def run_benchmark(catalogue, fasta, outdir, layouts=("SE150", "PE250"), depths=(10, 30, 100),
                  mixes=("100:0", "90:10", "70:30", "50:50"), stutter=True, orphan_fraction=0.1, seed=7,
                  params=None, noise=None, err_scale=1.0, only_ids=None, stutter_spec=None):
    """Returns (per_locus_rows, summary_rows); writes benchmark_per_locus.tsv, benchmark_summary.tsv and truth."""
    os.makedirs(outdir, exist_ok=True)
    params = params or Params()
    members = load_catalogue(catalogue) if isinstance(catalogue, str) else catalogue
    refs = read_fasta(fasta)
    attach_ref_tracts(members, refs)
    by_mid = {m.member_id: m for m in members}
    per, truths = [], []
    t0 = time.time()
    for li, layout in enumerate(layouts):
        for depth in depths:
            for mix in mixes:
                sid = "SIM_%s_%dx_%s_%s" % (layout, depth, mix.replace(":", "-"), "stut" if stutter else "nostut")
                spec = build_alleles(members, MIXES[mix], rot=li, only_ids=only_ids)
                pref = os.path.join(outdir, "reads_" + sid)
                s = seed + len(per) + depth
                files, truth = simulate_sample([x for x in members if x.member_id in spec], refs, spec, pref, layout, depth, s,
                                               (stutter_spec or StutterSpec()) if stutter else None,
                                               orphan_fraction if layout.startswith("PE") else 0.0,
                                               err_scale, sid, mix)
                ev = os.path.join(outdir, "evidence_%s.tsv.gz" % sid)
                extract_sample(members, sid, ev, files.get("se", []), files.get("pe", []),
                               files.get("orphan", []), params)
                rows = call_sample(members, ev, None, params, noise, refs)
                rmap = {r["member_id"]: r for r in rows}
                for t in truth:
                    if t['member_id'] not in spec: continue
                    m = by_mid[t["member_id"]]
                    r = rmap[t["member_id"]]
                    sp = spec.get(m.member_id)
                    if sp is None:
                        continue
                    mf = max(f for _, f in sp)
                    majors = [m.tract_len + d * m.motif_len for d, f in sp if f >= mf - 1e-9]
                    minor = [(d, f) for d, f in sp if f < 0.5 - 1e-9]
                    mixed_truth = bool(minor) and len(sp) > 1
                    tstate, _ = allele_state(m, refs[m.contig], make_tract(m, max(sp, key=lambda x: x[1])[0]), params)
                    per.append({"sample_id": sid, "layout": layout, "depth": depth, "mix": mix, "member_id": m.member_id,
                                "motif_len": m.motif_len, "context": m.context, "callable": r["callable"],
                                "uncallable_reason": r["uncallable_reason"], "true_major_bp": "/".join(str(x) for x in majors),
                                "called_dominant_bp": r["dominant_bp"], "dominant_correct":
                                "YES" if str(r["dominant_bp"]) in [str(x) for x in majors] else "NO",
                                "mixed_truth": "YES" if len(sp) > 1 else "NO", "mixed_flag": r["mixed_flag"],
                                "true_state_major": tstate, "called_state": r["functional_state"],
                                "confidence": r["confidence_class"], "spanning_reads": r["spanning_reads"]})
                truths.extend(truth)
                write_tsv(os.path.join(outdir, "truth_%s.tsv" % sid), truth, TRUTH_FIELDS)
                rmfiles = [f for f in os.listdir(outdir) if f.startswith("reads_" + sid)]
                for f in rmfiles:
                    os.remove(os.path.join(outdir, f))      # keep only evidence + truth (reads are regenerable by seed)
    fields = list(per[0].keys())
    write_tsv(os.path.join(outdir, "benchmark_per_locus.tsv"), per, fields)
    write_tsv(os.path.join(outdir, "benchmark_truth_all.tsv"), truths, TRUTH_FIELDS)
    summ = summarize(per)
    write_tsv(os.path.join(outdir, "benchmark_summary.tsv"), summ, list(summ[0].keys()))
    return per, summ


def summarize(per):
    g = defaultdict(list)
    for r in per:
        g[(r["layout"], r["depth"], r["mix"])].append(r)
    out = []
    for (lay, dep, mix), rs in sorted(g.items()):
        n = len(rs)
        call = [r for r in rs if r["callable"] == "YES"]
        correct = [r for r in call if r["dominant_correct"] == "YES"]
        mixed_truth = [r for r in call if r["mixed_truth"] == "YES"]
        pure = [r for r in call if r["mixed_truth"] == "NO"]
        st = [r for r in call if r["true_state_major"] in ("ON", "OFF")]
        out.append({"layout": lay, "depth": dep, "mix": mix, "loci": n, "callable": len(call),
                    "dominant_correct_of_callable": "%d/%d" % (len(correct), len(call)),
                    "minor_flag_MIXED": sum(1 for r in mixed_truth if r["mixed_flag"] == "MIXED"),
                    "minor_flag_CANDIDATE": sum(1 for r in mixed_truth if r["mixed_flag"] == "MIXTURE_CANDIDATE"),
                    "mixed_truth_callable": len(mixed_truth),
                    "false_MIXED_in_pure": sum(1 for r in pure if r["mixed_flag"] == "MIXED"),
                    "false_CANDIDATE_in_pure": sum(1 for r in pure if r["mixed_flag"] == "MIXTURE_CANDIDATE"),
                    "pure_callable": len(pure),
                    "ON_OFF_state_concordant_of_truth_ON_OFF": "%d/%d" % (
                        sum(1 for r in st if r["called_state"] == r["true_state_major"] or
                            (r["called_state"] in ("MIXED",) and mixed_truth)), len(st))})
    return out
