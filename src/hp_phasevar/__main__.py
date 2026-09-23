"""CLI: python -m hp_phasevar {extract,call,simulate,benchmark,estimate-noise,make-rules,make-test-catalogue}"""
import argparse
import json
import sys

from .common import load_params, load_catalogue, read_fasta, write_rules_tsv, read_tsv
from . import extract as ex, call as cl


def main(argv=None):
    ap = argparse.ArgumentParser(prog="hp_phasevar")
    sub = ap.add_subparsers(dest="cmd", required=True)

    p = sub.add_parser("extract", help="stage 1: reads -> evidence TSV (one sample = SE + PE + orphan files)")
    p.add_argument("--catalogue", required=True)
    p.add_argument("--sample", required=True)
    p.add_argument("--se", nargs="*", default=[])
    p.add_argument("--pe", nargs=2, action="append", default=[], metavar=("R1", "R2"))
    p.add_argument("--orphan", nargs="*", default=[])
    p.add_argument("--rules")
    p.add_argument("--out", required=True)
    p.add_argument("--max-fragments", type=int)

    p = sub.add_parser("call", help="stage 2: evidence -> per sample x locus report")
    p.add_argument("--catalogue", required=True)
    p.add_argument("--evidence", required=True)
    p.add_argument("--reference", required=True, help="reference FASTA (for coding-frame state)")
    p.add_argument("--rules")
    p.add_argument("--known-relations", help="optional TSV locus_id,copy_min,copy_max,state (FS05)")
    p.add_argument("--out", required=True)

    p = sub.add_parser("simulate", help="simulate reads + truth for one sample")
    p.add_argument("--catalogue", required=True)
    p.add_argument("--reference", required=True)
    p.add_argument("--prefix", required=True)
    p.add_argument("--layout", default="SE150")
    p.add_argument("--depth", type=float, default=30)
    p.add_argument("--mix", default="100:0")
    p.add_argument("--seed", type=int, default=1)
    p.add_argument("--stutter", action="store_true")
    p.add_argument("--orphan-fraction", type=float, default=0.0)

    p = sub.add_parser("benchmark", help="simulate -> extract -> call -> score")
    p.add_argument("--catalogue", required=True)
    p.add_argument("--reference", required=True)
    p.add_argument("--outdir", required=True)
    p.add_argument("--layouts", nargs="*", default=["SE150", "PE250"])
    p.add_argument("--depths", nargs="*", type=int, default=[10, 30, 100])
    p.add_argument("--mixes", nargs="*", default=["100:0", "90:10", "70:30", "50:50"])
    p.add_argument("--no-stutter", action="store_true")
    p.add_argument("--rules")
    p.add_argument("--noise-rules", help="TSV param_id,value overriding stutter parameters (from estimate-noise)")

    p = sub.add_parser("estimate-noise", help="estimate stutter/background parameters from pure-allele evidence")
    p.add_argument("--catalogue", required=True)
    p.add_argument("--evidence", nargs="+", required=True)
    p.add_argument("--rules")
    p.add_argument("--pure-members", help="file with one member_id per line known to be pure (synthetic truth)")
    p.add_argument("--min-dominant-fraction", type=float, default=0.95)
    p.add_argument("--out", required=True)

    p = sub.add_parser("make-rules", help="write default caller_rules_draft.tsv")
    p.add_argument("--out", required=True)

    p = sub.add_parser("make-test-catalogue", help="build the tiny TEST catalogue from a FASTA")
    p.add_argument("--reference", required=True)
    p.add_argument("--out", required=True)

    a = ap.parse_args(argv)
    if a.cmd == "make-rules":
        write_rules_tsv(a.out)
    elif a.cmd == "make-test-catalogue":
        from .testcat import build_tiny_catalogue
        rows = build_tiny_catalogue(a.reference, a.out)
        print("wrote %d loci" % len(rows))
    elif a.cmd == "extract":
        st = ex.extract_sample(a.catalogue, a.sample, a.out, a.se, [tuple(x) for x in a.pe], a.orphan,
                               load_params(a.rules), a.max_fragments)
        print(json.dumps(st))
    elif a.cmd == "call":
        kr = None
        if a.known_relations:
            kr = {}
            for r in read_tsv(a.known_relations):
                kr.setdefault(r["locus_id"], []).append((float(r["copy_min"]), float(r["copy_max"]), r["state"]))
        rows = cl.call_sample(a.catalogue, a.evidence, a.out, load_params(a.rules), None, a.reference, kr)
        from collections import Counter
        print(dict(Counter(r["functional_state"] for r in rows)))
    elif a.cmd == "simulate":
        from .simulate import simulate_sample, StutterSpec, TRUTH_FIELDS, attach_ref_tracts
        from .benchmark import build_alleles, MIXES
        from .common import write_tsv
        members = load_catalogue(a.catalogue)
        refs = read_fasta(a.reference)
        attach_ref_tracts(members, refs)
        files, truth = simulate_sample(members, refs, build_alleles(members, MIXES[a.mix]), a.prefix, a.layout,
                                       a.depth, a.seed, StutterSpec() if a.stutter else None,
                                       a.orphan_fraction, 1.0, "SIM", a.mix)
        write_tsv(a.prefix + ".truth.tsv", truth, TRUTH_FIELDS)
        print(json.dumps(files))
    elif a.cmd == "benchmark":
        from .benchmark import run_benchmark
        from .noise import NoiseModel
        params = load_params(a.rules)
        if a.noise_rules:
            for r in read_tsv(a.noise_rules):
                params.set(r["param_id"], r["value"])
        per, summ = run_benchmark(a.catalogue, a.reference, a.outdir, a.layouts, a.depths, a.mixes,
                                  not a.no_stutter, params=params)
        for r in summ:
            print("\t".join(str(v) for v in r.values()))
    elif a.cmd == "estimate-noise":
        from .noise import estimate_noise
        params = load_params(a.rules)
        members = load_catalogue(a.catalogue)
        pure = None
        if a.pure_members:
            pure = set(l.strip() for l in open(a.pure_members) if l.strip())
        obs = []
        for e in a.evidence:
            obs += cl.pure_observations(members, e, params, a.min_dominant_fraction, pure)
        est, report = estimate_noise(obs, params)
        with open(a.out, "w") as fh:
            fh.write("param_id\tvalue\tdescription\n")
            for k, v in est.items():
                fh.write("%s\t%s\testimated from %d pure-allele observations\n" % (k, v, len(obs)))
        print("\n".join(report))


if __name__ == "__main__":
    sys.exit(main())
