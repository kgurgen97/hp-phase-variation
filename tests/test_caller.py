"""Unit tests for the repeat-length caller.

Run with:
    python -m unittest discover -s tests -v

Set HP_PHASEVAR_TEST_FASTA to an H. pylori 26695 FASTA(.gz) to enable
the reference-dependent catalogue tests."""
import gzip
import os
import sys
import tempfile
import unittest
from pathlib import Path

HERE = Path(__file__).resolve()
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / "src"))

from hp_phasevar import common, extract, call, functional, noise, simulate, benchmark, testcat  # noqa: E402
from hp_phasevar.common import Params, Member, revcomp, canonical_motif  # noqa: E402

FASTA = os.environ.get("G5_TEST_FASTA") or str(
    ROOT / "data/reference/g4_identity/GCF_000008525.1_ASM852v1_genomic.fna.gz")
TMP = tempfile.mkdtemp(prefix="hp_phasevar_test_")
_CACHE = {}


def tiny():
    if "cat" not in _CACHE:
        path = os.path.join(TMP, "tiny.tsv")
        testcat.build_tiny_catalogue(FASTA, path)
        refs = common.read_fasta(FASTA)
        ms = common.load_catalogue(path)
        simulate.attach_ref_tracts(ms, refs)
        _CACHE["cat"] = (path, ms, refs)
    return _CACHE["cat"]


def write_fq(path, reads):
    with gzip.open(path, "wt") as fh:
        for i, (s, q) in enumerate(reads):
            fh.write("@r%d\n%s\n+\n%s\n" % (i, s, q if q else "I" * len(s)))


def row(m, tract, strand="+", q=35.0, lmm=0, rmm=0, var=".", src="single"):
    return {"kind": "SPAN", "member_id": m.member_id, "tract": tract, "strand": strand, "anchor_q": str(q),
            "l_mm": str(lmm), "r_mm": str(rmm), "variants": var, "source": src}


def mk_member(**kw):
    d = dict(locus_id="L1", member_id="L1.M1", contig="c", start=100, end=107, strand="+",
             motif_observed_plus="A", motif_canonical="A", motif_length=1, tract_length_bp=8, context="CDS",
             left_flank_40="C" * 40, right_flank_40="G" * 40, caller_eligible="YES", cds_start=1, cds_end=200,
             evidence_class="STRONG_PV_CANDIDATE", functional_interpretation_justified="YES",
             interruption_flag="NO")
    d.update(kw)
    return Member({k: str(v) for k, v in d.items()})


class TestBasics(unittest.TestCase):
    def test_revcomp(self):
        self.assertEqual(revcomp("ACGTN"), "NACGT")
        self.assertEqual(revcomp(revcomp("AACCGGTTA")), "AACCGGTTA")

    def test_motif_normalization(self):
        self.assertEqual(canonical_motif("CT"), "AG")
        self.assertEqual(canonical_motif("TC"), "AG")
        self.assertEqual(canonical_motif("GA"), "AG")
        self.assertEqual(canonical_motif("ATAT"), "AT")       # shortest period
        self.assertEqual(canonical_motif("TTT"), "A")
        self.assertEqual(canonical_motif("GGG"), "C")
        self.assertEqual(canonical_motif("TAG"), canonical_motif("CTA"))
        self.assertEqual(canonical_motif("TAG"), canonical_motif("CTA"[::-1].translate(str.maketrans("ACGT", "TGCA"))))

    def test_rules_tsv_matches_defaults(self):
        p = os.path.join(TMP, "rules.tsv")
        common.write_rules_tsv(p)
        self.assertEqual(dict(common.load_params(p)), dict(Params()))
        with self.assertRaises(KeyError):
            Params({"no_such_rule": 1})

    def test_binom_and_wilson(self):
        self.assertAlmostEqual(noise.binom_sf(1, 10, 0.5), 1 - 0.5 ** 10, places=9)
        self.assertAlmostEqual(noise.binom_sf(0, 10, 0.3), 1.0)
        self.assertLess(noise.binom_sf(20, 100, 0.02), 1e-10)
        lo, hi = noise.wilson(50, 100)
        self.assertTrue(0.39 < lo < 0.41 and 0.59 < hi < 0.61)


class TestRepeatCounting(unittest.TestCase):
    def setUp(self):
        self.m = mk_member()
        self.P = Params()

    def test_counts_and_filters(self):
        rows = ([row(self.m, "A" * 8)] * 5 + [row(self.m, "A" * 7)] * 3 + [row(self.m, "AAAATAAA")]      # 1 mm kept
                + [row(self.m, "AATATAAA")]           # 2 mm discarded as discordant
                + [row(self.m, "A" * 8, q=10.0)]      # low anchor quality
                + [row(self.m, "A" * 8, lmm=3, rmm=1)])  # too many anchor mismatches
        S = call.summarize_alleles(self.m, rows, self.P)
        self.assertEqual({L: len(v) for L, v in S["aligned"].items()}, {8: 6, 7: 3})  # 5 pure + 1 with 1 substitution
        self.assertEqual((S["lowq"], S["disc"]), (2, 1))

    def test_offunit_dinucleotide(self):
        m = mk_member(motif_observed_plus="CT", motif_canonical="AG", motif_length=2, tract_length_bp=10)
        rows = [row(m, "CT" * 5)] * 6 + [row(m, "CT" * 4)] * 2 + [row(m, "CT" * 4 + "C")] * 2
        S = call.summarize_alleles(m, rows, self.P)
        self.assertEqual({L: len(v) for L, v in S["aligned"].items()}, {10: 6, 8: 2})
        self.assertEqual(dict(S["off"]), {9: 2})


class TestAnchoring(unittest.TestCase):
    def setUp(self):
        self.P = Params()
        self.m = mk_member(left_flank_40="ACGTTGCAAGCTTGACCTAGGATCCGTAAGCTAGCTTGCA",
                           right_flank_40="GATCGGTACCTTAGCAAGCTTGCATCGATCGGCTAAGCTA", tract_length_bp=8)
        self.idx = extract.AnchorIndex([self.m], self.P)
        self.L, self.R = self.m.left40, self.m.right40

    def rd(self, tract, lpad="TTGCAGCATTCAGGC", rpad="ACGATTGGCATCCAG", lmut=(), rmut=()):
        l = list(self.L[-25:])
        for i in lmut:
            l[i] = "A" if l[i] != "A" else "C"
        r = list(self.R[:25])
        for i in rmut:
            r[i] = "A" if r[i] != "A" else "C"
        return lpad + "".join(l) + tract + "".join(r) + rpad

    def test_span_both_strands(self):
        s = self.rd("A" * 9)
        for seq, strand in ((s, "+"), (revcomp(s), "-")):
            res = extract.eval_read(self.idx, seq, "I" * len(seq))
            rec = res[0]
            self.assertEqual((rec[0], rec[1], rec[2]), ("SPAN", strand, "A" * 9))
            self.assertAlmostEqual(rec[3], 40.0)

    def test_mismatch_tolerance(self):
        ok = extract.eval_read(self.idx, self.rd("A" * 8, lmut=(2, 20), rmut=(5,)), "I" * 80)
        self.assertEqual(ok[0][0], "SPAN")
        self.assertEqual((ok[0][4], ok[0][5]), (2, 1))
        bad = extract.eval_read(self.idx, self.rd("A" * 8, lmut=(2, 10, 20)), "I" * 80)
        self.assertEqual(bad[0][0], "PART")
        self.assertEqual(bad[0][1], {"R"})

    def test_tract_reaching_read_end_is_partial(self):
        s = self.rd("A" * 8)
        cut = s[:s.index("A" * 8, 15 + 20) + 5]      # read ends inside the tract
        res = extract.eval_read(self.idx, cut, "I" * len(cut))
        self.assertEqual(res[0][0], "PART")
        self.assertEqual(res[0][1], {"L"})

    def test_unrelated_read(self):
        self.assertEqual(extract.eval_read(self.idx, "ACGT" * 30, "I" * 120), {})


class TestFunctionalState(unittest.TestCase):
    def build(self, strand, motif="A", tlen=6, motif_len=1):
        # CDS (gene orientation): ATG + AAA-tract in frame + blocks with stops in both shifted frames
        body = "ATG" + "GCT" * 3
        tract = motif * tlen
        pre = body
        post = ("GGTAGC" * 6) + ("CTAAGG" * 6) + "GCTGCT" + "TAA"
        gene = pre + tract + post
        assert len(gene) % 3 == 0
        ref_plus = "GGGGGGGGGG" + (gene if strand == "+" else revcomp(gene)) + "CCCCCCCCCCCCCCCCCCCC" * 40
        cs = 11
        ce = 10 + len(gene)
        if strand == "+":
            ts = cs + len(pre)
        else:
            ts = ce - len(pre) - len(tract) + 1
        te = ts + len(tract) - 1
        tplus = ref_plus[ts - 1:te]
        m = mk_member(contig="c", start=ts, end=te, strand=strand, cds_start=cs, cds_end=ce,
                      motif_observed_plus=tplus[:motif_len], motif_length=motif_len, tract_length_bp=len(tract),
                      left_flank_40=ref_plus[max(0, ts - 41):ts - 1].rjust(40, "N"), right_flank_40=ref_plus[te:te + 40])
        return m, ref_plus, tplus

    def test_frame_states_both_strands(self):
        P = Params()
        for strand in "+-":
            m, ref, tp = self.build(strand)
            base = tp
            unit = tp[0]
            self.assertEqual(functional.allele_state(m, ref, base, P)[0], "ON")
            self.assertEqual(functional.allele_state(m, ref, unit * (len(base) - 1), P)[0], "OFF")
            self.assertEqual(functional.allele_state(m, ref, unit * (len(base) + 1), P)[0], "OFF")
            self.assertEqual(functional.allele_state(m, ref, unit * (len(base) + 3), P)[0], "ON")
            self.assertEqual(functional.allele_state(m, ref, unit * (len(base) - 3), P)[0], "ON")
            st, basis = functional.allele_state(m, ref, base, P)
            self.assertIn("FS04", basis)
            self.assertIn("NOT assessed", basis)

    def test_motif_divisible_by_three_never_on_off(self):
        P = Params()
        m, ref, tp = self.build("+", motif="AAG", tlen=12, motif_len=3)
        for n in (2, 3, 5):
            st, basis = functional.allele_state(m, ref, tp[:3] * n, P)
            self.assertIn(st, ("ALLELE_STATE", "AMBIGUOUS"))
            self.assertNotIn(st, ("ON", "OFF"))

    def test_regulatory_repeat_only_rna_not_applicable(self):
        P = Params()
        self.assertEqual(functional.allele_state(mk_member(context="INTERGENIC_OTHER"), None, "A" * 9, P)[0], "REGULATORY_STATE")
        self.assertEqual(functional.allele_state(mk_member(context="PROMOTER_OR_INTERGENIC_UPSTREAM"), None, "A" * 9, P,
                                                 known_relations={"L1": [(10, 12, "ON")]})[0], "REGULATORY_STATE")
        self.assertEqual(functional.allele_state(mk_member(context="PROMOTER_OR_INTERGENIC_UPSTREAM"), None, "A" * 11, P,
                                                 known_relations={"L1": [(10, 12, "ON")]})[0], "ON")
        self.assertEqual(functional.allele_state(mk_member(evidence_class="REPEAT_ONLY"), None, "A" * 9, P)[0], "ALLELE_STATE")
        self.assertEqual(functional.allele_state(mk_member(context="RNA_GENE"), None, "A" * 9, P)[0], "ALLELE_STATE")
        self.assertEqual(functional.allele_state(mk_member(caller_eligible="NO reason=x"), None, "A" * 9, P)[0], "NOT_APPLICABLE")

    def test_disruptive_anchor_variant_gives_ambiguous(self):
        P = Params()
        m, ref, tp = self.build("+")
        gpos = m.end + 1            # first codon after the tract is GGT (in frame); edit it to TAA
        self.assertEqual(ref[gpos - 1:gpos + 2], "GGT")
        var = ["R0:G>T", "R1:G>A", "R2:T>A"]
        st, basis = functional.allele_state(m, ref, tp, P, variants=var)
        self.assertEqual(st, "AMBIGUOUS")
        self.assertIn("FS04", basis)


class TestCallerRules(unittest.TestCase):
    def setUp(self):
        self.P = Params()
        self.N = noise.NoiseModel(self.P)
        self.m = mk_member(context="INTERGENIC_OTHER", tract_length_bp=10)  # non-CDS: state logic aside
        self.m.tract_len = 10

    def go(self, allele_counts, strands=None, m=None, **kw):
        m = m or self.m
        rows = []
        for L, n in allele_counts.items():
            for i in range(n):
                s = "+" if (strands is None and i % 2 == 0) else ("-" if strands is None else strands(L, i))
                rows.append(row(m, "A" * L, strand=s))
        ev = {"span": rows, "PARTIAL_L": 4, "PARTIAL_R": 2, "CONFLICT": 0}
        return call.call_member(m, ev, "S", self.P, self.N, None, None)

    def test_pure_with_stutter_is_not_mixed(self):
        r = self.go({10: 190, 9: 5, 11: 2, 8: 1})
        self.assertEqual(r["callable"], "YES")
        self.assertEqual(r["dominant_bp"], 10)
        self.assertEqual(r["mixed_flag"], "NONE")
        self.assertEqual(r["confidence_class"], "HIGH")

    def test_never_mixed_at_or_below_noise(self):
        # 8% of reads one unit below on a 10x homopolymer: expected stutter (a+b*10=0.03)*safety 2=0.06 -> not MIXED
        r = self.go({10: 92, 9: 8})
        self.assertNotEqual(r["mixed_flag"], "MIXED")
        self.assertNotEqual(r["functional_state"], "MIXED")

    def test_mixed_alleles(self):
        r = self.go({10: 100, 13: 100})
        self.assertEqual(r["mixed_flag"], "MIXED")
        self.assertEqual(r["secondary_bp"], 13)
        r = self.go({10: 140, 12: 60})
        self.assertEqual(r["mixed_flag"], "MIXED")
        self.assertEqual(r["allele_counts"], "10bp/10c:140;12bp/12c:60")

    def test_minor_on_one_strand_only_is_candidate(self):
        r = self.go({10: 100, 12: 30}, strands=lambda L, i: "+" if L == 12 else ("+" if i % 2 else "-"))
        self.assertEqual(r["mixed_flag"], "MIXTURE_CANDIDATE")

    def test_low_depth_mixture_is_candidate_not_mixed(self):
        r = self.go({10: 9, 12: 4})
        self.assertEqual(r["callable"], "YES")
        self.assertEqual(r["mixed_flag"], "MIXTURE_CANDIDATE")

    def test_uncallable_cases(self):
        self.assertEqual(call.call_member(self.m, None, "S", self.P, self.N)["uncallable_reason"], "NO_SPANNING_READS")
        self.assertEqual(call.call_member(self.m, None, "S", self.P, self.N)["functional_state"], "UNCALLABLE")
        r = self.go({10: 5})
        self.assertTrue(r["uncallable_reason"].startswith("INSUFFICIENT_SPANNING_READS"))
        r = self.go({10: 5, 9: 4, 12: 4, 15: 5, 14: 4})
        self.assertTrue(r["uncallable_reason"].startswith(("NO_DOMINANT", "DIFFUSE")))
        ne = mk_member(caller_eligible="NO reason=not_unique", context="INTERGENIC_OTHER")
        r = call.call_member(ne, None, "S", self.P, self.N)
        self.assertEqual(r["functional_state"], "NOT_APPLICABLE")
        self.assertEqual(r["callable"], "NO")
        # non-unit-length reads dominate
        m2 = mk_member(motif_observed_plus="CT", motif_canonical="AG", motif_length=2, tract_length_bp=10,
                       context="INTERGENIC_OTHER")
        rows = [row(m2, "CT" * 5)] * 5 + [row(m2, "CT" * 4 + "C")] * 6
        r = call.call_member(m2, {"span": rows, "PARTIAL_L": 0, "PARTIAL_R": 0, "CONFLICT": 0}, "S", self.P, self.N)
        self.assertTrue(r["uncallable_reason"].startswith("EXCESS_NON_UNIT"))
        # low anchor quality reads are filtered out, leaving too few
        rows = [row(self.m, "A" * 10, q=12)] * 30
        r = call.call_member(self.m, {"span": rows, "PARTIAL_L": 0, "PARTIAL_R": 0, "CONFLICT": 0}, "S", self.P, self.N)
        self.assertEqual(r["uncallable_reason"], "NO_SPANNING_READS")

    def test_strand_bias_downgrades(self):
        r = self.go({10: 40}, strands=lambda L, i: "+")
        self.assertEqual(r["strand_flag"], "STRAND_BIASED")
        self.assertEqual(r["confidence_class"], "MEDIUM")
        p = Params({"strand_bias_action": "UNCALLABLE"})
        rows = [row(self.m, "A" * 10, strand="+")] * 40
        r = call.call_member(self.m, {"span": rows, "PARTIAL_L": 0, "PARTIAL_R": 0, "CONFLICT": 0}, "S", p, self.N)
        self.assertEqual(r["uncallable_reason"], "DOMINANT_ALLELE_STRAND_BIAS")

    def test_noise_estimation(self):
        obs = [{"mclass": "MONO", "copies": c, "n": 1000, "m1": int(1000 * (0.005 + 0.003 * c)), "p1": int(1000 * (0.005 + 0.003 * c) / 2),
                "m2": 2, "p2": 1, "far": 1, "offunit": 1, "n_all": 1002} for c in (8, 10, 12, 15)]
        est, rep = noise.estimate_noise(obs, Params())
        self.assertAlmostEqual(est["stutter_b_MONO"], 0.003, places=3)
        self.assertAlmostEqual(est["stutter_a_MONO"], 0.005, delta=0.001)
        self.assertAlmostEqual(est["stutter_up_ratio"], 0.5, places=1)
        self.assertNotIn("stutter_a_DI", est)     # too few observations: default retained


class TestExtractionInputs(unittest.TestCase):
    def setUp(self):
        self.P = Params()
        self.m = mk_member(left_flank_40="ACGTTGCAAGCTTGACCTAGGATCCGTAAGCTAGCTTGCA",
                           right_flank_40="GATCGGTACCTTAGCAAGCTTGCATCGATCGGCTAAGCTA", tract_length_bp=8)
        self.m.motif_len = 1

    def frag(self, tract, flen_extra=140):
        pad_l = "TTGCAGCATTCAGGCATCGATCAGGCTTACGATCGGATTTAGCCA"
        pad_r = "ACGATTGGCATCCAGGTTACGATTCGGATCCAGTTAGCTAGGCCT"
        return pad_l + self.m.left40[-25:] + tract + self.m.right40[:25] + pad_r

    def out_rows(self, path):
        with gzip.open(path, "rt") as fh:
            hdr = fh.readline().rstrip("\n").split("\t")
            return [dict(zip(hdr, l.rstrip("\n").split("\t"))) for l in fh]

    def test_se(self):
        f = self.frag("A" * 8)
        p = os.path.join(TMP, "se.fq.gz")
        write_fq(p, [(f, None), (revcomp(f), None)])
        out = os.path.join(TMP, "se.ev.gz")
        st = extract.extract_sample([self.m], "S1", out, se_files=[p], params=self.P)
        rows = self.out_rows(out)
        self.assertEqual(st["fragments"], 2)
        self.assertEqual(sorted(r["strand"] for r in rows if r["kind"] == "SPAN"), ["+", "-"])
        self.assertTrue(all(r["source"] == "single" for r in rows if r["kind"] == "SPAN"))

    def test_pe_overlapping_mates_not_double_counted_and_strand_convention(self):
        f = self.frag("A" * 8)                      # ~ 145 bp; mates of 120 bp overlap over the tract
        r1, r2 = f[:120], revcomp(f[-120:])
        p1, p2 = os.path.join(TMP, "r1.fq.gz"), os.path.join(TMP, "r2.fq.gz")
        write_fq(p1, [(r1, None)] * 3)
        write_fq(p2, [(r2, None)] * 3)
        out = os.path.join(TMP, "pe.ev.gz")
        extract.extract_sample([self.m], "S1", out, pe_pairs=[(p1, p2)], params=self.P)
        rows = [r for r in self.out_rows(out) if r["kind"] == "SPAN"]
        self.assertEqual(len(rows), 3)                            # one row per fragment, not two
        self.assertTrue(all(r["source"] == "mate1+mate2" and r["n"] == "2" for r in rows))
        self.assertTrue(all(r["strand"] == "+" for r in rows))

    def test_pe_single_mate_spans_strand_is_fragment_strand(self):
        f = self.frag("A" * 8)
        r1 = f[:40]                          # mate1 does not reach the tract (partial or nothing)
        r2 = revcomp(f[-125:])               # mate2 spans
        p1, p2 = os.path.join(TMP, "a1.fq.gz"), os.path.join(TMP, "a2.fq.gz")
        write_fq(p1, [(r1, None)])
        write_fq(p2, [(r2, None)])
        out = os.path.join(TMP, "pe2.ev.gz")
        extract.extract_sample([self.m], "S1", out, pe_pairs=[(p1, p2)], params=self.P)
        rows = [r for r in self.out_rows(out) if r["kind"] == "SPAN"]
        self.assertEqual(len(rows), 1)
        self.assertEqual((rows[0]["source"], rows[0]["strand"]), ("mate2", "+"))

    def test_pe_conflict_dropped(self):
        f1, f2 = self.frag("A" * 8), self.frag("A" * 9)
        r1, r2 = f1[:120], revcomp(f2[-120:])
        p1, p2 = os.path.join(TMP, "c1.fq.gz"), os.path.join(TMP, "c2.fq.gz")
        write_fq(p1, [(r1, None)])
        write_fq(p2, [(r2, None)])
        out = os.path.join(TMP, "pe3.ev.gz")
        st = extract.extract_sample([self.m], "S1", out, pe_pairs=[(p1, p2)], params=self.P)
        rows = self.out_rows(out)
        self.assertEqual(st["mate_conflict"], 1)
        self.assertEqual(len([r for r in rows if r["kind"] == "SPAN"]), 0)
        self.assertEqual([r["n"] for r in rows if r["kind"] == "CONFLICT"], ["1"])

    def test_orphans_add_to_same_sample(self):
        f = self.frag("A" * 8)
        se, o = os.path.join(TMP, "se2.fq.gz"), os.path.join(TMP, "orph.fq.gz")
        r1, r2 = f[:120], revcomp(f[-120:])
        p1, p2 = os.path.join(TMP, "d1.fq.gz"), os.path.join(TMP, "d2.fq.gz")
        write_fq(p1, [(r1, None)])
        write_fq(p2, [(r2, None)])
        write_fq(se, [(f, None)])
        write_fq(o, [(f, None), (self.frag("A" * 7), None)])
        out = os.path.join(TMP, "all.ev.gz")
        st = extract.extract_sample([self.m], "ONE", out, se_files=[se], pe_pairs=[(p1, p2)], orphan_files=[o], params=self.P)
        rows = [r for r in self.out_rows(out) if r["kind"] == "SPAN"]
        self.assertEqual({r["sample_id"] for r in rows}, {"ONE"})
        self.assertEqual(sorted(r["source"] for r in rows), ["mate1+mate2", "single", "unpaired", "unpaired"])
        ev, samples = call.read_evidence(out)
        self.assertEqual(samples, {"ONE"})
        self.assertEqual(len(ev[self.m.member_id]["span"]), 4)
        self.assertEqual(st["fragments"], 4)

    def test_partial_reads_counted(self):
        f = self.frag("A" * 8)
        left_only = f[:f.index("A" * 8) + 5]
        p = os.path.join(TMP, "part.fq.gz")
        write_fq(p, [(left_only, None)] * 3)
        out = os.path.join(TMP, "part.ev.gz")
        extract.extract_sample([self.m], "S1", out, se_files=[p], params=self.P)
        ev, _ = call.read_evidence(out)
        e = ev[self.m.member_id]
        self.assertEqual((len(e["span"]), e["PARTIAL_L"], e["PARTIAL_R"]), (0, 3, 0))


class TestOnTinyCatalogue(unittest.TestCase):
    @unittest.skipUnless(os.path.exists(FASTA), "reference FASTA not found")
    def test_catalogue_flanks_match_reference(self):
        path, ms, refs = tiny()
        self.assertGreaterEqual(len(ms), 8)
        for m in ms:
            g = refs[m.contig]
            self.assertEqual(g[m.start - 41:m.start - 1], m.left40)
            self.assertEqual(g[m.end:m.end + 40], m.right40)
            self.assertEqual(g[m.start - 1:m.end], m.ref_tract)

    @unittest.skipUnless(os.path.exists(FASTA), "reference FASTA not found")
    def test_smoke_roundtrip_se_pe_orphan(self):
        path, ms, refs = tiny()
        eligible = [m for m in ms if m.caller_eligible]
        for layout, orph in (("SE150", 0.0), ("PE250", 0.2)):
            spec = benchmark.build_alleles(ms, 0.0)
            pref = os.path.join(TMP, "sm_" + layout)
            files, truth = simulate.simulate_sample(ms, refs, spec, pref, layout, 60, 11, None, orph, 1.0, "T", "100:0")
            out = pref + ".ev.gz"
            extract.extract_sample(ms, "T", out, files.get("se", []), files.get("pe", []), files.get("orphan", []))
            rows = call.call_sample(ms, out, None, Params(), None, refs)
            tm = {t["member_id"]: t for t in truth}
            ok = 0
            for r in rows:
                m = next(x for x in ms if x.member_id == r["member_id"])
                if not m.caller_eligible:
                    self.assertEqual(r["functional_state"], "NOT_APPLICABLE")
                    continue
                self.assertEqual(r["callable"], "YES", r)
                self.assertEqual(str(r["dominant_bp"]), str(tm[r["member_id"]]["major_bp"]))
                self.assertEqual(r["mixed_flag"], "NONE")
                ok += 1
            self.assertEqual(ok, len(eligible))

    @unittest.skipUnless(os.path.exists(FASTA), "reference FASTA not found")
    def test_smoke_mixture_detected(self):
        path, ms, refs = tiny()
        spec = benchmark.build_alleles(ms, 0.5)
        pref = os.path.join(TMP, "mix")
        files, truth = simulate.simulate_sample(ms, refs, spec, pref, "SE150", 200, 5, None, 0.0, 1.0, "M", "50:50")
        out = pref + ".ev.gz"
        extract.extract_sample(ms, "M", out, files["se"])
        rows = call.call_sample(ms, out, None, Params(), None, refs)
        flags = [r["mixed_flag"] for r in rows if r["callable"] == "YES"]
        self.assertTrue(all(f == "MIXED" for f in flags), flags)


if __name__ == "__main__":
    unittest.main()
