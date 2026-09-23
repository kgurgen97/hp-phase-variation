"""Sequence utilities, catalogue/FASTA/FASTQ I/O and the rules-parameter table."""
import csv
import gzip
import hashlib
import io

_COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def revcomp(s):
    return s.translate(_COMP)[::-1]


def open_text(path, mode="rt"):
    if str(path).endswith(".gz"):
        return gzip.open(path, mode)
    return open(path, mode)


def rotations(s):
    return [s[i:] + s[:i] for i in range(len(s))]


def shortest_period(motif):
    n = len(motif)
    for p in range(1, n + 1):
        if n % p == 0 and motif[:p] * (n // p) == motif:
            return motif[:p]
    return motif


def canonical_motif(motif):
    """Lexicographically smallest rotation over both strands, reduced to shortest period (RD05)."""
    motif = shortest_period(motif.upper())
    cands = rotations(motif) + rotations(revcomp(motif))
    return min(cands)


def motif_class(motif_len):
    if motif_len == 1:
        return "MONO"
    if motif_len == 2:
        return "DI"
    return "POLY"


def sha256(s):
    return hashlib.sha256(s.encode()).hexdigest()


def read_fasta(path):
    seqs, name, buf = {}, None, []
    with open_text(path) as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                if name:
                    seqs[name] = "".join(buf).upper()
                name, buf = line[1:].split()[0], []
            else:
                buf.append(line)
    if name:
        seqs[name] = "".join(buf).upper()
    return seqs


def read_tsv(path):
    with open_text(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def write_tsv(path, rows, fields):
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "wt", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(fields)
        for r in rows:
            w.writerow([r.get(f, "") for f in fields])


def iter_fastq(path):
    """Yield (name, seq, qual_string). Phred+33."""
    with open_text(path) as fh:
        while True:
            h = fh.readline()
            if not h:
                return
            s = fh.readline().rstrip()
            fh.readline()
            q = fh.readline().rstrip()
            yield h[1:].split()[0], s.upper(), q


class Member:
    """One catalogue member with parsed numeric fields."""
    __slots__ = ("row", "locus_id", "member_id", "contig", "start", "end", "strand", "motif_plus",
                 "motif_canonical", "motif_len", "tract_len", "context", "left40", "right40",
                 "caller_eligible", "eligible_reason", "interrupted", "cds_start", "cds_end",
                 "evidence_class", "justified", "ref_tract", "gene_name", "locus_tag")

    def __init__(self, row):
        self.row = row
        g = row.get
        self.locus_id = g("locus_id")
        self.member_id = g("member_id")
        self.contig = g("contig")
        self.start = int(g("start"))
        self.end = int(g("end"))
        self.strand = g("strand") or "NA"
        self.motif_plus = (g("motif_observed_plus") or "").upper()
        self.motif_canonical = g("motif_canonical") or canonical_motif(self.motif_plus)
        self.motif_len = int(g("motif_length") or len(self.motif_plus))
        self.tract_len = int(g("tract_length_bp") or (self.end - self.start + 1))
        self.context = g("context") or ""
        self.left40 = (g("left_flank_40") or "").upper()
        self.right40 = (g("right_flank_40") or "").upper()
        ce = (g("caller_eligible") or "YES")
        self.caller_eligible = ce.split()[0].upper().startswith("YES") if ce else True
        self.eligible_reason = ce
        self.interrupted = (g("interruption_flag") or "NO").upper() in ("YES", "TRUE", "1")
        self.cds_start = int(g("cds_start")) if (g("cds_start") or "NA") not in ("NA", "") else None
        self.cds_end = int(g("cds_end")) if (g("cds_end") or "NA") not in ("NA", "") else None
        self.evidence_class = g("evidence_class") or ""
        self.justified = (g("functional_interpretation_justified") or "YES").upper()
        self.ref_tract = None
        self.gene_name = g("gene_name") or ""
        self.locus_tag = g("locus_tag") or ""

    @property
    def ref_copies(self):
        return self.tract_len / self.motif_len


def load_catalogue(path):
    members = [Member(r) for r in read_tsv(path)]
    ids = [m.member_id for m in members]
    if len(ids) != len(set(ids)):
        raise ValueError("duplicate member_id in catalogue")
    return members


# ---------------------------------------------------------------- parameters
# (param_id, default, description). Types follow the default. PROPOSED defaults, to be calibrated/frozen by main.
DEFAULT_PARAMS = [
    ("anchor_len", 25, "Anchor length (bp) taken from the flank next to the tract on each side (max 40)."),
    ("seed_k", 12, "k-mer length of the seed index built on anchors."),
    ("seed_stride", 4, "Reads are probed for seeds every this many bases (must be <= anchor_len-seed_k+1 for guaranteed exact-anchor hit)."),
    ("max_anchor_mismatch", 2, "Max mismatches allowed per anchor at extraction."),
    ("max_tract_scan_len", 150, "Longest read tract (bp) between anchors kept as evidence at extraction."),
    ("min_read_len", 40, "Reads shorter than this are ignored."),
    ("pe_conflict_policy", "DROP", "Overlapping mates spanning with different tracts: DROP (fragment discarded, counted) or BEST (higher anchor quality mate)."),
    ("max_tract_len_call", 100, "Tract lengths above this in the catalogue are not callable (RD07)."),
    ("min_anchor_mean_q", 20.0, "Minimum mean base quality over both anchors of a spanning read."),
    ("max_total_anchor_mismatch", 3, "Max summed left+right anchor mismatches for a read to count at call stage."),
    ("max_tract_mismatch", 1, "Unit-aligned tracts with more mismatches than this vs pure repeat are discarded as discordant (0 = strict)."),
    ("min_spanning_reads", 8, "Minimum spanning fragments (after filters) for a callable locus."),
    ("min_top2_fraction", 0.8, "Top two alleles must hold at least this fraction of unit-aligned reads, else UNCALLABLE (diffuse)."),
    ("min_dominant_fraction", 0.4, "Dominant allele must hold at least this fraction of unit-aligned reads, else UNCALLABLE."),
    ("max_offunit_fraction", 0.3, "Reads whose tract length is not a whole number of motif units (vs reference): if more than this fraction, UNCALLABLE."),
    ("strand_bias_min_reads", 10, "Strand-bias check needs at least this many dominant-allele fragments."),
    ("strand_bias_min_frac", 0.1, "Dominant allele fragments must have both strands at least this fraction, else strand bias."),
    ("strand_bias_action", "DOWNGRADE", "DOWNGRADE (confidence one class lower), UNCALLABLE, or IGNORE."),
    ("stutter_a_MONO", 0.005, "PLACEHOLDER stutter intercept for -1 unit stutter, mononucleotide (estimate from data)."),
    ("stutter_b_MONO", 0.0025, "PLACEHOLDER stutter slope per dominant copy, mononucleotide."),
    ("stutter_a_DI", 0.01, "PLACEHOLDER stutter intercept for -1 unit, dinucleotide."),
    ("stutter_b_DI", 0.004, "PLACEHOLDER stutter slope per copy, dinucleotide."),
    ("stutter_a_POLY", 0.003, "PLACEHOLDER stutter intercept for -1 unit, motif length 3-6."),
    ("stutter_b_POLY", 0.001, "PLACEHOLDER stutter slope per copy, motif length 3-6."),
    ("stutter_up_ratio", 0.5, "Rate of +1 unit stutter relative to -1 unit."),
    ("stutter_k2_ratio", 0.25, "Each additional unit of distance multiplies the stutter rate by this factor."),
    ("background_rate", 0.003, "Per-read rate of any other length allele (sequencing/mapping error floor), applies at >=3 units and to non-unit shifts."),
    ("noise_safety_factor", 2.0, "Multiplier on model noise rate used in the mixture test (conservative)."),
    ("mixed_alpha", 1e-6, "Binomial upper-tail p-value below which a minor allele exceeds noise (MIXED)."),
    ("mixed_candidate_alpha", 1e-3, "Looser p-value used for MIXTURE_CANDIDATE."),
    ("mixed_min_minor_fraction", 0.05, "Minimum minor allele fraction for MIXED."),
    ("mixed_min_minor_reads", 4, "Minimum minor allele fragments for MIXED."),
    ("mixed_min_total_reads", 20, "Minimum unit-aligned reads for MIXED."),
    ("mixed_min_minor_strand_frac", 0.1, "Minor allele reads must be on both strands, each at least this fraction, for MIXED."),
    ("candidate_min_fraction", 0.10, "Minor allele at or above this fraction but not meeting MIXED is MIXTURE_CANDIDATE."),
    ("candidate_min_reads", 3, "Minimum minor fragments for MIXTURE_CANDIDATE."),
    ("conf_high_min_reads", 20, "Spanning fragments for HIGH confidence."),
    ("conf_medium_min_reads", 10, "Spanning fragments for MEDIUM confidence."),
    ("conf_min_dominant_ci_low", 0.7, "Non-mixed call with Wilson lower bound of dominant fraction below this is downgraded one class."),
    ("cds_on_min_orf_fraction", 0.9, "ON requires ORF >= this fraction of annotated protein length (FS02)."),
    ("orf_downstream_extend", 600, "Reference bases appended past CDS end when translating a frameshifted CDS."),
    ("variant_min_fraction", 0.5, "Anchor mismatch at same site/base in >= this fraction of dominant reads counts as a sample variant (FS04)."),
    ("variant_min_reads", 3, "Minimum reads carrying an anchor variant to count as a sample variant."),
]


class Params(dict):
    def __init__(self, overrides=None):
        super().__init__({k: v for k, v, _ in DEFAULT_PARAMS})
        for k, v in (overrides or {}).items():
            self.set(k, v)

    def set(self, k, v):
        if k not in self:
            raise KeyError("unknown caller rule parameter: %s" % k)
        d = type(dict.__getitem__(self, k))
        dict.__setitem__(self, k, d(v) if d is not str else str(v))

    def __getattr__(self, k):
        try:
            return self[k]
        except KeyError:
            raise AttributeError(k)


def load_params(path=None):
    p = Params()
    if path:
        for r in read_tsv(path):
            p.set(r["param_id"], r["value"])
    return p


def write_rules_tsv(path):
    with open(path, "w", newline="") as fh:
        fh.write("param_id\tvalue\tdescription\n")
        for k, v, d in DEFAULT_PARAMS:
            fh.write("%s\t%s\t%s\n" % (k, v, d))
