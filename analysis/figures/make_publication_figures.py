#!/usr/bin/env python3
"""Render frozen publication Figures 1-4 from tracked Project A outputs.

Rendering only. This script must not change thresholds, disease labels, cohort
membership, statistical endpoints, or any frozen analytical result.
"""
from __future__ import annotations

import argparse
import csv
import json
from collections import Counter
from pathlib import Path

import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[2]
STYLE_PATH = ROOT / "analysis/figures/figure_style.json"
NUMBERS_PATH = ROOT / "results/manuscript/manuscript_numbers.tsv"


def read_tsv(path: Path):
    with path.open(newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def load_numbers():
    rows = read_tsv(NUMBERS_PATH)
    return {r["metric"]: r for r in rows}


def n_int(numbers, key):
    return int(numbers[key]["value"])


def n_float(numbers, key):
    return float(numbers[key]["value"])


def pct(num, den, digits=1):
    return 100.0 * num / den if den else float("nan")


def set_style(style):
    typ = style["typography"]
    plt.rcParams.update({
        "font.family": typ["font_family"],
        "font.size": typ["tick_size"],
        "axes.titlesize": typ["title_size"],
        "axes.labelsize": typ["label_size"],
        "xtick.labelsize": typ["tick_size"],
        "ytick.labelsize": typ["tick_size"],
        "text.color": style["canvas"]["foreground"],
        "axes.labelcolor": style["canvas"]["foreground"],
        "axes.edgecolor": style["canvas"]["panel_border"],
        "xtick.color": style["canvas"]["muted_text"],
        "ytick.color": style["canvas"]["muted_text"],
        "figure.facecolor": style["canvas"]["background"],
        "axes.facecolor": style["canvas"]["background"],
        "savefig.facecolor": style["canvas"]["background"],
        "axes.grid": False,
        "svg.fonttype": "none",
        "pdf.fonttype": 42,
    })


def clean_ax(ax, grid_axis=None):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    if grid_axis:
        ax.grid(axis=grid_axis, color="#dfe3e8", linewidth=0.6, alpha=0.8)
        ax.set_axisbelow(True)


def panel_label(ax, label, style):
    ax.text(
        -0.12, 1.08, label, transform=ax.transAxes,
        fontsize=style["typography"]["panel_label_size"],
        fontweight="bold", va="top", ha="left"
    )


def annotate_bar_values(ax, bars, values, fmt="{:,.0f}", pad=3):
    for bar, value in zip(bars, values):
        ax.annotate(
            fmt.format(value),
            (bar.get_x() + bar.get_width() / 2, bar.get_height()),
            xytext=(0, pad), textcoords="offset points",
            ha="center", va="bottom", fontsize=8
        )


def save_all(fig, outdir, stem, style):
    outdir.mkdir(parents=True, exist_ok=True)
    fig.savefig(outdir / f"{stem}.svg", bbox_inches="tight")
    fig.savefig(outdir / f"{stem}.pdf", bbox_inches="tight")
    fig.savefig(outdir / f"{stem}.png", dpi=style["export"]["dpi"], bbox_inches="tight")
    plt.close(fig)


def figure1(numbers, style, outdir):
    sem = style["semantic"]
    fig, axes = plt.subplots(2, 2, figsize=(7.2, 5.4))
    fig.subplots_adjust(wspace=0.36, hspace=0.52)

    ax = axes[0, 0]
    vals = [
        n_int(numbers, "catalogue_loci"),
        n_int(numbers, "catalogue_members"),
        n_int(numbers, "caller_eligible_loci"),
    ]
    labels = ["Repeat loci", "Catalogue members", "Caller-eligible loci"]
    colors = [sem["primary"], sem["neutral"], sem["eligible"]]
    bars = ax.bar(labels, vals, color=colors, width=0.62)
    annotate_bar_values(ax, bars, vals)
    ax.set_title("Catalogue scale")
    ax.set_ylabel("Count")
    ax.tick_params(axis="x", rotation=20)
    clean_ax(ax, "y")
    panel_label(ax, "A", style)

    ax = axes[0, 1]
    audit_labels = ["Boundary\nstable", "Boundary\nambiguous", "Orthology\nambiguous", "Insufficient\nto resolve"]
    audit_vals = [
        n_int(numbers, "boundary_stable_loci"),
        n_int(numbers, "boundary_ambiguous_loci"),
        n_int(numbers, "orthology_ambiguous_loci"),
        n_int(numbers, "insufficient_to_resolve_loci"),
    ]
    audit_colors = [sem["eligible"], sem["caution"], sem["failure"], sem["neutral"]]
    bars = ax.bar(audit_labels, audit_vals, color=audit_colors, width=0.68)
    annotate_bar_values(ax, bars, audit_vals)
    ax.set_title("Boundary / orthology audit")
    ax.set_ylabel("Loci")
    clean_ax(ax, "y")
    panel_label(ax, "B", style)

    ax = axes[1, 0]
    universe_labels = ["PRIMARY", "PROVISIONAL", "INELIGIBLE"]
    universe_vals = [
        n_int(numbers, "primary_technical_loci"),
        n_int(numbers, "provisional_technical_loci"),
        n_int(numbers, "primary_ineligible_loci"),
    ]
    universe_colors = [sem["primary"], sem["provisional"], sem["ineligible"]]
    bars = ax.bar(universe_labels, universe_vals, color=universe_colors, width=0.65)
    annotate_bar_values(ax, bars, universe_vals)
    ax.set_title("Frozen technical universe")
    ax.set_ylabel("Loci")
    clean_ax(ax, "y")
    panel_label(ax, "C", style)

    ax = axes[1, 1]
    evidence_labels = ["Repeat-only", "Strong PV\ncandidate", "Mixed", "Known PV\nin PRIMARY"]
    evidence_vals = [
        n_int(numbers, "primary_repeat_only_loci"),
        n_int(numbers, "primary_strong_pv_candidate_loci"),
        n_int(numbers, "primary_mixed_evidence_loci"),
        n_int(numbers, "primary_known_pv_loci"),
    ]
    evidence_colors = [
        sem["repeat_only"], sem["strong_pv_candidate"],
        sem["mixed_evidence"], sem["known_pv"]
    ]
    bars = ax.bar(evidence_labels, evidence_vals, color=evidence_colors, width=0.68)
    annotate_bar_values(ax, bars, evidence_vals)
    ax.set_title("Evidence within PRIMARY_TECHNICAL")
    ax.set_ylabel("Loci")
    clean_ax(ax, "y")
    panel_label(ax, "D", style)
    ax.text(
        0.01, -0.32,
        "Technical eligibility and phase-variation evidence are distinct annotations.",
        transform=ax.transAxes, fontsize=7.5, color=style["canvas"]["muted_text"]
    )

    fig.suptitle("Figure 1. Repeat-locus catalogue and frozen technical analysis universe", y=1.01, fontsize=12)
    save_all(fig, outdir, "fig1_catalogue_universe", style)


def figure2(numbers, style, outdir):
    sem = style["semantic"]
    fig, axes = plt.subplots(2, 2, figsize=(7.2, 5.5))
    fig.subplots_adjust(wspace=0.37, hspace=0.58)

    ax = axes[0, 0]
    empirical_total = n_int(numbers, "synthetic_empirical_total")
    callable_n = n_int(numbers, "synthetic_callable")
    exact_n = n_int(numbers, "synthetic_exact_dominant")
    confident_n = n_int(numbers, "synthetic_confident")
    false_n = n_int(numbers, "synthetic_false_confident_wrong")
    vals = [
        pct(callable_n, empirical_total),
        pct(exact_n, callable_n),
        pct(confident_n, callable_n),
        pct(false_n, confident_n),
    ]
    labels = ["Callable\n1780/2720", "Exact | callable\n1757/1780", "Confident | callable\n1595/1780", "Wrong | confident\n15/1595"]
    colors = [sem["eligible"], sem["primary"], sem["provisional"], sem["failure"]]
    bars = ax.bar(labels, vals, color=colors, width=0.68)
    annotate_bar_values(ax, bars, vals, fmt="{:.1f}%")
    ax.set_ylim(0, 105)
    ax.set_ylabel("Percent")
    ax.set_title("Synthetic empirical benchmark")
    clean_ax(ax, "y")
    panel_label(ax, "A", style)

    ax = axes[0, 1]
    truth_labels = ["High-confidence", "Ambiguous", "Uncallable"]
    truth_vals = [
        n_int(numbers, "hifi_truth_high_confidence"),
        n_int(numbers, "hifi_truth_ambiguous"),
        n_int(numbers, "hifi_truth_uncallable"),
    ]
    truth_colors = [sem["eligible"], sem["caution"], sem["ineligible"]]
    bars = ax.bar(truth_labels, truth_vals, color=truth_colors, width=0.68)
    annotate_bar_values(ax, bars, truth_vals)
    ax.set_title("HiFi truth availability")
    ax.set_ylabel("Isolate-locus pairs")
    clean_ax(ax, "y")
    panel_label(ax, "B", style)
    ax.text(
        0.99, 0.98, f"Total = {n_int(numbers, 'hifi_truth_total'):,}",
        transform=ax.transAxes, ha="right", va="top", fontsize=8,
        color=style["canvas"]["muted_text"]
    )

    ax = axes[1, 0]
    overall_n = n_int(numbers, "hifi_high_confidence_comparisons")
    overall_exact = n_int(numbers, "hifi_exact")
    nr_n = n_int(numbers, "hifi_nonreference_comparisons")
    nr_exact = n_int(numbers, "hifi_nonreference_exact")
    vals = [pct(overall_exact, overall_n), pct(nr_exact, nr_n)]
    labels = [f"All high-confidence\n{overall_exact}/{overall_n}", f"Non-reference\n{nr_exact}/{nr_n}"]
    bars = ax.bar(labels, vals, color=[sem["primary"], sem["provisional"]], width=0.58)
    annotate_bar_values(ax, bars, vals, fmt="{:.1f}%")
    ax.set_ylim(0, 102)
    ax.set_ylabel("Exact concordance (%)")
    ax.set_title("Short-read versus HiFi")
    clean_ax(ax, "y")
    panel_label(ax, "C", style)

    ax = axes[1, 1]
    rows = read_tsv(ROOT / "results/validation/hifi_nonreference_summary.tsv")
    rows = [r for r in rows if not r["stratum"].startswith("ALL_HIGH_CONFIDENCE")]
    strata = [r["stratum"] for r in rows]
    totals = [int(r["nonreference_comparisons"]) for r in rows]
    exact = [int(r["nonreference_exact"]) for r in rows]
    y = list(range(len(strata)))
    ax.barh(y, totals, color=sem["ineligible"], label="Observed non-reference truth")
    ax.barh(y, exact, color=sem["primary"], label="Exact")
    ax.set_yticks(y, strata)
    ax.invert_yaxis()
    ax.set_xlabel("Comparisons")
    ax.set_title("Non-reference validation coverage")
    clean_ax(ax, "x")
    panel_label(ax, "D", style)
    ax.legend(frameon=False, fontsize=7, loc="lower right")
    ax.text(
        0.0, -0.25,
        "Zero bars denote no observed non-reference truth; they do not validate that stratum.",
        transform=ax.transAxes, fontsize=7.3, color=style["canvas"]["muted_text"]
    )

    fig.suptitle("Figure 2. Technical validation of repeat-length calling", y=1.01, fontsize=12)
    save_all(fig, outdir, "fig2_technical_validation", style)


def figure3(numbers, style, outdir):
    sem = style["semantic"]
    fig, axes = plt.subplots(2, 2, figsize=(7.2, 5.6))
    fig.subplots_adjust(wspace=0.38, hspace=0.58)

    ax = axes[0, 0]
    labels = ["Stable", "Boundary\nambiguous", "Orthology\nambiguous", "Insufficient\nto resolve"]
    vals = [
        n_int(numbers, "boundary_stable_loci"),
        n_int(numbers, "boundary_ambiguous_loci"),
        n_int(numbers, "orthology_ambiguous_loci"),
        n_int(numbers, "insufficient_to_resolve_loci"),
    ]
    bars = ax.bar(labels, vals, color=[sem["eligible"], sem["caution"], sem["failure"], sem["neutral"]], width=0.68)
    annotate_bar_values(ax, bars, vals)
    ax.set_title("Catalogue-resolution failure modes")
    ax.set_ylabel("Loci")
    clean_ax(ax, "y")
    panel_label(ax, "A", style)

    ax = axes[0, 1]
    platform = read_tsv(ROOT / "results/callability/platform_transfer.tsv")
    counts = Counter(r["callability_class"] for r in platform)
    order = sorted(counts, key=lambda k: (-counts[k], k))
    vals = [counts[k] for k in order]
    labels = [k.replace("_", "\n") for k in order]
    colors = [sem["caution"] if k == "PLATFORM_CONFOUNDED_CALLABILITY" else sem["provisional"] for k in order]
    bars = ax.bar(labels, vals, color=colors, width=0.7)
    annotate_bar_values(ax, bars, vals)
    ax.set_title("Study/platform-associated callability")
    ax.set_ylabel("Loci")
    ax.tick_params(axis="x", labelsize=6.5)
    clean_ax(ax, "y")
    panel_label(ax, "B", style)
    ax.text(
        0.0, -0.31,
        f"{n_int(numbers, 'platform_study_confounded_callability_loci')} loci carry the frozen PLATFORM_CONFOUNDED_CALLABILITY label; descriptive, not causal.",
        transform=ax.transAxes, fontsize=7.2, color=style["canvas"]["muted_text"]
    )

    ax = axes[1, 0]
    universe = read_tsv(ROOT / "metadata/loci/analysis_universe.tsv")
    known = Counter(r["g7_locus_class"] for r in universe if r["known_pv"] == "YES")
    classes = ["PRIMARY_TECHNICAL", "PROVISIONAL_TECHNICAL", "PRIMARY_INELIGIBLE"]
    vals = [known[c] for c in classes]
    labels = ["PRIMARY", "PROVISIONAL", "INELIGIBLE"]
    bars = ax.bar(labels, vals, color=[sem["primary"], sem["provisional"], sem["known_pv"]], width=0.65)
    annotate_bar_values(ax, bars, vals)
    ax.set_title("Literature-anchored KNOWN_PV representation")
    ax.set_ylabel("Known-PV loci")
    clean_ax(ax, "y")
    panel_label(ax, "C", style)
    ax.text(
        0.0, -0.27,
        "No KNOWN_PV locus is in PRIMARY_TECHNICAL.",
        transform=ax.transAxes, fontsize=7.5, color=style["canvas"]["muted_text"]
    )

    ax = axes[1, 1]
    ax.axis("off")
    panel_label(ax, "D", style)
    ax.set_title("Reproducibility boundary")
    lines = [
        "Frozen generating source exported for G5-G8.",
        "Exact historical environment: partially recorded.",
        "Phase-1 call matrix: omitted large regenerable artifact;",
        "  archive blob SHA/size and regeneration route retained.",
        "G6 curated closure-table provenance: disclosed limitation.",
        "No new thresholds, cohorts, endpoints, or disease reclassification."
    ]
    y = 0.86
    for i, line in enumerate(lines):
        weight = "bold" if i in (0, 5) else "normal"
        ax.text(0.02, y, line, transform=ax.transAxes, va="top", fontsize=8.2, fontweight=weight)
        y -= 0.13 if i != 2 else 0.10

    fig.suptitle("Figure 3. Technical failure modes and interpretation boundary", y=1.01, fontsize=12)
    save_all(fig, outdir, "fig3_failure_modes", style)


def figure4(numbers, style, outdir):
    sem = style["semantic"]
    fig, axes = plt.subplots(2, 2, figsize=(7.2, 5.7))
    fig.subplots_adjust(wspace=0.38, hspace=0.60)

    ax = axes[0, 0]
    vals = [
        n_int(numbers, "confounding_single_country_studies"),
        n_int(numbers, "confounding_multicountry_studies"),
    ]
    labels = ["Single-country", "Multicountry"]
    bars = ax.bar(labels, vals, color=[sem["neutral"], sem["caution"]], width=0.58)
    annotate_bar_values(ax, bars, vals)
    ax.set_title("Public-cohort geography structure")
    ax.set_ylabel("Studies")
    clean_ax(ax, "y")
    panel_label(ax, "A", style)
    ax.text(
        0.02, 0.96,
        f"{n_int(numbers, 'confounding_correa_studies')} studies; {n_int(numbers, 'confounding_correa_biosamples')} Correa-labelled BioSamples",
        transform=ax.transAxes, va="top", fontsize=7.8, color=style["canvas"]["muted_text"]
    )

    ax = axes[0, 1]
    cohort_labels = ["PRJNA360417\nNAG", "PRJNA360417\nIM", "PRJNA678459\nAG", "PRJNA678459\nGC", "PRJNA1103397\nGC"]
    cohort_vals = [
        n_int(numbers, "prjna360417_nag"),
        n_int(numbers, "prjna360417_im"),
        n_int(numbers, "prjna678459_ag"),
        n_int(numbers, "prjna678459_gc"),
        n_int(numbers, "prjna1103397_gc"),
    ]
    bars = ax.bar(cohort_labels, cohort_vals, color=[sem["provisional"], sem["primary"], sem["provisional"], sem["primary"], sem["neutral"]], width=0.68)
    annotate_bar_values(ax, bars, cohort_vals)
    ax.set_title("Pilot cohort sample sizes")
    ax.set_ylabel("Patients / selected isolates")
    ax.tick_params(axis="x", labelsize=6.7)
    clean_ax(ax, "y")
    panel_label(ax, "B", style)

    ax = axes[1, 0]
    labels = ["PRJNA360417", "PRJNA678459", "PRJNA1103397"]
    primary_vals = [
        n_int(numbers, "prjna360417_primary_screen_pass"),
        n_int(numbers, "prjna678459_primary_screen_pass"),
        n_int(numbers, "prjna1103397_primary_screen_pass"),
    ]
    provisional_vals = [0, n_int(numbers, "prjna678459_provisional_screen_pass"), 0]
    x = list(range(len(labels)))
    b1 = ax.bar(x, primary_vals, color=sem["primary"], width=0.58, label="PRIMARY passed screen")
    b2 = ax.bar(x, provisional_vals, bottom=primary_vals, color=sem["provisional"], width=0.58, label="PROVISIONAL passed screen")
    annotate_bar_values(ax, b1, primary_vals)
    for bar, base, value in zip(b2, primary_vals, provisional_vals):
        if value:
            ax.annotate(
                str(value), (bar.get_x() + bar.get_width()/2, base + value),
                xytext=(0, 3), textcoords="offset points", ha="center", va="bottom", fontsize=8
            )
    ax.set_xticks(x, labels)
    ax.set_title("Locus-screen attrition")
    ax.set_ylabel("Loci passing frozen screen")
    clean_ax(ax, "y")
    panel_label(ax, "C", style)
    ax.legend(frameon=False, fontsize=7)
    ax.text(
        0.0, -0.27, "PRIMARY denominator = 261 loci; PROVISIONAL denominator = 39 loci.",
        transform=ax.transAxes, fontsize=7.3, color=style["canvas"]["muted_text"]
    )

    ax = axes[1, 1]
    ax.axis("off")
    panel_label(ax, "D", style)
    ax.set_title("Strict endpoint status")
    lines = [
        ("PRJNA360417 NAG vs IM", "NOT_COMPUTABLE", n_float(numbers, "prjna360417_diagnostic_exact_p")),
        ("PRJNA678459 AG vs GC", "NOT_COMPUTABLE", n_float(numbers, "prjna678459_diagnostic_exact_p")),
    ]
    y = 0.83
    for cohort, status, pval in lines:
        ax.text(0.02, y, cohort, transform=ax.transAxes, fontsize=8.5, fontweight="bold", va="top")
        ax.text(0.02, y-0.11, f"Strict: {status}", transform=ax.transAxes, fontsize=8.2, color=sem["failure"], va="top")
        ax.text(0.02, y-0.21, f"Below-floor diagnostic exact p = {pval:.6f}", transform=ax.transAxes, fontsize=7.8, color=style["canvas"]["muted_text"], va="top")
        y -= 0.38
    ax.text(
        0.02, 0.06,
        "Pre-specified floor: >=10 jointly callable loci/sample-pair.\n"
        "PRJNA678459 structure: screening caution only "
        f"(R²={n_float(numbers, 'structure_prjna678459_r2'):.4f}, p={n_float(numbers, 'structure_prjna678459_exact_p'):.6f}).\n"
        f"PRJNA360417 structure: {numbers['structure_prjna360417']['value']}.",
        transform=ax.transAxes, fontsize=7.2, color=style["canvas"]["muted_text"], va="bottom"
    )

    fig.suptitle("Figure 4. Public-cohort stress test under frozen rules", y=1.01, fontsize=12)
    save_all(fig, outdir, "fig4_public_cohort_stress_test", style)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", default="figures/generated")
    args = parser.parse_args()

    style = json.loads(STYLE_PATH.read_text())
    set_style(style)
    numbers = load_numbers()
    outdir = ROOT / args.output_dir

    figure1(numbers, style, outdir)
    figure2(numbers, style, outdir)
    figure3(numbers, style, outdir)
    figure4(numbers, style, outdir)
    print(f"Rendered publication Figures 1-4 to {outdir}")


if __name__ == "__main__":
    main()
