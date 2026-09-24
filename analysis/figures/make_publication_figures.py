#!/usr/bin/env python3
"""Render frozen publication Figures 1-4 from tracked Project A outputs.

Rendering only. This script must not change thresholds, disease labels, cohort
membership, statistical endpoints, or any frozen analytical result.
"""
from __future__ import annotations

import argparse
import csv
import json
import textwrap
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
    return {r["metric"]: r for r in read_tsv(NUMBERS_PATH)}


def n_int(numbers, key):
    return int(numbers[key]["value"])


def n_float(numbers, key):
    return float(numbers[key]["value"])


def pct(num, den):
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
        -0.10, 1.05, label, transform=ax.transAxes,
        fontsize=style["typography"]["panel_label_size"],
        fontweight="bold", va="top", ha="left"
    )


def annotate_vertical(ax, bars, values, fmt="{:,.0f}", pad=3):
    for bar, value in zip(bars, values):
        ax.annotate(
            fmt.format(value),
            (bar.get_x() + bar.get_width() / 2, bar.get_height()),
            xytext=(0, pad), textcoords="offset points",
            ha="center", va="bottom", fontsize=8
        )


def horizontal_bars(ax, labels, values, colors, title, xlabel, value_fmt="{:,.0f}", xmax=None):
    y = list(range(len(labels)))
    bars = ax.barh(y, values, color=colors, height=0.62)
    ax.set_yticks(y, labels)
    ax.invert_yaxis()
    ax.set_title(title, pad=6)
    ax.set_xlabel(xlabel)
    clean_ax(ax, "x")
    max_v = max(values) if values else 1
    if xmax is None:
        xmax = max_v * 1.22 if max_v > 0 else 1
    ax.set_xlim(0, xmax)
    for bar, value in zip(bars, values):
        ax.annotate(
            value_fmt.format(value),
            (bar.get_width(), bar.get_y() + bar.get_height() / 2),
            xytext=(4, 0), textcoords="offset points",
            ha="left", va="center", fontsize=8
        )
    return bars


def save_all(fig, outdir, stem, style):
    outdir.mkdir(parents=True, exist_ok=True)
    fig.savefig(outdir / f"{stem}.svg", bbox_inches="tight")
    fig.savefig(outdir / f"{stem}.pdf", bbox_inches="tight")
    fig.savefig(outdir / f"{stem}.png", dpi=style["export"]["dpi"], bbox_inches="tight")
    plt.close(fig)


def figure1(numbers, style, outdir):
    sem = style["semantic"]
    fig, axes = plt.subplots(2, 2, figsize=(7.4, 5.7), layout="constrained")

    ax = axes[0, 0]
    horizontal_bars(
        ax,
        ["Repeat loci", "Catalogue members", "Caller-eligible loci"],
        [
            n_int(numbers, "catalogue_loci"),
            n_int(numbers, "catalogue_members"),
            n_int(numbers, "caller_eligible_loci"),
        ],
        [sem["primary"], sem["neutral"], sem["eligible"]],
        "Catalogue scale",
        "Count",
    )
    panel_label(ax, "A", style)

    ax = axes[0, 1]
    horizontal_bars(
        ax,
        ["Boundary stable", "Boundary ambiguous", "Orthology ambiguous", "Insufficient to resolve"],
        [
            n_int(numbers, "boundary_stable_loci"),
            n_int(numbers, "boundary_ambiguous_loci"),
            n_int(numbers, "orthology_ambiguous_loci"),
            n_int(numbers, "insufficient_to_resolve_loci"),
        ],
        [sem["eligible"], sem["caution"], sem["failure"], sem["neutral"]],
        "Boundary / orthology audit",
        "Loci",
    )
    panel_label(ax, "B", style)

    ax = axes[1, 0]
    horizontal_bars(
        ax,
        ["PRIMARY", "PROVISIONAL", "INELIGIBLE"],
        [
            n_int(numbers, "primary_technical_loci"),
            n_int(numbers, "provisional_technical_loci"),
            n_int(numbers, "primary_ineligible_loci"),
        ],
        [sem["primary"], sem["provisional"], sem["ineligible"]],
        "Frozen technical universe",
        "Loci",
    )
    panel_label(ax, "C", style)

    ax = axes[1, 1]
    horizontal_bars(
        ax,
        ["Repeat-only", "Strong PV candidate", "Mixed evidence", "Known PV in PRIMARY"],
        [
            n_int(numbers, "primary_repeat_only_loci"),
            n_int(numbers, "primary_strong_pv_candidate_loci"),
            n_int(numbers, "primary_mixed_evidence_loci"),
            n_int(numbers, "primary_known_pv_loci"),
        ],
        [sem["repeat_only"], sem["strong_pv_candidate"], sem["mixed_evidence"], sem["known_pv"]],
        "Evidence within PRIMARY_TECHNICAL",
        "Loci",
        xmax=285,
    )
    ax.text(
        0.98, 0.03,
        "Technical eligibility ≠ phase-variation evidence.",
        transform=ax.transAxes, ha="right", va="bottom",
        fontsize=7.2, color=style["canvas"]["muted_text"]
    )
    panel_label(ax, "D", style)

    fig.suptitle("Figure 1. Repeat-locus catalogue and frozen technical analysis universe", fontsize=12)
    save_all(fig, outdir, "fig1_catalogue_universe", style)


def figure2(numbers, style, outdir):
    sem = style["semantic"]
    fig, axes = plt.subplots(2, 2, figsize=(7.4, 5.9), layout="constrained")

    ax = axes[0, 0]
    empirical_total = n_int(numbers, "synthetic_empirical_total")
    callable_n = n_int(numbers, "synthetic_callable")
    exact_n = n_int(numbers, "synthetic_exact_dominant")
    confident_n = n_int(numbers, "synthetic_confident")
    false_n = n_int(numbers, "synthetic_false_confident_wrong")
    rate_labels = [
        "Callable  [1780/2720]",
        "Exact among callable  [1757/1780]",
        "Confident among callable  [1595/1780]",
        "Wrong among confident  [15/1595]",
    ]
    rate_vals = [
        pct(callable_n, empirical_total),
        pct(exact_n, callable_n),
        pct(confident_n, callable_n),
        pct(false_n, confident_n),
    ]
    horizontal_bars(
        ax, rate_labels, rate_vals,
        [sem["eligible"], sem["primary"], sem["provisional"], sem["failure"]],
        "Synthetic empirical benchmark", "Percent", value_fmt="{:.1f}%", xmax=108
    )
    panel_label(ax, "A", style)

    ax = axes[0, 1]
    truth_labels = ["High-confidence", "Ambiguous", "Uncallable"]
    truth_vals = [
        n_int(numbers, "hifi_truth_high_confidence"),
        n_int(numbers, "hifi_truth_ambiguous"),
        n_int(numbers, "hifi_truth_uncallable"),
    ]
    horizontal_bars(
        ax, truth_labels, truth_vals,
        [sem["eligible"], sem["caution"], sem["ineligible"]],
        "HiFi truth availability", "Isolate-locus pairs"
    )
    ax.text(
        0.98, 0.05, f"Total = {n_int(numbers, 'hifi_truth_total'):,}",
        transform=ax.transAxes, ha="right", va="bottom",
        fontsize=7.5, color=style["canvas"]["muted_text"]
    )
    panel_label(ax, "B", style)

    ax = axes[1, 0]
    overall_n = n_int(numbers, "hifi_high_confidence_comparisons")
    overall_exact = n_int(numbers, "hifi_exact")
    nr_n = n_int(numbers, "hifi_nonreference_comparisons")
    nr_exact = n_int(numbers, "hifi_nonreference_exact")
    vals = [pct(overall_exact, overall_n), pct(nr_exact, nr_n)]
    labels = ["All high-confidence", "Non-reference"]
    bars = ax.bar(labels, vals, color=[sem["primary"], sem["provisional"]], width=0.58)
    annotate_vertical(ax, bars, vals, fmt="{:.1f}%")
    ax.set_ylim(0, 104)
    ax.set_ylabel("Exact concordance (%)")
    ax.set_title("Short-read versus HiFi", pad=6)
    clean_ax(ax, "y")
    ax.text(0, 2, f"{overall_exact}/{overall_n}", ha="center", va="bottom", fontsize=7, color=style["canvas"]["muted_text"])
    ax.text(1, 2, f"{nr_exact}/{nr_n}", ha="center", va="bottom", fontsize=7, color=style["canvas"]["muted_text"])
    panel_label(ax, "C", style)

    ax = axes[1, 1]
    rows = read_tsv(ROOT / "results/validation/hifi_nonreference_summary.tsv")
    rows = [r for r in rows if not r["stratum"].startswith("ALL_HIGH_CONFIDENCE")]
    strata = [r["stratum"] for r in rows]
    totals = [int(r["nonreference_comparisons"]) for r in rows]
    exact = [int(r["nonreference_exact"]) for r in rows]
    y = list(range(len(strata)))
    ax.barh(y, totals, color=sem["ineligible"], height=0.62, label="Observed non-reference truth")
    ax.barh(y, exact, color=sem["primary"], height=0.62, label="Exact")
    ax.set_yticks(y, strata)
    ax.invert_yaxis()
    ax.set_xlabel("Comparisons")
    ax.set_title("Non-reference validation coverage", pad=6)
    clean_ax(ax, "x")
    ax.legend(frameon=False, fontsize=6.8, loc="lower right")
    for yi, total in enumerate(totals):
        if total == 0:
            ax.text(0.8, yi, "0", va="center", ha="left", fontsize=7.2,
                    color=style["canvas"]["muted_text"])
    panel_label(ax, "D", style)

    fig.suptitle("Figure 2. Technical validation of repeat-length calling", fontsize=12)
    save_all(fig, outdir, "fig2_technical_validation", style)


def figure3(numbers, style, outdir):
    sem = style["semantic"]
    fig, axes = plt.subplots(2, 2, figsize=(8.6, 5.9), layout="constrained")

    ax = axes[0, 0]
    horizontal_bars(
        ax,
        ["Stable", "Boundary ambig.", "Orthology ambig.", "Insufficient"],
        [
            n_int(numbers, "boundary_stable_loci"),
            n_int(numbers, "boundary_ambiguous_loci"),
            n_int(numbers, "orthology_ambiguous_loci"),
            n_int(numbers, "insufficient_to_resolve_loci"),
        ],
        [sem["eligible"], sem["caution"], sem["failure"], sem["neutral"]],
        "Catalogue-resolution states",
        "Loci",
    )
    panel_label(ax, "A", style)

    ax = axes[0, 1]
    platform = read_tsv(ROOT / "results/callability/platform_transfer.tsv")
    counts = Counter(r["callability_class"] for r in platform)
    preferred = [
        "CALLABLE_BOTH", "PARTIAL", "UNCALLABLE_BOTH",
        "PLATFORM_CONFOUNDED_CALLABILITY", "NA"
    ]
    order = [k for k in preferred if k in counts] + [k for k in sorted(counts) if k not in preferred]
    pretty = {
        "CALLABLE_BOTH": "Callable both",
        "PARTIAL": "Partial",
        "UNCALLABLE_BOTH": "Uncallable both",
        "PLATFORM_CONFOUNDED_CALLABILITY": "Confounded label",
        "NA": "NA",
    }
    labels = [pretty.get(k, k.replace("_", " ").title()) for k in order]
    vals = [counts[k] for k in order]
    colors = []
    for k in order:
        if k == "PLATFORM_CONFOUNDED_CALLABILITY":
            colors.append(sem["caution"])
        elif k == "UNCALLABLE_BOTH":
            colors.append(sem["neutral"])
        else:
            colors.append(sem["provisional"])
    horizontal_bars(ax, labels, vals, colors, "Study/platform-associated callability", "Loci")
    panel_label(ax, "B", style)

    ax = axes[1, 0]
    universe = read_tsv(ROOT / "metadata/loci/analysis_universe.tsv")
    known = Counter(r["g7_locus_class"] for r in universe if r["known_pv"] == "YES")
    classes = ["PRIMARY_TECHNICAL", "PROVISIONAL_TECHNICAL", "PRIMARY_INELIGIBLE"]
    vals = [known[c] for c in classes]
    horizontal_bars(
        ax,
        ["PRIMARY", "PROVISIONAL", "INELIGIBLE"],
        vals,
        [sem["primary"], sem["provisional"], sem["known_pv"]],
        "Literature-anchored KNOWN_PV",
        "Known-PV loci",
        xmax=max(max(vals + [1]) * 1.25, 10),
    )
    panel_label(ax, "C", style)

    ax = axes[1, 1]
    ax.axis("off")
    ax.set_title("Reproducibility boundary", pad=6)
    panel_label(ax, "D", style)
    blocks = [
        ("Frozen generating source exported for G5-G8.", True),
        ("Exact historical environment was only partially recorded.", False),
        ("Phase-1 call matrix is an omitted large regenerable artifact; archive blob SHA/size and the regeneration route are retained.", False),
        ("G6 curated closure-table provenance remains a disclosed limitation.", False),
        ("No new thresholds, cohorts, endpoints, or disease reclassification.", True),
    ]
    y = 0.88
    for text_value, bold in blocks:
        wrapped = textwrap.fill(text_value, width=48)
        ax.text(
            0.02, y, wrapped, transform=ax.transAxes, va="top",
            fontsize=8.0, fontweight="bold" if bold else "normal"
        )
        y -= 0.16 + 0.06 * (wrapped.count("\n"))

    fig.suptitle("Figure 3. Technical failure modes and interpretation boundary", fontsize=12)
    save_all(fig, outdir, "fig3_failure_modes", style)


def figure4(numbers, style, outdir):
    sem = style["semantic"]
    fig, axes = plt.subplots(2, 2, figsize=(7.5, 6.1), layout="constrained")

    ax = axes[0, 0]
    vals = [
        n_int(numbers, "confounding_single_country_studies"),
        n_int(numbers, "confounding_multicountry_studies"),
    ]
    labels = ["Single-country", "Multicountry"]
    bars = ax.bar(labels, vals, color=[sem["neutral"], sem["caution"]], width=0.58)
    annotate_vertical(ax, bars, vals)
    ax.set_title("Public-cohort geography structure", pad=6)
    ax.set_ylabel("Studies")
    ax.set_ylim(0, 46)
    clean_ax(ax, "y")
    panel_label(ax, "A", style)

    ax = axes[0, 1]
    cohort_labels = ["360417 NAG", "360417 IM", "678459 AG", "678459 GC", "1103397 GC"]
    cohort_vals = [
        n_int(numbers, "prjna360417_nag"),
        n_int(numbers, "prjna360417_im"),
        n_int(numbers, "prjna678459_ag"),
        n_int(numbers, "prjna678459_gc"),
        n_int(numbers, "prjna1103397_gc"),
    ]
    horizontal_bars(
        ax, cohort_labels, cohort_vals,
        [sem["provisional"], sem["primary"], sem["provisional"], sem["primary"], sem["neutral"]],
        "Pilot cohort sample sizes", "Patients / selected isolates", xmax=7
    )
    panel_label(ax, "B", style)

    ax = axes[1, 0]
    labels = ["PRJNA360417", "PRJNA678459", "PRJNA1103397"]
    primary_vals = [
        n_int(numbers, "prjna360417_primary_screen_pass"),
        n_int(numbers, "prjna678459_primary_screen_pass"),
        n_int(numbers, "prjna1103397_primary_screen_pass"),
    ]
    provisional_vals = [0, n_int(numbers, "prjna678459_provisional_screen_pass"), 0]
    y = list(range(len(labels)))
    ax.barh(y, primary_vals, color=sem["primary"], height=0.62, label="PRIMARY passed screen")
    ax.barh(y, provisional_vals, left=primary_vals, color=sem["provisional"], height=0.62, label="PROVISIONAL passed screen")
    ax.set_yticks(y, labels)
    ax.invert_yaxis()
    ax.set_xlim(0, 8)
    ax.set_xlabel("Loci passing frozen screen")
    ax.set_title("Locus-screen attrition", pad=6)
    clean_ax(ax, "x")
    for i, v in enumerate(primary_vals):
        ax.text(v + 0.10, i, str(v), va="center", ha="left", fontsize=8)
    if provisional_vals[1]:
        ax.text(primary_vals[1] + provisional_vals[1] + 0.10, 1, str(provisional_vals[1]), va="center", ha="left", fontsize=8)
    ax.legend(frameon=False, fontsize=6.8, loc="lower right")
    panel_label(ax, "C", style)

    ax = axes[1, 1]
    ax.axis("off")
    ax.set_title("Strict endpoint status", pad=6)
    panel_label(ax, "D", style)

    ax.text(0.02, 0.88, "PRJNA360417 NAG vs IM", transform=ax.transAxes, fontsize=8.4, fontweight="bold", va="top")
    ax.text(0.02, 0.76, "Strict: NOT_COMPUTABLE", transform=ax.transAxes, fontsize=8.2, color=sem["failure"], va="top")
    ax.text(
        0.02, 0.66,
        f"Below-floor diagnostic exact p = {n_float(numbers, 'prjna360417_diagnostic_exact_p'):.6f}",
        transform=ax.transAxes, fontsize=7.5, color=style["canvas"]["muted_text"], va="top"
    )

    ax.text(0.02, 0.49, "PRJNA678459 AG vs GC", transform=ax.transAxes, fontsize=8.4, fontweight="bold", va="top")
    ax.text(0.02, 0.37, "Strict: NOT_COMPUTABLE", transform=ax.transAxes, fontsize=8.2, color=sem["failure"], va="top")
    ax.text(
        0.02, 0.27,
        f"Below-floor diagnostic exact p = {n_float(numbers, 'prjna678459_diagnostic_exact_p'):.6f}",
        transform=ax.transAxes, fontsize=7.5, color=style["canvas"]["muted_text"], va="top"
    )

    footer = (
        "Coverage floor: >=10 jointly callable loci/sample-pair.\n"
        f"PRJNA678459 structure screen: R²={n_float(numbers, 'structure_prjna678459_r2'):.4f}, "
        f"p={n_float(numbers, 'structure_prjna678459_exact_p'):.6f} — confounding caution only.\n"
        f"PRJNA360417 structure: {numbers['structure_prjna360417']['value']}."
    )
    ax.text(0.02, 0.04, footer, transform=ax.transAxes, fontsize=6.8, color=style["canvas"]["muted_text"], va="bottom")

    fig.suptitle("Figure 4. Public-cohort stress test under frozen rules", fontsize=12)
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
