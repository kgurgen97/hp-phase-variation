#!/usr/bin/env python3
from __future__ import annotations

import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch

from io import manuscript_numbers, n_int, write_tsv
from style import P, apply_style, clean_axis, panel_label, save_figure, source_path

apply_style()
N = manuscript_numbers()

catalogue = n_int(N, "catalogue_loci")
eligible = n_int(N, "caller_eligible_loci")
stable = n_int(N, "boundary_stable_loci")
primary = n_int(N, "primary_technical_loci")
provisional = n_int(N, "provisional_technical_loci")
ineligible = n_int(N, "primary_ineligible_loci")

audit = [
    ("Boundary stable", n_int(N, "boundary_stable_loci"), P["primary_blue"]),
    ("Boundary ambiguous", n_int(N, "boundary_ambiguous_loci"), P["warning_coral"]),
    ("Orthology ambiguous", n_int(N, "orthology_ambiguous_loci"), P["ambiguity_purple"]),
    ("Insufficient to resolve", n_int(N, "insufficient_to_resolve_loci"), P["neutral_gray"]),
]
primary_evidence = [
    ("REPEAT_ONLY", n_int(N, "primary_repeat_only_loci"), P["neutral_gray"]),
    ("STRONG_PV_CANDIDATE", n_int(N, "primary_strong_pv_candidate_loci"), P["candidate_magenta"]),
    ("Mixed member-level evidence", n_int(N, "primary_mixed_evidence_loci"), P["ambiguity_purple"]),
]

source_rows = []
for stage, value in [
    ("catalogue_loci", catalogue),
    ("caller_eligible_loci", eligible),
    ("boundary_stable_loci", stable),
    ("primary_technical_loci", primary),
    ("provisional_technical_loci", provisional),
    ("primary_ineligible_loci", ineligible),
]:
    source_rows.append({"panel": "A/C", "category": stage, "value": value})
for label, value, _ in audit:
    source_rows.append({"panel": "B", "category": label, "value": value})
for label, value, _ in primary_evidence:
    source_rows.append({"panel": "D", "category": label, "value": value})
source_rows.append({"panel": "D", "category": "KNOWN_PV_in_PRIMARY", "value": n_int(N, "primary_known_pv_loci")})
source_rows.append({"panel": "D", "category": "KNOWN_PV_catalogue_total", "value": 8})
write_tsv(source_path("figure1_source.tsv"), source_rows, ["panel", "category", "value"])

fig = plt.figure(figsize=(7.2, 5.4))
gs = fig.add_gridspec(2, 2, height_ratios=[1.0, 1.05], hspace=0.42, wspace=0.38)

# A: process flow
ax = fig.add_subplot(gs[0, 0])
panel_label(ax, "A")
ax.set_axis_off()
boxes = [
    (0.04, 0.58, 0.24, 0.25, "Repeat catalogue", f"{catalogue:,} loci\n4,009 members", P["light_blue"]),
    (0.38, 0.58, 0.24, 0.25, "Caller-eligible", f"{eligible:,} loci", P["primary_blue"]),
    (0.72, 0.58, 0.24, 0.25, "Boundary stable", f"{stable:,} loci", P["support_green"]),
    (0.38, 0.08, 0.24, 0.25, "PRIMARY", f"{primary} loci", P["primary_blue"]),
    (0.72, 0.08, 0.24, 0.25, "PROVISIONAL", f"{provisional} loci", P["provisional_teal"]),
]
for x, y, w, h, title, value, color in boxes:
    patch = FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.012,rounding_size=0.02",
                           linewidth=0.7, edgecolor=P["axis"], facecolor=color, alpha=0.22,
                           transform=ax.transAxes)
    ax.add_patch(patch)
    ax.text(x+w/2, y+h*0.66, title, ha="center", va="center", transform=ax.transAxes,
            fontweight="bold", fontsize=8)
    ax.text(x+w/2, y+h*0.32, value, ha="center", va="center", transform=ax.transAxes, fontsize=8)
arrow = dict(arrowstyle="-|>", lw=0.8, color=P["axis"], mutation_scale=9)
ax.annotate("", xy=(0.37,0.705), xytext=(0.29,0.705), xycoords=ax.transAxes, arrowprops=arrow)
ax.annotate("", xy=(0.71,0.705), xytext=(0.63,0.705), xycoords=ax.transAxes, arrowprops=arrow)
ax.annotate("", xy=(0.50,0.35), xytext=(0.78,0.56), xycoords=ax.transAxes, arrowprops=arrow)
ax.annotate("", xy=(0.78,0.35), xytext=(0.82,0.56), xycoords=ax.transAxes, arrowprops=arrow)
ax.text(0.5, 0.94, "Technical filtering, not disease association", ha="center",
        va="center", transform=ax.transAxes, fontsize=8.5)

# B: audit
ax = fig.add_subplot(gs[0, 1])
panel_label(ax, "B")
labels = [x[0] for x in audit][::-1]
vals = [x[1] for x in audit][::-1]
cols = [x[2] for x in audit][::-1]
y = range(len(labels))
ax.barh(list(y), vals, color=cols, edgecolor="none", height=0.68)
for i, v in enumerate(vals):
    ax.text(v + 18, i, f"{v}", va="center", fontsize=7.5)
ax.set_yticks(list(y), labels)
ax.set_xlim(0, 570)
ax.set_xlabel("Loci")
ax.set_title("Boundary / orthology audit", loc="left", fontweight="bold")
clean_axis(ax, "x")

# C: final universe
ax = fig.add_subplot(gs[1, 0])
panel_label(ax, "C")
cats = ["PRIMARY_TECHNICAL", "PROVISIONAL_TECHNICAL", "PRIMARY_INELIGIBLE"]
vals = [primary, provisional, ineligible]
cols = [P["primary_blue"], P["provisional_teal"], P["ineligible_gray"]]
ax.barh([2,1,0], vals, color=cols, height=0.62, edgecolor="none")
for yy, v in zip([2,1,0], vals):
    ax.text(v+18, yy, f"{v:,}", va="center", fontsize=8)
ax.set_yticks([2,1,0], cats)
ax.set_xlim(0, 1220)
ax.set_xlabel("Loci")
ax.set_title("Final technical universe", loc="left", fontweight="bold")
clean_axis(ax, "x")

# D: primary evidence composition
ax = fig.add_subplot(gs[1, 1])
panel_label(ax, "D")
labs = [x[0] for x in primary_evidence]
vals = [x[1] for x in primary_evidence]
cols = [x[2] for x in primary_evidence]
left = 0
for lab, val, col in zip(labs, vals, cols):
    ax.barh([0], [val], left=left, color=col, edgecolor="white", linewidth=0.8, height=0.42, label=lab)
    if val >= 10:
        ax.text(left + val/2, 0, f"{val}", ha="center", va="center",
                fontsize=8, color="white" if col != P["neutral_gray"] else P["text"], fontweight="bold")
    left += val
ax.set_xlim(0, primary)
ax.set_yticks([])
ax.set_xlabel("PRIMARY_TECHNICAL loci")
ax.set_title("Evidence class within PRIMARY", loc="left", fontweight="bold")
ax.legend(frameon=False, loc="upper center", bbox_to_anchor=(0.5, -0.20), ncol=1)
ax.text(0.98, 0.78, "KNOWN_PV in PRIMARY\n0 / 8",
        transform=ax.transAxes, ha="right", va="center", fontsize=10, fontweight="bold",
        color=P["known_pv_gold"])
ax.text(0.98, 0.55, "Evidence class ≠ technical eligibility",
        transform=ax.transAxes, ha="right", va="center", fontsize=7.5)
clean_axis(ax)

fig.suptitle("Repeat-locus discovery and conservative technical eligibility", x=0.05, ha="left",
             y=1.01, fontsize=11, fontweight="bold")
save_figure(fig, "figure1_catalogue_universe")
plt.close(fig)
