from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[2]
STYLE = json.loads((ROOT / "config/figure_style.json").read_text())
RAW = STYLE["palette"]

# Compatibility aliases used by the figure scripts. These are deterministic
# semantic mappings onto the literal, reference-derived palette stored in the
# JSON source of truth; no new colors are invented here.
P = {
    **RAW,
    "off_white": RAW["background"],
    "axis": RAW["muted_text"],
    "grid": RAW["neutral_light"],
    "primary_blue": RAW["primary_indigo"],
    "primary_blue_dark": RAW["text"],
    "provisional_teal": RAW["teal"],
    "candidate_magenta": RAW["mauve"],
    "known_pv_gold": RAW["gold"],
    "ambiguity_purple": RAW["plum"],
    "warning_coral": RAW["orange"],
    "support_green": RAW["green"],
    "neutral_gray": RAW["neutral"],
    "ineligible_gray": RAW["neutral_light"],
    "light_blue": RAW["blue"],
    "light_pink": RAW["salmon"],
}


def apply_style():
    s = STYLE["svg"]
    plt.rcParams.update({
        "font.family": "Arial",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
        "font.size": 8.0,
        "axes.labelsize": 8.0,
        "axes.titlesize": 9.0,
        "xtick.labelsize": 7.0,
        "ytick.labelsize": 7.0,
        "legend.fontsize": 7.0,
        "axes.linewidth": 0.6,
        "xtick.major.width": 0.6,
        "ytick.major.width": 0.6,
        "axes.edgecolor": P["axis"],
        "axes.labelcolor": P["text"],
        "text.color": P["text"],
        "xtick.color": P["axis"],
        "ytick.color": P["axis"],
        "figure.facecolor": P["background"],
        "axes.facecolor": P["background"],
        "savefig.facecolor": P["background"],
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "svg.fonttype": "none",
    })


def clean_axis(ax, grid_axis=None):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_linewidth(0.6)
    ax.spines["bottom"].set_linewidth(0.6)
    if grid_axis:
        ax.grid(axis=grid_axis, color=P["grid"], linewidth=0.4, zorder=0)
    ax.set_axisbelow(True)


def panel_label(ax, label):
    ax.text(-0.12, 1.07, label, transform=ax.transAxes, fontsize=11,
            fontweight="bold", va="top", ha="left", color=P["text"])


def save_figure(fig, stem):
    out = ROOT / "figures" / "main"
    out.mkdir(parents=True, exist_ok=True)
    for ext in ("png", "svg", "pdf"):
        kwargs = {"bbox_inches": "tight"}
        if ext == "png":
            kwargs["dpi"] = 400
        fig.savefig(out / f"{stem}.{ext}", **kwargs)


def source_path(name):
    out = ROOT / "figures" / "source_data"
    out.mkdir(parents=True, exist_ok=True)
    return out / name
