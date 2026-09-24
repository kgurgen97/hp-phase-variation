from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[2]
STYLE = json.loads((ROOT / "config/figure_style.json").read_text())
P = STYLE["palette"]


def apply_style():
    m = STYLE["matplotlib"]
    plt.rcParams.update({
        "font.family": m["font_family"],
        "font.size": m["font_size"],
        "axes.labelsize": m["axes_labelsize"],
        "axes.titlesize": m["axes_titlesize"],
        "xtick.labelsize": m["xtick_labelsize"],
        "ytick.labelsize": m["ytick_labelsize"],
        "legend.fontsize": m["legend_fontsize"],
        "axes.linewidth": m["axes_linewidth"],
        "xtick.major.width": m["tick_width"],
        "ytick.major.width": m["tick_width"],
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
