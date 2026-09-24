#!/usr/bin/env python3
from __future__ import annotations

from collections import defaultdict

import matplotlib.pyplot as plt

from figure_io import manuscript_numbers, n_float, n_int, read_tsv, write_tsv
from style import P, apply_style, clean_axis, panel_label, save_figure, source_path

apply_style()
N = manuscript_numbers()

rows = [
    {"panel":"A","metric":"synthetic_exact","numerator":1757,"denominator":1780,"fraction":1757/1780},
    {"panel":"B","metric":"synthetic_false_confident_wrong","numerator":15,"denominator":1595,"fraction":15/1595},
    {"panel":"D","metric":"hifi_exact","numerator":911,"denominator":913,"fraction":911/913},
    {"panel":"E","metric":"hifi_reference_exact","numerator":829,"denominator":829,"fraction":1.0},
    {"panel":"E","metric":"hifi_nonreference_exact","numerator":82,"denominator":84,"fraction":82/84},
    {"panel":"F","metric":"hifi_high_confidence_truth","numerator":1034,"denominator":5980,"fraction":1034/5980},
    {"panel":"F","metric":"hifi_ambiguous_truth","numerator":1453,"denominator":5980,"fraction":1453/5980},
    {"panel":"F","metric":"hifi_uncallable_truth","numerator":3493,"denominator":5980,"fraction":3493/5980},
]
write_tsv(source_path("figure2_source.tsv"), rows, ["panel","metric","numerator","denominator","fraction"])

callability = [r for r in read_tsv("results/callability/phase1_callability.tsv") if r["record_type"] == "SAMPLE"]
by_study = defaultdict(list)
for r in callability:
    by_study[r["bioproject"]].append(float(r["callable_fraction"]))

fig = plt.figure(figsize=(7.2, 7.0))
gs = fig.add_gridspec(3, 2, hspace=0.55, wspace=0.42)

def accuracy_panel(ax, label, num, den, title, color, error_color=None):
    panel_label(ax, label)
    frac = num/den
    ax.barh([0], [frac], color=color, height=0.36)
    if error_color:
        ax.barh([0], [1-frac], left=[frac], color=error_color, height=0.36)
    ax.set_xlim(0, 1)
    ax.set_yticks([])
    ax.set_xticks([0, .25, .5, .75, 1], ["0", "25", "50", "75", "100"])
    ax.set_xlabel("Percent")
    ax.set_title(title, loc="left", fontweight="bold")
    ax.text(frac/2, 0, f"{num:,}/{den:,}\n{100*frac:.1f}%", ha="center", va="center",
            fontsize=10, fontweight="bold", color="white" if frac > 0.4 else P["text"])
    clean_axis(ax)

accuracy_panel(fig.add_subplot(gs[0,0]), "A", 1757, 1780,
               "Synthetic exact dominant allele", P["primary_blue"], P["ineligible_gray"])

ax=fig.add_subplot(gs[0,1]); panel_label(ax,"B")
wrong=15/1595
ax.barh([0], [wrong], color=P["warning_coral"], height=0.36)
ax.barh([0], [1-wrong], left=[wrong], color=P["ineligible_gray"], height=0.36)
ax.set_xlim(0,0.03); ax.set_yticks([])
ax.set_xticks([0,.01,.02,.03],["0","1","2","3"])
ax.set_xlabel("Percent")
ax.set_title("False-confident wrong dominant allele", loc="left", fontweight="bold")
ax.text(0.029,0,f"15/1,595 = {100*wrong:.2f}%",ha="right",va="center",fontsize=9,fontweight="bold")
clean_axis(ax)

ax=fig.add_subplot(gs[1,0]); panel_label(ax,"C")
order=sorted(by_study)
data=[by_study[k] for k in order]
bp=ax.boxplot(data, vert=True, patch_artist=True, widths=0.55, showfliers=False,
              medianprops={"color":P["text"],"linewidth":1.0},
              boxprops={"edgecolor":P["axis"],"linewidth":0.7},
              whiskerprops={"color":P["axis"],"linewidth":0.7},
              capprops={"color":P["axis"],"linewidth":0.7})
for patch,col in zip(bp["boxes"],[P["light_blue"],P["primary_blue"],P["provisional_teal"]]):
    patch.set_facecolor(col); patch.set_alpha(0.65)
for i,k in enumerate(order,1):
    vals=by_study[k]
    offsets=[(j-(len(vals)-1)/2)*0.012 for j in range(len(vals))]
    ax.scatter([i+o for o in offsets], vals, s=7, color=P["text"], alpha=0.45, zorder=3)
ax.set_xticks(range(1,len(order)+1),order,rotation=20,ha="right")
ax.set_ylim(0,1); ax.set_ylabel("Callable fraction")
ax.set_title("Phase-1 callability under frozen rules",loc="left",fontweight="bold")
clean_axis(ax,"y")

accuracy_panel(fig.add_subplot(gs[1,1]), "D", 911, 913,
               "Short-read vs HiFi exact concordance", P["primary_blue_dark"], P["ineligible_gray"])

ax=fig.add_subplot(gs[2,0]); panel_label(ax,"E")
labs=["Reference length","Non-reference"]
nums=[829,82]; dens=[829,84]; fracs=[1.0,82/84]
bars=ax.barh([1,0], fracs, color=[P["light_blue"],P["candidate_magenta"]], height=0.52)
for y,n,d,f in zip([1,0],nums,dens,fracs):
    ax.text(min(f,0.98)-0.02,y,f"{n}/{d}\n{100*f:.1f}%",ha="right",va="center",
            fontsize=8,fontweight="bold")
ax.set_yticks([1,0],labs); ax.set_xlim(0.90,1.005)
ax.set_xlabel("Exact concordance")
ax.set_title("HiFi reference vs non-reference alleles",loc="left",fontweight="bold")
clean_axis(ax,"x")

ax=fig.add_subplot(gs[2,1]); panel_label(ax,"F")
truth=[("HIGH_CONFIDENCE",1034,P["primary_blue"]),
       ("AMBIGUOUS",1453,P["ambiguity_purple"]),
       ("UNCALLABLE",3493,P["ineligible_gray"])]
left=0
for lab,val,col in truth:
    ax.barh([0],[val],left=left,color=col,height=0.42,edgecolor="white",linewidth=0.8,label=lab)
    if val>700:
        ax.text(left+val/2,0,f"{val:,}",ha="center",va="center",fontsize=8,fontweight="bold",
                color="white" if col!=P["ineligible_gray"] else P["text"])
    left+=val
ax.set_xlim(0,5980); ax.set_yticks([]); ax.set_xlabel("Isolate × locus pairs")
ax.set_title("HiFi truth availability (total 5,980)",loc="left",fontweight="bold")
ax.legend(frameon=False,loc="upper center",bbox_to_anchor=(0.5,-0.23),ncol=1)
clean_axis(ax)

fig.text(0.05,0.005,
         "HiFi values are descriptive per-comparison concordance, not 913 independent biological replicates. "
         "Most comparisons are reference-length; KNOWN_PV and mixture validation remain insufficient.",
         fontsize=6.6,ha="left",va="bottom")
fig.suptitle("Technical validation of dominant repeat-length allele calling",x=0.05,ha="left",
             y=1.01,fontsize=11,fontweight="bold")
save_figure(fig,"figure2_validation")
plt.close(fig)
