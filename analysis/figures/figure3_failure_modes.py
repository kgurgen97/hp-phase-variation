#!/usr/bin/env python3
from __future__ import annotations

from collections import Counter

import matplotlib.pyplot as plt

from io import manuscript_numbers, n_int, read_tsv, write_tsv
from style import P, apply_style, clean_axis, panel_label, save_figure, source_path

apply_style()
N = manuscript_numbers()

universe = read_tsv("metadata/loci/analysis_universe.tsv")
ineligible = [r for r in universe if r["g7_locus_class"] == "PRIMARY_INELIGIBLE"]

reason_counts = Counter()
for r in ineligible:
    reasons = r["reasons"]
    if reasons and reasons != "NA":
        for token in reasons.split(";"):
            reason_counts[token] += 1

def pretty_reason(x):
    mapping = {
        "BOUNDARY_AUDIT_ORTHOLOGY_AMBIGUOUS": "Orthology ambiguous",
        "BOUNDARY_AUDIT_INSUFFICIENT_TO_RESOLVE": "Insufficient to resolve",
        "BOUNDARY_AUDIT_BOUNDARY_AMBIGUOUS": "Boundary ambiguous",
        "NOT_CALLABLE_OR_NOT_CALLER_ELIGIBLE_UNDER_G5": "Not caller-eligible / callable",
        "G5_PLATFORM_CONFOUNDED_CALLABILITY": "Platform/study-confounded callability",
        "KNOWN_PV_INSUFFICIENT_EMPIRICAL_VALIDATION": "KNOWN_PV insufficient validation",
        "MOTIF_CLASS_NOT_PRIMARY_VALIDATED": "Motif class not primary-validated",
    }
    return mapping.get(x, x.replace("_", " ").title())

top_reasons = reason_counts.most_common(6)
motif = Counter(r["motif_class"] for r in universe if r["g7_locus_class"] == "PRIMARY_TECHNICAL")
motif_order = sorted(motif.items(), key=lambda x: (-x[1], x[0]))

platform_total = n_int(N, "caller_eligible_loci")
platform_conf = n_int(N, "platform_study_confounded_callability_loci")
known_total = 8
known_primary = n_int(N, "primary_known_pv_loci")

mstatus = {r["motif_class"]: r for r in read_tsv("metadata/validation/motif_validation_status.tsv")}
known_status = mstatus["KNOWN_PV"]["status_after_rescue"]

src=[]
for status, metric in [
    ("Boundary stable","boundary_stable_loci"),
    ("Boundary ambiguous","boundary_ambiguous_loci"),
    ("Orthology ambiguous","orthology_ambiguous_loci"),
    ("Insufficient to resolve","insufficient_to_resolve_loci")]:
    src.append({"panel":"A","category":status,"value":n_int(N,metric)})
src += [
    {"panel":"B","category":"platform_study_confounded","value":platform_conf},
    {"panel":"B","category":"other_caller_eligible","value":platform_total-platform_conf},
    {"panel":"D","category":"KNOWN_PV_catalogue","value":known_total},
    {"panel":"D","category":"KNOWN_PV_primary","value":known_primary},
]
for reason,count in top_reasons:
    src.append({"panel":"C","category":reason,"value":count})
for key,val in motif_order:
    src.append({"panel":"E","category":key,"value":val})
write_tsv(source_path("figure3_source.tsv"),src,["panel","category","value"])

fig=plt.figure(figsize=(7.2,7.2))
gs=fig.add_gridspec(3,2,hspace=0.55,wspace=0.46)

ax=fig.add_subplot(gs[0,0]); panel_label(ax,"A")
cats=["BOUNDARY_STABLE","BOUNDARY_AMBIGUOUS","ORTHOLOGY_AMBIGUOUS","INSUFFICIENT_TO_RESOLVE"]
vals=[n_int(N,"boundary_stable_loci"),n_int(N,"boundary_ambiguous_loci"),
      n_int(N,"orthology_ambiguous_loci"),n_int(N,"insufficient_to_resolve_loci")]
cols=[P["primary_blue"],P["warning_coral"],P["ambiguity_purple"],P["neutral_gray"]]
ax.barh([3,2,1,0],vals,color=cols,height=0.62)
for y,v in zip([3,2,1,0],vals): ax.text(v+15,y,str(v),va="center",fontsize=7.5)
ax.set_yticks([3,2,1,0],["Stable","Boundary\nambiguous","Orthology\nambiguous","Insufficient\nto resolve"])
ax.set_xlim(0,550);ax.set_xlabel("Loci");ax.set_title("Catalogue boundary / orthology failure modes",loc="left",fontweight="bold")
clean_axis(ax,"x")

ax=fig.add_subplot(gs[0,1]); panel_label(ax,"B")
ax.barh([0],[platform_conf],color=P["warning_coral"],height=0.42)
ax.barh([0],[platform_total-platform_conf],left=[platform_conf],color=P["ineligible_gray"],height=0.42)
ax.set_xlim(0,platform_total);ax.set_yticks([]);ax.set_xlabel("Caller-eligible loci")
ax.set_title("Platform/study-confounded callability",loc="left",fontweight="bold")
ax.text(platform_conf/2,0,f"{platform_conf}\nconfounded",ha="center",va="center",fontsize=8,fontweight="bold",color="white")
ax.text(platform_conf+(platform_total-platform_conf)/2,0,f"{platform_total-platform_conf}\nother",ha="center",va="center",fontsize=8)
ax.text(0.5,-0.42,"Association/confounding label; not a causal platform effect",transform=ax.transAxes,ha="center",fontsize=6.8)
clean_axis(ax)

ax=fig.add_subplot(gs[1,:]); panel_label(ax,"C")
labs=[pretty_reason(k) for k,_ in top_reasons][::-1]
vals=[v for _,v in top_reasons][::-1]
ax.barh(range(len(labs)),vals,color=P["neutral_gray"],height=0.62)
for i,v in enumerate(vals): ax.text(v+10,i,str(v),va="center",fontsize=7)
ax.set_yticks(range(len(labs)),labs)
ax.set_xlabel("PRIMARY_INELIGIBLE loci carrying reason (reasons can overlap)")
ax.set_title("Most frequent technical exclusion reasons",loc="left",fontweight="bold")
clean_axis(ax,"x")

ax=fig.add_subplot(gs[2,0]); panel_label(ax,"D")
ax.set_axis_off()
ax.text(0.02,0.87,"Literature-anchored KNOWN_PV",fontsize=9,fontweight="bold",transform=ax.transAxes)
ax.text(0.02,0.61,f"{known_primary} / {known_total}",fontsize=24,fontweight="bold",
        color=P["known_pv_gold"],transform=ax.transAxes)
ax.text(0.02,0.43,"entered PRIMARY_TECHNICAL",fontsize=8,transform=ax.transAxes)
ax.text(0.02,0.20,f"HiFi motif-validation status:\n{known_status}",fontsize=7.2,transform=ax.transAxes)
ax.text(0.98,0.16,"Biological evidence did not override\ntechnical exclusion.",
        ha="right",va="bottom",fontsize=7.2,transform=ax.transAxes,color=P["axis"])

ax=fig.add_subplot(gs[2,1]); panel_label(ax,"E")
labels=[x[0] for x in motif_order][::-1]; vals=[x[1] for x in motif_order][::-1]
ax.barh(range(len(labels)),vals,color=P["primary_blue"],height=0.58)
for i,v in enumerate(vals): ax.text(v+2,i,str(v),va="center",fontsize=7)
ax.set_yticks(range(len(labels)),labels)
ax.set_xlabel("PRIMARY_TECHNICAL loci")
ax.set_title("Motif composition of PRIMARY",loc="left",fontweight="bold")
clean_axis(ax,"x")

fig.suptitle("Failure modes define the usable repeat-locus universe",x=0.05,ha="left",y=1.005,fontsize=11,fontweight="bold")
save_figure(fig,"figure3_failure_modes")
plt.close(fig)
