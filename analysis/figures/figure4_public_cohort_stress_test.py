#!/usr/bin/env python3
from __future__ import annotations

from collections import Counter

import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, PowerNorm

from io import manuscript_numbers, n_float, n_int, read_tsv, write_tsv
from style import P, apply_style, clean_axis, panel_label, save_figure, source_path

apply_style()
N=manuscript_numbers()
stages=["NAG","AG","IM","DYS","GC"]
study=read_tsv("results/confounding/stage_by_study.tsv")
study=sorted(study,key=lambda r:(-int(r["total_correa"]),r["bioproject"]))
matrix=[[int(r[s]) for s in stages] for r in study]
study_names=[r["bioproject"] for r in study]

platform=read_tsv("results/confounding/stage_by_platform.tsv")
platform_totals=[(r["platform_classes"],int(r["total_correa"])) for r in platform]
platform_totals=sorted(platform_totals,key=lambda x:(-x[1],x[0]))

src=[]
for r in study:
    for s in stages:
        src.append({"panel":"A","category":r["bioproject"],"subgroup":s,"value":r[s]})
for k,v in platform_totals:
    src.append({"panel":"B","category":k,"subgroup":"all","value":v})
src += [
    {"panel":"B","category":"single_country_studies","subgroup":"all","value":n_int(N,"confounding_single_country_studies")},
    {"panel":"B","category":"multicountry_studies","subgroup":"all","value":n_int(N,"confounding_multicountry_studies")},
    {"panel":"C","category":"PRJNA360417_primary_pass","subgroup":"4/261","value":4},
    {"panel":"C","category":"PRJNA678459_primary_pass","subgroup":"6/261","value":6},
    {"panel":"C","category":"PRJNA678459_provisional_pass","subgroup":"1/39","value":1},
    {"panel":"C","category":"PRJNA1103397_primary_pass","subgroup":"0/261","value":0},
    {"panel":"D","category":"PRJNA360417_diagnostic_p","subgroup":"diagnostic_only","value":N["prjna360417_diagnostic_exact_p"]["value"]},
    {"panel":"D","category":"PRJNA678459_diagnostic_p","subgroup":"diagnostic_only","value":N["prjna678459_diagnostic_exact_p"]["value"]},
    {"panel":"E","category":"PRJNA678459_structure_R2","subgroup":"SCREENING_STRUCTURE_ONLY","value":N["structure_prjna678459_r2"]["value"]},
    {"panel":"E","category":"PRJNA678459_structure_p","subgroup":"SCREENING_STRUCTURE_ONLY","value":N["structure_prjna678459_exact_p"]["value"]},
]
write_tsv(source_path("figure4_source.tsv"),src,["panel","category","subgroup","value"])

fig=plt.figure(figsize=(7.4,9.1))
gs=fig.add_gridspec(3,2,height_ratios=[2.15,1.0,1.0],hspace=0.55,wspace=0.42)

ax=fig.add_subplot(gs[0,:]); panel_label(ax,"A")
cmap=LinearSegmentedColormap.from_list("bg_blues",[P["off_white"],P["light_blue"],P["primary_blue_dark"]])
im=ax.imshow(matrix,aspect="auto",cmap=cmap,norm=PowerNorm(gamma=0.55,vmin=0,vmax=max(max(x) for x in matrix)))
ax.set_xticks(range(len(stages)),stages)
ax.set_yticks(range(len(study_names)),study_names,fontsize=4.8)
ax.tick_params(length=0)
ax.set_title("Disease-stage labels are concentrated within studies",loc="left",fontweight="bold")
cb=fig.colorbar(im,ax=ax,fraction=0.018,pad=0.015)
cb.set_label("BioSamples",fontsize=7); cb.ax.tick_params(labelsize=6,width=0.4)
ax.text(1.0,1.035,"45 studies · 793 BioSamples · 41 single-country studies",
        transform=ax.transAxes,ha="right",va="bottom",fontsize=8,fontweight="bold")

ax=fig.add_subplot(gs[1,0]); panel_label(ax,"B")
top=platform_totals
labels=[x[0] for x in top][::-1];vals=[x[1] for x in top][::-1]
ax.barh(range(len(labels)),vals,color=P["neutral_gray"],height=0.6)
for i,v in enumerate(vals): ax.text(v+8,i,str(v),va="center",fontsize=7)
ax.set_yticks(range(len(labels)),labels)
ax.set_xlabel("Correa-labelled BioSamples")
ax.set_title("Sequencing-configuration composition",loc="left",fontweight="bold")
clean_axis(ax,"x")
ax.text(0.98,0.04,"Country structure: 41/45 single-country\n4/45 multi-country",
        transform=ax.transAxes,ha="right",va="bottom",fontsize=7.2)

ax=fig.add_subplot(gs[1,1]); panel_label(ax,"C")
cohorts=["PRJNA360417\nNAG vs IM\n5 vs 6","PRJNA678459\nAG vs GC\n5 vs 5","PRJNA1103397\nGC only\n6 isolates"]
passed=[4,6,0]
den=[261,261,261]
cols=[P["primary_blue"],P["primary_blue"],P["ineligible_gray"]]
ax.bar(range(3),passed,color=cols,width=0.55)
for i,(v,d) in enumerate(zip(passed,den)):
    ax.text(i,max(v,0)+0.25,f"{v}/{d}",ha="center",va="bottom",fontsize=8,fontweight="bold")
ax.text(1,7.2,"+ 1/39 PROVISIONAL",ha="center",fontsize=7,color=P["provisional_teal"],fontweight="bold")
ax.set_xticks(range(3),cohorts,fontsize=6.8)
ax.set_ylim(0,9);ax.set_ylabel("PRIMARY loci passing screen")
ax.set_title("Disease-blind screen attrition",loc="left",fontweight="bold")
clean_axis(ax,"y")

ax=fig.add_subplot(gs[2,0]); panel_label(ax,"D")
ax.set_axis_off()
ax.text(0.02,0.90,"Strict cohort endpoint",fontsize=9,fontweight="bold",transform=ax.transAxes)
for y,cohort,pv in [
    (0.62,"PRJNA360417",N["prjna360417_diagnostic_exact_p"]["value"]),
    (0.28,"PRJNA678459",N["prjna678459_diagnostic_exact_p"]["value"])]:
    ax.text(0.02,y,cohort,fontsize=8,fontweight="bold",transform=ax.transAxes)
    ax.text(0.42,y,"NOT_COMPUTABLE",fontsize=9,fontweight="bold",color=P["neutral_gray"],transform=ax.transAxes)
    ax.text(0.42,y-0.13,f"diagnostic exact p = {pv}",fontsize=6.8,color=P["axis"],transform=ax.transAxes)
ax.text(0.02,0.02,"Pre-specified floor: ≥10 jointly callable screened loci per sample pair",fontsize=6.8,transform=ax.transAxes)

ax=fig.add_subplot(gs[2,1]); panel_label(ax,"E")
ax.set_axis_off()
ax.text(0.02,0.90,"Population-structure caution",fontsize=9,fontweight="bold",transform=ax.transAxes)
ax.text(0.02,0.66,"PRJNA678459",fontsize=8,fontweight="bold",transform=ax.transAxes)
ax.text(0.02,0.48,f"pseudo-F {N['structure_prjna678459_pseudo_f']['value']}   R² {N['structure_prjna678459_r2']['value']}",
        fontsize=8,transform=ax.transAxes)
ax.text(0.02,0.31,f"exact p {N['structure_prjna678459_exact_p']['value']}   252 permutations",
        fontsize=8,transform=ax.transAxes)
ax.text(0.02,0.14,"SCREENING_STRUCTURE_ONLY\nconfounding caution; not a repeat-disease finding",
        fontsize=7,color=P["warning_coral"],fontweight="bold",transform=ax.transAxes)
ax.text(0.98,0.14,"PRJNA360417\nNOT_PERFORMED",ha="right",fontsize=7.5,color=P["neutral_gray"],fontweight="bold",transform=ax.transAxes)

fig.suptitle("Public gastric-disease cohorts stress-test the frozen framework",x=0.05,ha="left",y=0.995,fontsize=11,fontweight="bold")
save_figure(fig,"figure4_public_cohort_stress_test")
plt.close(fig)
