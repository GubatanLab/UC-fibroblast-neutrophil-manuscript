from pathlib import Path
import json
import re
import shutil
import tempfile
import zipfile

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Circle
from scipy.stats import spearmanr
import seaborn as sns
from PIL import Image

Image.MAX_IMAGE_PIXELS = None

ROOT = Path(r"input_data/codex")
RESULTS = ROOT / "giotto_codex_results"
SOURCE = RESULTS / "fap_a5b1_neutrophil_reprogramming" / "tables"
COEX = RESULTS / "a5b1_FAP_high_low_UC_fibroblasts" / "tables"
CELL_CSV = RESULTS / "additional_analyses" / "codex_extended_cells.csv"
RAW_QPTIFF = Path(r"input_data/codex_images\20240420_Pos1_IBD_Scan1-001.qptiff")
CONTROL_ARCHIVE = Path(r"input_data/codex_images\Control-20240731T174939Z-001.zip")
CONTROL_MEMBER = "Control/20240420_Pos2_control_Scan1.qptiff"
OUT = Path(r"input_data/mouse\CODEX_Figure_2")
OUT.mkdir(parents=True, exist_ok=True)
HIGH_TABLES = Path(r"input_data/mouse\CODEX_high_yield_analyses\tables")

GROUPS = ["Control", "UC_Noninflamed", "UC_Inflamed"]
GROUP_LABEL = {"Control":"Control", "UC_Noninflamed":"UC noninflamed", "UC_Inflamed":"UC inflamed"}
GROUP_COLOR = {"Control":"#4C78A8", "UC_Noninflamed":"#C9A227", "UC_Inflamed":"#D95F5F"}
STATES = ["FAP+/a5B1-", "FAP-/a5B1+", "FAP+/a5B1+", "FAP-/a5B1-"]
STATE_LABEL = {"FAP+/a5B1-":"FAP+ / α5β1−", "FAP-/a5B1+":"FAP− / α5β1+",
               "FAP+/a5B1+":"FAP+ / α5β1+", "FAP-/a5B1-":"FAP− / α5β1−"}
STATE_COLOR = {"FAP+/a5B1-":"#E45756", "FAP-/a5B1+":"#4C78A8",
               "FAP+/a5B1+":"#7A5195", "FAP-/a5B1-":"#B8B8B8"}
MARKERS = ["CXCR4", "OSM", "CD16", "CD11b"]
REPRESENTATIVE = {"Control":"P06", "UC_Noninflamed":"P17", "UC_Inflamed":"P12"}
CHANNEL_PAGES_QUARTER = {"DAPI":111, "CD66b":118, "FAP":142, "a5B1":151}
CHANNEL_COLOR = {"DAPI":np.array([0.13,0.25,1.00]), "CD66b":np.array([0.00,0.95,1.00]),
                 "FAP":np.array([1.00,0.00,0.52]), "a5B1":np.array([0.15,1.00,0.20])}
BROAD_ORDER = ["Epithelial Cell","Enteroendocrine Cell","Endothelial Cell","Fibroblast",
               "Macrophage","Dendritic Cell","Neutrophil","CD4 T","CD8 T","TReg",
               "B Cell","Plasma B Cell"]
BROAD_LABEL = {"Epithelial Cell":"Epithelial","Enteroendocrine Cell":"Enteroendocrine",
               "Endothelial Cell":"Endothelial","Fibroblast":"Fibroblast","Macrophage":"Macrophage",
               "Dendritic Cell":"Dendritic","Neutrophil":"Neutrophil","CD4 T":"CD4 T",
               "CD8 T":"CD8 T","TReg":"Treg","B Cell":"B cell","Plasma B Cell":"Plasma B"}
BROAD_COLOR = dict(zip(BROAD_ORDER,["#4E79A7","#A0CBE8","#59A14F","#8CD17D","#E15759","#FF9D9A",
                                      "#F28E2B","#B6992D","#EDC948","#B07AA1","#76B7B2","#9C755F"]))
NEUT_ORDER = ["Neutrophil","Neutrophil_PADI4","Neutrophil_MX1","Neutrophil_CXCR4","Neutrophil_OSM"]
NEUT_LABEL = {"Neutrophil":"Unpolarized","Neutrophil_PADI4":"PADI4+","Neutrophil_MX1":"MX1+",
              "Neutrophil_CXCR4":"CXCR4+","Neutrophil_OSM":"OSM+"}
NEUT_COLOR = {"Neutrophil":"#D9D9D9","Neutrophil_PADI4":"#A56CC1","Neutrophil_MX1":"#28B6A8",
              "Neutrophil_CXCR4":"#56A0D3","Neutrophil_OSM":"#FF5A5F"}

mpl.rcParams.update({
    "font.family":"DejaVu Sans", "font.size":8.5, "axes.titlesize":10,
    "axes.labelsize":8.5, "xtick.labelsize":7.5, "ytick.labelsize":7.5,
    "axes.linewidth":0.8, "pdf.fonttype":42, "ps.fonttype":42,
    "svg.fonttype":"none"
})
sns.set_style("ticks")


def panel_label(ax, letter, x=-0.10, y=1.07):
    ax.text(x, y, letter, transform=ax.transAxes, fontsize=16, fontweight="bold",
            va="top", ha="left", clip_on=False)


def clean(ax):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


def node(ax, xy, wh, text, fc, ec="#333333", fontsize=8.5, weight="normal"):
    x, y = xy; w, h = wh
    patch = FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.018,rounding_size=0.025",
                           facecolor=fc, edgecolor=ec, linewidth=1.0)
    ax.add_patch(patch)
    ax.text(x+w/2, y+h/2, text, ha="center", va="center", fontsize=fontsize,
            fontweight=weight, color="#202020")
    return patch


def arrow(ax, start, end, color="#4A4A4A", lw=1.4, style="-|>", rad=0.0):
    a = FancyArrowPatch(start, end, arrowstyle=style, mutation_scale=10,
                        linewidth=lw, color=color, connectionstyle=f"arc3,rad={rad}")
    ax.add_patch(a)
    return a


def parse_scale_factor(description):
    match = re.search(r"<ScaleFactor>(.*?)</ScaleFactor>", description or "")
    if not match:
        raise ValueError("QPTIFF ScaleFactor was not found")
    return float(match.group(1))


def normalize_plane(array, low, high, gamma=.72):
    x = np.clip((array.astype(np.float32)-low)/max(high-low,1e-6),0,1)
    return np.power(x,gamma)


def composite_planes(planes, scaling):
    rgb=np.zeros((*next(iter(planes.values())).shape,3),dtype=np.float32)
    for channel,array in planes.items():
        signal=normalize_plane(array,*scaling[channel])
        rgb += signal[...,None]*CHANNEL_COLOR[channel]
    return 1-np.exp(-1.10*rgb)


# Load patient-level source tables.
abundance = pd.read_csv(SOURCE / "01_fibroblast_state_abundance_by_patient.csv")
spatial = pd.read_csv(SOURCE / "03_spatial_enrichment_by_patient.csv")
markers = pd.read_csv(SOURCE / "08_neutrophil_marker_proximity_effects_by_patient.csv")
multiscale = pd.read_csv(HIGH_TABLES / "analysis_3_multiscale_patient_effects.csv")
multiscale_stats = pd.read_csv(HIGH_TABLES / "analysis_3_multiscale_patient_statistics.csv")
distance_bins = pd.read_csv(HIGH_TABLES / "analysis_4_patient_distance_bin_marker_gradients.csv")
gradient_stats = pd.read_csv(HIGH_TABLES / "analysis_4_patient_gradient_statistics.csv")
coex = pd.read_csv(COEX / "01_patient_a5B1_summary_median_split.csv")
correlations = pd.read_csv(COEX / "05_patient_FAP_a5B1_correlations.csv")
cohort = pd.read_csv(r"input_data/mouse\CODEX_Control_UC_FAP_neutrophil_interactions\tables\01_patient_cell_counts.csv")
thresholds = pd.read_csv(SOURCE / "00_fibroblast_state_thresholds.csv").set_index("marker")
fap_cut = float(thresholds.loc["FAP", "primary_threshold_clr"])
a5_cut = float(thresholds.loc["a5B1", "primary_threshold_clr"])

# Cell data used for representative maps and continuous co-expression display.
usecols = ["cell_ID","PatientID","Diagnosis2","cell_type","FAPa_cell_clr","a5B1_cell_clr",
           "a5B1_cell_raw","x","y"]
cells = pd.read_csv(CELL_CSV, usecols=usecols, low_memory=False)
cells["is_fibroblast"] = cells.cell_type.eq("Fibroblast")
cells["is_neutrophil"] = cells.cell_type.str.startswith("Neutrophil", na=False)
cells["broad_type"] = np.where(cells.is_neutrophil, "Neutrophil", cells.cell_type)
cells["fibroblast_state"] = "Other"
fap = cells.FAPa_cell_clr.ge(fap_cut); a5 = cells.a5B1_cell_clr.ge(a5_cut)
cells.loc[cells.is_fibroblast & fap & ~a5, "fibroblast_state"] = "FAP+/a5B1-"
cells.loc[cells.is_fibroblast & ~fap & a5, "fibroblast_state"] = "FAP-/a5B1+"
cells.loc[cells.is_fibroblast & fap & a5, "fibroblast_state"] = "FAP+/a5B1+"
cells.loc[cells.is_fibroblast & ~fap & ~a5, "fibroblast_state"] = "FAP-/a5B1-"

fig = plt.figure(figsize=(22, 32.0), facecolor="white")
outer = GridSpec(4, 2, figure=fig, height_ratios=[0.92, 3.58, 1.12, 1.55],
                 width_ratios=[1, 1], hspace=0.27, wspace=0.24,
                 left=0.042, right=0.995, top=0.985, bottom=0.026)

# A — neutrophil abundance.
ax = fig.add_subplot(outer[0,0]); panel_label(ax,"A",-0.10,1.03)
plot = cohort.copy(); plot["Group"] = plot.Diagnosis2.map(GROUP_LABEL)
order = [GROUP_LABEL[g] for g in GROUPS]
sns.boxplot(data=plot, x="Group", y="fraction_all_cells_neutrophils", order=order,
            palette=[GROUP_COLOR[g] for g in GROUPS], width=0.58, showfliers=False,
            linewidth=1.0, ax=ax)
sns.stripplot(data=plot, x="Group", y="fraction_all_cells_neutrophils", order=order,
              color="#202020", size=4.7, jitter=0.13, ax=ax)
ax.set_title("Neutrophil expansion in active UC", loc="left", fontweight="bold")
ax.set_xlabel(""); ax.set_ylabel("Neutrophils (% of all cells)")
ax.yaxis.set_major_formatter(mpl.ticker.PercentFormatter(1, decimals=0)); clean(ax)
meds = plot.groupby("Diagnosis2").fraction_all_cells_neutrophils.median()
for i,g in enumerate(GROUPS):
    ax.text(i, meds[g] + 0.010, f"{100*meds[g]:.1f}%", ha="center", va="bottom", fontweight="bold")
ax.text(0.35,0.97,"Kruskal–Wallis FDR = 0.004\nInflamed vs control FDR = 0.007\nInflamed vs noninflamed FDR = 0.016",
        transform=ax.transAxes, ha="left", va="top", fontsize=7.3,
        bbox=dict(boxstyle="round,pad=0.3", facecolor="white", edgecolor="#BBBBBB", alpha=.9))

# B — co-expression and patient-level ratio.
sub = GridSpecFromSubplotSpec(1,2,subplot_spec=outer[0,1],width_ratios=[1.38,0.72],wspace=0.27)
axc1 = fig.add_subplot(sub[0,0]); panel_label(axc1,"B",-0.14,1.04)
fib_uc = cells[cells.is_fibroblast & cells.Diagnosis2.isin(["UC_Noninflamed","UC_Inflamed"])].copy()
hb = axc1.hexbin(fib_uc.FAPa_cell_clr, fib_uc.a5B1_cell_clr, gridsize=52, mincnt=1,
                 cmap="magma", bins="log", linewidths=0, rasterized=True)
axc1.axvline(fap_cut,color="#E45756",lw=1.2,ls="--"); axc1.axhline(a5_cut,color="#4C78A8",lw=1.2,ls="--")
axc1.set_xlabel("FAP CLR intensity"); axc1.set_ylabel("α5β1 CLR intensity")
axc1.set_title("FAP–α5β1 coupling in UC fibroblasts", loc="left", fontweight="bold")
clean(axc1)
median_rho = correlations.spearman_rho_FAP_a5B1_raw.median()
axc1.text(0.03,0.97,f"35,078 fibroblasts\nMedian patient Spearman ρ = {median_rho:.2f}",transform=axc1.transAxes,
          va="top",ha="left",fontsize=7.4,color="white",
          bbox=dict(boxstyle="round,pad=.25",facecolor="#222222",edgecolor="none",alpha=.78))
cb=fig.colorbar(hb,ax=axc1,fraction=.045,pad=.02); cb.set_label("Cell density (log)",fontsize=7); cb.ax.tick_params(labelsize=6)

axc2 = fig.add_subplot(sub[0,1])
coex_plot = coex.copy(); coex_plot["Group"] = coex_plot.Diagnosis2.map(GROUP_LABEL)
sns.boxplot(data=coex_plot,x="Group",y="log2_raw_a5B1_high_low",order=[GROUP_LABEL[g] for g in GROUPS[1:]],
            palette=[GROUP_COLOR[g] for g in GROUPS[1:]],width=.55,showfliers=False,linewidth=1,ax=axc2)
sns.stripplot(data=coex_plot,x="Group",y="log2_raw_a5B1_high_low",order=[GROUP_LABEL[g] for g in GROUPS[1:]],
              color="#222222",size=4,jitter=.12,ax=axc2)
axc2.axhline(0,color="#555555",lw=.9,ls="--")
axc2.set_title("α5β1 in FAP-high\nvs FAP-low",fontweight="bold")
axc2.set_xlabel(""); axc2.set_ylabel("log2 median ratio")
axc2.tick_params(axis="x",rotation=18); clean(axc2)
axc2.text(.5,.96,"20% higher overall\npaired FDR = 0.022",transform=axc2.transAxes,ha="center",va="top",fontsize=7.3,
          bbox=dict(boxstyle="round,pad=.25",facecolor="white",edgecolor="#BBBBBB"))

# C — representative maps registered to fluorescence and whole-slide image.
control_temp_root=Path(tempfile.mkdtemp(prefix="codex_control_qptiff_"))
with zipfile.ZipFile(CONTROL_ARCHIVE) as archive:
    archive.extract(CONTROL_MEMBER,control_temp_root)
control_qptiff=control_temp_root/CONTROL_MEMBER
uc_raw=Image.open(RAW_QPTIFF); control_raw=Image.open(control_qptiff)
uc_raw.seek(0); full_scale=parse_scale_factor(str(uc_raw.tag_v2.get(270,"")))
control_raw.seek(0); control_scale=parse_scale_factor(str(control_raw.tag_v2.get(270,"")))
quarter_scale=full_scale/4.0
crop_boxes={}; raw_crops={}
for g in GROUPS:
    z=cells[cells.PatientID.eq(REPRESENTATIVE[g])]
    pad=90.0
    xmin=max(0,float(z.x.min()-pad)); xmax=float(z.x.max()+pad)
    ymin=max(0,float(z.y.min()-pad)); ymax=float(z.y.max()+pad)
    box=(int(xmin*quarter_scale),int(ymin*quarter_scale),int(xmax*quarter_scale),int(ymax*quarter_scale))
    crop_boxes[g]=(xmin,xmax,ymin,ymax,box)
    raw=control_raw if g=="Control" else uc_raw
    for channel,page in CHANNEL_PAGES_QUARTER.items():
        raw.seek(page); raw_crops[(g,channel)]=np.asarray(raw.crop(box),dtype=np.uint8)
scaling={}
for channel in CHANNEL_PAGES_QUARTER:
    sampled=np.concatenate([raw_crops[(g,channel)][::8,::8].ravel() for g in GROUPS])
    scaling[channel]=(float(np.percentile(sampled,1.0)),float(np.percentile(sampled,99.7)))
uc_raw.seek(276); uc_slide=np.asarray(uc_raw,dtype=np.uint8)
control_raw.seek(221); control_slide=np.asarray(control_raw,dtype=np.uint8)
uc_raw.close(); control_raw.close(); shutil.rmtree(control_temp_root,ignore_errors=True)
uc_slide_norm=normalize_plane(uc_slide,float(np.percentile(uc_slide,1)),float(np.percentile(uc_slide,99.7)),gamma=.68)
control_slide_norm=normalize_plane(control_slide,float(np.percentile(control_slide,1)),float(np.percentile(control_slide,99.7)),gamma=.68)

panel_c_aspects=[]
for g in GROUPS:
    xmin,xmax,ymin,ymax,_=crop_boxes[g]
    panel_c_aspects.append((xmax-xmin)/(ymax-ymin))
sub = GridSpecFromSubplotSpec(4,4,subplot_spec=outer[1,:],
                              width_ratios=[.30]+panel_c_aspects,
                              height_ratios=[1,1,1,.105],wspace=-0.035,hspace=0.006)
locator=fig.add_subplot(sub[:3,0]); locator.set_xticks([]); locator.set_yticks([]); locator.axis("off")
locator.set_title("Whole-slide locators",fontsize=8.5,fontweight="bold",pad=1)
control_loc=locator.inset_axes([.04,.55,.92,.40]); uc_loc=locator.inset_axes([.04,.04,.92,.46])
control_loc.imshow(control_slide_norm,cmap="gray",vmin=0,vmax=1); uc_loc.imshow(uc_slide_norm,cmap="gray",vmin=0,vmax=1)
for a,title in [(control_loc,"Control scan"),(uc_loc,"UC scan")]:
    a.set_xticks([]); a.set_yticks([]); a.set_title(title,fontsize=7.5,fontweight="bold",pad=2)
    for spine in a.spines.values(): spine.set_color("#777777"); spine.set_linewidth(.6)
for g,target_ax,slide_width,full_width,scan_scale in [
    ("Control",control_loc,control_slide.shape[1],13440,control_scale),
    ("UC_Noninflamed",uc_loc,uc_slide.shape[1],27840,full_scale),
    ("UC_Inflamed",uc_loc,uc_slide.shape[1],27840,full_scale),
]:
    xmin,xmax,ymin,ymax,_=crop_boxes[g]; thumb_scale=slide_width/full_width
    rect=mpl.patches.Rectangle((xmin*scan_scale*thumb_scale,ymin*scan_scale*thumb_scale),
        (xmax-xmin)*scan_scale*thumb_scale,(ymax-ymin)*scan_scale*thumb_scale,
        fill=False,edgecolor=GROUP_COLOR[g],linewidth=2.0)
    target_ax.add_patch(rect)
    target_ax.text(xmin*scan_scale*thumb_scale,ymin*scan_scale*thumb_scale-3,REPRESENTATIVE[g],
                   color=GROUP_COLOR[g],fontsize=7.0,fontweight="bold")

raw_axes=[]; annotation_axes=[]; neutrophil_axes=[]
for j,g in enumerate(GROUPS):
    z=cells[cells.PatientID.eq(REPRESENTATIVE[g])]
    xmin,xmax,ymin,ymax,_=crop_boxes[g]
    planes={ch:raw_crops[(g,ch)] for ch in CHANNEL_PAGES_QUARTER}
    rgb=composite_planes(planes,scaling)

    # Fluorescence plus the original mechanistic fibroblast-state overlay.
    ax=fig.add_subplot(sub[0,j+1]); raw_axes.append(ax)
    ax.imshow(rgb,extent=[xmin,xmax,ymax,ymin],interpolation="nearest",aspect="equal",rasterized=True)
    for state in STATES:
        q=z[z.fibroblast_state.eq(state)]
        ax.scatter(q.x,q.y,s=3.8,facecolors="none",edgecolors=STATE_COLOR[state],linewidths=.38,alpha=.90,rasterized=True)
    q=z[z.is_neutrophil]
    ax.scatter(q.x,q.y,s=3.8,facecolors="#FFFFFF",edgecolors="#161616",linewidths=.30,alpha=.92,rasterized=True)
    ax.set_xlim(xmin,xmax); ax.set_ylim(ymax,ymin); ax.set_xticks([]); ax.set_yticks([])
    ax.set_title(f"{GROUP_LABEL[g]}  |  {REPRESENTATIVE[g]}",fontsize=9,fontweight="bold")
    for spine in ax.spines.values(): spine.set_color(GROUP_COLOR[g]); spine.set_linewidth(1.1)
    bar=500
    bx=xmin+.06*(xmax-xmin); by=ymax-.06*(ymax-ymin)
    ax.plot([bx,bx+bar],[by,by],color="white",lw=3,solid_capstyle="butt")
    ax.text(bx+bar/2,by-.025*(ymax-ymin),"500 µm",color="white",ha="center",va="top",fontsize=6.8)

    # Full segmented-cell annotations at identical coordinates.
    ax=fig.add_subplot(sub[1,j+1]); annotation_axes.append(ax)
    ax.imshow(rgb*.27,extent=[xmin,xmax,ymax,ymin],interpolation="nearest",aspect="equal",rasterized=True)
    for celltype in BROAD_ORDER:
        q=z[z.broad_type.eq(celltype)]
        ax.scatter(q.x,q.y,s=3.1,c=BROAD_COLOR[celltype],alpha=.76,linewidths=0,rasterized=True)
    ax.set_xlim(xmin,xmax); ax.set_ylim(ymax,ymin); ax.set_xticks([]); ax.set_yticks([])
    for spine in ax.spines.values(): spine.set_color(GROUP_COLOR[g]); spine.set_linewidth(1.0)

    # Neutrophil subtype map at identical coordinates.
    ax=fig.add_subplot(sub[2,j+1]); neutrophil_axes.append(ax)
    gray=np.clip(rgb.mean(axis=2),0,1)
    ax.imshow(gray,cmap="gray",vmin=0,vmax=1,extent=[xmin,xmax,ymax,ymin],interpolation="nearest",
              aspect="equal",alpha=.38,rasterized=True)
    for subtype in NEUT_ORDER:
        q=z[z.cell_type.eq(subtype)]
        ax.scatter(q.x,q.y,s=8.5,c=NEUT_COLOR[subtype],alpha=.93,edgecolors="#202020",
                   linewidths=.16,rasterized=True)
    ax.set_xlim(xmin,xmax); ax.set_ylim(ymax,ymin); ax.set_xticks([]); ax.set_yticks([])
    for spine in ax.spines.values(): spine.set_color(GROUP_COLOR[g]); spine.set_linewidth(1.0)
    ax.text(.975,.955,f"n={z.is_neutrophil.sum():,}",transform=ax.transAxes,ha="right",va="top",fontsize=6.2,
            bbox=dict(boxstyle="round,pad=.18",facecolor="white",edgecolor="none",alpha=.82))
panel_label(locator,"C",-0.18,1.06)
locator.text(0,1.09,"Registered CODEX microscopy and spatial annotations",transform=locator.transAxes,
             ha="left",va="bottom",fontsize=10,fontweight="bold")
for ax,label in zip([raw_axes[0],annotation_axes[0],neutrophil_axes[0]],
                    ["Fluorescence + FAP/α5β1 states","Complete cell annotations","Neutrophil subtypes"]):
    ax.text(-.035,.5,label,transform=ax.transAxes,rotation=90,ha="right",va="center",
            fontsize=7.2,fontweight="bold",clip_on=False)
channel_handles=[mpl.lines.Line2D([],[],marker='o',ls='',color=CHANNEL_COLOR[ch],label=ch,markersize=4) for ch in CHANNEL_PAGES_QUARTER]
state_handles=[mpl.lines.Line2D([],[],marker='o',ls='',markerfacecolor='none',markeredgecolor=STATE_COLOR[s],label=STATE_LABEL[s],markersize=4) for s in STATES]
state_handles.append(mpl.lines.Line2D([],[],marker='o',ls='',markerfacecolor='white',markeredgecolor='#161616',label='Neutrophil',markersize=4))
broad_handles=[mpl.lines.Line2D([],[],marker='o',ls='',color=BROAD_COLOR[x],label=BROAD_LABEL[x],markersize=3.8) for x in BROAD_ORDER]
neut_handles=[mpl.lines.Line2D([],[],marker='o',ls='',markerfacecolor=NEUT_COLOR[x],markeredgecolor="#202020",
                              markeredgewidth=.2,label=NEUT_LABEL[x],markersize=4.2) for x in NEUT_ORDER]
legend_axes=[]
for j in range(3):
    la=fig.add_subplot(sub[3,j+1]); la.axis("off"); legend_axes.append(la)
for dx,trio,la in zip([.022,0,-.022],zip(raw_axes,annotation_axes,neutrophil_axes),legend_axes):
    for a in (*trio,la):
        p=a.get_position()
        a.set_position([p.x0+dx,p.y0,p.width,p.height])
legend_axes[0].legend(handles=channel_handles+state_handles,loc="center",ncol=5,
                      frameon=False,fontsize=5.55,handletextpad=.22,columnspacing=.48)
legend_axes[1].legend(handles=broad_handles,loc="center",ncol=6,
                      frameon=False,fontsize=5.25,handletextpad=.18,columnspacing=.42)
legend_axes[2].legend(handles=neut_handles,loc="center",ncol=5,
                      frameon=False,fontsize=5.65,handletextpad=.20,columnspacing=.50)

# D — spatial topology heatmap.
ax = fig.add_subplot(outer[2,0]); panel_label(ax,"D",-0.10,1.04)
e=spatial[spatial.radius_um.eq(50)].groupby(["Diagnosis2","fibroblast_state"]).log2_observed_expected.median().unstack()
e=e.reindex(index=GROUPS,columns=STATES)
sns.heatmap(e,ax=ax,cmap="vlag",center=0,annot=True,fmt=".2f",linewidths=.7,linecolor="white",
            cbar_kws={"label":"Median log2 observed / expected edges","shrink":.78})
ax.set_title("Abundance-adjusted neutrophil topology (50 µm)",loc="left",fontweight="bold")
ax.set_xlabel(""); ax.set_ylabel("")
ax.set_yticklabels([GROUP_LABEL[g] for g in GROUPS],rotation=0)
ax.set_xticklabels([STATE_LABEL[s] for s in STATES],rotation=20,ha="right")
ax.text(0,-.25,"Within-patient label permutations preserve fibroblast abundance and tissue geometry",
        transform=ax.transAxes,fontsize=7.2,color="#4A4A4A")

# E — multiscale spatial interaction curves.
sub = GridSpecFromSubplotSpec(1,2,subplot_spec=outer[2,1],wspace=0.18)
e_axes=[]
for j,state in enumerate(["FAP+/a5B1-", "FAP+/a5B1+"]):
    ax=fig.add_subplot(sub[0,j]); e_axes.append(ax)
    q=multiscale[multiscale.fibroblast_state.eq(state)]
    for g in GROUPS:
        gg=q[q.Diagnosis2.eq(g)]
        agg=gg.groupby("radius_um").log2_observed_expected.agg(
            median="median", q1=lambda x: x.quantile(.25), q3=lambda x: x.quantile(.75)
        ).reset_index()
        ax.fill_between(agg.radius_um,agg.q1,agg.q3,color=GROUP_COLOR[g],alpha=.12,lw=0)
        ax.plot(agg.radius_um,agg["median"],marker="o",ms=4.0,lw=1.5,
                color=GROUP_COLOR[g],label=GROUP_LABEL[g])
    ax.axhline(0,color="#777777",ls="--",lw=.8)
    ax.axvline(50,color="#BBBBBB",ls=":",lw=.9)
    ax.set_xticks([15,25,50,100,150])
    ax.set_xlabel("Neighborhood radius (µm)")
    ax.set_title(STATE_LABEL[state],fontweight="bold",color=STATE_COLOR[state])
    if j==0:
        ax.set_ylabel("Patient-median log2 observed / expected edges")
        ax.legend(frameon=False,fontsize=6.6,loc="best")
    else:
        ax.set_ylabel("")
    clean(ax)
panel_label(e_axes[0],"E",-0.24,1.08)
e_axes[0].text(0,1.16,"Multiscale neutrophil–fibroblast topology",transform=e_axes[0].transAxes,
               ha="left",va="bottom",fontsize=10,fontweight="bold")
e_axes[0].text(0,-.27,"Lines show patient medians; bands show interquartile ranges. Dotted line marks 50 µm.",
               transform=e_axes[0].transAxes,fontsize=6.8,color="#4A4A4A")

# F — compact continuous phenotype gradients, with G placed to the right.
bottom = GridSpecFromSubplotSpec(1, 2, subplot_spec=outer[3, :],
                                 width_ratios=[1.58, 0.70], wspace=0.27)
sub = GridSpecFromSubplotSpec(2, 2, subplot_spec=bottom[0, 0],
                              wspace=0.22, hspace=0.36)
f_axes=[]
bin_order=["0–25","25–50","50–100",">100"]
primary_bins=distance_bins[distance_bins.fibroblast_state.eq("FAP+/a5B1-")].copy()
for j,marker in enumerate(MARKERS):
    row, col = divmod(j, 2)
    ax=fig.add_subplot(sub[row,col],sharey=f_axes[0] if f_axes else None); f_axes.append(ax)
    for g in GROUPS:
        q=primary_bins[(primary_bins.marker.eq(marker)) & primary_bins.Diagnosis2.eq(g)]
        agg=q.groupby("distance_bin").median_marker_robust_z.agg(
            median="median",q1=lambda x:x.quantile(.25),q3=lambda x:x.quantile(.75)
        ).reindex(bin_order)
        xvals=np.arange(len(bin_order))
        ax.fill_between(xvals,agg.q1,agg.q3,color=GROUP_COLOR[g],alpha=.12,lw=0)
        ax.plot(xvals,agg["median"],marker="o",ms=4.0,lw=1.5,color=GROUP_COLOR[g],label=GROUP_LABEL[g])
    ax.axhline(0,color="#888888",ls="--",lw=.8)
    ax.set_xticks(range(4),bin_order,rotation=22,ha="right")
    ax.set_title(marker,fontweight="bold")
    if row == 1:
        ax.set_xlabel("Distance to nearest FAP+ / α5β1− fibroblast (µm)")
    else:
        ax.set_xlabel("")
    if col==0:
        ax.set_ylabel("Patient-median marker intensity\n(within-patient robust z)")
        if j == 0:
            ax.legend(frameon=False,fontsize=6.4,loc="best")
    else:
        ax.tick_params(labelleft=False)
    inflamed=gradient_stats[(gradient_stats.Diagnosis2.eq("UC_Inflamed")) &
                            (gradient_stats.fibroblast_state.eq("FAP+/a5B1-")) &
                            (gradient_stats.marker.eq(marker))]
    if len(inflamed) and float(inflamed.iloc[0].fdr)<.10:
        ax.text(.98,.96,f"UC inflamed q={float(inflamed.iloc[0].fdr):.3f}",transform=ax.transAxes,
                ha="right",va="top",fontsize=6.6,color=GROUP_COLOR["UC_Inflamed"],fontweight="bold")
    clean(ax)
panel_label(f_axes[0],"F",-0.16,1.12)
f_axes[0].text(0,1.20,"Neutrophil phenotype gradients around FAP+ / α5β1− fibroblasts",
               transform=f_axes[0].transAxes,ha="left",va="bottom",fontsize=10,fontweight="bold")

# G — patient-level bridge association.
ax=fig.add_subplot(bottom[0,1]); panel_label(ax,"G",-0.15,1.04)
s=spatial[(spatial.fibroblast_state.eq("FAP+/a5B1-")) & spatial.radius_um.eq(50)]
m=markers[(markers.fibroblast_state.eq("FAP+/a5B1-")) & markers.feature.eq("NF-niche program") &
          markers.method.eq("Proximal-vs-distant median difference")]
gdat=s.merge(m,on=["PatientID","Diagnosis2","fibroblast_state"],suffixes=("_spatial","_program"))
rho,p=spearmanr(gdat.log2_observed_expected,gdat.effect)
for group in GROUPS:
    q=gdat[gdat.Diagnosis2.eq(group)]
    ax.scatter(q.log2_observed_expected,q.effect,s=42,c=GROUP_COLOR[group],edgecolor="white",linewidth=.6,
               label=GROUP_LABEL[group],alpha=.92)
x=gdat.log2_observed_expected.to_numpy(); y=gdat.effect.to_numpy(); coef=np.polyfit(x,y,1)
xx=np.linspace(x.min(),x.max(),100); ax.plot(xx,np.polyval(coef,xx),color="#333333",lw=1.3)
ax.axhline(0,color="#AAAAAA",lw=.8); ax.axvline(0,color="#AAAAAA",lw=.8)
ax.set_xlabel("FAP+ / α5β1− niche enrichment\n(log2 observed / expected edges)")
ax.set_ylabel("Proximal-neutrophil niche program\n(near − far effect)")
ax.set_title("Patient-level spatial–phenotypic coupling",loc="left",fontweight="bold")
ax.text(.03,.96,f"Spearman ρ = {rho:.2f}\nP = {p:.3f} (exploratory)",transform=ax.transAxes,va="top",ha="left",fontsize=7.5,
        bbox=dict(boxstyle="round,pad=.25",facecolor="white",edgecolor="#BBBBBB"))
ax.legend(frameon=False,loc="lower right",fontsize=7); clean(ax)

for ext in ["png","pdf","svg"]:
    path=OUT/f"Figure_2_CODEX_FAP_a5B1_NAMPT.{ext}"
    if ext=="png": fig.savefig(path,dpi=300,facecolor="white")
    else: fig.savefig(path,facecolor="white")
plt.close(fig)

manifest={"source_object":"UCCODEX1_Annotated_repaired_withUMAP.rds",
          "source_uc_qptiff":str(RAW_QPTIFF),"source_control_qptiff_archive":str(CONTROL_ARCHIVE),
          "source_control_qptiff_member":CONTROL_MEMBER,
          "representative_samples":REPRESENTATIVE,"groups":GROUPS,
          "FAP_CLR_threshold":fap_cut,"a5B1_CLR_threshold":a5_cut,
          "panel_C_qptiff_quarter_resolution_pages":CHANNEL_PAGES_QUARTER,
          "panel_C_coordinate_scale_pixels_per_um":quarter_scale,
          "panel_C_patient_bounds_um":{g:list(crop_boxes[g][:4]) for g in GROUPS},
          "panel_E_multiscale_radii_um":[15,25,35,50,75,100,150],
          "panel_F_distance_bins_um":["0-25","25-50","50-100",">100"],
          "panel_G_spearman_rho":float(rho),"panel_G_p_value":float(p),
          "note":"Panel H was removed; Figure 2 now ends with the patient-level spatial-phenotypic association in panel G."}
(OUT/"Figure_2_manifest.json").write_text(json.dumps(manifest,indent=2),encoding="utf-8")
print(OUT / "Figure_2_CODEX_FAP_a5B1_NAMPT.png")
