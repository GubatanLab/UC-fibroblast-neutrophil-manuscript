from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from PIL import Image

ROOT=Path("giotto_codex_results")
FIG=ROOT/"figures"
NF=ROOT/"neutrophil_spatial_analyses/figures"

def crop_white(path,pad=20):
    im=(path if isinstance(path,Image.Image) else Image.open(path)).convert("RGB")
    a=np.asarray(im)
    mask=np.any(a<248,axis=2)
    yy,xx=np.where(mask)
    if not len(xx): return im
    x0=max(0,xx.min()-pad); x1=min(a.shape[1],xx.max()+pad+1)
    y0=max(0,yy.min()-pad); y1=min(a.shape[0],yy.max()+pad+1)
    return im.crop((x0,y0,x1,y1))

core=crop_white(FIG/"Figure_22_CODEX_core_panels_compact.png",pad=28)
main_panels=[
    ("F",NF/"Figure_N04_neutrophil_subtypes_by_neighborhood.png"),
    ("G",NF/"Figure_N06_neutrophil_compartment_localization.png"),
    ("H",NF/"Figure_N07_vascular_to_epithelial_gradient.png"),
    ("I",NF/"Figure_N08_fibroblast_proximity_marker_states.png"),
]
supp_panels=[
    ("A",NF/"Figure_N01_neutrophil_fibroblast_proximity.png"),
    ("B",NF/"Figure_N05_neutrophil_spatial_aggregation.png"),
    ("C",NF/"Figure_N02_neutrophil_fibroblast_state_enrichment.png"),
    ("D",NF/"Figure_N03_distance_dependent_fibroblast_association.png"),
]

w,h=core.size
ab=crop_white(core.crop((0,0,w,int(h*.282))),pad=8)
crow=crop_white(FIG/"Figure_22_panel_C_compact.png",pad=6)
drow=crop_white(FIG/"Figure_22_panel_D_compact_CODEX.png",pad=6)
erow=crop_white(FIG/"Figure_22_panel_E_compact_CODEX.png",pad=6)

fig=plt.figure(figsize=(30,27),facecolor="white")
ax=fig.add_axes([.025,.735,.95,.25]); ax.imshow(ab); ax.axis("off")
ax=fig.add_axes([.025,.535,.95,.18]); ax.imshow(crow); ax.axis("off")
ax=fig.add_axes([.025,.285,.64,.225]); ax.imshow(drow); ax.axis("off")
ax=fig.add_axes([.025,.035,.64,.225]); ax.imshow(erow); ax.axis("off")
right_y=[.405,.280,.155,.030]
for (lab,path),y in zip(main_panels,right_y):
    ax=fig.add_axes([.685,y,.30,.112])
    ax.imshow(crop_white(path,pad=6)); ax.axis("off")
    ax.text(.006,.985,lab,transform=ax.transAxes,va="top",ha="left",fontsize=25,
            fontweight="bold",color="black",bbox=dict(facecolor="white",edgecolor="none",pad=1.5),zorder=10)

out=FIG/"Figure_22_CODEX_manuscript_multipanel"
fig.savefig(out.with_suffix(".png"),dpi=300,bbox_inches="tight",facecolor="white")
fig.savefig(out.with_suffix(".pdf"),bbox_inches="tight",facecolor="white")
plt.close(fig)

fig=plt.figure(figsize=(30,16),facecolor="white")
gs=fig.add_gridspec(2,2,hspace=.08,wspace=.05,left=.025,right=.985,top=.975,bottom=.025)
for i,(lab,path) in enumerate(supp_panels):
    ax=fig.add_subplot(gs[i//2,i%2]); ax.imshow(crop_white(path,pad=8)); ax.axis("off")
    ax.text(.006,.985,lab,transform=ax.transAxes,va="top",ha="left",fontsize=26,
            fontweight="bold",color="black",bbox=dict(facecolor="white",edgecolor="none",pad=1.5),zorder=10)
out=FIG/"Figure_S_CODEX_neutrophil_spatial_analyses"
fig.savefig(out.with_suffix(".png"),dpi=300,bbox_inches="tight",facecolor="white")
fig.savefig(out.with_suffix(".pdf"),bbox_inches="tight",facecolor="white")
plt.close(fig)
print("Saved main Figure 22 (A-I) and supplementary neutrophil figure (A-D)")
