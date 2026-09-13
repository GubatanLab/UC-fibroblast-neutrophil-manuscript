"""Matched, compact Figure 4/5 layouts; source values and embedding geometry retained."""
from pathlib import Path
import math

import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.lines import Line2D
from matplotlib.ticker import MaxNLocator
import numpy as np
import pandas as pd
from PIL import Image

HERE = Path(__file__).resolve().parent
DATA = HERE / "source_data"
W, H = 183., 170.
INK, MUTED, GRID = "#25282C", "#5F666D", "#E1E5E8"
BLUE, RED, PURPLE, TEAL, ORANGE = "#3B6FB6", "#C84C4C", "#8E24AA", "#168C88", "#D97724"
COMPARTMENTS = ["Epithelial", "Immune", "Stromal"]
UMAP_FILES = {4: "FAP_ablation_UMAP_coordinates_and_annotations.csv.gz", 5: "figure5_umap_display_cells.csv.gz"}
CONFIG = {
    4: dict(title="FAP ablation reduces DSS-colitis pathology", image="FAPTK DSS.png",
            boxes=[(0,0,826,704),(832,0,1658,704),(1667,0,2493,704)],
            labels=["FAP-TK + PBS", "FAP-TK + PBS + DSS", "FAP-TK + GCV + DSS"],
            ticks=["PBS", "PBS +\nDSS", "GCV +\nDSS"],
            values=[[[0,0,0,0,0,1],[7,7,7,7,8,8],[1,2,5,5,5,5]],
                    [[0,0,1,2,3,3],[7.2,7.4,8.8,10,10.2,10.4],[3,3.1,3.3,3.5,4.2,5]],
                    [[1.5,2.4,2.8,3,3.5,3.8],[6.3,6.5,11,12,12.5,14],[2,2.8,3,3.5,4,4]]],
            stars=[["****","*","ns"],["****","****","*"],["**","**","ns"]],
            contrasts=["DSS_induction","FAP_ablation"],
            contrast_titles=["DSS induction: DSS versus PBS", "FAP ablation: GCV + DSS versus DSS"]),
    5: dict(title="α5β1 blockade reduces DSS-colitis pathology", image="a5B1i DSS.png",
            boxes=[(15,0,847,704),(874,0,1706,704),(1751,0,2583,704)],
            labels=["WT + PBS", "WT + PBS + DSS", "WT + α5β1i + DSS"],
            ticks=["PBS", "DSS", "DSS +\nα5β1i"],
            values=[[[0,0,0,0,1,1],[7,8,8,8,10,10],[1,1,1,2,2,2,2,2,2,4,5,5]],
                    [[1,1.5,2,2,2.5,3],[9,10,11,14,16,18],[.5,.7,.8,1,1,1.2,1.3,1.5,1.6,1.8,2,2.2]],
                    [[1.5,2,2.5,3,3.5,4],[11,12,12.5,13,14,16.5],[1.4,1.7,2,2.2,2.4,2.6,2.8,3,3.3,3.5,4,4.2]]],
            stars=[["****","**","*"],["****","****","ns"],["*","***","ns"]],
            contrasts=["Severe DSS colitis - Control", "Blockade - severe DSS colitis"],
            contrast_titles=["DSS induction: DSS versus control", "α5β1 blockade: blockade versus DSS"]),
}

SHORT = {
 "Mature absorptive colonocyte":"Mature absorptive", "Absorptive-secretory transitional epithelial":"Abs.-secr. transitional",
 "Secretory-absorptive transitional epithelial":"Secr.-abs. transitional", "IFN/MHC-II-responsive epithelial":"IFN/MHC-II-responsive",
 "Reg4+ deep-crypt secretory epithelial":"Reg4+ deep-crypt secr.", "Colonocyte-like epithelial cell":"Colonocyte-like",
 "Satb2+ epithelial-like transitional cell":"Satb2+ transitional", "Zg16+ secretory epithelial":"Zg16+ secretory",
 "Inflammatory Zg16+ secretory epithelial":"Inflam. Zg16+ secretory", "Enteroendocrine-secretory transitional cell":"Enteroendocrine-secr. trans.",
 "Injury/inflammatory epithelial":"Injury/inflammatory", "IFN-responsive absorptive colonocyte":"IFN-resp. colonocyte",
 "Enteroendocrine cell":"Enteroendocrine", "Mature B cell":"Mature B", "Conventional B cell":"Conventional B",
 "Inflammatory monocyte/macrophage":"Inflam. mono./macro.", "Folr2+ resident macrophage":"Folr2+ res. macro.",
 "Gamma-delta T/NK-like cell":"γδ T/NK-like", "Activated B cell":"Activated B", "Naive conventional T cell":"Naive conv. T",
 "Foxp3+ regulatory T cell":"Foxp3+ Treg", "Cxcr2+ neutrophil":"Cxcr2+ neutrophil", "Conventional T cell":"Conventional T",
 "Conventional dendritic cell":"Conventional DC", "Plasma cell":"Plasma", "Inflammatory neutrophil":"Inflam. neutrophil",
 "Follicular B cell":"Follicular B", "IL-5/IL-13+ ILC2":"IL5/IL13+ ILC2", "Cycling immune cell":"Cycling immune",
 "C1q+ macrophage":"C1q+ macrophage", "Activated dendritic cell":"Activated DC", "IL-17+ gamma-delta T cell":"IL17+ γδ T",
 "IL-22+ ILC3-like cell":"IL22+ ILC3-like", "Gamma-delta T/NK cell":"γδ T/NK", "Naive T cell":"Naive T",
 "Mast cell/basophil":"Mast/basophil", "Ccr3+ eosinophil":"Ccr3+ eosinophil", "Neutrophil MX1":"MX1 neutrophil",
 "Neutrophil PADI4":"PADI4 neutrophil", "Neutrophil OSM":"OSM neutrophil", "Neutrophil CXCR4":"CXCR4 neutrophil",
 "ECM-rich fibroblast":"ECM-rich fibroblast", "Il1r1+ fibroblast":"Il1r1+ fibroblast",
 "Ccl7+ Cxcl1+ inflammatory fibroblast":"Ccl7+ Cxcl1+ inflam.", "Clec3b+ Has1+ matrix fibroblast":"Clec3b+ Has1+ matrix",
 "Ackr4+ Pi16-like fibroblast":"Ackr4+ Pi16-like", "Cxcl13+ fibroblast":"Cxcl13+ fibroblast",
 "Apod+ Wnt5a+ trophic fibroblast":"Apod+ Wnt5a+ trophic", "Fibroblast (subtype unresolved)":"Unresolved fibroblast",
 "Wif1+ BMP-niche fibroblast":"Wif1+ BMP-niche", "Pdgfra+ Wnt5a+ fibroblast":"Pdgfra+ Wnt5a+",
}


def text(fig,x,y,s,size=5.3,weight="normal",color=INK,**kw):
    return fig.text(x/W,1-y/H,s,fontsize=size,fontweight=weight,color=color,va=kw.pop("va","top"),**kw)


def axis(fig,x,y,w,h,label):
    return fig.add_axes([x/W,1-(y+h)/H,w/W,h/H],label=label)


def heading(fig,x,y,letter,s):
    text(fig,x,y-.2,letter,7,"bold")
    text(fig,x+4.5,y,s,6.5,"bold")


def line(fig,x1,y1,x2,y2,color=GRID,lw=.5,**kw):
    fig.add_artist(Line2D([x1/W,x2/W],[1-y1/H,1-y2/H],transform=fig.transFigure,color=color,linewidth=lw,**kw))


def clean(ax,grid="y"):
    ax.spines[["top","right"]].set_visible(False)
    ax.spines[["bottom","left"]].set_linewidth(.5)
    ax.tick_params(labelsize=5,length=1.8,width=.5,pad=1.2,colors=INK)
    ax.set_axisbelow(True)
    ax.grid(axis=grid,color=GRID,linewidth=.35)


def draw_histology(fig,number):
    cfg=CONFIG[number]
    heading(fig,3,3,"a",cfg["title"])
    image=Image.open(DATA/f"figure{number}"/cfg["image"]).convert("RGBA")
    source=Image.new("RGB",image.size,"white");source.paste(image,mask=image.getchannel("A"))
    for i,box in enumerate(cfg["boxes"]):
        x=8+i*27
        crop=source.crop(box)
        h=25.5*crop.height/crop.width
        ax=axis(fig,x,9,25.5,h,f"histology_{i}")
        ax.imshow(crop,interpolation="none");ax.axis("off")
        # Preserve the existing provisional scale fraction, with black text and no backing box.
        ax.plot([.68,.91],[.075,.075],transform=ax.transAxes,color="white",lw=2.5,solid_capstyle="butt")
        ax.plot([.68,.91],[.075,.075],transform=ax.transAxes,color="black",lw=1.3,solid_capstyle="butt")
        ax.text(.795,.17,"100 µm",transform=ax.transAxes,fontsize=5,color="black",ha="center",va="center",
                path_effects=[pe.withStroke(linewidth=.9,foreground="white")])
        text(fig,x+12.75,32,cfg["labels"][i],5.0,"bold",ha="center")
    titles=["Histologic\ninflammation","Colonic\nneutrophils","Total\nfibroblasts"]
    ylabels=["Histology score","Immune cells (%)","Total cells (%)"]
    for i,groups in enumerate(cfg["values"]):
        left=11+i*27
        text(fig,left+9.5,36.4,titles[i],5.5,"bold",ha="center",linespacing=1.05)
        ax=axis(fig,left,42,19.4,26.5,f"quant_{i}")
        rng=np.random.default_rng(20260901)
        ymax=max(max(g) for g in groups)
        group_colors=[INK,RED,PURPLE] if number==4 else [BLUE,ORANGE,TEAL]
        for j,(vals,color) in enumerate(zip(groups,group_colors)):
            vals=np.asarray(vals,float)
            jitter=rng.uniform(-.13,.13,len(vals))
            ax.scatter(j+jitter,vals,s=9,color=color,edgecolors="white",linewidths=.25,zorder=4)
            mean,sem=vals.mean(),vals.std(ddof=1)/np.sqrt(len(vals))
            ax.errorbar(j,mean,yerr=sem,fmt="none",ecolor=color,elinewidth=.7,capsize=1.8,zorder=5)
            ax.plot([j-.19,j+.19],[mean,mean],color=color,lw=1.1,zorder=5)
        for k,(first,second) in enumerate([(0,1),(1,2),(0,2)]):
            y=ymax*(1.12+k*.14);dh=ymax*.025
            ax.plot([first,first,second,second],[y,y+dh,y+dh,y],color=INK,lw=.55)
            ax.text((first+second)/2,y+dh*1.3,cfg["stars"][i][k],ha="center",va="bottom",fontsize=5)
        ax.set_xlim(-.5,2.5);ax.set_ylim(0,ymax*1.61)
        ax.set_xticks(range(3),cfg["ticks"])
        ax.yaxis.set_major_locator(MaxNLocator(nbins=4,steps=[1,2,5,10]))
        ax.set_ylabel(ylabels[i],fontsize=5,labelpad=1.2)
        clean(ax)


def load_umaps():
    return {n:pd.read_csv(DATA/f"figure{n}"/name) for n,name in UMAP_FILES.items()}


def palettes(all_data):
    colors=list(mpl.colormaps["tab20"].colors)+list(mpl.colormaps["Dark2"].colors)
    out={}
    for comp in COMPARTMENTS:
        states=sorted(set().union(*(set(d.loc[d.compartment.eq(comp),"level3"]) for d in all_data.values())))
        assert len(states)<=len(colors)
        out[comp]={s:colors[i] for i,s in enumerate(states)}
    return out


def umap_glyph(ax):
    for end in [(.22,.03),(.03,.24)]:
        ax.annotate("",xy=end,xytext=(.03,.03),xycoords="axes fraction",
                    arrowprops=dict(arrowstyle="-|>",lw=.5,color=MUTED,mutation_scale=4,shrinkA=0,shrinkB=0))
    ax.text(.25,-.03,"UMAP1",transform=ax.transAxes,fontsize=5,color=MUTED,va="bottom")
    ax.text(-.025,.28,"UMAP2",transform=ax.transAxes,fontsize=5,color=MUTED,rotation=90,va="bottom")


def draw_umaps(fig,number,all_data):
    heading(fig,92,3,"b","Compartment-specific UMAPs")
    data=all_data[number];colors=palettes(all_data)
    enlarged=number==5
    map_w=38.0 if enlarged else 26.0
    map_h=19.2*map_w/26.0
    row_step=36.0 if enlarged else 23.5
    for row,comp in enumerate(COMPARTMENTS):
        top=8.5+row*row_step
        z=data.loc[data.compartment.eq(comp)]
        order=z.level3.value_counts().index.tolist()
        text(fig,94,top,f"{comp} | {len(z):,} cells",5.8 if enlarged else 5.5,"bold")
        text(fig,179,top,f"{len(order)} Level 3 states",5.3 if enlarged else 5,ha="right",color=MUTED)
        ax=axis(fig,93,top+3.4,map_w,map_h,f"umap_{comp}")
        x0,x1=z.UMAP_1.min(),z.UMAP_1.max();y0,y1=z.UMAP_2.min(),z.UMAP_2.max()
        cx,cy=(x0+x1)/2,(y0+y1)/2
        span_y=max((y1-y0)*1.15,(x1-x0)*1.15/(26/19.2));span_x=span_y*26/19.2
        cx-=span_x*.02;cy-=span_y*.02
        z=z.sample(frac=1,random_state=45)
        ax.scatter(z.UMAP_1,z.UMAP_2,c=[colors[comp][s] for s in z.level3],s=.6 if enlarged else .4,alpha=.82,linewidths=0,rasterized=True)
        ax.set_xlim(cx-span_x/2,cx+span_x/2);ax.set_ylim(cy-span_y/2,cy+span_y/2)
        ax.set_aspect("equal");ax.axis("off");umap_glyph(ax)
        cols=(2 if comp=="Immune" else 1) if enlarged else (3 if comp=="Immune" else 2)
        legend_left,legend_width,legend_height=(133,46.2,27.6) if enlarged else (121,59.2,19)
        rows=math.ceil(len(order)/cols);step=legend_height/rows;colw=legend_width/cols
        for i,state in enumerate(order):
            col,r=divmod(i,rows);x=legend_left+col*colw;y=top+4.2+r*step
            line(fig,x,y+.7,x,y+.7,color=colors[comp][state],marker="o",markersize=3 if enlarged else 2.5)
            text(fig,x+1.65,y,SHORT[state],5.3 if enlarged and comp!="Immune" else 5.0)


def draw_milo(fig,number):
    cfg=CONFIG[number]
    data=pd.read_csv(DATA/f"figure{number}"/f"Figure{number}C_miloR_summary.csv")
    stacked=number==5
    heading(fig,3,78,"c","MiloR neighborhood fractions" if stacked else "MiloR neighborhood remodeling (FDR < 0.05)")
    legend_positions=[(50,BLUE,"Depleted"),(71,RED,"Enriched")] if stacked else [(143,BLUE,"Depleted"),(164,RED,"Enriched")]
    for x,color,label in legend_positions:
        line(fig,x,79.1,x+3,79.1,color=color,lw=2)
        text(fig,x+4.1,79.1,label,5,va="center")
    for i,contrast in enumerate(cfg["contrasts"]):
        x=24 if stacked else 24+i*94
        title_y=83+i*18 if stacked else 82.7
        plot_y=86.5+i*18 if stacked else 87
        text(fig,x,title_y,cfg["contrast_titles"][i],5.4,"bold")
        ax=axis(fig,x,plot_y,60.5,8.5 if stacked else 11.5,f"milo_{i}")
        p=data.loc[data.contrast.eq(contrast)].set_index("compartment").loc[COMPARTMENTS]
        ax.barh(np.arange(3),-p.depleted/p.neighborhoods,color=BLUE,height=.58,label="Depleted")
        ax.barh(np.arange(3),p.enriched/p.neighborhoods,color=RED,height=.58,label="Enriched")
        ax.axvline(0,color=INK,lw=.6);ax.set_xlim(-.85,.85);ax.set_ylim(2.55,-.55)
        ax.set_yticks(range(3),COMPARTMENTS);ax.set_xticks([-.8,-.4,0,.4,.8])
        clean(ax,grid="x");ax.spines["left"].set_visible(False);ax.tick_params(axis="y",length=0)
        if not stacked:ax.set_xlabel("Fraction of tested neighborhoods",fontsize=5,labelpad=1)
    if stacked:
        text(fig,8,83,"FDR < 0.05",5,color=MUTED)
        text(fig,54.25,117.3,"Fraction of tested neighborhoods",5,ha="center")


LR_ROWS=[
 ("Fibroblast -> Neutrophil","Il1b","Il1r2"),("Fibroblast -> Neutrophil","Il1b","Il1rap"),
 ("Fibroblast -> Neutrophil","Cxcl5","Cxcr2"),("Fibroblast -> Neutrophil","Cxcl2","Cxcr2"),
 ("Fibroblast -> Neutrophil","H2.DMb1","Cd74"),("Fibroblast -> Neutrophil","Vcam1","Itgb2"),
 ("Fibroblast -> Neutrophil","Il33","Il1rap"),("Fibroblast -> Neutrophil","Cxcl1","Cxcr2"),
 ("Fibroblast -> Neutrophil","Tgfbi","Itgb1"),("Neutrophil -> Fibroblast","Thbs1","Itga4"),
 ("Neutrophil -> Fibroblast","Il1a","Il1r1"),("Neutrophil -> Fibroblast","Tgm2","Itga4"),
 ("Neutrophil -> Fibroblast","Tnf","Tnfrsf1b"),("Neutrophil -> Fibroblast","Nampt","Itga5"),
 ("Neutrophil -> Fibroblast","Il10","Il10ra"),("Neutrophil -> Fibroblast","S100a9","Alcam"),
 ("Neutrophil -> Fibroblast","Osm","Osmr"),("Neutrophil -> Fibroblast","Osm","Il6st")]


def draw_signalling(fig):
    heading(fig,3,107,"d","Reciprocal signalling routes")
    for x,color,marker,label in [(8,RED,"o","DSS versus PBS"),(44,PURPLE,"^","DSS versus ablation")]:
        line(fig,x,111.8,x,111.8,color=color,marker=marker,markersize=3)
        text(fig,x+2,111.8,label,5.0,va="center")
    data=pd.read_csv(DATA/"figure4"/"both_comparisons_nonambient_interactions.csv")
    for block,direction in enumerate(["Fibroblast -> Neutrophil","Neutrophil -> Fibroblast"]):
        y=116.5+block*21.5
        text(fig,8,y-3,direction.replace(" -> "," → "),5.4,"bold")
        ax=axis(fig,32,y,55.2,17.6,f"signalling_{block}")
        selected=[r for r in LR_ROWS if r[0]==direction]
        for row,(_,ligand,receptor) in enumerate(selected):
            q=data.loc[data.direction.eq(direction)&data.ligand.eq(ligand)&data.receptor.eq(receptor)].set_index("comparison")
            assert len(q)==2
            values=[q.loc[c,"score_delta"] for c in ["PBS_DSS_vs_PBS","PBS_DSS_vs_GCV_DSS"]]
            ax.plot(values,[row,row],color="#B8BEC4",lw=.8,zorder=1)
            ax.scatter(values[0],row,s=10,color=RED,marker="o",zorder=2)
            ax.scatter(values[1],row,s=14,color=PURPLE,marker="^",zorder=3)
        ax.set_ylim(8.6,-.6);ax.set_xlim(0,.88)
        ax.set_yticks(range(9),[f"{l} → {r}" for _,l,r in selected]);ax.set_xticks([0,.2,.4,.6,.8])
        clean(ax,grid="x");ax.tick_params(axis="y",length=0,pad=1.5);ax.spines["left"].set_visible(False)
        for t in ax.get_yticklabels():t.set_fontstyle("italic")
        if block==0:ax.tick_params(axis="x",labelbottom=False,length=0);ax.spines["bottom"].set_visible(False)
        else:ax.set_xlabel("DSS-associated priority-score increase",fontsize=5,labelpad=1)


PATH_NAMES={"Tnfa Signaling via Nfkb":"TNF/NF-κB", "Apoptosis":"Apoptosis", "Hypoxia":"Hypoxia",
            "Inflammatory Response":"Inflammation", "Interferon Gamma Response":"IFN-γ response",
            "Epithelial Mesenchymal Transition":"EMT", "Il6 Jak Stat3 Signaling":"IL6/JAK/STAT3"}


def draw_pathways(fig):
    heading(fig,96,107,"e","Predicted receiver-cell pathways")
    data=pd.read_csv(DATA/"figure4"/"Figure_4F_downstream_target_pathways.csv")
    cmap=LinearSegmentedColormap.from_list("pathways",["#F4E8D4","#A96FAC","#772B8C"])
    norm=Normalize(0,6)
    text(fig,100.5,111,"Target genes",5,color=MUTED)
    for x,value in [(117,10),(126,30),(135,60)]:
        line(fig,x,112.4,x,112.4,color="#777777",marker="o",markersize=np.sqrt(value*.95),markerfacecolor="white",markeredgewidth=.5)
        text(fig,x+1.8,112.4,str(value),5,va="center")
    cax=axis(fig,150,111.5,27,1.4,"pathway_colorbar")
    cb=fig.colorbar(mpl.cm.ScalarMappable(norm=norm,cmap=cmap),cax=cax,orientation="horizontal",ticks=[0,3,6])
    cb.outline.set_visible(False);cb.solids.set_rasterized(False);cb.solids.set_edgecolor("face")
    cb.ax.tick_params(labelsize=5,length=1,pad=.5,width=.4)
    text(fig,163.5,117.0,"−log10 FDR (capped at 6)",5,ha="center",color=MUTED)
    configs=[
      ("Fibroblast -> Neutrophil",121,10.2,["Tnfa Signaling via Nfkb","Apoptosis","Hypoxia"],
       ["IL1/IL33","CXCL-CXCR2","TGF-beta/ECM","VCAM adhesion","MHC-II/CD74"],
       ["IL1/IL33","CXCL/\nCXCR2","TGFβ/\nECM","VCAM\nadhesion","MHC-II/\nCD74"]),
      ("Neutrophil -> Fibroblast",143.4,13.0,["Tnfa Signaling via Nfkb","Inflammatory Response","Interferon Gamma Response","Epithelial Mesenchymal Transition","Il6 Jak Stat3 Signaling"],
       ["Integrin/ECM","S100A9-ALCAM","OSM","TNF","IL1","IL10"],
       ["Integrin/\nECM","S100A9/\nALCAM","OSM","TNF","IL1","IL10"])]
    for block,(direction,y,h,pathways,families,labels) in enumerate(configs):
        text(fig,100.5,y-2.6,direction.replace(" -> "," → "),5.2,"bold")
        ax=axis(fig,131,y,47.5,h,f"pathways_{block}")
        q=data.loc[data.direction.eq(direction)].copy()
        x=q.lr_pathway.map({v:i for i,v in enumerate(families)});yy=q.downstream_pathway.map({v:i for i,v in enumerate(pathways)})
        assert x.notna().all() and yy.notna().all()
        ax.scatter(x,yy,s=q.overlap*.95,c=q.enrichment_strength,cmap=cmap,norm=norm,edgecolors="#695B6C",linewidths=.35,zorder=3)
        ax.set_xlim(-.5,len(families)-.5);ax.set_ylim(len(pathways)-.5,-.5)
        ax.set_xticks(range(len(families)),labels);ax.set_yticks(range(len(pathways)),[PATH_NAMES[p] for p in pathways])
        clean(ax,grid="both");ax.tick_params(length=0,pad=1.4)
        ax.spines[["left","bottom"]].set_visible(False)


def draw_programs(fig):
    heading(fig,3,121,"d","Neutrophil programs along pseudotime")
    text(fig,7.5,125.4,"Source mouse-by-bin means ± s.e.m.; observed bins only",5,color=MUTED)
    conditions=[("Control",BLUE,"Control"),("Severe DSS colitis",ORANGE,"DSS colitis"),("DSS + alpha5beta1 blockade",TEAL,"DSS + α5β1i")]
    for x,(_,color,label) in zip([105,128,155],conditions):
        line(fig,x,122.2,x+3.5,122.2,color=color,lw=1)
        text(fig,x+4.6,122.2,label,5.1,va="center")
    data=pd.read_csv(DATA/"figure5"/"neutrophil_program_scores_across_pseudotime.csv")
    for i,gene in enumerate(["OSM","CXCR4","PADI4","MX1"]):
        x=12+i*43
        text(fig,x+18,131.0,gene,6.2,"bold",ha="center")
        ax=axis(fig,x,136,35.8,21,f"trajectory_{gene}")
        lower,upper=np.inf,-np.inf
        for condition,color,label in conditions:
            p=data.loc[data.score_program.eq("Neutrophil "+gene)&data.condition.eq(condition)].sort_values("bin_midpoint")
            xx,yy,se=p.bin_midpoint.to_numpy(),p["mean"].to_numpy(),p["sem"].to_numpy()
            ax.plot(xx,yy,color=color,lw=1,marker="o",markersize=1.4,markeredgewidth=0,label=label)
            ax.fill_between(xx,yy-se,yy+se,color=color,alpha=.13,lw=0)
            lower=min(lower,np.nanmin(yy-se));upper=max(upper,np.nanmax(yy+se))
        margin=(upper-lower)*.05
        ax.set_ylim(lower-margin,upper+margin);ax.set_xlim(0,1);ax.set_xticks([0,.25,.5,.75,1],["0",".25",".5",".75","1"])
        ax.axhline(0,color="#A9AFB6",lw=.5);clean(ax)
        ax.yaxis.set_major_locator(MaxNLocator(nbins=5))
        ax.set_xlabel("Pseudotime",fontsize=5.3,labelpad=1.3)
        if i==0:ax.set_ylabel("Module score (z)",fontsize=5.3,labelpad=1.5)


def build_figure(number,save):
    mpl.rcParams.update({"font.family":"Arial","font.size":5.3,"pdf.fonttype":42,"ps.fonttype":42,"savefig.facecolor":"white"})
    fig=plt.figure(figsize=(W/25.4,H/25.4),facecolor="white")
    draw_histology(fig,number);draw_umaps(fig,number,load_umaps());draw_milo(fig,number)
    if number==4:draw_signalling(fig);draw_pathways(fig)
    else:draw_programs(fig)
    save(fig,number)
