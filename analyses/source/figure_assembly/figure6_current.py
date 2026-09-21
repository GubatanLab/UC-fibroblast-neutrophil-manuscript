"""Rebuild current Figure 6 from deposited tables; canonical PDFs are untouched.

Requires NumPy, pandas, matplotlib, PyMuPDF, Pillow, and installed Arial fonts.
Output is written under results/figure6_current. This is a portable path
adaptation of the approved six-panel assembly; estimates and layout are unchanged.
"""
from pathlib import Path
import sys, json, hashlib, subprocess
ROOT = Path(__file__).resolve().parents[3]
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle
import pymupdf as fitz
from PIL import Image

OUT = ROOT/'results/figure6_current'
QA = OUT/'qa'
OUT.mkdir(parents=True, exist_ok=True); QA.mkdir(parents=True, exist_ok=True)
S = ROOT/'source_data/figure_06'
P = S/'epithelial_programs_2026-09-07'
DATA = S/'epithelial_panels_2026-09-21'
sha = lambda p: hashlib.sha256(p.read_bytes()).hexdigest()
protected = list(S.rglob('*')) + list((ROOT/'figures/main').glob('*'))
before = {str(p.relative_to(ROOT)): sha(p) for p in protected if p.is_file()}
ef = pd.read_csv(P/'epithelial_program_effects.csv')
programs = list(json.loads((P/'prespecified_programs.json').read_text()))
all_ep = ef[ef.scope.eq('All epithelial')]
within = ef[ef.scope.eq('Level 3')]
composition = pd.read_csv(S/'Figure_6_Level3_display_source.csv')
compartment = pd.read_csv(S/'Figure_6_compartment_recovery_source.csv')
z = pd.read_csv(DATA/'All_39_genes_row_zscores.csv', index_col=0)
expression = pd.read_csv(DATA/'All_39_genes_stored_log2_CPM.csv', index_col=0)
inventory = pd.read_csv(DATA/'Panel_E_mouse_inventory.csv')
profiles = pd.read_csv(DATA/'Panel_F_existing_mouse_program_scores.csv')
centroids = pd.read_csv(DATA/'Panel_F_display_group_means.csv', index_col=0)
coverage = pd.read_csv(DATA/'Gene_detection_by_group.csv').set_index('gene')
keys = inventory.key.tolist()
genes = ['Cxcl1','Cxcl2','Cxcl5','Ccl20','Muc2','Tff3','Fcgbp','Tjp1','Ocln','Cldn7','Slc26a3','Aqp8','Alpi','Clu','Ly6a','Tacstd2']
assert len(all_ep)==12 and len(within)==30 and len(composition)==18 and len(compartment)==6
assert len(inventory)==len(profiles)==22
assert inventory.groupby('group').size().to_dict()=={'DSS':6,'FAP ablation':4,'alpha5beta1 blockade':12}
assert profiles.groupby('group').size().to_dict()==inventory.groupby('group').size().to_dict()
assert int(coverage.loc['Tacstd2','mice_detected'])==4
W,H = 180,170
INK='#25282C'; MUTED='#525960'; PUR='#7651A6'; TEAL='#168C88'; GRAY='#7F8A96'
groups=[('FAP ablation',PUR,'o'),('alpha5beta1 blockade',TEAL,'^')]
allgroups=[('DSS',GRAY,'s')]+groups
plt.rcParams.update({'font.family':'Arial','font.size':5.4,'axes.labelsize':5.5,'xtick.labelsize':5.2,'ytick.labelsize':5.4,'axes.linewidth':.5,'pdf.fonttype':42,'ps.fonttype':42,'svg.fonttype':'none','savefig.facecolor':'white','axes.unicode_minus':False})
fig=plt.figure(figsize=(W/25.4,H/25.4),facecolor='white')
def ax(x,y,w,h): return fig.add_axes([x/W,1-(y+h)/H,w/W,h/H])
def txt(x,y,s,size=5.5,bold=False,**kw): return fig.text(x/W,1-y/H,s,fontsize=size,weight='bold' if bold else 'normal',va='top',color=INK,**kw)
def head(x,y,letter,s): txt(x,y,letter,8,True);txt(x+4,y,s,7,True)
def clean(a):
    a.spines[['top','right']].set_visible(False)
    a.tick_params(length=1.7,width=.5,pad=1.2)
def cells(a,mat,limit):
    cm=plt.get_cmap('RdBu_r').copy();cm.set_bad('#E8EAED')
    n,m=mat.shape
    im=a.pcolormesh(np.arange(m+1)-.5,np.arange(n+1)-.5,np.ma.masked_invalid(mat),cmap=cm,vmin=-limit,vmax=limit,shading='flat',rasterized=False)
    a.set_xlim(-.5,m-.5);a.set_ylim(n-.5,-.5);a.spines[:].set_visible(False)
    a.tick_params(length=0,pad=1.3)
    return im
def cb_style(cb):
    cb.outline.set_visible(False);cb.ax.tick_params(length=1,pad=.8,labelsize=5)
    cb.solids.set_rasterized(False)
    cb.solids.set_edgecolor('face')

# a: all 12 existing effects, confidence intervals, and q values.
head(2,2,'a','Epithelial RNA programs')
txt(6,6.9,'Whole epithelium | DSS n=6; ablation n=4; blockade n=12',5.1)
legend=[Line2D([0],[0],color=c,marker=m,lw=.8,markersize=3.1,label='FAP ablation' if i==0 else 'α5β1 blockade') for i,(_,c,m) in enumerate(groups)]
fig.legend(handles=legend,loc='upper left',bbox_to_anchor=(6/W,1-10.3/H),borderaxespad=0,frameon=False,ncol=2,fontsize=5.3,handlelength=1.2,columnspacing=1)
a=ax(37,17,39,29.5)
for j,(g,col,marker) in enumerate(groups):
    data=all_ep[all_ep.group.eq(g)].set_index('program').loc[programs]
    yy=np.arange(6)+(j-.5)*.33
    a.errorbar(data.effect,yy,xerr=[data.effect-data.low,data.high-data.effect],fmt=marker,ms=3.1,lw=.8,c=col,capsize=1.3,capthick=.6)
    for yy1,r in zip(yy,data.itertuples()):
        a.text(1.13,yy1,'<.001' if r.q<.001 else f'{r.q:.3f}',transform=a.get_yaxis_transform(),va='center',fontsize=5.1,color=INK)
        a.plot(1.06,yy1,marker=marker,color=col,ms=2.1,transform=a.get_yaxis_transform(),clip_on=False)
a.set_yticks(range(6),['Inflammatory\nchemokines','Interferon response','Barrier / junctions','Mucus / secretory','Injury-associated\nregeneration','Absorptive\ndifferentiation'])
a.set_ylim(5.55,-.55);a.set_xlim(-7,4.2);a.set_xticks([-6,-3,0,3])
a.axvline(0,c='#8E949A',lw=.55);a.set_xlabel('Intervention - DSS (log2 CPM)',labelpad=1.4)
a.text(1.13,1.06,'q',transform=a.transAxes,fontsize=5.2)
clean(a);a.tick_params(axis='y',length=0,pad=3,labelsize=5.4)
txt(6,52.2,'Mean differences; pointwise 95% mouse-bootstrap CIs',5)

# b: two compact heatmaps preserve all states and all 30 supported entries.
head(94,2,'b','Within-state epithelial programs')
txt(98,6.9,'≥20 cells per mouse; n = DSS / intervention',5.1)
states=list(within.state.unique());limit=float(abs(within.effect).max())
state_labels=['Satb2+\ntransitional','Mature absorptive\ncolonocyte','IFN/MHC-II-\nresponsive']
short=['Chemo.','IFN','Barrier','Mucus','Injury','Absorp.']
for j,(g,col,_) in enumerate(groups):
    y0=15+j*18.7
    txt(129,y0-3.5,'FAP ablation' if j==0 else 'α5β1 blockade',5.7,True)
    fig.lines.append(Line2D([126/W,128/W],[1-(y0-2.5)/H]*2,transform=fig.transFigure,lw=2,color=col))
    a=ax(127,y0,50,11.8)
    m=within[within.group.eq(g)].pivot(index='state',columns='program',values='effect').reindex(index=states,columns=programs)
    im=cells(a,m.to_numpy(),limit)
    ns=['6/3','6/4','6/4'] if j==0 else ['NA','6/12','6/5']
    a.set_yticks(range(3),[f'{label} ({n})' for label,n in zip(state_labels,ns)]);a.set_xticks(range(6),short);a.tick_params(labelsize=5,pad=1.2)
    for yy in [.5,1.5]:a.axhline(yy,c='white',lw=.4)
    for xx in np.arange(.5,5,1):a.axvline(xx,c='white',lw=.4)
    for i,state in enumerate(states):
        for k,prog in enumerate(programs):
            rr=within[(within.group==g)&(within.state==state)&(within.program==prog)]
            if len(rr) and rr.iloc[0].q<.05:
                a.text(k,i,'*',ha='center',va='center',fontsize=7,color='white' if abs(rr.iloc[0].effect)>.55*limit else INK)
cb=fig.colorbar(im,cax=ax(127,50,21,1.2),orientation='horizontal',ticks=[-4,0,4]);cb_style(cb)
txt(150,49.6,'Δ log2 CPM',5)
txt(98,51.8,'*q<0.05',5)
fig.patches.append(Rectangle((151/W,1-54.4/H),1.7/W,1.7/H,transform=fig.transFigure,facecolor='#E8EAED',edgecolor='none'))
txt(154,52.5,'Insufficient support',5)

# c: all 18 previously displayed states and original statistical calls.
head(2,57,'c','Cell-state composition')
mat=composition[['effect_percentage_points_ablation','effect_percentage_points_blockade']].to_numpy()
comp_limit=float(abs(mat).max());a=ax(45,66,47,39.6);im=cells(a,mat,comp_limit)
labels=['Absorptive-secretory transitional','Mature absorptive colonocyte','Zg16+ secretory epithelial','Injury/inflammatory epithelial','Satb2+ transitional','IFN/MHC-II-responsive','Gamma-delta T/NK-like','Ccr3+ eosinophil','Cxcr2+ neutrophil','Inflammatory mono./macro.','Conventional B cell','Mature B cell','Ackr4+ Pi16-like fibroblast','Cxcl13+ fibroblast','Ccl7+ Cxcl1+ inflam. fibroblast','Clec3b+ Has1+ matrix fibroblast','ECM-rich fibroblast','Il1r1+ fibroblast']
a.set_yticks(range(18),labels);a.set_xticks([0,1],['FAP ablation','α5β1 blockade'])
a.tick_params(length=0,pad=1.1,top=True,labeltop=True,bottom=False,labelbottom=False,labelsize=5.1)
for i,row in composition.iterrows():
    for j,suffix in enumerate(['ablation','blockade']):
        q=row['fdr_global_scope_'+suffix];mark='**' if q<.01 else '*' if q<.05 else ''
        if mark:a.text(j,i,mark,ha='center',va='center',fontsize=5.6,color='white' if abs(mat[i,j])>.35*comp_limit else INK)
for yy in [5.5,11.5]:a.axhline(yy,c='white',lw=1.1)
cb=fig.colorbar(im,cax=ax(45,108,45,1.2),orientation='horizontal',ticks=[-60,-30,0,30,60]);cb_style(cb)
txt(5,108,'Within-compartment change\n(percentage points)',5,linespacing=1.0)

# d: same six effects and intervals for recovered proportions.
head(98,57,'d','Recovered-cell proportions')
txt(102,61.9,'DSS n=6; ablation n=5; blockade n=12',5.1)
a=ax(121,70,56,24)
comps=['Epithelial','Immune','Stromal']
for j,(g,col,m) in enumerate(groups):
    v=compartment[compartment.intervention==g].set_index('compartment').loc[comps]
    yy=np.arange(3)+(j-.5)*.24
    a.errorbar(v.effect_percentage_points,yy,xerr=[v.effect_percentage_points-v.ci_low,v.ci_high-v.effect_percentage_points],fmt=m,ms=3.5,lw=.8,c=col,capsize=1.4)
a.set_yticks(range(3),comps);a.set_ylim(2.45,-.45);a.set_xlim(-75,90);a.set_xticks([-60,-30,0,30,60,90])
a.axvline(0,c='#8E949A',lw=.55);a.set_xlabel('Intervention - DSS (percentage points)',labelpad=2)
clean(a);a.tick_params(axis='y',length=0,pad=3)
txt(102,102,'C: *q<0.05; **q<0.01',5)
txt(102,106,'Shared DSS reference; separate interventions.',5)
txt(102,110,'RNA programs do not measure epithelial repair.',5)

# e: unchanged 16-gene z-scores, all 22 individual mice, vector cells.
head(2,116,'e','Epithelial gene expression')
txt(6,120.6,'Mature absorptive colonocytes | 16 selected genes',5.1)
a=ax(22,130,70,29.7);im=cells(a,z.loc[genes,keys].to_numpy(),2.5)
a.set_yticks(range(16),[g+('†' if g=='Tacstd2' else '') for g in genes],fontsize=5.2,fontstyle='italic')
a.set_xticks(range(22),inventory.mouse_id.tolist(),rotation=90,fontsize=5)
for edge in [3.5,6.5,9.5,12.5]:a.axhline(edge,c='white',lw=.7)
for edge in [5.5,9.5]:a.axvline(edge,c='white',lw=1.2)
for start,n,(g,col,m),label in zip([0,6,10],[6,4,12],allgroups,['DSS (6)','Ablation (4)','α5β1 blockade (12)']):
    left=22+70*start/22;gw=70*n/22
    bar=ax(left,128.3,gw,.9);bar.set_facecolor(col);bar.set_xticks([]);bar.set_yticks([]);bar.spines[:].set_visible(False)
    txt(left+gw/2,124.8,label,5,ha='center')
cb=fig.colorbar(im,cax=ax(22,166.1,20,1),orientation='horizontal',ticks=[-2.5,0,2.5]);cb_style(cb)
txt(46,165.1,'Row z-score (log2 CPM)\n†Tacstd2 detected: 4/22 mice',5,linespacing=1.05)

# f: two facets use exact previously plotted mouse and group-mean coordinates.
head(98,116,'f','Individual-mouse epithelial profiles')
txt(102,120.6,'Whole epithelium | n=6 / 4 / 12 mice',5.1)
handles=[Line2D([0],[0],marker=m,linestyle='none',color=c,markersize=2.8,label={'DSS':'DSS','FAP ablation':'Ablation','alpha5beta1 blockade':'α5β1 blockade'}[g]) for g,c,m in allgroups]
fig.legend(handles=handles,loc='upper left',bbox_to_anchor=(102/W,1-124.3/H),frameon=False,ncol=3,fontsize=5.1,handletextpad=.25,columnspacing=.65,borderaxespad=0)
for j,(program,title) in enumerate([('Absorptive differentiation','Absorptive differentiation'),('Mucus / secretory','Mucus / secretory')]):
    a=ax(110+j*39,135,28,25)
    for g,c,m in allgroups:
        r=profiles[profiles.group==g]
        a.scatter(r['Inflammatory chemokines'],r[program],s=12,color=c,marker=m,edgecolor='white',linewidth=.3,alpha=.9,zorder=3)
        a.scatter(centroids.loc[g,'Inflammatory chemokines'],centroids.loc[g,program],s=44,c=c,marker='+',linewidths=1.15,zorder=5)
    a.set_title(title,fontsize=5.2,pad=3.5)
    a.set_xlim(-2.25,5.8);a.set_xticks([-2,0,2,4])
    if j==0:a.set_ylim(4.1,7.6);a.set_yticks([4.5,5.5,6.5,7.5])
    else:a.set_ylim(4.4,10.5);a.set_yticks([5,7,9])
    a.set_ylabel('Program score (log2 CPM)',fontsize=5.1,labelpad=1.6)
    clean(a);a.tick_params(labelsize=5.1)
txt(144,163.4,'Inflammatory chemokine score (log2 CPM)',5.1,ha='center')
txt(102,167.8,'Small symbols, mice; large crosses, group means.',5)

name='Figure_6_Nature_180x170_Rebuilt'
PDF=OUT/f'{name}.pdf'
fig.savefig(PDF,dpi=450,metadata={'Title':'Figure 6 - epithelial responses, Nature dimensions (rebuild)','Subject':'180 x 170 mm; deposited data; six panels; rebuilt copy'})
fig.savefig(OUT/f'{name}.svg')
fig.savefig(OUT/f'{name}_450dpi.png',dpi=450)
fig.savefig(OUT/f'{name}_preview.png',dpi=200)
plt.close(fig)
with Image.open(OUT/f'{name}_450dpi.png') as im:
    im.convert('RGB').save(OUT/f'{name}_450dpi_LZW.tiff',compression='tiff_lzw',dpi=(450,450))

# Independent checks of text size, page bounds, vector output, and file protection.
with fitz.open(PDF) as d:
    p=d[0];spans=[s for b in p.get_text('dict')['blocks'] if 'lines' in b for l in b['lines'] for s in l['spans']]
    assert abs(p.rect.width*25.4/72-W)<.01 and abs(p.rect.height*25.4/72-H)<.01
    assert min(s['size'] for s in spans)>=4.99
    assert all(s['size']<=7.01 or s['text'] in 'abcdef' for s in spans)
    out_of_bounds=[s['text'] for s in spans if not p.rect.contains(fitz.Rect(s['bbox']))]
    assert not out_of_bounds, out_of_bounds
    assert len(p.get_images())==0
    font_sizes=sorted(set(round(s['size'],2) for s in spans))
    font_names=sorted(set(s['font'] for s in spans))
    assert all('Arial' in f for f in font_names)
assert all(sha(ROOT/p)==h for p,h in before.items())
audit={'status':'preview only; canonical unchanged','dimensions_mm':[W,H],'font_sizes_pt':font_sizes,'font_names':font_names,'raster_images_in_pdf':0,'panel_labels':'8 pt bold lowercase','body_text':'5-7 pt Arial, embedded TrueType; SVG text editable','new_tests_or_scores':False,'panel_A_effects_CIs_q_values':12,'panel_B_supported_effects_q_values':30,'panel_B_unsupported_cells':6,'panel_C_state_effects':36,'panel_D_compartment_effects_CIs':6,'panel_E_genes':genes,'panel_E_mice':22,'panel_F_mice_per_facet':22,'protected_files_unchanged':before,'visual_review':'pending','format_sources':['https://www.nature.com/nature/for-authors/initial-submission','https://research-figure-guide.nature.com/figures/preparing-figures-our-specifications/']}
(OUT/'verification.json').write_text(json.dumps(audit,indent=2),encoding='utf-8')
print(json.dumps({'pdf':str(PDF),'dimensions_mm':[W,H],'font_sizes':font_sizes,'fonts':font_names,'canonical_unchanged':True,'visual_review':'pending'},indent=2))
