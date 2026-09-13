from pathlib import Path
import sys,json,hashlib,shutil
sys.stdout.reconfigure(encoding='utf-8')
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'tmp/canonical_nature_20260903/packages'))
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.ticker import MaxNLocator
import pymupdf as fitz
from PIL import Image

OUT=ROOT/'output/pdf/Extended_Data_1_Nature_2026-09-12';OUT.mkdir(parents=True,exist_ok=True)
SRC=ROOT/'Figure 1/source_data/program_panels_2026-09-12'
OLD=ROOT/'output/Canonical_Nature_Revision_2026-09-03/source_data/figure1'
INK='#2F3337'; MUTED='#666D75'; GOLD='#C9A227'; RED='#D44B50'
CMAP=LinearSegmentedColormap.from_list('program',['#3B75B9','#F7F7F7','#B53A45'])
plt.rcParams.update({'font.family':'Arial','font.size':6,'axes.labelsize':6,'axes.titlesize':6.2,
 'xtick.labelsize':5.5,'ytick.labelsize':5.5,'pdf.fonttype':42,'ps.fonttype':42,'svg.fonttype':'none',
 'text.color':INK,'axes.labelcolor':INK,'axes.edgecolor':INK,'xtick.color':INK,'ytick.color':INK,
 'axes.spines.top':False,'axes.spines.right':False,'axes.linewidth':.5,
 'xtick.major.width':.5,'ytick.major.width':.5,'xtick.major.size':2,'ytick.major.size':2})
fig=plt.figure(figsize=(180/25.4,170/25.4))
def txt(x,y,s,**kw):fig.text(x/180,1-y/170,s,**kw)
def axes(x,y,w,h):return fig.add_axes([x/180,1-(y+h)/170,w/180,h/170])
def heading(letter,title,y,x=2):
 txt(x,y,letter,fontsize=8,fontweight='bold',va='top')
 txt(x+5,y+.2,title,fontsize=6.5,fontweight='bold',va='top')
heading('a','Paired stromal expression of fibroblast-associated programs',2)
txt(7,8,'Matched noninflamed and inflamed UC tissue; 14 patients',fontsize=5.5,va='top',color=MUTED)
fib=pd.read_csv(SRC/'panel_D_paired_fibroblast_scores.csv')
tests=pd.read_csv(SRC/'panel_D_all_8_program_statistics.csv').set_index('program')
progs=['FAP_inflammatory','Neutrophil_recruitment','OSM_response','ECM_remodeling',
       'Alpha5Beta1_adhesion','TGFb_response','NFkB_AP1','YAP_mechanotransduction']
labels=['FAP-inflammatory','Neutrophil recruitment','OSM response','Matrix remodeling',
        'α5β1 adhesion','TGFβ response','NFκB / AP1','YAP mechanotransduction']
ids=sorted(fib.PatientID.unique());assert len(ids)==14
counts=fib.pivot(index='PatientID',columns='condition',values='n_cells').reindex(ids)
for i,(pr,label) in enumerate(zip(progs,labels)):
 ax=axes(13+(i%4)*43,20+(i//4)*36,32,22)
 vals=fib.pivot(index='PatientID',columns='condition',values=pr).reindex(ids)
 jit=np.linspace(-.04,.04,14)
 for j,(_,r) in enumerate(vals.iterrows()):ax.plot(np.array([0,1])+jit[j],[r['Uninflamed UC'],r['Inflamed UC']],color='#BCC1C7',lw=.5,zorder=1)
 for cond,pos,col in [('Uninflamed UC',0,GOLD),('Inflamed UC',1,RED)]:
  ok=counts[cond].to_numpy()>=10
  ax.scatter(pos+jit[ok],vals.loc[ok,cond],s=8,color=col,edgecolor='white',lw=.25,zorder=3)
  ax.scatter(pos+jit[~ok],vals.loc[~ok,cond],s=8,facecolor='white',edgecolor=col,lw=.65,zorder=3)
  ax.plot([pos-.14,pos+.14],[vals[cond].median()]*2,color=INK,lw=.9,zorder=4)
 ax.set_xlim(-.25,1.25);ax.set_ylim(-vals.max().max()*.025,vals.max().max()*1.10)
 ax.set_xticks([0,1],['Noninflamed','Inflamed']);ax.tick_params(axis='x',length=0,pad=3)
 ax.set_title(label,fontweight='bold',pad=13)
 ax.text(.5,1.04,'q = '+format(tests.loc[pr,'FDR'],'.3g'),transform=ax.transAxes,ha='center',fontsize=5.5)
 ax.yaxis.set_major_locator(MaxNLocator(nbins=3));ax.grid(axis='y',lw=.35,color='#E5E7EB');ax.set_axisbelow(True)
 if i%4==0:ax.set_ylabel('Mean expression',labelpad=3)
txt(7,86,'Paired Wilcoxon q (BH, 8 programs). Lines: patients; bars: medians; open points: <10 recovered cells.',fontsize=5.3,color=MUTED)

heading('b','Pooled RNA programs and annotation controls',94)
txt(7,100,'All conditions; original gene modules',fontsize=5.3,color=MUTED,va='top')
table=pd.read_csv(SRC/'22_neutrophil_state_program_summary.tsv',sep='\t').set_index('neut_state')
states=['Erythrocyte','Monocyte-derived macrophage','Neutrophil CXCR4','Neutrophil MX1','Neutrophil OSM','Neutrophil PADI4','Neutrophil lowRNA']
assert set(states)==set(table.index)
nprogs=['Recruitment_migration','Retention_aging','Degranulation','NET_associated','Oxidative_burst','Inflammatory','Interferon','Survival_immaturity','Tissue_injury']
v=table.loc[states,nprogs];z=((v-v.mean())/v.std(ddof=0)).T
ax=axes(31,108,52,42)
im=ax.pcolormesh(np.arange(8)-.5,np.arange(10)-.5,z.to_numpy(),cmap=CMAP,vmin=-2,vmax=2,shading='flat')
ax.set_ylim(8.5,-.5);ax.set_xlim(-.5,6.5)
ax.set_yticks(range(9),['Recruitment / migration','Retention / aging','Granule / protease','NET-associated','Oxidative burst','Inflammatory','Interferon','Survival / immaturity','Tissue injury'])
ax.set_xticks(range(7),['RBC','Mac.','CXCR4','MX1','OSM','PADI4','Low RNA'],rotation=45,ha='right')
ax.tick_params(length=0,pad=2,labelsize=5.2)
ax.set_xticks(np.arange(-.5,7),minor=True);ax.set_yticks(np.arange(-.5,9),minor=True)
ax.grid(which='minor',color='white',lw=.6);ax.tick_params(which='minor',length=0)
for sp in ax.spines.values():sp.set_visible(False)
ca=axes(34,161,32,2)
cb=fig.colorbar(im,cax=ca,orientation='horizontal',ticks=[-2,0,2]);cb.solids.set_rasterized(False);cb.solids.set_edgecolor('face');cb.ax.tick_params(size=1.5,pad=1,labelsize=5)
txt(69,163,'Program z-score',fontsize=5)
z.to_csv(OUT/'ED1b_display_z_scores.csv')

heading('c','Paired ligand–receptor coexpression',94,x=93)
txt(98,100,'9 matched UC patient pairs',fontsize=5.3,color=MUTED,va='top')
co=pd.read_csv(OLD/'paired_lr_interaction_scores.csv')
ct=pd.read_csv(OLD/'paired_lr_interaction_tests.csv')
ct.to_csv(OUT/'coexpression_all_five_original_tests.csv',index=False)
routes=['CXCL1 -> CXCR2','NAMPT -> ITGA5+ITGB1','OSM -> OSMR']
for i,(route,label) in enumerate(zip(routes,['CXCL1–CXCR2','NAMPT–α5β1*','OSM–OSMR'])):
 ax=axes(100+i*27,117,20,29)
 vals=co[co.interaction==route].pivot(index='patient',columns='group',values='score');assert len(vals)==9
 for _,r in vals.iterrows():ax.plot([0,1],[r['Uninflamed.UC'],r['Inflamed.UC']],color='#BCC1C7',lw=.55,zorder=1)
 for cond,pos,col in [('Uninflamed.UC',0,GOLD),('Inflamed.UC',1,RED)]:
  ax.scatter([pos]*len(vals),vals[cond],s=9,color=col,edgecolor='white',lw=.3,zorder=3)
  ax.plot([pos-.12,pos+.12],[vals[cond].median()]*2,color=INK,lw=.9,zorder=4)
 ax.set_xlim(-.25,1.25);ax.set_ylim(0,vals.max().max()*1.13)
 ax.set_xticks([0,1],['Noninf.','Inf.']);ax.tick_params(axis='x',length=0,pad=3)
 ax.yaxis.set_major_locator(MaxNLocator(nbins=3));ax.grid(axis='y',color='#E5E7EB',lw=.35);ax.set_axisbelow(True)
 ax.set_title(label,fontweight='bold',fontsize=5.7,pad=23)
 test=ct.set_index('interaction').loc[route]
 delta=test.median_delta;q=test.q_value
 ax.text(.5,1.08,f'Δ = +{delta:.2f}\nq = {q:.3f}',ha='center',transform=ax.transAxes,fontsize=5.2,linespacing=1.4)
 if i==0:ax.set_ylabel('Coexpression score',labelpad=2)
txt(98,157,'Paired Wilcoxon q (BH, 5 original routes).',fontsize=5.1,color=MUTED)
txt(98,162,'*Receptor summary: min(ITGA5, ITGB1).',fontsize=5.1,color=MUTED)
txt(98,167,'Δ: median paired change; lines: patients.',fontsize=5.1,color=MUTED)
pdf=OUT/'Extended_Data_Figure_01.pdf'
fig.savefig(pdf);fig.savefig(OUT/'Extended_Data_Figure_01.svg');plt.close(fig)
doc=fitz.open(pdf);p=doc[0]
pix=p.get_pixmap(dpi=300,alpha=False);pix.save(OUT/'Extended_Data_Figure_01.png')
Image.frombytes('RGB',(pix.width,pix.height),pix.samples).save(OUT/'Extended_Data_Figure_01.tiff',compression='tiff_lzw',dpi=(300,300))
p.get_pixmap(dpi=210,alpha=False).save(OUT/'preview.png')
spans=[s for b in p.get_text('dict')['blocks'] for l in b.get('lines',[]) for s in l['spans']]
assert min(s['size'] for s in spans)>=4.99
assert not p.get_images()
for s in spans:assert p.rect.contains(fitz.Rect(s['bbox'])),s
print(json.dumps({'output':str(OUT),'dimensions_mm':[p.rect.width*25.4/72,p.rect.height*25.4/72],
 'min_text_pt':min(s['size'] for s in spans),'raster_objects':len(p.get_images()),'filesize':pdf.stat().st_size},indent=2))
