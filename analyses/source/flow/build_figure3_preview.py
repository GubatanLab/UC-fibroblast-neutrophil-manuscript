from pathlib import Path
import sys,json,hashlib,copy,shutil,zipfile
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'analysis/packages'))
import numpy as np,pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
from pypdf import PdfReader,PdfWriter,Transformation
from pypdf.generic import RectangleObject
OUT=ROOT/'outputs/07-figure3-preview';TMP=ROOT/'analysis/regating/figure3_preview';TMP.mkdir(exist_ok=True)
CAN=Path('./Figure 3')
source=CAN/'Figure 3.pdf';original_hash=hashlib.sha256(source.read_bytes()).hexdigest()
T=pd.read_csv(OUT/'FACS_context_interactions.csv');D=pd.read_csv(OUT/'FACS_coculture_effects_figure3.csv');R=pd.read_csv(OUT/'FACS_sample_expression.csv');MARKERS=['OSM','CXCR4','MX1','PADI4','MPO']
plt.rcParams.update({'font.family':'Arial','font.size':7,'axes.labelsize':7,'axes.titlesize':8,'xtick.labelsize':6,'ytick.labelsize':7,'axes.linewidth':.6,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42,'svg.fonttype':'none'})
W,H=518.74,710
fig=plt.figure(figsize=(W/72,H/72));fig.patch.set_facecolor('white')
def text(x,y,s,**kw):fig.text(x/W,1-y/H,s,va='top',**kw)
def axes(x,y,w,h):return fig.add_axes([x/W,1-(y+h)/H,w/W,h/H])
text(9,6,'FIGURE 3 PREVIEW  |  Existing RNA + new FACS analyses',fontsize=8,color='#626D79')
text(9,29,'a   Culture design and readouts',fontweight='bold',fontsize=10)
ax=axes(10,49,238,104);ax.axis('off')
entries=[('Neutrophils only','Untreated / FK-866 / ATN-161'),('Control fibroblast co-culture','Untreated'),('UC fibroblast co-culture','Untreated / FK-866 / ATN-161 / dual')]
for i,(left,right) in enumerate(entries):
    y=.97-i*.27;ax.add_patch(FancyBboxPatch((0,y-.22),.99,.22,boxstyle='round,pad=.006',facecolor=['#F0F3F6','#EBF1F6','#F8EEEB'][i],edgecolor='#B8C2CC',lw=.6));ax.text(.025,y-.045,left,fontsize=7,fontweight='bold',va='top');ax.text(.025,y-.135,right,fontsize=6.4,va='top')
ax.text(0,.08,'RNA: recorded donor matching; independent preculture reference',fontsize=5.8,va='top');ax.text(0,-.03,'FACS: all neutrophils; 47 samples; n = 6/arm except control FB n = 5',fontsize=5.8,va='top')
text(302,174,'d   Protein response to ATN-161',fontweight='bold',fontsize=9)
text(302,190,'UC co-culture: ATN-161 minus untreated',fontsize=6.5,color='#626D79')
def forest(ax,data,estimate,lo,hi,pcol,style='primary'):
    for i,m in enumerate(MARKERS):
        r=data[data.marker==m].iloc[0];x=r[estimate];a=r[lo];b=r[hi]
        ax.errorbar(x,i,xerr=[[x-a],[b-x]],fmt='o',ms=4,color='#0072B2' if x<0 else '#D55E00',elinewidth=.9,capsize=2)
        ax.text(1.02,i,f'{r[pcol]:.4f}',transform=ax.get_yaxis_transform(),fontsize=6.5,va='center')
    ax.axvline(0,color='#8A929C',ls='--',lw=.7);ax.set_yticks(range(5),MARKERS);ax.set_ylim(4.7,-.7);ax.set_xlabel('Effect / untreated reference SD');ax.text(1.02,1.03,'Adj. P',transform=ax.transAxes,fontsize=6.5);ax.spines['left'].set_visible(False);ax.tick_params(axis='y',length=0)
ax=axes(331,216,128,166);forest(ax,D,'standardized_difference','standardized_ci95_low','standardized_ci95_high','p_holm_5_figure3')
text(302,413,'Median antibody fluorescence; 6 samples/group\n95% CIs; exact tests; Holm across 5 markers',fontsize=6.2,color='#626D79')
# e has a preserved exposure-only elastase plot and a new clearly separate MPO measurement.
text(9,441,'e   Release assay and intracellular MPO',fontweight='bold',fontsize=9)
text(11,461,'Existing assay: exposure only',fontsize=6.5,color='#626D79')
text(163,461,'New FACS: UC co-culture',fontsize=6.5,color='#626D79')
ax=axes(171,485,102,122)
for i,g in enumerate(['N+IF','N+IF+ATN']):
    x=R[(R.marker=='MPO')&(R.group==g)].sort_values('sample').median_fi.to_numpy();ax.scatter(i+np.linspace(-.08,.08,len(x)),x,s=13,color=['#0072B2','#D55E00'][i],zorder=3);ax.errorbar(i,x.mean(),yerr=x.std(ddof=1),fmt='_',color='black',capsize=2,markersize=10,lw=.7)
ax.set_xticks([0,1],['Untreated','ATN-161']);ax.set_ylabel('Median MPO fluorescence (a.u.)');ax.set_ylim(-350,5300);ax.set_xlim(-.4,1.4)
p=D[D.marker=='MPO'].iloc[0].p_holm_5_figure3;ax.text(.5,5050,f'Adj. P = {p:.4f}',ha='center',fontsize=6.5)
text(10,630,'Elastase panel: source observations/statistical labels retained.\nNo blockade arm; assay metadata still require confirmation.\nMPO: antibody abundance; points and mean ± SD.\nThese measurements do not establish NET formation.',fontsize=6.2,color='#626D79')
text(302,441,'f   Fibroblast-context interaction',fontweight='bold',fontsize=9)
text(302,459,'ATN effect in co-culture minus its effect alone',fontsize=6.1,color='#626D79')
ax=axes(331,485,128,151);forest(ax,T[T.metric=='median_fi'],'standardized_interaction','standardized_ci_low','standardized_ci_high','p_holm_5')
text(302,660,'Unpaired four-group contrast; 6 samples/arm\nWelch-Satterthwaite inference; Holm across 5 markers',fontsize=6,color='#626D79')
text(9,692,'Preview only. FCS marker labels used. Protein intensity, RNA-state abundance and release are distinct endpoints.',fontsize=6,color='#626D79')
base=TMP/'base.pdf';fig.savefig(base);plt.close(fig)
reader=PdfReader(source);src=reader.pages[0];sw=float(src.mediabox.width);sh=float(src.mediabox.height)
dest=PdfReader(base).pages[0]
placements=[]
def place(name,clip,target):
    x0,y0,x1,y1=clip;tx,ty,tw,th=target;s=min(tw/(x1-x0),th/(y1-y0));p=copy.deepcopy(src);p.cropbox=RectangleObject([x0,sh-y1,x1,sh-y0])
    dx=tx-x0*s;dy=H-ty-(y1-y0)*s-(sh-y1)*s
    dest.merge_transformed_page(p,Transformation().scale(s).translate(dx,dy));placements.append(dict(panel=name,source_top_left_clip=clip,target_top_left_box=target,scale=s))
place('b preserved RNA state interactions',[265,0,518.74,137],[265,27,249,134])
place('c preserved RNA program heatmap',[0,134,317,420],[5,174,285,257])
# Crop the plot only, excluding its old heading and panel label.
place('e preserved elastase exposure assay',[10,459,122,593],[10,480,142,143])
writer=PdfWriter();writer.add_page(dest)
# Supporting all-sample FACS plots and mean-versus-median interaction sensitivity.
fig,axs=plt.subplots(3,2,figsize=(7.2,8.0));groups=['N','N+NF','N+IF','N+IF+FK','N+IF+ATN','N+IF+ATN+FK','N+FK','N+ATN'];labels=['N','N+NF','N+IF','N+IF\n+FK','N+IF\n+ATN','N+IF\n+ATN+FK','N+FK','N+ATN']
for ax,m in zip(axs.flat,MARKERS):
    for i,g in enumerate(groups):
        x=R[(R.marker==m)&(R.group==g)].sort_values('sample').median_fi.to_numpy();ax.scatter(i+np.linspace(-.1,.1,len(x)),x,s=12,color='#D55E00' if 'ATN' in g else '#0072B2',zorder=3);ax.errorbar(i,x.mean(),yerr=x.std(ddof=1),fmt='_',color='black',capsize=2,markersize=8,lw=.6)
    ax.set_title(m,loc='left',fontweight='bold');ax.set_xticks(range(8),labels,fontsize=5.5);ax.set_ylabel('Median fluorescence (a.u.)');ax.ticklabel_format(axis='y',style='plain',useOffset=False)
ax=axs.flat[-1];ax.axis('off');lines=['Interaction sensitivity','Holm-adjusted P (5 markers)','Marker      Median       Mean']
for m in MARKERS:
    a=T[(T.marker==m)&(T.metric=='median_fi')].iloc[0];b=T[(T.marker==m)&(T.metric=='mean_fi')].iloc[0];lines.append(f'{m:<9} {a.p_holm_5:.4f}       {b.p_holm_5:.4f}')
ax.text(0,1,'\n\n'.join(lines),va='top',fontsize=8,fontfamily='DejaVu Sans Mono')
fig.suptitle('Supporting FACS data: all biological samples and sensitivity',x=.08,ha='left',fontsize=11,fontweight='bold');fig.text(.08,.025,'n = 6 per condition except N+NF (n = 5). Original neutrophil parent; no CXCR4 restriction.\nMX1/PADI4 use FCS assignments. Missing viability and treatment/plate blocks limit inference.',fontsize=7)
fig.subplots_adjust(left=.1,right=.97,top=.93,bottom=.12,hspace=.55,wspace=.34);supp=TMP/'support.pdf';fig.savefig(supp);plt.close(fig);writer.add_page(PdfReader(supp).pages[0])
with (OUT/'Figure_3_recommended_preview.pdf').open('wb') as f:writer.write(f)
assert hashlib.sha256(source.read_bytes()).hexdigest()==original_hash
(OUT/'source_panel_preservation.json').write_text(json.dumps({'source':str(source),'sha256':original_hash,'canonical_source_unchanged':True,'placements':placements},indent=2),encoding='utf-8')
legend='''# Recommended Figure 3 preview legend

**Figure 3. Fibroblast exposure and ATN-161 treatment are associated with distinct RNA and protein responses in neutrophils.** (a) Recorded culture conditions and assay readouts. FACS observations are independent biological samples; matching across treatment conditions or to RNA donors is not established. (b) Existing RNA-state CLR interaction panel, preserved from the canonical figure. Its source corrected interactions are nonsignificant and have not been recalculated. (c) Existing 14-program RNA heatmap, preserved with its source pairing, eligible-donor counts, estimates and contrast-specific BH correction. Gene-program scores describe RNA, not measured function. Full source definitions remain in the included original legend. (d) ATN-161-associated protein changes in UC fibroblast co-culture, using median compensated fluorescence in all neutrophils. Effects and unadjusted Welch 95% CIs are divided by the marker-specific pooled within-group SD of untreated N and N+IF samples solely for display. P values are exact two-sided unpaired permutation tests, Holm-adjusted across five markers, n = 6/group. (e) Left: the existing exposure-only elastase-release assay, preserved without reanalysis, including source statistical labels and error bars. Exact assay n, test, correction and error-bar definition require source confirmation; there is no blockade arm. Right: MPO antibody fluorescence in UC co-culture without/with ATN-161, one point per biological sample, mean ± SD, n = 6/group; the adjusted P value is from the same five-marker family as d. Neither measurement alone establishes NET formation. (f) Protein context interactions: [N+IF+ATN - N+IF] - [N+ATN - N], n = 6 independent samples in each of four groups. Points and unadjusted 95% CIs use a heteroskedastic contrast with Welch-Satterthwaite inference and Holm correction across five markers. Display units follow d and are not calibrated protein units or Hedges g intervals.

**Supporting page.** Individual sample medians for all five markers and all eight culture conditions, with mean ± SD, and interaction-test sensitivity using sample mean fluorescence. n = 6 except N+NF (n = 5). Median and mean results can differ; CXCR4 interaction evidence is summary-dependent.

The original RNA contrast scatterplot and predicted ligand-receptor panel are proposed for Extended Data, not deleted from the canonical files. The original canonical legend is included for full RNA and elastase provenance. New FACS families are exploratory and separate from the source RNA testing families and earlier FACS analyses. The FCS data were not assumed paired merely because RNA donors were recorded as paired. No new functional, contact, signaling, neutralization or migration experiment has been performed.

Markers follow FCS labels: OSM R670-A, CXCR4 V610-A, MX1 B710-A, PADI4 B515-A, MPO YG610-A. The protocol's reverse B515/B710 annotation remains unresolved. A viability channel is absent in 44/47 FCS files. Compensation, parent gates and source outcomes were retained. Treatment/plate blocks remain potential confounding. Receptor-specific or fibroblast-specific causality is not established by these observations. This is a proposed layout and analysis preview; the canonical manuscript and figures are unchanged.
'''
(OUT/'Recommended_Figure3_legend.md').write_text(legend,encoding='utf-8');shutil.copy2(CAN/'Figure 3 legend.md',OUT/'Original_Figure3_legend.md');shutil.copy2(__file__,OUT/Path(__file__).name)
print('Created two-page Figure 3 preview; canonical panels preserved.')
