from pathlib import Path
import sys,json,re
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'analysis/packages'))
import numpy as np,pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
OUT=ROOT/'analysis/regating';DEST=ROOT/'outputs/02-gating-reanalysis';FIG=DEST/'figures';FIG.mkdir(exist_ok=True)
QA=OUT/'pdf_qa';QA.mkdir(exist_ok=True)
plt.rcParams.update({'font.family':'Arial','font.size':10,'svg.fonttype':'none','pdf.fonttype':42,'axes.spines.top':False,'axes.spines.right':False})
meta=json.loads((OUT/'event_cache_index.json').read_text());R=pd.read_csv(DEST/'revised_marker_measurements.csv');R=R[R['quantile']==.99];T=pd.read_csv(DEST/'gate_thresholds.csv');T=T[T['quantile']==.99]
COLORS={'N':'#596673','N+NF':'#3E8196','N+IF':'#BE5E35','N+IF+FK':'#D39252','N+IF+ATN':'#327E69','N+IF+ATN+FK':'#735890','N+FK':'#959FAA','N+ATN':'#99BFA9','WD':'#BD633E','FD':'#267A78','FC':'#667A91','FG':'#8D699E'}
CO=['N','N+NF','N+IF','N+IF+FK','N+IF+ATN','N+IF+ATN+FK','N+FK','N+ATN']
LABELS={'N':'N','N+NF':'N + NF','N+IF':'N + IF','N+IF+FK':'N + IF\n+ FK','N+IF+ATN':'N + IF\n+ ATN','N+IF+ATN+FK':'N + IF\n+ ATN + FK','N+FK':'N + FK','N+ATN':'N + ATN'}
def save(fig,name):
    for ext in ['png','svg']:fig.savefig(FIG/f'{name}.{ext}',dpi=180,bbox_inches='tight')
    plt.close(fig)
for kind in ['co','mouse']:
    for metric,title in [('percent','Above-control fractions'),('median_fi','Median fluorescence per sample')]:
        fig,axes=plt.subplots(2,2,figsize=(12,8));order=CO if kind=='co' else ['WD','FD','FC','FG']
        for ax,marker in zip(axes.flat,T[T.experiment==kind].marker):
            q=R[(R.experiment==kind)&(R.marker==marker)]
            for i,g in enumerate(order):
                values=q[q.group==g].sort_values('sample')[metric].dropna().to_numpy();jitter=np.linspace(-.12,.12,len(values))
                ax.scatter(i+jitter,values,color=COLORS[g],s=32,edgecolors='white',lw=.5,zorder=3)
                ax.errorbar(i,np.mean(values),yerr=np.std(values,ddof=1),fmt='_',color='black',capsize=4,markersize=18,lw=1)
            ax.set_title(marker+(' signal; annotation unresolved' if marker in ['B515','B710'] else ''),loc='left',fontweight='bold',fontsize=11)
            ax.set_xticks(range(len(order)),[LABELS.get(g,g) for g in order],fontsize=8)
            ax.set_ylabel('Percent of original neutrophil parent' if metric=='percent' else 'Median compensated fluorescence (a.u.)',fontsize=9)
            ax.grid(axis='y',color='#E4E8EC');ax.set_axisbelow(True)
            if metric=='percent':ax.set_ylim(0,105)
            else:ax.axhline(0,color='#8B949D',lw=.6)
        fig.suptitle(('Co-culture' if kind=='co' else 'DSS mice')+': '+title+' - review analysis',x=.06,ha='left',fontsize=16,fontweight='bold')
        foot='All available biological samples; points and mean +/- SD. Control 99th-percentile thresholds; original neutrophil parents retained.\n'
        foot+=('Co-culture: above-unstained sensitivity analysis; no matched FMOs or uniform viability data. B515/B710 identities unresolved.\nN=neutrophils; NF/IF=non-inflamed/inflamed fibroblasts; FK=FK-866; ATN=ATN-161.' if kind=='co' else 'DSS: matched marker FMOs have 41-173 neutrophils. R710 is provisionally the PADI4 FMO.\nWD=WT+DSS; FD=FAP+GCV+DSS; FC=FAP+PBS+water; FG=FAP+GCV+water. Frequency is not tissue cell count.')
        fig.text(.06,.018,foot,fontsize=8.5,color='#515B66');fig.subplots_adjust(top=.9,bottom=.17,hspace=.38,wspace=.3);save(fig,f'{kind}_{metric}')

def natural(s):return [int(x) if x.isdigit() else x for x in re.split(r'(\d+)',s)]
samples=sorted([m for m in meta if m['biological']],key=lambda m:(m['experiment'],natural(m['sample'])))
with PdfPages(DEST/'Neutrophil_gating_review.pdf') as pdf:
    for page,m in enumerate(samples,1):
        kind=m['experiment'];z=np.load(OUT/'events'/(m['key']+'.npz'));ev=z['events'];parent=z['old_3' if kind=='co' else 'old_5'];params=T[T.experiment==kind].set_index('marker');n=int(parent.sum());markers=list(params.index)
        fig,axes=plt.subplots(2,3,figsize=(14,9));axes=axes.ravel()
        # Back-gating onto lineage channels: gray singlets and all selected neutrophils.
        xch,ych=('V710-A','UV586-A') if kind=='co' else ('V610-A','U670-A');xx=ev[:,m['channels'].index(xch)];yy=ev[:,m['channels'].index(ych)]
        display=np.linspace(0,len(ev)-1,min(len(ev),4500),dtype=int)
        axes[0].scatter(np.arcsinh(xx[display]/150),np.arcsinh(yy[display]/150),s=2,color='#A9B0B5',alpha=.25,rasterized=True)
        axes[0].scatter(np.arcsinh(xx[parent]/150),np.arcsinh(yy[parent]/150),s=3,color='#23677F',alpha=.5,rasterized=True)
        axes[0].set_xlabel(('CD16' if kind=='co' else 'Ly6G')+' - '+xch);axes[0].set_ylabel('CD11b - '+ych);axes[0].set_title('Original neutrophil parent: back-gating',loc='left',fontsize=10,fontweight='bold')
        for ax,marker in zip(axes[1:5],markers):
            row=params.loc[marker];control_name=Path(row.control_file).name.removesuffix('.fcs');cm=next(x for x in meta if x['experiment']==kind and x['sample']==control_name)
            cz=np.load(OUT/'events'/(cm['key']+'.npz'));cp=np.ones(len(cz['events']),bool) if kind=='co' else cz['old_5'];cv=cz['events'][cp,cm['channels'].index(row.channel)]
            vals=ev[parent,m['channels'].index(row.channel)];tx=np.arcsinh(vals/150);tc=np.arcsinh(cv/150);edges=np.linspace(min(tx.min(),tc.min())-.1,max(tx.max(),tc.max())+.1,75)
            ax.hist(tc,bins=edges,density=True,histtype='stepfilled',color='#B66D3E',alpha=.3,label=f'Control (n={len(cv)})')
            ax.hist(tx,bins=edges,density=True,histtype='step',color='#23677F',lw=1.5,label=f'{m["sample"]} (n={n})')
            ax.axvline(np.arcsinh(row.threshold/150),color='#272F35',ls='--',lw=1,label='Review cutoff')
            count=int(np.sum(vals>=row.threshold));ax.set_title(f'{marker}: {count}/{n} ({100*count/n:.2f}%)',loc='left',fontsize=10,fontweight='bold');ax.set_xlabel(row.channel+' (compensated)');ax.set_ylabel('Density');ax.legend(fontsize=7)
        xr=params.loc['CXCR4'];yr=params.loc['OSM'];xx=ev[parent,m['channels'].index(xr.channel)];yy=ev[parent,m['channels'].index(yr.channel)];joint=(xx>=xr.threshold)&(yy>=yr.threshold)
        axes[5].scatter(np.arcsinh(xx/150),np.arcsinh(yy/150),c=np.where(joint,'#327E69','#8299AA'),s=4,alpha=.5,rasterized=True)
        axes[5].axvline(np.arcsinh(xr.threshold/150),color='black',ls='--',lw=1);axes[5].axhline(np.arcsinh(yr.threshold/150),color='black',ls='--',lw=1);axes[5].set_xlabel('CXCR4 - '+xr.channel);axes[5].set_ylabel('OSM - '+yr.channel);axes[5].set_title(f'Consistent CXCR4 AND OSM: {joint.sum()}/{n}',loc='left',fontsize=10,fontweight='bold')
        fig.suptitle(f'{"Co-culture" if kind=="co" else "DSS mice"} | {m["sample"]} | {m["group"]} | {n:,} neutrophil events',x=.06,ha='left',fontsize=17,fontweight='bold')
        note='Co-culture gates show signal above the available unstained control; matched FMOs and uniform viability data are missing.' if kind=='co' else 'Mouse marker gates use matching FMOs. Small FMO parent counts limit tail precision; R710/PADI4 FMO identity is provisional.'
        fig.text(.06,.045,note+'\nAll axes show asinh(fluorescence / 150). All events contribute to counts; gray parent background alone is capped for display.\nOriginal scatter/singlet/lineage gates retained. B515/B710 co-culture marker annotations remain unresolved.',fontsize=8.5,color='#515B66')
        fig.text(.94,.035,str(page),ha='right',fontsize=9);fig.subplots_adjust(left=.06,right=.96,top=.9,bottom=.16,hspace=.4,wspace=.33)
        pdf.savefig(fig)
        if m['sample'] in ['P1-1','P3-1','P4-1','P5-1','WD1','FD5','FC3','FG1']:fig.savefig(QA/f'page_{page:02d}_{m["sample"]}.png',dpi=100)
        plt.close(fig)
        if page%10==0:print('Plotted',page,'of',len(samples),flush=True)
print('Created review PDF with',len(samples),'pages')
