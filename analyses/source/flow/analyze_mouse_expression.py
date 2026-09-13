"""Mouse follow-up: all-neutrophil expression, fixed contrasts and FMO sensitivity."""
from pathlib import Path
import sys,json,hashlib,itertools,shutil,zipfile,html
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'analysis/mouse_packages'))
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

OUT=ROOT/'outputs/06-mouse-neutrophils';OUT.mkdir(exist_ok=True)
CACHE=ROOT/'analysis/regating/events'
MARKERS={'OSM':'B515-A','CXCR4':'U820-A','MX1':'Y670-A','PADI4':'R710-A'}
ORDER=['WD','FD','FC','FG'];PAIRS=[('WD','FD'),('FC','FG')]
GROUPS={'WD':'WT + DSS','FD':'FAP + GCV + DSS','FC':'FAP + PBS + water','FG':'FAP + GCV + water'}
COLORS={'WD':'#34699a','FD':'#d47934','FC':'#427f66','FG':'#9566a4'}
prior=pd.read_csv(ROOT/'outputs/02-gating-reanalysis/revised_marker_measurements.csv')
prior=prior[(prior.experiment=='mouse')&(prior['quantile']==.99)]
inventory=pd.read_csv(ROOT/'analysis/results/raw_fcs_inventory.csv').set_index('file')
rows=[];events={};nfiles=0
for mp in sorted(CACHE.glob('mouse_*.json')):
    m=json.loads(mp.read_text())
    if not m['biological']:continue
    digest=hashlib.sha256((ROOT/m['file']).read_bytes()).hexdigest()
    assert digest==inventory.loc[m['file'],'sha256']
    nfiles+=1
    with np.load(CACHE/(m['key']+'.npz')) as z:
        parent=z['old_5'];n=int(parent.sum())
        for marker,ch in MARKERS.items():
            assert m['markers'][m['channels'].index(ch)].upper()==marker.upper()
            v=z['events'][parent,m['channels'].index(ch)].copy()
            assert np.isfinite(v).all() and ch in m['matrix_info']['detectors']
            old=prior[(prior['sample']==m['sample'])&(prior.channel==ch)].iloc[0]
            assert n==old.neutrophil_events and np.isclose(np.median(v),old.median_fi,atol=1e-8,rtol=0)
            events[m['sample'],marker]=v
            rows.append(dict(sample=m['sample'],group=m['group'],marker=marker,channel=ch,neutrophil_events=n,low_parent_n=n<100,median_fi=np.median(v),mean_fi=v.mean(),q25=np.quantile(v,.25),q75=np.quantile(v,.75),fmo_percent=old.percent,fmo_positive_events=old.positive_events,source_file=m['file'],source_sha256=digest,cache_sha256=hashlib.sha256((CACHE/(m['key']+'.npz')).read_bytes()).hexdigest(),compensation_sha256=m['matrix_info']['matrix_sha256']))
R=pd.DataFrame(rows);assert len(R)==92 and nfiles==23
assert R.groupby('group')['sample'].nunique().to_dict()=={'FC':6,'FD':5,'FG':6,'WD':6}
R.to_csv(OUT/'Mouse_sample_measurements.csv',index=False)
S=R.groupby(['marker','channel','group'],sort=False).agg(n=('sample','size'),mean_sample_median=('median_fi','mean'),sd_sample_median=('median_fi','std'),mean_sample_mean=('mean_fi','mean'),mean_fmo_percent=('fmo_percent','mean'),sd_fmo_percent=('fmo_percent','std'),min_neutrophil_events=('neutrophil_events','min'),max_neutrophil_events=('neutrophil_events','max')).reset_index()
S.to_csv(OUT/'Mouse_group_summary.csv',index=False)
tests=[];checked=0
for metric in ['median_fi','mean_fi','fmo_percent']:
    for marker in MARKERS:
        for ref,trt in PAIRS:
            x=R[(R.marker==marker)&(R.group==trt)][metric].to_numpy()
            y=R[(R.marker==marker)&(R.group==ref)][metric].to_numpy()
            pooled=np.r_[x,y];diff=x.mean()-y.mean()
            alloc=np.array(list(itertools.combinations(range(len(pooled)),len(x))))
            sums=pooled[alloc].sum(axis=1);null=sums/len(x)-(pooled.sum()-sums)/len(y)
            p=float(np.mean(abs(null)>=abs(diff)-max(1e-9,abs(diff)*1e-12)))
            check=stats.permutation_test((x,y),lambda a,b:abs(a.mean()-b.mean()),vectorized=False,n_resamples=np.inf,alternative='greater').pvalue
            assert np.isclose(p,check,rtol=0,atol=1e-12),(marker,metric,p,check)
            checked+=1
            vx=x.var(ddof=1)/len(x);vy=y.var(ddof=1)/len(y)
            df=(vx+vy)**2/(vx**2/(len(x)-1)+vy**2/(len(y)-1)) if vx+vy else np.nan
            half=stats.t.ppf(.975,df)*np.sqrt(vx+vy)
            tests.append(dict(marker=marker,metric=metric,reference=ref,comparison=trt,n_reference=len(y),n_comparison=len(x),mean_reference=y.mean(),mean_comparison=x.mean(),difference=diff,ci95_low=diff-half,ci95_high=diff+half,p_exact=p,allocations=len(null)))
T=pd.DataFrame(tests)
for metric,ix in T.groupby('metric').groups.items():
    ids=T.loc[ix].sort_values('p_exact',kind='stable').index
    T.loc[ids,'p_holm_8']=np.minimum(1,np.maximum.accumulate(T.loc[ids,'p_exact'].to_numpy()*np.arange(8,0,-1)))
T.to_csv(OUT/'Mouse_comparisons.csv',index=False)
primary=T[T.metric=='median_fi']
freq=pd.read_csv(ROOT/'outputs/02-gating-reanalysis/revised_marker_measurements.csv')
freq[freq.experiment=='mouse'].to_csv(OUT/'Mouse_FMO_threshold_sensitivity.csv',index=False)
thr=pd.read_csv(ROOT/'outputs/02-gating-reanalysis/gate_thresholds.csv')
thr[thr.experiment=='mouse'].to_csv(OUT/'Mouse_FMO_thresholds.csv',index=False)
for filename in ['dss_corrected_sample_key.csv','dss_preserved_exclusions.csv']:
    shutil.copy2(ROOT/'analysis/results'/filename,OUT/filename)

plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'axes.spines.top':False,'axes.spines.right':False,'svg.fonttype':'none','pdf.fonttype':42})
reps=[]
with PdfPages(OUT/'Mouse_neutrophil_figures.pdf') as pdf:
    for marker,ch in MARKERS.items():
        fig,axes=plt.subplots(2,2,figsize=(10,7.8))
        fig.suptitle(f'{marker} | Mouse neutrophils',fontsize=18,x=.08,ha='left',y=.975)
        for ax,metric,ylabel in [(axes[0,0],'median_fi','Median compensated fluorescence (a.u.)'),(axes[0,1],'fmo_percent','Above FMO threshold (% of neutrophils)')]:
            for i,g in enumerate(ORDER):
                rr=R[(R.marker==marker)&(R.group==g)].sort_values('sample');v=rr[metric].to_numpy()
                xx=i+np.linspace(-.10,.10,len(v))
                ax.scatter(xx,v,c=COLORS[g],s=30,zorder=3)
                low=rr.low_parent_n.to_numpy()
                ax.scatter(xx[low],v[low],facecolors='none',edgecolors='black',s=75,zorder=4)
                ax.errorbar(i,v.mean(),yerr=v.std(ddof=1),fmt='_',color='black',capsize=4,markersize=15)
            ax.set_xticks(range(4),[f'{g}\nn = {6 if g!="FD" else 5}' for g in ORDER]);ax.set_ylabel(ylabel)
            tt=T[(T.marker==marker)&(T.metric==metric)]
            ax.set_title(' | '.join(f'{r.comparison} vs {r.reference}: P = {r.p_holm_8:.3g}' for r in tt.itertuples()),fontsize=9,pad=12)
        pooled=np.concatenate([np.arcsinh(v/150) for (s,mk),v in events.items() if mk==marker])
        bins=np.linspace(pooled.min()-.05,pooled.max()+.05,75)
        for ax,(ref,trt) in zip(axes[1],PAIRS):
            for g in [ref,trt]:
                rr=R[(R.marker==marker)&(R.group==g)].copy()
                rr['distance']=abs(rr.median_fi-rr.median_fi.median())
                rep=rr.sort_values(['distance','sample']).iloc[0]['sample']
                reps.append(dict(marker=marker,group=g,sample=rep,rule='Closest sample median to group median of sample medians; sample ID breaks ties'))
                for sample in sorted(rr['sample']):
                    counts,_=np.histogram(np.arcsinh(events[sample,marker]/150),bins)
                    ax.stairs(counts/counts.max()*100,bins,color=COLORS[g],alpha=.19,lw=.7)
                counts,_=np.histogram(np.arcsinh(events[rep,marker]/150),bins)
                ax.stairs(counts/counts.max()*100,bins,color=COLORS[g],lw=1.8,label=f'{g}: {rep}')
            ticks=[v for v in [-10000,-1000,0,1000,10000,100000,1000000] if bins[0]<=np.arcsinh(v/150)<=bins[-1]]
            ax.set_xticks(np.arcsinh(np.array(ticks)/150),[f'{v:,}' for v in ticks]);ax.tick_params(axis='x',labelsize=7)
            ax.set_xlabel(f'{marker} fluorescence (a.u.; asinh scale)');ax.set_ylabel('Normalized to mode (%)');ax.legend(frameon=False,fontsize=8)
            ax.set_title(f'{GROUPS[trt]} vs {GROUPS[ref]}',fontsize=10)
        footer='WD: WT + DSS   |   FD: FAP + GCV + DSS   |   FC: FAP + PBS + water   |   FG: FAP + GCV + water\nPoints: individual mice; bars: mean +/- SD. Black rings: <100 neutrophils. Thick histograms: representatives; thin: all mice.\nP values: exact unpaired tests, Holm-adjusted across 8 comparisons separately for intensity and FMO frequency.'
        if marker=='PADI4':footer+='\nPADI4 detector identity matches mouse FCS labels; the R710 FMO omission identity remains provisional.'
        fig.text(.08,.022,footer,fontsize=8,va='bottom',linespacing=1.5)
        fig.subplots_adjust(left=.09,right=.97,top=.88,bottom=.23,hspace=.52,wspace=.30)
        pdf.savefig(fig);fig.savefig(OUT/f'Mouse_{marker}.png',dpi=220);fig.savefig(OUT/f'Mouse_{marker}.svg');plt.close(fig)
pd.DataFrame(reps).to_csv(OUT/'Representative_samples.csv',index=False)

def table(headers,rr):return '| '+' | '.join(headers)+' |\n|'+'|'.join(['---']*len(headers))+'|\n'+'\n'.join('| '+' | '.join(map(str,r))+' |' for r in rr)
report='''# Mouse neutrophil FACS expression analysis

All 23 available biological mouse samples were analyzed using the previously corrected compensated event cache. Primary endpoints are per-mouse median OSM, CXCR4, MX1 and PADI4 fluorescence in the entire retained live-neutrophil parent, including CXCR4-negative cells. No marker-positive restriction was applied to expression measurements. The mouse panel does not contain MPO, so an MPO analysis cannot be reproduced in these mice.

## Design and sample inclusion

WD = WT+DSS (n=6); FD = FAP+GCV+DSS (n=5); FC = FAP+PBS+water (n=6); FG = FAP+GCV+water (n=6). These group descriptions are inherited from the corrected dataset analysis; the exact FAP model construct is not established here. FD3 is absent. WD3 retains its corrected identity. Historical undocumented Prism exclusion flags are preserved in the supplied audit but do not remove mice from this analysis. All 23 mice contribute one value per marker; events are not treated as independent replicates.

The two contrasts follow the existing mouse analysis: FD minus WD and FG minus FC. FD versus WD differs in genotype and treatment, and cannot isolate GCV treatment or fibroblast depletion. FG versus FC compares treatment within the recorded FAP/water condition. There is no complete matched genotype-by-DSS-by-treatment design. Pairing, randomization and experimental blocks are unconfirmed; tests assume independent mice.

## Primary expression results

Only OSM fluorescence in FG versus FC meets the adjusted P < 0.05 threshold (Holm P = 0.03463). The difference is -12,009.32 a.u.; group means of mouse medians are 769.40 (FC) and -11,239.92 (FG). This shift includes strongly negative compensated signals and requires staining/compensation review before interpretation as reduced OSM protein. No primary marker meets the adjusted threshold in FD versus WD. The other water-condition markers do not meet it either; lack of significance is not evidence of equivalence. No FMO-frequency comparison meets the adjusted threshold.

Exact two-sided unpaired label-permutation tests use absolute differences in group means of mouse medians. Every allocation is enumerated: 462 for FD versus WD and 924 for FG versus FC. Holm adjustment covers all eight primary tests (four markers, two contrasts). Secondary arithmetic-mean fluorescence and FMO frequency each have a separate eight-test family. This is exploratory follow-up, without adjustment jointly across historical analyses. Welch 95% confidence intervals are unadjusted and descriptive.

'''+table(['Marker','Contrast','Reference mean','Comparison mean','Difference (95% CI), a.u.','Holm P'],[[r.marker,f'{r.comparison} - {r.reference}',f'{r.mean_reference:.2f}',f'{r.mean_comparison:.2f}',f'{r.difference:.2f} ({r.ci95_low:.2f}, {r.ci95_high:.2f})',f'{r.p_holm_8:.5f}'] for r in primary.itertuples()])+'''

## Gating, quality and interpretation

The original CELLS/s1/s2/cell1/cd11b+/neutrophil parent retains scatter, singlets, viability and CD45/CD11b/Ly6G lineage selection. Original compensation coefficients were preserved using the corrected physical-detector aliases. Parent gates were retained, not independently redesigned or certified. Sample parent sizes range from 68 to 1,691 events. FC3 (68) and FC4 (72) are flagged but retained; sparse samples limit precision, especially for rare-positive estimates. Native FlowJo recalculation has not been performed. Prior corrected-import comparisons against saved marker gates still had sparse-gate discrepancies, up to 9.411 percentage points for FD5 CXCR4; agreement with the preceding corrected cache does not resolve that limitation.

Mouse marker labels map OSM to B515-A, CXCR4 to U820-A, MX1 to Y670-A and PADI4 to R710-A. The co-culture MX1/PADI4 annotation conflict does not apply to these mouse detector labels. The R710 FMO is only provisionally identified as the PADI4-omission control. OSM, MX1 and R710 FMO parents have 55, 41 and 95 events; CXCR4 has 173. Primary FMO thresholds are immediately above the conservative empirical 99th percentile, with 98th and 99.5th percentile results supplied for sensitivity. Broad control distributions and sparse FMO parents prevent confident rare-positive assignments. FMO uncertainty is not included in P values. Zero above-threshold events do not establish biological absence.

Fluorescence values are compensated arbitrary units, not calibrated protein quantities or enzyme activity. Negative values are retained. No pooled-event inference, background subtraction, imputation or outcome-driven exclusion was used. Frequencies describe the gated neutrophil parent and are not tissue cell counts. Differences do not establish a cell-specific mechanism.

## Figure legend and verification

Each marker page shows all individual mouse medians and FMO frequencies with mean +/- SD, plus both contrast histograms. Black rings identify low-count mice. Thin lines show all mice; thick lines show the sample closest to the group median of sample medians, breaking ties by sample ID. Representatives need not be the same mouse across markers and are not paired. Histograms use unsmoothed compensated signals, asinh cofactor 150, identical bin edges within each marker, and independent normalization to peak bin height. Histogram height does not measure cell abundance. FMO gates are supporting review results.

All 23 raw-file SHA-256 hashes match the previous inventory; all 92 marker medians and parent counts match the corrected reanalysis. All 24 exact P values were checked independently with SciPy permutation enumeration using the same absolute-difference statistic. Source hashes, cache hashes and compensation hashes are recorded per sample. These checks establish computational consistency, not biological gate validity. Original source files and previous analyses are unchanged.

Reproduce by running analysis/analyze_mouse_expression.py from this project with the corrected event cache and numpy, pandas, scipy and matplotlib installed. The ZIP contains source tables, threshold sensitivity, sample identity/exclusion audits, figures and this script; large raw files and event caches remain in the project.
'''
(OUT/'Mouse_neutrophil_report.md').write_text(report,encoding='utf-8')
body='<h1>Mouse neutrophil FACS analysis</h1><p>23 mice | Four markers | All retained live neutrophils</p>'+''.join(f'<h2>{m}</h2><img src="Mouse_{m}.png" style="width:100%;max-width:1100px">' for m in MARKERS)+'<h2>Full report and methods</h2><pre style="white-space:pre-wrap;font:15px/1.6 system-ui">'+html.escape(report)+'</pre>'
(OUT/'Mouse_neutrophil_report.html').write_text('<!doctype html><meta charset="utf-8"><title>Mouse neutrophil analysis</title><body style="max-width:1150px;margin:40px auto;padding:20px;font-family:system-ui;color:#173342">'+body,encoding='utf-8')
(OUT/'verification.json').write_text(json.dumps({'raw_hashes_checked':nfiles,'medians_and_parent_counts_checked':len(R),'exact_tests_checked':checked,'passed':True,'versions':{'numpy':np.__version__,'pandas':pd.__version__,'scipy':__import__('scipy').__version__,'matplotlib':matplotlib.__version__}},indent=2))
shutil.copy2(__file__,OUT/Path(__file__).name)
manifest=[dict(file=p.name,sha256=hashlib.sha256(p.read_bytes()).hexdigest(),bytes=p.stat().st_size) for p in sorted(OUT.iterdir()) if p.is_file() and p.name!='manifest.json']
(OUT/'manifest.json').write_text(json.dumps(manifest,indent=2))
with zipfile.ZipFile(ROOT/'outputs/Mouse_neutrophil_analysis.zip','w',zipfile.ZIP_DEFLATED) as z:
    for p in OUT.iterdir():
        if p.is_file():z.write(p,'Mouse_neutrophils/'+p.name)
print(primary.round(5).to_string(index=False))
print('Verified',nfiles,'mouse files;',len(R),'marker measurements;',checked,'exact tests')
