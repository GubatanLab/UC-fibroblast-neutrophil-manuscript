from pathlib import Path
import sys, json, shutil, hashlib
sys.stdout.reconfigure(encoding='utf-8')
ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/'tmp/canonical_nature_20260903/packages'))
import numpy as np
import pandas as pd
from scipy.stats import friedmanchisquare, wilcoxon, rankdata
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.ticker import MaxNLocator
import pymupdf as fitz
from PIL import Image

OUT = ROOT/'output/pdf/Figure_1_Additional_Panels_2026-09-12'
SRC = OUT/'source_data'
QA = OUT/'quality_checks'
QA.mkdir(exist_ok=True)
HUM = ROOT/'output/Nature_Extended_Data_Consolidated_2026-09-03/source_data/human_discovery'
STATES = ['Neutrophil OSM','Neutrophil CXCR4','Neutrophil PADI4','Neutrophil MX1']
PROGS = ['Inflammatory','Recruitment_migration','Retention_aging','Degranulation','Oxidative_burst','Interferon']
LABELS = ['Inflammatory','Recruitment /\nmigration','Retention /\naging','Granule /\nprotease','Oxidative\nburst','Interferon']
FIBPROGS = ['FAP_inflammatory','Neutrophil_recruitment','OSM_response','ECM_remodeling','Alpha5Beta1_adhesion']
FIBLABELS = ['FAP-inflammatory','Neutrophil\nrecruitment','OSM response','Matrix remodeling','α5β1 adhesion']
COLORS = ['#D44B50','#6F4C9B','#E6862B','#3B75B9']
BLUE, RED, GOLD, INK, MUTED = '#3B75B9','#D44B50','#C9A227','#2F3337','#666D75'
CMAP = LinearSegmentedColormap.from_list('study_diverging',[BLUE,'#F7F7F7','#B53A45'])
plt.rcParams.update({'font.family':'Arial','font.size':7,'axes.titlesize':7,'axes.labelsize':6.8,
    'xtick.labelsize':6.3,'ytick.labelsize':6.3,'text.color':INK,'axes.labelcolor':INK,
    'axes.edgecolor':INK,'xtick.color':INK,'ytick.color':INK,'pdf.fonttype':42,'ps.fonttype':42,
    'svg.fonttype':'none','axes.spines.top':False,'axes.spines.right':False,'axes.linewidth':.55,
    'xtick.major.width':.5,'ytick.major.width':.5,'xtick.major.size':2,'ytick.major.size':2})

def bh(p):
    p = np.array(p,dtype=float)
    order = np.argsort(p)
    result = np.empty_like(p)
    result[order] = np.minimum(1,np.minimum.accumulate((p[order]*len(p)/np.arange(1,len(p)+1))[::-1])[::-1])
    return result

def zscore(frame):
    return (frame-frame.mean(axis=0))/frame.std(axis=0,ddof=0)

def aggregate(cells,threshold=20):
    use = cells[cells.condition.eq('Inflamed UC') & cells.neut_state.isin(STATES)].copy()
    ct = use.groupby(['PatientID','neut_state']).size().unstack(fill_value=0).reindex(columns=STATES)
    ids = ct.index[ct.ge(threshold).all(axis=1)].tolist()
    cols = cells.columns[5:].tolist()
    sample = use[use.PatientID.isin(ids)].groupby(['PatientID','neut_state'])[cols].mean()
    means = sample.groupby('neut_state').mean().reindex(STATES)
    return ids, ct, sample, means

def friedman_table(sample):
    rows=[]
    for col in sample.columns:
        x=sample[col].unstack().reindex(columns=STATES)
        res=friedmanchisquare(*[x[s].to_numpy() for s in STATES])
        rows.append({'program':col,'n_patients':len(x),'statistic':res.statistic,
                     'p_value_asymptotic':res.pvalue,'kendall_W':res.statistic/(len(x)*3)})
    tab=pd.DataFrame(rows)
    tab['q_BH_all_9_programs']=bh(tab.p_value_asymptotic)
    return tab

for name in ['09_fibroblast_biopsy_programs.tsv','11_paired_program_tests.tsv','22_neutrophil_state_program_summary.tsv']:
    shutil.copy2(HUM/name,SRC/name)
fib=pd.read_csv(SRC/'09_fibroblast_biopsy_programs.tsv',sep='\t')
ftests=pd.read_csv(SRC/'11_paired_program_tests.tsv',sep='\t')
ftests=ftests[ftests.compartment.eq('Fibroblast')].set_index('program')
fib_ids=sorted(set(fib.loc[fib.condition.eq('Uninflamed UC'),'PatientID']) & set(fib.loc[fib.condition.eq('Inflamed UC'),'PatientID']))
assert len(fib_ids)==14
fibpaired=fib[fib.PatientID.isin(fib_ids) & fib.condition.isin(['Uninflamed UC','Inflamed UC'])].copy()
assert not fibpaired.duplicated(['PatientID','condition']).any()
fibpaired.to_csv(SRC/'panel_D_paired_fibroblast_scores.csv',index=False)
fcheck=[]
for key in ftests.index:
    x=fibpaired.pivot(index='PatientID',columns='condition',values=key)
    p=wilcoxon(x['Inflamed UC'],x['Uninflamed UC'],method='approx',correction=True).pvalue
    assert np.isclose(p,ftests.loc[key,'p_value'],atol=1e-12)
    assert np.isclose(np.median(x['Inflamed UC']-x['Uninflamed UC']),ftests.loc[key,'median_change'])
    fcheck.append(p)
assert np.allclose(bh(fcheck),ftests.FDR.to_numpy(),atol=1e-12)
ftests.to_csv(SRC/'panel_D_all_8_program_statistics.csv')
fib_counts=fibpaired.pivot(index='PatientID',columns='condition',values='n_cells')
fsens=[]
for threshold in [1,5,10,20]:
    eligible=fib_counts.index[fib_counts.ge(threshold).all(axis=1)]
    rows=[]
    for prog in ftests.index:
        x=fibpaired[fibpaired.PatientID.isin(eligible)].pivot(index='PatientID',columns='condition',values=prog)
        diff=x['Inflamed UC']-x['Uninflamed UC']
        rows.append({'min_cells_per_patient_condition':threshold,'program':prog,'n_pairs':len(x),
            'median_change':float(diff.median()),
            'p_asymptotic_continuity_corrected':wilcoxon(diff,method='approx',correction=True).pvalue,
            'p_exact':wilcoxon(diff,method='exact').pvalue})
    t=pd.DataFrame(rows)
    t['q_BH_8_asymptotic']=bh(t.p_asymptotic_continuity_corrected)
    t['q_BH_8_exact']=bh(t.p_exact)
    fsens.append(t)
pd.concat(fsens).to_csv(SRC/'fibroblast_minimum_cell_count_sensitivity.csv',index=False)

original=pd.read_csv(SRC/'neutrophil_original_cell_scores.csv')
clean=pd.read_csv(SRC/'neutrophil_marker_excluded_cell_scores.csv')
pooled=original.groupby('neut_state')[original.columns[5:]].mean()
oldpooled=pd.read_csv(SRC/'22_neutrophil_state_program_summary.tsv',sep='\t').set_index('neut_state')
pool_error=float(np.max(np.abs(pooled.reindex(oldpooled.index).to_numpy()-oldpooled.drop(columns='n').to_numpy())))
assert pool_error < 1e-10
ids, counts, sample, means=aggregate(clean)
_,_,orig_sample,orig_means=aggregate(original)
assert len(ids)==11
stats=friedman_table(sample)
stats.to_csv(SRC/'panel_E_all_9_program_state_tests.csv',index=False)
counts.to_csv(SRC/'inflamed_UC_all_patient_state_cell_counts.csv')
counts.loc[ids].to_csv(SRC/'panel_E_included_patient_state_cell_counts.csv')
sample.to_csv(SRC/'panel_E_marker_excluded_patient_state_means.csv')
means.to_csv(SRC/'panel_E_marker_excluded_equal_patient_means.csv')
zscore(means).to_csv(SRC/'panel_E_marker_excluded_z_scores.csv')
orig_sample.to_csv(SRC/'neutrophil_original_patient_state_means.csv')
orig_means.to_csv(SRC/'neutrophil_original_equal_patient_means.csv')
allstates=original.groupby(['condition','neut_state']).agg(n_cells=('cell','size'),n_patients=('PatientID','nunique'))
allstates.to_csv(SRC/'all_conditions_annotation_control_inventory.csv')
sens=[]
for variant,cells in [('Original',original),('Marker-excluded',clean)]:
    for threshold in [10,20,50]:
        eligible,ct,agg,m=aggregate(cells,threshold)
        stat=friedman_table(agg).set_index('program')
        for prog in m.columns:
            sens.append({'variant':variant,'min_cells_per_state':threshold,'n_complete_patients':len(eligible),
                'program':prog,'highest_mean_state':m[prog].idxmax(),'q_BH_all_9_programs':stat.loc[prog,'q_BH_all_9_programs'],
                'correlation_with_primary_4_state_profile':float(m[prog].corr(means[prog]))})
pd.DataFrame(sens).to_csv(SRC/'neutrophil_threshold_and_marker_sensitivity.csv',index=False)

def save_panel(fig,stem):
    fig.savefig(OUT/(stem+'.pdf'),metadata={'Title':stem.replace('_',' '),'Subject':'Proposed Figure 1 addition; source-based exploratory analysis','Author':'UC fibroblast-neutrophil study'})
    fig.savefig(OUT/(stem+'.svg'))
    fig.savefig(OUT/(stem+'.png'),dpi=450,facecolor='white')
    with Image.open(OUT/(stem+'.png')) as im:
        im.convert('RGB').save(OUT/(stem+'.tiff'),compression='tiff_lzw',dpi=(450,450))
    plt.close(fig)

# Paired observations retain the source scores and the original eight-test correction.
fig=plt.figure(figsize=(183/25.4,72/25.4))
fig.text(.012,.955,'D',fontsize=10,fontweight='bold',va='top')
fig.text(.043,.953,'Paired fibroblast program changes',fontsize=8.2,fontweight='bold',va='top')
fig.text(.043,.874,'Matched noninflamed and inflamed UC tissue | 14 patients',fontsize=6.8,color=MUTED)
for i,(prog,label) in enumerate(zip(FIBPROGS,FIBLABELS)):
    left=.072+i*.185
    ax=fig.add_axes([left,.25,.142,.48])
    x=fibpaired.pivot(index='PatientID',columns='condition',values=prog).reindex(fib_ids)
    # Use the same patient offset at both visits, preserving pairing without concealing tied points.
    offsets=np.linspace(-.046,.046,len(x))
    for j,(_,row) in enumerate(x.iterrows()):
        ax.plot(np.array([0,1])+offsets[j],[row['Uninflamed UC'],row['Inflamed UC']],lw=.62,color='#BCC1C7',zorder=1)
    for condition,pos,color in [('Uninflamed UC',0,GOLD),('Inflamed UC',1,RED)]:
        adequate=fib_counts.reindex(fib_ids)[condition].to_numpy()>=10
        ax.scatter(pos+offsets[adequate],x.loc[adequate,condition],s=11,color=color,edgecolor='white',linewidth=.3,zorder=3)
        ax.scatter(pos+offsets[~adequate],x.loc[~adequate,condition],s=11,facecolor='white',edgecolor=color,linewidth=.75,zorder=3)
        ax.plot([pos-.12,pos+.12],[x[condition].median()]*2,color=INK,lw=1.1,zorder=4)
    ax.set_xlim(-.25,1.25)
    ax.set_ylim(bottom=-x.max().max()*.027,top=x.max().max()*1.13)
    ax.set_xticks([0,1],['Noninflamed','Inflamed'],rotation=25,ha='right')
    ax.tick_params(axis='x',length=0,pad=3)
    ax.yaxis.set_major_locator(MaxNLocator(nbins=4))
    ax.grid(axis='y',color='#E5E7EB',lw=.45,zorder=0)
    ax.set_axisbelow(True)
    ax.set_title(label,fontsize=6.9,fontweight='bold',pad=13)
    ax.text(.5,1.025,'q = '+format(ftests.loc[prog,'FDR'],'.3g'),transform=ax.transAxes,ha='center',fontsize=6.2)
    if i==0: ax.set_ylabel('Mean module expression',labelpad=3)
fig.text(.072,.084,'Lines connect patients; black bars show medians. Open points: fewer than 10 recovered cells.',fontsize=6.2,color=MUTED)
fig.text(.072,.032,'q: paired Wilcoxon test, BH correction across 8 programs; original 14-patient analysis retained.',fontsize=6.2,color=MUTED)
save_panel(fig,'Proposed_Figure_1D_Fibroblast_Programs')

# Primary heatmap: matching patients across all states, existing scores minus label markers.
fig=plt.figure(figsize=(183/25.4,78/25.4))
fig.text(.012,.958,'E',fontsize=10,fontweight='bold',va='top')
fig.text(.043,.956,'Neutrophil state-associated RNA programs',fontsize=8.2,fontweight='bold',va='top')
fig.text(.043,.875,'Inflamed UC | 11 matched patients across all 4 states | ≥20 cells per patient-state',fontsize=6.8,color=MUTED)
ax=fig.add_axes([.192,.304,.69,.43])
mat=zscore(means)[PROGS].to_numpy()
im=ax.pcolormesh(np.arange(7)-.5,np.arange(5)-.5,mat,cmap=CMAP,vmin=-1.8,vmax=1.8,shading='flat',rasterized=False)
ax.set_xlim(-.5,5.5)
ax.set_ylim(3.5,-.5)
ax.set_xticks(range(6),LABELS)
ax.set_yticks(range(4),['OSM-associated','CXCR4-associated','PADI4-associated','MX1-associated'])
ax.tick_params(length=0,pad=5)
for tick,color in zip(ax.get_yticklabels(),COLORS):
    tick.set_color(color); tick.set_fontweight('bold')
for spine in ax.spines.values(): spine.set_visible(False)
ax.set_xticks(np.arange(-.5,6,1),minor=True)
ax.set_yticks(np.arange(-.5,4,1),minor=True)
ax.grid(which='minor',color='white',linewidth=1.1)
ax.tick_params(which='minor',length=0)
q=stats.set_index('program')['q_BH_all_9_programs']
for col,prog in enumerate(PROGS):
    ax.text(col,-.69,'q < 0.001' if q[prog]<.001 else f'q = {q[prog]:.3f}',ha='center',va='bottom',fontsize=6.1,clip_on=False)
    for row in range(4):
        val=mat[row,col]
        ax.text(col,row,f'{val:+.1f}',ha='center',va='center',fontsize=7.1,color='white' if abs(val)>1.1 else INK)
cax=fig.add_axes([.909,.337,.013,.36])
cb=fig.colorbar(im,cax=cax,ticks=[-1.5,0,1.5])
cb.solids.set_rasterized(False)
cb.solids.set_edgecolor('face')
cb.ax.tick_params(labelsize=6.1)
cb.set_label('Program z-score',fontsize=6.3,labelpad=4)
fig.text(.043,.116,'OSM, CXCR4, PADI4 and MX1 excluded from module scoring. Equal patient weights; z-scores across states within each program.',fontsize=6.1,color=MUTED)
fig.text(.043,.061,'q: exploratory within-patient Friedman test, BH correction across 9 programs. RNA programs do not measure effector activity.',fontsize=6.1,color=MUTED)
save_panel(fig,'Proposed_Figure_1E_Neutrophil_Programs')

# One-page composition retains vector artwork at its original width.
combined=fitz.open()
mm=72/25.4
page=combined.new_page(width=183*mm,height=153*mm)
for stem,y,h in [('Proposed_Figure_1D_Fibroblast_Programs',0,72),('Proposed_Figure_1E_Neutrophil_Programs',75,78)]:
    with fitz.open(OUT/(stem+'.pdf')) as src:
        page.show_pdf_page(fitz.Rect(0,y*mm,183*mm,(y+h)*mm),src,0)
combined.set_metadata({'title':'Proposed additional Figure 1 panels D and E','subject':'Review copy; canonical figures and manuscript have not been relabeled'})
combined.save(OUT/'Figure_1_Two_Additional_Panels.pdf',garbage=4,deflate=True)
page.get_pixmap(matrix=fitz.Matrix(250/72,250/72),alpha=False).save(QA/'combined_preview.png')
combined.close()

# Supporting checks are figures for QA, not additional proposed main panels.
fig,axes=plt.subplots(2,3,figsize=(10.8,5.4))
for r,(variant,cells) in enumerate([('Original scores',original),('Marker-excluded scores',clean)]):
    for c,threshold in enumerate([10,20,50]):
        eligible,_,agg,m=aggregate(cells,threshold)
        ax=axes[r,c]
        ax.imshow(zscore(m)[PROGS],cmap=CMAP,vmin=-1.8,vmax=1.8,aspect='auto')
        ax.set_yticks(range(4),[s.replace('Neutrophil ','') for s in STATES])
        ax.set_xticks(range(6),['Inflam.','Recruit.','Retention','Granule','Oxidative','IFN'],rotation=40,ha='right')
        ax.set_title(f'{variant}\n≥{threshold} cells/state; {len(eligible)} complete patients',fontsize=8)
        ax.tick_params(length=0)
fig.suptitle('Sensitivity to state markers and the minimum cell count',fontweight='bold',fontsize=11)
fig.tight_layout(rect=[0,0,1,.94])
fig.savefig(QA/'marker_and_cell_count_sensitivity.png',dpi=220)
plt.close(fig)
fig,axes=plt.subplots(2,3,figsize=(10.8,5.8))
for prog,label,ax in zip(PROGS,LABELS,axes.flat):
    m=sample[prog].unstack().reindex(columns=STATES)
    for patient,row in m.iterrows():
        ax.plot(range(4),row,color='#C3C7CD',lw=.6,alpha=.75)
    for j,color in enumerate(COLORS):
        ax.scatter(np.full(len(m),j),m.iloc[:,j],s=12,color=color,edgecolor='white',lw=.25)
    ax.set_xticks(range(4),[s.replace('Neutrophil ','') for s in STATES])
    ax.set_title(label.replace('\n',' '),fontsize=9,fontweight='bold')
    ax.set_ylabel('Mean module expression')
fig.suptitle('Marker-excluded programs in the same 11 patients',fontweight='bold',fontsize=11)
fig.tight_layout(rect=[0,0,1,.95])
fig.savefig(QA/'individual_patient_profiles.png',dpi=220)
plt.close(fig)

legend=f'''# Proposed additional Figure 1 panels

**D. Paired fibroblast program changes.** Mean gene-module expression in matched noninflamed and inflamed UC tissue from 14 patients. Points represent the mean across recovered stromal/fibroblast-compartment cells within each patient-condition group; lines connect the same patient, and black bars indicate group medians. Open points flag patient-condition groups with fewer than 10 recovered cells. Existing scores are the arithmetic mean of normalized RNA expression across genes present in the original RNA matrix and curated module definitions. Five of the eight originally tested fibroblast programs are shown. Two-sided paired Wilcoxon signed-rank tests used the asymptotic approximation with continuity correction. The displayed q values preserve Benjamini-Hochberg correction across all eight original fibroblast-program tests. Scores, gene definitions and all eight test results are provided in the source data. The α5β1-adhesion program did not meet the q < 0.05 threshold. This analysis compares tissue inflammation conditions across the broader recovered compartment and is distinct from current Figure 1C's inflammatory-versus-other-fibroblast comparison. Cell recovery is sparse in several patient-condition groups: five pairs have fewer than five cells in at least one condition, including three with one noninflamed cell. Requiring at least five cells per condition retains nine pairs and preserves the FAP-inflammatory and recruitment findings at q < 0.05, but not the OSM-response and matrix-remodeling findings; requiring at least 10 or 20 cells retains only six or five pairs, respectively, with no program passing BH correction. All thresholds, effect estimates, original-method and exact-test sensitivity results are included. These limitations should accompany interpretation of the 14-patient result.

**E. Neutrophil state-associated RNA programs.** Heatmap of six selected curated RNA programs across OSM-, CXCR4-, PADI4- and MX1-associated neutrophil states in inflamed UC tissue. The same 11 patients contribute to every state, with at least 20 recovered cells per patient-state; included cell counts are OSM, {counts.loc[ids,STATES[0]].sum():,}; CXCR4, {counts.loc[ids,STATES[1]].sum():,}; PADI4, {counts.loc[ids,STATES[2]].sum():,}; and MX1, {counts.loc[ids,STATES[3]].sum():,}. Existing state assignments were retained. To reduce direct dependence on state-naming genes, OSM, CXCR4, PADI4 and MX1 were removed wherever present in the original module definitions before scores were recomputed from the existing normalized RNA matrix. Cell scores were averaged within each patient-state, then averaged equally across patients. Colors and printed values are z-scores of these equal-patient means across the four states within each program (population SD); they indicate relative profiles within a column and do not compare absolute expression across programs. The original Degranulation module is labeled Granule/protease, and its gene membership is unchanged. Exploratory Friedman tests compare all four states within patients using unscaled patient-state scores; asymptotic P values were BH-adjusted across all nine original neutrophil programs, including the three not displayed. These omnibus tests do not identify individual state pairs as significantly different. Original-score, ≥10-cell and ≥50-cell sensitivity analyses are supplied separately. Healthy/noninflamed samples, the low-RNA state and non-neutrophil annotation controls were excluded from this focused main display; their source summaries and counts are retained. RNA programs are descriptive molecular features and do not establish measured secretion, oxidative activity, NET formation or state transitions. Discovery and coculture modules are separate definitions and analyses.

## Proposed Results text

In the original 14-patient paired analysis, fibroblast FAP-inflammatory, neutrophil-recruitment, OSM-response and extracellular-matrix-remodeling scores increased in inflamed tissue, whereas the broader α5β1-adhesion score did not pass multiple-testing correction (proposed Figure 1D). Sparse cell recovery in several noninflamed samples and loss of statistical support under stricter cell-count thresholds limit these contrasts. A complementary analysis of neutrophil states within the same 11 inflamed UC patients identified distinct RNA-program profiles after exclusion of the state-naming genes: OSM-associated cells had the highest mean inflammatory score, PADI4-associated cells the highest granule/protease and recruitment/migration scores, and MX1-associated cells the highest interferon score (proposed Figure 1E). These complementary stromal and neutrophil profiles motivated spatial assessment of fibroblast-neutrophil neighborhoods and experimental testing of fibroblast effects on neutrophil protein phenotypes and effector-associated release.

## Placement

D and E are proposed positions following current Figure 1C. These files contain only the two requested additions. Existing main figures and the September 12 manuscript have not been changed or renumbered by this operation. The proposed Results text should be inserted only when the new panel numbering is integrated.
'''
(OUT/'Figure_1_Additional_Panel_Legends_and_Results.md').write_text(legend,encoding='utf-8')
(OUT/'Figure_1_Additional_Panel_Legends.txt').write_text(legend,encoding='utf-8')

defs=pd.read_csv(SRC/'original_program_gene_definitions.csv')
coverage=pd.read_csv(SRC/'neutrophil_program_gene_coverage.csv')
assert coverage.present.all()
audit=pd.read_csv(SRC/'expression_reproduction_audit.csv')
assert audit.max_abs_difference.max()<1e-7
report={'fibroblast_paired_patients':14,'fibroblast_original_q_values_reproduced':True,
    'fibroblast_pairs_with_under_5_cells_in_either_condition':int(fib_counts.lt(5).any(axis=1).sum()),
    'fibroblast_pairs_eligible_at_min_10_cells':int(fib_counts.ge(10).all(axis=1).sum()),
    'neutrophil_included_patients':ids,'neutrophil_included_cells':int(counts.loc[ids].to_numpy().sum()),
    'primary_condition':'Inflamed UC','primary_min_cells_per_patient_state':20,
    'equal_patient_weights':True,'state_markers_excluded':['OSM','CXCR4','PADI4','MX1'],
    'original_pooled_summary_max_abs_error':pool_error,
    'original_cell_score_max_abs_reproduction_error':float(audit.max_abs_difference.max()),
    'inference':'exploratory Friedman asymptotic omnibus tests; BH across all 9 neutrophil programs',
    'no_pairwise_state_inference':True,'existing_canonical_assets_modified':False,
    'important_sensitivity':'CXCR4-associated cells are not the highest retention/aging state; this program includes genes beyond CXCR4, and its profile changes upon marker exclusion.',
    'sources':[]}
for path in [Path('input_data/human/run_uc_fibro_neutrophil_analysis.R'),
             Path('input_data/human/UC_fibroblast_neutrophil_analysis/objects/neutrophil_cell_scores.rds'),
             HUM/'09_fibroblast_biopsy_programs.tsv',HUM/'11_paired_program_tests.tsv',HUM/'22_neutrophil_state_program_summary.tsv']:
    report['sources'].append({'path':str(path),'sha256':hashlib.sha256(path.read_bytes()).hexdigest()})
report['raw_expression_object']=str(Path('input_data/mouse/UCGNE_fibroblast_neutrophil_subclustering/objects/UCGNE_colon_neutrophils_subclustered.rds'))
(OUT/'analysis_verification.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
(OUT/'README.md').write_text('''# Additional Figure 1 panels - September 12, 2026

Two proposed main panels, supplied individually as vector PDF and SVG, plus 450-dpi PNG/TIFF and a combined one-page vector PDF. Panel letters D and E are proposed insertion positions. Canonical figures and manuscript are unchanged.

- D: paired fibroblast programs, 14 patients; original eight-test correction preserved.
- E: neutrophil state programs, 11 matched inflamed UC patients; equal patient weights; state-naming genes excluded; exploratory repeated-patient omnibus tests with correction across all nine programs.
- `source_data`: original and derived scores, exact gene lists, full statistics, patient-state cell counts and sensitivity results.
- `quality_checks`: original-versus-marker-excluded heatmaps at multiple cell-count thresholds and individual patient profiles.
- `Figure_1_Additional_Panel_Legends_and_Results.md`: full legends, analysis definitions and proposed Results text.
- `analysis_verification.json`: reproduction checks and source provenance.

Interpretation: the new heatmap compares states within inflamed UC. It is not evidence for significant whole-neutrophil program changes between inflamed and noninflamed tissue. Granule/protease, retention/aging and oxidative-burst labels refer to gene lists, not directly measured functions. Inclusion depends on cell recovery and the complete-state criterion, so the 11-patient analysis is a selected subset of the tissue cohort.
''',encoding='utf-8')
print(json.dumps({'output':str(OUT),'patients':ids,'cell_counts':counts.loc[ids].sum().to_dict(),'neutrophil_tests':stats.to_dict('records'),'pdfs':[p.name for p in OUT.glob('*.pdf')]},indent=2))
