# Academic Figure Skill Asset Confirmation (verified against assets/figures/)
# (a,b) threshold curves -> LineTrend/plot_trend.py -> cross-type param inherit
# (supp a-d) paired-error curves -> LineTrend/plot_trend.py -> cross-type param inherit
# Inspected trend_by_month.png. Cumulative publication areas are semantically
# different from paired temperature-threshold effects, so use baseline colours,
# open axes and line/band structure; draw all source thresholds and annual curves.
# RULE: native run = load pre-rendered PNG; param inherit = draw from source data.
# Academic Figure Skill Typography Baseline — COPY VERBATIM, place at TOP of script
import matplotlib as mpl
mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "Liberation Sans"],
    "font.size": 8,
    "axes.titlesize": 8,
    "axes.labelsize": 8,
    "xtick.labelsize": 7,
    "ytick.labelsize": 7,
    "legend.fontsize": 8,
    "figure.titlesize": 9,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.linewidth": 0.6,
    "xtick.direction": "out",
    "ytick.direction": "out",
    "xtick.major.width": 0.6,
    "ytick.major.width": 0.6,
    "legend.frameon": False,
})

# Academic Figure Skill Nature/Cell/Science Color Palette -- COPY VERBATIM
CATEGORICAL = ["#2166AC", "#B2182B", "#1B7837", "#F1A340", "#762A83", "#666666"]
CATEGORICAL_EXTENDED = [
    "#2166AC", "#B2182B", "#1B7837", "#F1A340", "#762A83", "#666666",
    "#4393C3", "#D6604D", "#5AAE61", "#B35806", "#9970AB", "#999999",
]
DIVERGING   = ["#2166AC", "#F7F7F7", "#B2182B"]
SEQUENTIAL  = ["#F7FBFF", "#6BAED6", "#08306B"]
ACCENT_RED  = "#B2182B"
GREY        = "#999999"
BLACK       = "#222222"

# Academic Figure Skill Export Baseline — COPY VERBATIM
mpl.rcParams.update({
    "pdf.fonttype": 42,         # TrueType font embedding
    "svg.fonttype": "none",     # editable text in SVG
    "savefig.bbox": "tight",    # trim whitespace
    "savefig.dpi": 300,
})

def save_cns_figure(fig, filename):
    """Standard Academic Figure Skill export: vector PDF + 300dpi PNG preview."""
    fig.savefig(f"{filename}.pdf", bbox_inches="tight", dpi=300)
    fig.savefig(f"{filename}.png", bbox_inches="tight", dpi=300)
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from common import OUT,config,save_json
FIG=OUT/'figures';FIG.mkdir(exist_ok=True)
GEN=OUT/'generated';GEN.mkdir(exist_ok=True)
annual=pd.read_csv(OUT/'annual_footprints.csv');summary=pd.read_csv(OUT/'coverage_summary.csv')
fig,axes=plt.subplots(1,2,figsize=(183/25.4,78/25.4),sharex=True,sharey=True)
for ax,outcome,title,color,letter in zip(axes,['A','C'],['Simultaneous exceedance area','Largest connected exceedance area'],[CATEGORICAL[0],CATEGORICAL[1]],['a','b']):
 for year,d in annual[annual.outcome==outcome].groupby('year'):
  d=d.sort_values('threshold_c');ax.plot(d.threshold_c,100*d.difference,color=GREY,alpha=.22,lw=.45,zorder=1)
 d=summary[summary.outcome==outcome].sort_values('threshold_c')
 ax.fill_between(d.threshold_c,100*d.lower,100*d.upper,color=color,alpha=.2,zorder=2)
 ax.plot(d.threshold_c,100*d.difference,color=color,lw=1.4,zorder=3)
 ax.axhline(0,color=BLACK,lw=.6,linestyle='--',alpha=.5)
 ax.set_title(title,loc='left',pad=10);ax.text(-.14,1.08,letter,transform=ax.transAxes,fontweight='bold',fontsize=10)
 ax.set_xlim(d.threshold_c.min(),d.threshold_c.max());ax.set_xticks([10,15,20,25,28.5]);ax.set_xlabel('Wet-bulb threshold (°C)')
axes[0].set_ylabel('High − middle difference (percentage points)')
fig.subplots_adjust(left=.10,right=.98,bottom=.20,top=.84,wspace=.18)
save_cns_figure(fig,FIG/'application_coverage');plt.close(fig)

paired=pd.read_csv(OUT/'threshold_paired_improvement.csv')
fig,axes=plt.subplots(2,2,figsize=(183/25.4,127/25.4),sharex=True)
for i,endpoint in enumerate(['A','C']):
 for j,subset in enumerate(['all','high']):
  ax=axes[i,j]
  for ref,color,style,name in [('simple',CATEGORICAL[0],'-','Simple − graph'),('variogram',CATEGORICAL[3],'--','Variogram − graph')]:
   d=paired[(paired.outcome==endpoint)&(paired.subset==subset)&(paired.reference==ref)].sort_values('threshold_c')
   ax.fill_between(d.threshold_c,d.lower_pp,d.upper_pp,color=color,alpha=.17)
   ax.plot(d.threshold_c,d.improvement_pp,color=color,linestyle=style,lw=1.2,label=name)
  ax.axhline(0,color=BLACK,lw=.6)
  ax.set_title(('Exceedance area' if i==0 else 'Largest connected area')+' · '+('all days' if j==0 else 'high days'),loc='left')
  ax.text(-.14,1.08,'abcd'[i*2+j],transform=ax.transAxes,fontweight='bold',fontsize=10)
  ax.set_xlim(9.5,28.5);ax.set_xticks([10,15,20,25,28.5])
  if j==0:ax.set_ylabel('MAE reduction (percentage points)')
  if i==1:ax.set_xlabel('Wet-bulb threshold (°C)')
handles,labels=axes[0,0].get_legend_handles_labels();fig.legend(handles,labels,loc='lower center',bbox_to_anchor=(.5,-.015),ncol=2)
fig.subplots_adjust(left=.11,right=.98,bottom=.16,top=.94,hspace=.40,wspace=.25)
save_cns_figure(fig,FIG/'application_increment');plt.close(fig)

stats=pd.read_csv(OUT/'comparison_summary.csv')
methods=[('mean','Mean and season'),('simple','Strong simple baseline'),('variogram','Simple + variogram'),('graph','Simple + graph'),('interpolation','Spatial interpolation')]
rows=[]
for method,label in methods:
 d=stats[stats.method==method].set_index(['outcome','subset'])
 vals=[d.loc[(o,g),'mae_pp'] for o in ['A','C','A_new344'] for g in ['all','high']]
 rows.append(label+' & '+' & '.join(f'{v:.3f}' for v in vals)+r' \\')
tex=r'''\begin{table}[!t]
\centering\small
\caption{Held-out mean absolute error in area percentage points. Errors are
averaged over all 77 thresholds (9.5--28.5$^\circ$C at 0.25$^\circ$C steps),
then equally over months and 33 evaluation summers. All-day errors use 3,036
June--August daily fields in 1991--2025 excluding 2015 and 2022; high-day errors
use the 792 original upper-quartile days. The last two columns use the 344
additional sites. All inputs come from the original 121 sites.}
\label{tab:application-comparison}
\begin{tabular}{@{}lrrrrrr@{}}
\toprule
& \multicolumn{2}{c}{Exceedance area} & \multicolumn{2}{c}{Largest component} & \multicolumn{2}{c}{Additional sites}\\
\cmidrule(lr){2-3}\cmidrule(lr){4-5}\cmidrule(lr){6-7}
Method & All & High & All & High & All & High\\
\midrule
'''+ '\n'.join(rows)+r'''
\bottomrule
\end{tabular}
\end{table}
'''
(GEN/'application_comparison.tex').write_text(tex)
# Full quantitative table of paired increments (positive means smaller graph MAE).
pairs=pd.read_csv(OUT/'paired_improvement_summary.csv');rows=[]
for r in pairs.itertuples():
 label={'A':'Area','C':'Connected area','A_new344':'Additional sites'}[r.outcome]
 rows.append(f'{label} & {r.subset} & {r.reference} & ${r.improvement_pp:.4f}$ & $[{r.lower_pp:.4f},{r.upper_pp:.4f}]$ & {r.improved_years}/33'+r' \\')
(GEN/'application_paired.tex').write_text(r'''\begin{table}[H]
\centering\small
\caption{Paired reduction in held-out mean absolute error after adding graph
information, in percentage points. Positive values favour the graph model.
Intervals use five-year calendar-block resampling of the annual paired errors,
conditional on the fitted cross-validation predictions.}
\label{tab:application-paired}
\begin{tabular}{@{}lllrrr@{}}
\toprule
Outcome & Days & Reference & Reduction & 95\% interval & Improved summers\\
\midrule
'''+ '\n'.join(rows)+r'''
\bottomrule
\end{tabular}
\end{table}
''')
# Consolidate the existing event-definition analyses without rerunning them.
event=pd.read_csv(OUT.parents[0]/'results/event_partition/summary.csv')
labels={'Quartiles':'Quartiles','Tertiles':'Tertiles','20--80':'20th/80th','30--70':'30th/70th','25--70':'25th/70th','25--80':'25th/80th'}
rows=[]
for key,name in labels.items():
 d=event[event.partition==key].set_index('scale')
 rows.append(name+' & '+ ' & '.join(f'${100*d.loc[k,"estimate"]:.2f}$' for k in [0,1,5])+r' \\')
cont=pd.read_csv(OUT.parents[0]/'results/extension_continuous_profile_summary.csv')
print('Figures/tables generated. Continuous summary columns:',cont.columns.tolist())
save_json('figure_statistics.json',dict(archetype='quantitative_grid',primary_panels=2,supplement_panels=4,all_thresholds=77,annual_replicates=33,centre='equal-month and equal-summer mean',interval='pointwise percentile bands from 1999 five-year calendar moving-block resamples',multiplicity='no simultaneous significance claims',prediction_intervals='conditional on fitted nested blocked predictions',data_sources=['coverage_summary.csv','annual_footprints.csv','threshold_paired_improvement.csv','comparison_summary.csv','paired_improvement_summary.csv'],font='Arial',vector_pdf=True,png_dpi=300,source_rows_discarded='none; predeclared evaluation excludes development years',asset='LineTrend/plot_trend.py inspected; cross-type parameter inheritance because cumulative counts differ from threshold contrasts'))
