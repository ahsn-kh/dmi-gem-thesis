#!/usr/bin/env python3
"""Generate standardized thesis figures from compact frozen source tables.

Run:
    python code/generate_thesis_figures.py

Outputs PNG (600 dpi), vector PDF, and SVG into figures/{png,pdf,svg}.
"""
from pathlib import Path
import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch, FancyBboxPatch

BASE = Path(__file__).resolve().parents[1]
DATA = BASE / "source_tables"
PNG = BASE / "figures" / "png"
PDF = BASE / "figures" / "pdf"
SVG = BASE / "figures" / "svg"
for d in (PNG, PDF, SVG): d.mkdir(parents=True, exist_ok=True)

BLUE = "#2F6B8A"
ORANGE = "#C56B2D"
DARK = "#252525"
MID = "#737373"
LIGHT = "#D9D9D9"
PALE_BLUE = "#DCEAF1"
PALE_ORANGE = "#F3E1D2"
PALE_GRAY = "#F2F2F2"
GRID = "#E5E5E5"
FONT = "Liberation Sans"
THS = np.array([5,10,15,20,25,30,35,40,45])

plt.rcParams.update({
    "font.family": FONT,
    "font.size": 9,
    "axes.titlesize": 10,
    "axes.titleweight": "bold",
    "axes.labelsize": 9,
    "xtick.labelsize": 8,
    "ytick.labelsize": 8,
    "legend.fontsize": 8,
    "figure.titlesize": 11,
    "figure.titleweight": "bold",
    "axes.spines.top": False,
    "axes.spines.right": False,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "savefig.bbox": "tight",
    "savefig.pad_inches": 0.08,
})

def clean_ax(ax, grid="y"):
    ax.tick_params(length=3, width=.8, color=MID)
    ax.spines["left"].set_color(MID); ax.spines["bottom"].set_color(MID)
    if grid:
        ax.grid(axis=grid, color=GRID, lw=.7, zorder=0)
    ax.set_axisbelow(True)

def panel(ax, letter, x=-0.12, y=1.08):
    ax.text(x, y, letter, transform=ax.transAxes, fontsize=12, fontweight="bold", va="top", ha="left")

def save(fig, stem):
    fig.savefig(PNG / f"{stem}.png", dpi=600)
    fig.savefig(PDF / f"{stem}.pdf")
    fig.savefig(SVG / f"{stem}.svg")
    plt.close(fig)

def tissue_label(x):
    return {"TibialisAnt":"Tibialis anterior"}.get(x,x)

# ---------------- Figure 2.1 ----------------
def fig_2_1():
    fig, ax = plt.subplots(figsize=(8.2, 3.25))
    ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.axis("off")
    items = [
        (0.04, "1", "Transcriptomic data", "GSE17576\nGSE305719"),
        (0.28, "2", "Expression evidence", "Cutoff → gene\nconfidence → GPR"),
        (0.52, "3", "CORDA reconstruction", "Context-specific\nmuscle model"),
        (0.76, "4", "DMI-GEM analysis", "Fixed DMI targets,\nATP demand, FVA + sampling"),
    ]
    for i,(x,n,title,body) in enumerate(items):
        edge = ORANGE if i==1 else BLUE
        face = PALE_ORANGE if i==1 else PALE_BLUE
        box=FancyBboxPatch((x,.33),.19,.38,boxstyle="round,pad=0.012,rounding_size=0.018",fc=face,ec=edge,lw=1.2)
        ax.add_patch(box)
        ax.text(x+.095,.78,n,ha="center",va="center",fontweight="bold",color=edge,fontsize=9,
                bbox=dict(boxstyle="circle,pad=.25",fc="white",ec=edge,lw=1))
        ax.text(x+.095,.64,title,ha="center",va="center",fontweight="bold",fontsize=9)
        ax.text(x+.095,.47,body,ha="center",va="center",fontsize=8,color=DARK,linespacing=1.2)
        if i<3:
            ax.annotate("",xy=(x+.235,.52),xytext=(x+.205,.52),arrowprops=dict(arrowstyle="->",color=MID,lw=1.2))
    ax.plot([.28,.95],[.19,.19],color=MID,lw=1)
    ax.plot([.28,.28],[.19,.25],color=MID,lw=1); ax.plot([.95,.95],[.19,.25],color=MID,lw=1)
    ax.text(.615,.10,"Only transcriptomic source and expression cutoff varied; downstream DMI-GEM settings were fixed.",ha="center",va="center",fontsize=8,color=DARK)
    ax.set_title("Study workflow", pad=8)
    save(fig,"Figure_2_1_revised")

# ---------------- Figure 2.2: keep former C only ----------------
def fig_2_2():
    d=pd.read_csv(DATA/'gse17576_per_gene_mean_expression.csv')
    c=pd.read_csv(DATA/'gse17576_gene_confidence.csv').set_index('threshold')
    vals=d.per_gene_mean_expression.to_numpy()
    fig,ax=plt.subplots(figsize=(7.5,4.25))
    ax.hist(vals,bins=55,color=BLUE,alpha=.88,edgecolor="white",linewidth=.3,zorder=2)
    specs=[('P5',MID,'--','right',-4,1.105),('P25',ORANGE,'-','left',4,1.035),('P45',MID,'--','right',-4,1.105),('P50',MID,':','left',4,1.035)]
    ymax_data=ax.get_ylim()[1]
    ax.set_ylim(0,ymax_data*1.16)  # reserve a two-row label band above the histogram
    for lab,col,ls,ha,dx,yfac in specs:
        x=float(c.loc[lab,'expression_cutoff'])
        ax.axvline(x,color=col,ls=ls,lw=1.5,zorder=3)
        ax.annotate(f"{lab}\n{x:.2f}",(x,ymax_data*yfac),xytext=(dx,0),textcoords='offset points',
                    ha=ha,va='center',fontsize=8,color=col,fontweight='bold' if lab=='P25' else 'normal')
    ax.set_xlabel("Per-gene mean expression across the three conditions (log$_2$)")
    ax.set_ylabel("Genes")
    ax.set_title("GSE17576 expression cutoffs used for reconstruction")
    clean_ax(ax,'y')
    ax.text(.99,.92,"1,619 of 1,865 iMM1865 genes\nmatched the processed dataset",transform=ax.transAxes,ha='right',va='top',fontsize=8,color=DARK)
    save(fig,"Figure_2_2_revised_C_only")

# ---------------- Figure 2.3: keep former A only ----------------
def fig_2_3():
    d=pd.read_csv(DATA/'gse305719_library_sizes.csv')
    order_codes=['gas','sol','bat','ewat','iwat','liver','panc','gut','hypo','pit','neo','hipp','olf']
    labels=['Gastrocnemius','Soleus','Brown adipose tissue','Epididymal WAT','Inguinal WAT','Liver','Pancreas','Small intestine','Hypothalamus','Pituitary','Neocortex','Hippocampus','Olfactory bulb']
    data=[d.loc[d.organ_code==c,'library_size_million'].dropna().to_numpy() for c in order_codes]
    fig,ax=plt.subplots(figsize=(7.6,5.35))
    pos=np.arange(len(data))+1
    bp=ax.boxplot(data,vert=False,positions=pos,widths=.55,patch_artist=True,showfliers=False,
                  medianprops=dict(color=DARK,lw=1.1),whiskerprops=dict(color=MID,lw=.8),capprops=dict(color=MID,lw=.8))
    for i,b in enumerate(bp['boxes']):
        if i<2: b.set(facecolor=BLUE,alpha=.85,edgecolor=BLUE)
        else: b.set(facecolor=LIGHT,alpha=.75,edgecolor=MID)
    rng=np.random.default_rng(42)
    for i,arr in enumerate(data):
        y=np.full(arr.size,pos[i])+rng.normal(0,.045,arr.size)
        ax.scatter(arr,y,s=7,color=BLUE if i<2 else MID,alpha=.28,zorder=3,linewidths=0)
    ax.set_yticks(pos,labels); ax.invert_yaxis()
    for i,t in enumerate(ax.get_yticklabels()):
        if i<2: t.set_color(BLUE); t.set_fontweight('bold')
    ax.set_xlabel("Raw library size (million reads)"); ax.set_title("GSE305719 raw library sizes across organs")
    clean_ax(ax,'x')
    ax.text(.99,.03,"Gastrocnemius and Soleus were used for the reconstruction subset.",transform=ax.transAxes,ha='right',va='bottom',fontsize=8,color=DARK)
    save(fig,"Figure_2_3_revised_A_only")

# ---------------- Figure 2.4 ----------------
def fig_2_4():
    d=pd.read_csv(DATA/'dmi_targets.csv')
    order=['Gastrocnemius','Soleus','TibialisAnt']; labels=['Gastrocnemius','Soleus','Tibialis\nanterior']
    fig=plt.figure(figsize=(8.4,6.3)); gs=fig.add_gridspec(2,2,height_ratios=[1,1.05],hspace=.55,wspace=.28)
    ax1=fig.add_subplot(gs[0,0]);ax2=fig.add_subplot(gs[0,1]);ax3=fig.add_subplot(gs[1,:]);
    offsets={45:-.19,48:-.11,55:-.03,56:.05,49:.13,51:.21}
    markers={'dia':'o','norm':'s'}; colors={'dia':ORANGE,'norm':BLUE}; names={'dia':'HFD','norm':'Control'}
    for ax,col,title,ylabel in [(ax1,'HEX1_target','HEX1 targets from MR$_{glc}$','Target flux (mmol gDW$^{-1}$ h$^{-1}$)'),(ax2,'PDHm_target','PDHm targets from V$_{ox}$','Target flux (mmol gDW$^{-1}$ h$^{-1}$)')]:
        for xi,tis in enumerate(order):
            z=d[d.tissue==tis]
            for _,r in z.iterrows():
                x=xi+offsets[int(r.mouseID)]
                ax.scatter(x,r[col],s=38,marker=markers[r.group],facecolor=colors[r.group],edgecolor='white',linewidth=.6,zorder=4)
            for grp in ['dia','norm']:
                zz=z[z.group==grp][col]
                if len(zz):
                    x0=xi+(-.08 if grp=='dia' else .17)
                    ax.plot([x0-.06,x0+.06],[zz.median(),zz.median()],color=colors[grp],lw=2.0,zorder=5)
        ax.set_xticks(range(3),labels);ax.set_ylabel(ylabel);ax.set_title(title);clean_ax(ax,'y')
    panel(ax1,'A');panel(ax2,'B')
    handles=[Line2D([0],[0],marker='o',ls='',markerfacecolor=ORANGE,markeredgecolor='white',markersize=7,label='HFD'),Line2D([0],[0],marker='s',ls='',markerfacecolor=BLUE,markeredgecolor='white',markersize=7,label='Control')]
    fig.legend(handles=handles,loc='center',bbox_to_anchor=(.50,.505),ncol=2,frameon=False,columnspacing=1.4,handletextpad=.5)
    ax3.set_xlim(0,1);ax3.set_ylim(0,1);ax3.axis('off');panel(ax3,'C',x=-.025,y=1.06);ax3.set_title("Constraint mapping",pad=4)
    boxes=[(.03,.56,.17,.24,PALE_BLUE,BLUE,"DMI rates","MR$_{glc}$, V$_{ox}$"),(.27,.56,.17,.24,PALE_GRAY,MID,"Global scaling","k = 0.0301891"),(.51,.56,.18,.24,PALE_BLUE,BLUE,"Reaction targets","MR$_{glc}$ → HEX1\nV$_{ox}$ → PDHm"),(.77,.56,.20,.24,PALE_ORANGE,ORANGE,"Soft intervals","HEX1 ±20%, weight 1.00\nPDHm ±50%, weight 0.75")]
    for x,y,w,h,fc,ec,head,body in boxes:
        p=FancyBboxPatch((x,y),w,h,boxstyle='round,pad=.012,rounding_size=.014',fc=fc,ec=ec,lw=1.1);ax3.add_patch(p)
        ax3.text(x+w/2,y+h*.68,head,ha='center',va='center',fontweight='bold',fontsize=8.5);ax3.text(x+w/2,y+h*.33,body,ha='center',va='center',fontsize=8)
    for x in [.20,.44,.69]: ax3.annotate('',xy=(x+.065,.68),xytext=(x+.015,.68),arrowprops=dict(arrowstyle='->',lw=1.1,color=MID))
    ax3.text(.5,.26,"All reconstructions: fixed ATP demand = 3.5 mmol gDW$^{-1}$ h$^{-1}$ and identical nutrient/reaction-bound policy",ha='center',va='center',fontsize=8.5,color=DARK)
    fig.suptitle("DMI-derived targets and model constraints",y=.995)
    save(fig,"Figure_2_4_revised")

# ---------------- Figure 3.1 ----------------
def fig_3_1():
    conf=pd.read_csv(DATA/'gse17576_gene_confidence.csv'); conf=conf[conf['threshold'].isin([f'P{x}' for x in THS])].copy(); st=pd.read_csv(DATA/'gse17576_structure_atp.csv');chg=pd.read_csv(DATA/'gse17576_adjacent_changes.csv')
    fig,axs=plt.subplots(2,2,figsize=(8.5,6.5));ax=axs[0,0]
    x=np.arange(len(conf)); bottoms=np.zeros(len(conf)); cols=[BLUE,'#9FBFD0',ORANGE,LIGHT]; labs=['HC','MC','LC','NC']
    for lab,col in zip(labs,cols):
        ax.bar(x,conf[lab],bottom=bottoms,color=col,width=.72,label=lab,edgecolor='white',linewidth=.25);bottoms+=conf[lab].to_numpy()
    ax.set_xticks(x,[f"P{t}" for t in THS]);ax.set_ylabel("Model genes");ax.set_ylim(0,float(bottoms.max())*1.13);ax.set_title("Gene confidence classes");clean_ax(ax,'y');panel(ax,'A');ax.legend(ncol=4,loc='upper center',bbox_to_anchor=(.5,.995),frameon=False,columnspacing=1.1,handlelength=1.1)
    ax=axs[0,1];ax.plot(st.threshold,st.reaction_count,marker='o',color=BLUE,lw=1.8,ms=4.5);ax.set_xticks(THS,[f'P{x}' for x in THS]);ax.set_ylabel("Reactions");ax.set_title("Reaction count");clean_ax(ax,'y');panel(ax,'B');
    for th in [5,25,45]:
        r=st[st.threshold==th].iloc[0];ax.annotate(f"{int(r.reaction_count):,}",(th,r.reaction_count),xytext=(0,7),textcoords='offset points',ha='center',fontsize=8,color=DARK)
    ax=axs[1,0];y=np.arange(len(chg));ax.barh(y,-chg.lost,color=BLUE,alpha=.85,label='Lost');ax.barh(y,chg.gained,color=ORANGE,alpha=.9,label='Gained');ax.axvline(0,color=DARK,lw=.8);ax.set_yticks(y,chg.transition);ax.set_ylim(len(chg)-.5,-1.05);ax.set_xlabel("Reaction count change");ax.set_title("Reaction gains and losses");clean_ax(ax,'x');panel(ax,'C');ax.legend(ncol=2,frameon=False,loc='upper center',bbox_to_anchor=(.5,.995))
    ax=axs[1,1];xx=np.arange(len(chg));ax.plot(xx,chg.jaccard,marker='o',color=BLUE,lw=1.8,ms=4.5);imin=chg.jaccard.idxmin();i=list(chg.index).index(imin);ax.scatter(i,chg.loc[imin,'jaccard'],s=42,color=ORANGE,zorder=5);ax.annotate(f"lowest: {chg.loc[imin,'jaccard']:.3f}",(i,chg.loc[imin,'jaccard']),xytext=(8,-16),textcoords='offset points',fontsize=8,color=ORANGE);ax.set_xticks(xx,chg.transition,rotation=35,ha='right');ax.set_ylim(.91,.97);ax.set_ylabel("Jaccard index");ax.set_title("Similarity of adjacent reaction sets");clean_ax(ax,'y');panel(ax,'D')
    fig.suptitle("Effect of expression cutoff on GSE17576 model structure",y=.995);fig.subplots_adjust(hspace=.48,wspace=.32)
    save(fig,"Figure_3_1_revised")

# ---------------- Figure 3.2 ----------------
def fig_3_2():
    st=pd.read_csv(DATA/'gse17576_structure_atp.csv');su=pd.read_csv(DATA/'p35_succact_intervention.csv')
    fig,axs=plt.subplots(1,2,figsize=(8.4,3.75),gridspec_kw={'wspace':.34})
    ax=axs[0];ax.plot(st.threshold,st.max_ATP,color=BLUE,marker='o',lw=1.8,ms=4.5);q=st[st.threshold==35].iloc[0];ax.scatter(35,q.max_ATP,color=ORANGE,s=45,zorder=4);ax.axhline(3.5,color=DARK,ls='--',lw=1,label='Fixed ATP demand');ax.set_xticks(THS,[f'P{x}' for x in THS]);ax.set_ylabel("Maximum ATP flux (mmol gDW$^{-1}$ h$^{-1}$)");ax.set_title("ATP capacity across expression cutoffs");clean_ax(ax,'y');panel(ax,'A');ax.legend(frameon=False,loc='upper right')
    ax=axs[1];x=np.arange(len(su));cols=[BLUE,ORANGE,DARK,BLUE];ax.scatter(x,su.max_ATP,s=55,c=cols,zorder=4);ax.vlines(x,3.5,su.max_ATP,colors=cols,lw=1.3,alpha=.75);ax.axhline(3.5,color=DARK,ls='--',lw=1);ax.set_xticks(x,su.condition);ax.set_ylabel("Maximum ATP flux (mmol gDW$^{-1}$ h$^{-1}$)");ax.set_title("Effect of SUCCACT inactivation on P35");clean_ax(ax,'y');panel(ax,'B')
    for i,v in enumerate(su.max_ATP):
        if v>4.8:
            ax.text(i,v-.12,f"{v:.3f}",ha='center',va='top',fontsize=8)
        else:
            ax.text(i,v+.10,f"{v:.3f}",ha='center',va='bottom',fontsize=8)
    fig.suptitle("ATP capacity and the P35 SUCCACT intervention",y=.995)
    save(fig,"Figure_3_2_revised")

# ---------------- Figure 3.3 ----------------
def fig_3_3():
    d=pd.read_csv(DATA/'anchor_flux_summary.csv');g=d[d.source=='GSE17576'].sort_values('threshold');f=pd.read_csv(DATA/'gse17576_pdhm_fva_summary.csv')
    fig=plt.figure(figsize=(8.5,6.25));gs=fig.add_gridspec(2,2,height_ratios=[1,1.05],hspace=.55,wspace=.30);a=fig.add_subplot(gs[0,0]);b=fig.add_subplot(gs[0,1]);c=fig.add_subplot(gs[1,:])
    for ax,prefix,col,title,region,ylim in [(a,'HEX1',BLUE,'HEX1',(0.8,1.2),(0.70,1.25)),(b,'PDHm',ORANGE,'PDHm',(0.5,1.5),(0.30,1.75))]:
        ax.axhspan(region[0],region[1],color=LIGHT,alpha=.35,zorder=0);ax.axhline(1,color=DARK,ls='--',lw=1)
        ax.fill_between(g.threshold,g[f'{prefix}_q25_rel'],g[f'{prefix}_q75_rel'],color=col,alpha=.20,lw=0)
        ax.plot(g.threshold,g[f'{prefix}_median_rel'],color=col,marker='o',lw=1.8,ms=4.5)
        ax.set_xticks(THS,[f'P{x}' for x in THS]);ax.set_ylim(*ylim);ax.set_ylabel("Sampled median / DMI target");ax.set_title(title);clean_ax(ax,'y')
    panel(a,'A');panel(b,'B')
    c.axhline(1,color=DARK,ls='--',lw=1);c.axhspan(.5,1.5,color=LIGHT,alpha=.25,zorder=0)
    c.vlines(f.threshold,f.fva_min_rel,f.fva_max_rel,color=BLUE,lw=2.0,alpha=.75,label='PDHm FVA range');c.scatter(g.threshold,g.PDHm_median_rel,color=ORANGE,s=36,zorder=3,label='Sampled median')
    c.set_xticks(THS,[f'P{x}' for x in THS]);c.set_ylim(.3,1.75);c.set_ylabel("PDHm flux / DMI target");c.set_title("PDHm FVA range and sampled median");clean_ax(c,'y');panel(c,'C',x=-.055);c.legend(frameon=False,ncol=2,loc='upper center',bbox_to_anchor=(.5,-.20))
    fig.suptitle("HEX1 and PDHm across GSE17576 expression cutoffs",y=.995)
    save(fig,"Figure_3_3_revised")

# ---------------- Figure 3.4 ----------------
def fig_3_4():
    d=pd.read_csv(DATA/'pathway_kendall_w.csv');order=['TCA cycle','Redox shuttles','Pyruvate oxidation','Oxidative phosphorylation','Lactate metabolism','Ketone metabolism','Glycolysis','Glutamate/aKG metabolism','Glucose uptake/phosphorylation','Fatty acid uptake & beta-oxidation','Creatine phosphate buffering','Glutamine metabolism']
    fig,axs=plt.subplots(1,2,figsize=(8.6,5.35),sharey=True,gridspec_kw={'wspace':.10})
    for ax,src,col,letter in zip(axs,['GSE17576','GSE305719'],[BLUE,ORANGE],['A','B']):
        z=d[d.source==src].set_index('pathway').reindex(order).reset_index();y=np.arange(len(z))
        ax.scatter(z.W,y,s=38,facecolors=[col if s else 'white' for s in z.significant_q_lt_0_05],edgecolors=col,lw=1.2,zorder=3)
        for yi,r in z.iterrows():
            if r.significant_q_lt_0_05:
                ax.text(r.W+.025,yi,f"W={r.W:.3f}",va='center',fontsize=7.5,color=col)
        ax.set_xlim(0,.78);ax.set_xlabel("Kendall's W");ax.set_title(src);clean_ax(ax,'x');panel(ax,letter)
    axs[0].set_yticks(np.arange(len(order)),order);axs[0].invert_yaxis()
    legend=[Line2D([0],[0],marker='o',ls='',mfc=DARK,mec=DARK,markersize=6,label='q < 0.05'),Line2D([0],[0],marker='o',ls='',mfc='white',mec=DARK,markersize=6,label='q ≥ 0.05')]
    fig.legend(handles=legend,loc='upper center',bbox_to_anchor=(.5,.94),ncol=2,frameon=False)
    fig.suptitle("Effect of expression cutoff on pathway flux variability",y=.995)
    save(fig,"Figure_3_4_revised")

# ---------------- Figure 3.5 ----------------
def fig_3_5():
    d=pd.read_csv(DATA/'pathway_hfd_effect_ranges.csv');paths=['Fatty acid uptake & beta-oxidation','Glutamate/aKG metabolism','Glycolysis','Lactate metabolism','Oxidative phosphorylation','Pyruvate oxidation','TCA cycle'];tissues=['Gastrocnemius','Soleus','TibialisAnt']
    fig,axs=plt.subplots(3,1,figsize=(8.5,7.25),sharex=True,gridspec_kw={'hspace':.28});offs={'GSE17576':-.13,'GSE305719':.13};cols={'GSE17576':BLUE,'GSE305719':ORANGE}
    for ax,tis,letter in zip(axs,tissues,['A','B','C']):
        y=np.arange(len(paths));ax.axvline(0,color=DARK,lw=1)
        for src in ['GSE17576','GSE305719']:
            z=d[(d.source==src)&(d.tissue==tis)].set_index('pathway').reindex(paths).reset_index();yy=y+offs[src]
            ax.hlines(yy,z['min'],z['max'],color=cols[src],lw=1.4,alpha=.72);ax.scatter(z['median'],yy,color=cols[src],s=28,zorder=3)
        ax.set_yticks(y,paths);ax.invert_yaxis();ax.set_title(tissue_label(tis),loc='left');clean_ax(ax,'x');panel(ax,letter,x=-.13,y=1.04)
    axs[-1].set_xlabel("log$_2$(HFD median / Control median)");axs[-1].set_xlim(-1.2,1.9)
    handles=[Line2D([0],[0],marker='o',color=BLUE,lw=1.4,label='GSE17576'),Line2D([0],[0],marker='o',color=ORANGE,lw=1.4,label='GSE305719')]
    fig.legend(handles=handles,ncol=2,loc='upper center',bbox_to_anchor=(.55,.96),frameon=False)
    fig.suptitle("HFD–Control pathway effects across expression cutoffs",y=.995)
    save(fig,"Figure_3_5_revised")

# ---------------- Figure 3.6 ----------------
def fig_3_6():
    s175=pd.read_csv(DATA/'gse17576_structure_atp.csv');s305=pd.read_csv(DATA/'gse305719_structure_atp.csv');cross=pd.read_csv(DATA/'cross_source_structure.csv');anc=pd.read_csv(DATA/'anchor_flux_summary.csv');sat=pd.read_csv(DATA/'pdhm_vs_atp_saturation.csv')
    fig,axs=plt.subplots(2,3,figsize=(9.0,7.0));
    # A reaction count
    ax=axs[0,0];ax.plot(s175.threshold,s175.reaction_count,color=BLUE,marker='o',ms=4,lw=1.7,label='GSE17576');ax.plot(s305.threshold,s305.reaction_count,color=ORANGE,marker='o',ms=4,lw=1.7,label='GSE305719');ax.set_title('Reaction count');ax.set_ylabel('Reactions');ax.set_xticks(THS,[f'P{x}' for x in THS]);clean_ax(ax,'y');panel(ax,'A')
    # B Jaccard
    ax=axs[0,1];ax.plot(cross.threshold_num,cross.jaccard,color=DARK,marker='o',ms=4,lw=1.7);ax.set_title('Reaction-set similarity');ax.set_ylabel('Cross-source Jaccard index');ax.set_xticks(THS,[f'P{x}' for x in THS]);ax.set_ylim(.86,.95);clean_ax(ax,'y');panel(ax,'B')
    # C ATP
    ax=axs[0,2];ax.plot(s175.threshold,s175.max_ATP,color=BLUE,marker='o',ms=4,lw=1.7);ax.plot(s305.threshold,s305.max_ATP,color=ORANGE,marker='o',ms=4,lw=1.7);ax.axhline(3.5,color=DARK,ls='--',lw=1);ax.set_title('Maximum ATP capacity');ax.set_ylabel('mmol gDW$^{-1}$ h$^{-1}$');ax.set_xticks(THS,[f'P{x}' for x in THS]);clean_ax(ax,'y');panel(ax,'C')
    # D/E anchors
    for ax,prefix,title,letter,ylim in [(axs[1,0],'HEX1','HEX1 sampled flux','D',(.72,.88)),(axs[1,1],'PDHm','PDHm sampled flux','E',(.35,.88))]:
        ax.axhline(1,color=DARK,ls='--',lw=1)
        for src,col in [('GSE17576',BLUE),('GSE305719',ORANGE)]:
            z=anc[anc.source==src].sort_values('threshold');ax.fill_between(z.threshold,z[f'{prefix}_q25_rel'],z[f'{prefix}_q75_rel'],color=col,alpha=.13,lw=0);ax.plot(z.threshold,z[f'{prefix}_median_rel'],color=col,marker='o',ms=3.8,lw=1.6)
        ax.set_title(title);ax.set_ylabel('Median / DMI target');ax.set_xticks(THS,[f'P{x}' for x in THS]);ax.set_ylim(*ylim);clean_ax(ax,'y');panel(ax,letter)
    # F saturation relation
    ax=axs[1,2]
    for src,col in [('GSE17576',BLUE),('GSE305719',ORANGE)]:
        z=sat[sat.source==src];ax.scatter(z.ATP_demand_fraction,z.PDHm_median_rel,color=col,s=35,alpha=.9,label=src)
    ax.set_xlabel('ATP demand / maximum ATP capacity');ax.set_ylabel('PDHm median / DMI target');ax.set_title('PDHm and ATP-demand fraction');clean_ax(ax,'both');panel(ax,'F',x=-.17);ax.text(.03,.96,'Pearson r = 0.940\nSpearman ρ = 0.670',transform=ax.transAxes,ha='left',va='top',fontsize=8,color=DARK)
    fig.legend(handles=[Line2D([0],[0],color=BLUE,marker='o',label='GSE17576'),Line2D([0],[0],color=ORANGE,marker='o',label='GSE305719')],ncol=2,loc='upper center',bbox_to_anchor=(.5,.965),frameon=False)
    fig.suptitle("Effect of transcriptomic dataset on model structure and DMI-constrained fluxes",y=.995);fig.subplots_adjust(hspace=.42,wspace=.34)
    save(fig,"Figure_3_6_revised")

# ---------------- Figure 3.7: combine old C/D ----------------
def fig_3_7():
    a=pd.read_csv(DATA/'cross_source_pathway_agreement.csv');h=pd.read_csv(DATA/'hex1_hfd_control_differences.csv');td=pd.read_csv(DATA/'hex1_target_hfd_control_differences.csv').set_index('tissue')
    fig,axs=plt.subplots(1,2,figsize=(8.6,4.2),gridspec_kw={'width_ratios':[.8,1.35],'wspace':.34})
    ax=axs[0];ax.plot(a.threshold,100*a.direction_agreement_fraction,color=BLUE,marker='o',ms=4.5,lw=1.7);mean=100*a.direction_agreement_fraction.mean();ax.axhline(mean,color=DARK,ls='--',lw=1);ax.text(44,mean+2,f'mean {mean:.1f}%',ha='right',va='bottom',fontsize=8,color=DARK);ax.set_xticks(THS,[f'P{x}' for x in THS]);ax.set_ylabel('Same HFD–Control direction (%)');ax.set_ylim(30,90);ax.set_title('Cross-source pathway agreement');clean_ax(ax,'y');panel(ax,'A')
    ax=axs[1]
    styles={'TibialisAnt':'-','Soleus':'--'}
    for src,col in [('GSE17576',BLUE),('GSE305719',ORANGE)]:
        for tis in ['TibialisAnt','Soleus']:
            z=h[(h.source==src)&(h.tissue==tis)].sort_values('threshold');ax.plot(z.threshold,z.hex1_hfd_minus_control,color=col,ls=styles[tis],marker='o' if tis=='TibialisAnt' else 's',ms=4,lw=1.6)
    for tis in ['TibialisAnt','Soleus']:
        y=float(td.loc[tis,'target_hfd_minus_control']);ax.axhline(y,color=MID,ls=styles[tis],lw=1,alpha=.8);ax.text(45.7,y,f'{tissue_label(tis)} DMI target',fontsize=7.5,color=MID,va='center',clip_on=False)
    ax.set_xlim(4,51);ax.set_xticks(THS,[f'P{x}' for x in THS]);ax.set_ylabel('HEX1 HFD − Control (mmol gDW$^{-1}$ h$^{-1}$)');ax.set_title('HEX1 HFD–Control difference');clean_ax(ax,'y');panel(ax,'B')
    source_legend=[Line2D([0],[0],color=BLUE,lw=2,label='GSE17576'),Line2D([0],[0],color=ORANGE,lw=2,label='GSE305719')]
    tissue_legend=[Line2D([0],[0],color=DARK,ls='-',marker='o',label='Tibialis anterior'),Line2D([0],[0],color=DARK,ls='--',marker='s',label='Soleus')]
    fig.legend(handles=source_legend+tissue_legend,frameon=False,loc='lower center',bbox_to_anchor=(.70,.015),ncol=2,columnspacing=1.6,handlelength=2.2)
    fig.subplots_adjust(bottom=.22)
    fig.suptitle('Diet-related results across reconstruction choices',y=.995)
    save(fig,"Figure_3_7_revised")

# Optional Methods QC: old Figure 2.2 A/B moved out of the main text.
def supplementary_g175_qc():
    expr_path=DATA/'normalized_expression.csv.xz'; ann_path=DATA/'sample_annotation.csv.xz'
    expr=pd.read_csv(expr_path,index_col=0) if expr_path.exists() else None
    ann=pd.read_csv(ann_path,index_col=0) if ann_path.exists() else None
    if expr is None or ann is None: return
    from sklearn.decomposition import PCA
    fig,axs=plt.subplots(1,2,figsize=(8.3,3.8),gridspec_kw={'wspace':.34})
    ax=axs[0];conds=['LFD_day0','LFD_day56','HFD_day56'];cols=[BLUE,'#7EA6BB',ORANGE];x=0
    for cond,col in zip(conds,cols):
        ids=ann.index[ann.condition==cond]
        vals=[expr[s].to_numpy() for s in ids]
        parts=ax.violinplot(vals,positions=np.arange(x,x+len(ids)),showmedians=True,showextrema=False,widths=.8)
        for b in parts['bodies']: b.set_facecolor(col);b.set_edgecolor(col);b.set_alpha(.65)
        parts['cmedians'].set_color(DARK);x+=len(ids)
    ax.axvline(9.5,color=LIGHT,lw=1);ax.axvline(19.5,color=LIGHT,lw=1);ax.set_xticks([4.5,14.5,24.5],['LFD day 0','LFD day 56','HFD day 56']);ax.set_ylabel('Expression (log$_2$)');ax.set_title('Sample expression distributions');clean_ax(ax,'y');panel(ax,'A')
    var=expr.var(axis=1).nlargest(1000).index;X=expr.loc[var].T;X=(X-X.mean())/X.std(ddof=0);pc=PCA(n_components=2).fit_transform(X)
    ax=axs[1]
    for cond,col in zip(conds,cols):
        mask=ann.loc[X.index,'condition'].eq(cond).to_numpy();ax.scatter(pc[mask,0],pc[mask,1],s=30,color=col,label=cond.replace('_',' '),alpha=.85)
    ax.set_xlabel('PC1');ax.set_ylabel('PC2');ax.set_title('PCA of 1,000 most variable genes');clean_ax(ax,'both');panel(ax,'B');ax.legend(frameon=False,loc='upper left',bbox_to_anchor=(1.01,1))
    fig.suptitle('GSE17576 preprocessing checks',y=.995);save(fig,'Supplementary_Methods_GSE17576_QC')

if __name__ == '__main__':
    fig_2_1();fig_2_2();fig_2_3();fig_2_4();fig_3_1();fig_3_2();fig_3_3();fig_3_4();fig_3_5();fig_3_6();fig_3_7();supplementary_g175_qc()
    print('Wrote figures to', BASE/'figures')
