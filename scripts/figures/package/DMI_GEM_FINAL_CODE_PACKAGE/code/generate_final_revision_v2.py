#!/usr/bin/env python3
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import FancyBboxPatch
from PIL import Image

BASE = Path('/mnt/data/thesis_figure_revision_20260902')
DATA = BASE / 'source_tables'
OUT = Path('/mnt/data/final_thesis_build/final_v2_figures')
EQOUT = Path('/mnt/data/final_thesis_build/final_v2_equations')
OUT.mkdir(parents=True, exist_ok=True)
EQOUT.mkdir(parents=True, exist_ok=True)

BLUE = '#2F6B8A'; ORANGE = '#C56B2D'; DARK = '#252525'; MID = '#737373'; LIGHT = '#D9D9D9'; GRID='#E5E5E5'
THS=np.array([5,10,15,20,25,30,35,40,45])
FONT='Liberation Sans'

plt.rcParams.update({
    'font.family': FONT, 'font.size': 10.5, 'axes.titlesize': 11.5, 'axes.titleweight':'bold',
    'axes.labelsize':10.5,'xtick.labelsize':9.5,'ytick.labelsize':9.5,'legend.fontsize':9.5,
    'figure.titlesize':12.5,'figure.titleweight':'bold','axes.spines.top':False,'axes.spines.right':False,
    'pdf.fonttype':42,'ps.fonttype':42,'savefig.bbox':'tight','savefig.pad_inches':0.10,
})

def clean_ax(ax, grid='y'):
    ax.tick_params(length=3,width=.8,color=MID)
    ax.spines['left'].set_color(MID); ax.spines['bottom'].set_color(MID)
    if grid=='both': ax.grid(axis='both',color=GRID,lw=.7,zorder=0)
    elif grid: ax.grid(axis=grid,color=GRID,lw=.7,zorder=0)
    ax.set_axisbelow(True)

def panel(ax, letter, x=-0.10, y=1.03):
    ax.text(x,y,letter,transform=ax.transAxes,fontsize=14,fontweight='bold',va='bottom',ha='left')

def tissue_label(x): return {'TibialisAnt':'Tibialis anterior'}.get(x,x)

def save(fig, stem):
    fig.savefig(OUT/f'{stem}.png',dpi=600)
    fig.savefig(OUT/f'{stem}.pdf')
    plt.close(fig)

# 2.1: deliberately simple, black-and-white, contained text

def fig_2_1():
    fig,ax=plt.subplots(figsize=(9.2,3.55))
    ax.set_xlim(0,1); ax.set_ylim(0,1); ax.axis('off')
    items=[
        (.025,'1','Transcriptomic data','GSE17576\nGSE305719'),
        (.275,'2','Expression evidence','Cutoff → gene\nconfidence → GPR'),
        (.525,'3','CORDA\nreconstruction','Context-specific\nmuscle model'),
        (.775,'4','DMI-GEM analysis','Fixed DMI targets,\nATP demand,\nFVA + sampling'),
    ]
    w=.19; h=.37
    for i,(x,n,title,body) in enumerate(items):
        box=FancyBboxPatch((x,.33),w,h,boxstyle='round,pad=0.012,rounding_size=0.017',fc='white',ec='black',lw=1.25)
        ax.add_patch(box)
        ax.text(x+w/2,.79,n,ha='center',va='center',fontweight='bold',fontsize=10,
                bbox=dict(boxstyle='circle,pad=.23',fc='white',ec='black',lw=1.1))
        ax.text(x+w/2,.61,title,ha='center',va='center',fontweight='bold',fontsize=9.2,linespacing=1.15)
        ax.text(x+w/2,.465,body,ha='center',va='center',fontsize=8.8,color=DARK,linespacing=1.15)
        if i<3:
            ax.annotate('',xy=(x+w+.055,.515),xytext=(x+w+.010,.515),arrowprops=dict(arrowstyle='->',color=DARK,lw=1.2))
    ax.plot([.275,.965],[.19,.19],color=MID,lw=1)
    ax.plot([.275,.275],[.19,.245],color=MID,lw=1); ax.plot([.965,.965],[.19,.245],color=MID,lw=1)
    ax.text(.62,.105,'Only transcriptomic source and expression cutoff varied; all downstream DMI-GEM settings were fixed.',
            ha='center',va='center',fontsize=8.8,color=DARK)
    ax.set_title('Study workflow',pad=8)
    save(fig,'Figure_2_1_FINAL')

# 3.1: A/B across top; C full width; D full width. Eliminates C overlap entirely.

def fig_3_1():
    conf=pd.read_csv(DATA/'gse17576_gene_confidence.csv')
    conf=conf[conf['threshold'].isin([f'P{x}' for x in THS])].copy()
    st=pd.read_csv(DATA/'gse17576_structure_atp.csv')
    chg=pd.read_csv(DATA/'gse17576_adjacent_changes.csv')
    fig=plt.figure(figsize=(10.0,9.4))
    gs=fig.add_gridspec(3,2,height_ratios=[1.0,1.05,1.0],hspace=.55,wspace=.28)
    a=fig.add_subplot(gs[0,0]); b=fig.add_subplot(gs[0,1]); c=fig.add_subplot(gs[1,:]); d=fig.add_subplot(gs[2,:])
    # A
    x=np.arange(len(conf)); bottoms=np.zeros(len(conf)); cols=[BLUE,'#9FBFD0',ORANGE,LIGHT]; labs=['HC','MC','LC','NC']
    for lab,col in zip(labs,cols):
        a.bar(x,conf[lab],bottom=bottoms,color=col,width=.70,label=lab,edgecolor='white',linewidth=.3)
        bottoms += conf[lab].to_numpy()
    a.set_xticks(x,[f'P{t}' for t in THS]); a.set_ylabel('Model genes'); a.set_ylim(0,float(bottoms.max())*1.16)
    a.set_title('Gene confidence classes',pad=8); clean_ax(a,'y'); panel(a,'A')
    a.legend(ncol=4,loc='upper center',bbox_to_anchor=(.5,1.02),frameon=False,columnspacing=1.1,handlelength=1.1)
    # B
    b.plot(st.threshold,st.reaction_count,marker='o',color=BLUE,lw=2.0,ms=5)
    b.set_xticks(THS,[f'P{x}' for x in THS]); b.set_ylabel('Reactions'); b.set_title('Reaction count',pad=8); clean_ax(b,'y'); panel(b,'B')
    for th in [5,25,45]:
        r=st[st.threshold==th].iloc[0]; b.annotate(f"{int(r.reaction_count):,}",(th,r.reaction_count),xytext=(0,8),textcoords='offset points',ha='center',fontsize=9,color=DARK)
    # C full-width
    y=np.arange(len(chg)); c.barh(y,-chg.lost,color=BLUE,alpha=.85,label='Lost',height=.72); c.barh(y,chg.gained,color=ORANGE,alpha=.9,label='Gained',height=.72)
    c.axvline(0,color=DARK,lw=.9); c.set_yticks(y,chg.transition); c.invert_yaxis(); c.set_xlabel('Reaction count change'); c.set_title('Reaction gains and losses',pad=8)
    clean_ax(c,'x'); panel(c,'C',x=-.06,y=1.03)
    c.legend(ncol=2,frameon=False,loc='upper right',bbox_to_anchor=(.99,1.00))
    c.margins(y=.08)
    # D full-width
    xx=np.arange(len(chg)); d.plot(xx,chg.jaccard,marker='o',color=BLUE,lw=2.0,ms=5)
    imin=chg.jaccard.idxmin(); i=list(chg.index).index(imin); d.scatter(i,chg.loc[imin,'jaccard'],s=65,color=ORANGE,zorder=5)
    d.annotate(f"lowest: {chg.loc[imin,'jaccard']:.3f}",(i,chg.loc[imin,'jaccard']),xytext=(10,-19),textcoords='offset points',fontsize=9,color=ORANGE)
    d.set_xticks(xx,chg.transition,rotation=25,ha='right'); d.set_ylim(.91,.97); d.set_ylabel('Jaccard index'); d.set_title('Similarity of adjacent reaction sets',pad=8)
    clean_ax(d,'y'); panel(d,'D',x=-.06,y=1.03)
    fig.suptitle('Effect of expression cutoff on GSE17576 model structure',y=.985)
    fig.subplots_adjust(top=.91,bottom=.075,left=.09,right=.98)
    save(fig,'Figure_3_1_FINAL')

# 3.4: explicit title/legend band and panel letters separated from labels

def fig_3_4():
    d=pd.read_csv(DATA/'pathway_kendall_w.csv')
    order=['TCA cycle','Redox shuttles','Pyruvate oxidation','Oxidative phosphorylation','Lactate metabolism','Ketone metabolism','Glycolysis','Glutamate/aKG metabolism','Glucose uptake/phosphorylation','Fatty acid uptake & beta-oxidation','Creatine phosphate buffering','Glutamine metabolism']
    fig,axs=plt.subplots(1,2,figsize=(10.4,6.5),sharey=True,gridspec_kw={'wspace':.17})
    for ax,src,col,letter in zip(axs,['GSE17576','GSE305719'],[BLUE,ORANGE],['A','B']):
        z=d[d.source==src].set_index('pathway').reindex(order).reset_index(); y=np.arange(len(z))
        ax.scatter(z.W,y,s=70,facecolors=[col if s else 'white' for s in z.significant_q_lt_0_05],edgecolors=col,lw=1.5,zorder=3)
        for yi,r in z.iterrows():
            if r.significant_q_lt_0_05: ax.text(r.W+.03,yi,f"W={r.W:.3f}",va='center',fontsize=8.8,color=col)
        ax.set_xlim(0,.84); ax.set_xlabel("Kendall's W"); ax.set_title(src,pad=12); clean_ax(ax,'x'); panel(ax,letter,x=-.07,y=1.025)
    axs[0].set_yticks(np.arange(len(order)),order); axs[0].invert_yaxis()
    legend=[Line2D([0],[0],marker='o',ls='',mfc=DARK,mec=DARK,markersize=7,label='q < 0.05'),Line2D([0],[0],marker='o',ls='',mfc='white',mec=DARK,markersize=7,label='q ≥ 0.05')]
    fig.legend(handles=legend,loc='upper center',bbox_to_anchor=(.5,.915),ncol=2,frameon=False,columnspacing=2.0)
    fig.suptitle('Effect of expression cutoff on pathway flux variability',y=.985)
    fig.subplots_adjust(top=.84,bottom=.085,left=.20,right=.98)
    save(fig,'Figure_3_4_FINAL')

# 3.5: taller, larger panels and generous inter-panel gaps

def fig_3_5():
    d=pd.read_csv(DATA/'pathway_hfd_effect_ranges.csv')
    paths=['Fatty acid uptake & beta-oxidation','Glutamate/aKG metabolism','Glycolysis','Lactate metabolism','Oxidative phosphorylation','Pyruvate oxidation','TCA cycle']
    tissues=['Gastrocnemius','Soleus','TibialisAnt']; offs={'GSE17576':-.13,'GSE305719':.13}; cols={'GSE17576':BLUE,'GSE305719':ORANGE}
    fig,axs=plt.subplots(3,1,figsize=(10.4,9.7),sharex=True,gridspec_kw={'hspace':.36})
    for ax,tis,letter in zip(axs,tissues,['A','B','C']):
        y=np.arange(len(paths)); ax.axvline(0,color=DARK,lw=1)
        for src in ['GSE17576','GSE305719']:
            z=d[(d.source==src)&(d.tissue==tis)].set_index('pathway').reindex(paths).reset_index(); yy=y+offs[src]
            ax.hlines(yy,z['min'],z['max'],color=cols[src],lw=1.8,alpha=.72); ax.scatter(z['median'],yy,color=cols[src],s=46,zorder=3)
        ax.set_yticks(y,paths); ax.invert_yaxis(); ax.set_title(tissue_label(tis),loc='left',pad=10); clean_ax(ax,'x'); panel(ax,letter,x=-.12,y=1.045)
        ax.margins(y=.10)
    axs[-1].set_xlabel('log$_2$(HFD median / Control median)'); axs[-1].set_xlim(-1.2,1.9)
    handles=[Line2D([0],[0],marker='o',color=BLUE,lw=1.8,label='GSE17576'),Line2D([0],[0],marker='o',color=ORANGE,lw=1.8,label='GSE305719')]
    fig.legend(handles=handles,ncol=2,loc='upper center',bbox_to_anchor=(.56,.94),frameon=False)
    fig.suptitle('HFD–Control pathway effects across expression cutoffs',y=.985)
    fig.subplots_adjust(top=.88,bottom=.07,left=.24,right=.98)
    save(fig,'Figure_3_5_FINAL')

# 3.6: bigger six-panel figure, with no compressed titles/labels

def fig_3_6():
    s175=pd.read_csv(DATA/'gse17576_structure_atp.csv'); s305=pd.read_csv(DATA/'gse305719_structure_atp.csv'); cross=pd.read_csv(DATA/'cross_source_structure.csv')
    anc=pd.read_csv(DATA/'anchor_flux_summary.csv'); sat=pd.read_csv(DATA/'pdhm_vs_atp_saturation.csv')
    fig,axs=plt.subplots(2,3,figsize=(11.0,9.3),gridspec_kw={'hspace':.40,'wspace':.34})
    ax=axs[0,0]; ax.plot(s175.threshold,s175.reaction_count,color=BLUE,marker='o',ms=5,lw=2,label='GSE17576'); ax.plot(s305.threshold,s305.reaction_count,color=ORANGE,marker='o',ms=5,lw=2,label='GSE305719'); ax.set_title('Reaction count',pad=9); ax.set_ylabel('Reactions'); ax.set_xticks(THS,[f'P{x}' for x in THS]); clean_ax(ax,'y'); panel(ax,'A',x=-.11)
    ax=axs[0,1]; ax.plot(cross.threshold_num,cross.jaccard,color=DARK,marker='o',ms=5,lw=2); ax.set_title('Reaction-set similarity',pad=9); ax.set_ylabel('Cross-source Jaccard index'); ax.set_xticks(THS,[f'P{x}' for x in THS]); ax.set_ylim(.86,.95); clean_ax(ax,'y'); panel(ax,'B',x=-.11)
    ax=axs[0,2]; ax.plot(s175.threshold,s175.max_ATP,color=BLUE,marker='o',ms=5,lw=2); ax.plot(s305.threshold,s305.max_ATP,color=ORANGE,marker='o',ms=5,lw=2); ax.axhline(3.5,color=DARK,ls='--',lw=1); ax.set_title('Maximum ATP capacity',pad=9); ax.set_ylabel('mmol gDW$^{-1}$ h$^{-1}$'); ax.set_xticks(THS,[f'P{x}' for x in THS]); clean_ax(ax,'y'); panel(ax,'C',x=-.11)
    for ax,prefix,title,letter,ylim in [(axs[1,0],'HEX1','HEX1 sampled flux','D',(.72,.88)),(axs[1,1],'PDHm','PDHm sampled flux','E',(.35,.88))]:
        for src,col in [('GSE17576',BLUE),('GSE305719',ORANGE)]:
            z=anc[anc.source==src].sort_values('threshold'); ax.fill_between(z.threshold,z[f'{prefix}_q25_rel'],z[f'{prefix}_q75_rel'],color=col,alpha=.13,lw=0); ax.plot(z.threshold,z[f'{prefix}_median_rel'],color=col,marker='o',ms=4.6,lw=1.8)
        ax.set_title(title,pad=9); ax.set_ylabel('Median / DMI target'); ax.set_xticks(THS,[f'P{x}' for x in THS]); ax.set_ylim(*ylim); clean_ax(ax,'y'); panel(ax,letter,x=-.11)
    ax=axs[1,2]
    for src,col in [('GSE17576',BLUE),('GSE305719',ORANGE)]:
        z=sat[sat.source==src]; ax.scatter(z.ATP_demand_fraction,z.PDHm_median_rel,color=col,s=54,alpha=.9,label=src)
    ax.set_xlabel('ATP demand / maximum ATP capacity'); ax.set_ylabel('PDHm median / DMI target'); ax.set_title('PDHm and ATP-demand fraction',pad=9); clean_ax(ax,'both'); panel(ax,'F',x=-.11)
    ax.text(.04,.95,'Pearson r = 0.940\nSpearman ρ = 0.670',transform=ax.transAxes,ha='left',va='top',fontsize=9.2,color=DARK,bbox=dict(facecolor='white',edgecolor='none',alpha=.85,pad=1.5))
    fig.legend(handles=[Line2D([0],[0],color=BLUE,marker='o',label='GSE17576'),Line2D([0],[0],color=ORANGE,marker='o',label='GSE305719')],ncol=2,loc='upper center',bbox_to_anchor=(.5,.94),frameon=False)
    fig.suptitle('Effect of transcriptomic dataset on model structure and DMI-constrained fluxes',y=.985)
    fig.subplots_adjust(top=.87,bottom=.07,left=.08,right=.98)
    save(fig,'Figure_3_6_FINAL')

# Consistent math: same STIX font size and same 300-dpi physical scale.

def render_eq(name, math):
    with plt.rc_context({'font.family':'STIXGeneral','mathtext.fontset':'stix'}):
        fig=plt.figure(figsize=(12,1.3),dpi=300)
        fig.patch.set_alpha(0)
        ax=fig.add_axes([0,0,1,1]); ax.axis('off')
        ax.text(.5,.5,math,ha='center',va='center',fontsize=16,color='black')
        out=EQOUT/f'{name}.png'
        fig.savefig(out,dpi=300,transparent=True,bbox_inches='tight',pad_inches=0.025)
        plt.close(fig)
        return out

def equations():
    eqs={
      'eq_2_1': r'$\bar{x}_{g}=\frac{1}{3}\sum_{c=1}^{3}\bar{x}_{g,c},\qquad T_q=Q_q(\bar{x}_{g})$',
      'eq_2_2': r'$C_R^{\mathrm{OR}}=\max(C_{g_1},\ldots,C_{g_m}),\qquad C_R^{\mathrm{AND}}=\min(C_{g_1},\ldots,C_{g_n})$',
      'eq_2_3': r'$P(v)=\sum_i\frac{w_i}{\sigma_i}(s_{i,L}+s_{i,U}),\qquad s_{i,L}=\max(0,L_i-v_i),\qquad s_{i,U}=\max(0,v_i-U_i)$',
      'eq_2_3b': r'$P(v)\leq 1.5P^{*}+0.1$',
      'eq_2_4': r'$v_j^{\min}=\min_v v_j,\qquad v_j^{\max}=\max_v v_j\quad\mathrm{subject\ to}\quad Sv=0,\;l\leq v\leq u,\;P(v)\leq 1.5P^{*}+0.1$',
      'eq_2_5': r'$J(A,B)=\frac{|A\cap B|}{|A\cup B|}$',
      'eq_2_6a': r'$A_P^{(s)}=\frac{1}{|R_P|}\sum_{r\in R_P}|v_r^{(s)}|$',
      'eq_2_6b': r'$\mathrm{IQR}_P=Q_{0.75}(A_P)-Q_{0.25}(A_P)$',
      'eq_2_7': r'$E_{p,m,q}=\log_2\!\left(\frac{\operatorname{median}(A_{p,m,q}^{\mathrm{HFD}})}{\operatorname{median}(A_{p,m,q}^{\mathrm{Control}})}\right)$',
    }
    for k,v in eqs.items(): render_eq(k,v)

if __name__=='__main__':
    fig_2_1(); fig_3_1(); fig_3_4(); fig_3_5(); fig_3_6(); equations()
    print('Generated final v2 figures and equations')
