import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'src'))
from plot_common import *
from scipy.stats import gaussian_kde

def main():
    d=read('adni_scores.csv');base=d[d.JSD_BASELINE_ROW].copy()
    tr=read('adni_trajectory_included.csv');tr=tr[tr.ANALYSIS.eq('TRAJECTORY_ALL_AVAILABLE_L0')]
    effects=read('adni_trajectory_effects.csv');effects=effects[effects.MODEL.eq('L0')]
    curves=read('adni_adjusted_curves.csv');fig=figure(170)
    panel(fig,.025,.975,'A','Baseline score distributions')
    panel(fig,.545,.975,'B','Adjusted annual score change')
    ax=style(fig.add_axes([.115,.585,.355,.325]),True)
    rng=np.random.default_rng(20261102)
    for i,g in enumerate(['CN','MCI','AD']):
        vals=base.loc[base.DX_STD.eq(g),'JSD_SCORE_Z'].to_numpy();col=COLORS[g]
        yy=np.linspace(vals.min(),vals.max(),256);density=gaussian_kde(vals)(yy);dd=.28*density/density.max()
        ax.fill_betweenx(yy,i-dd,i,color=col,alpha=.2,lw=0);ax.plot(i-dd,yy,color=col,lw=.8)
        ax.scatter(i+.05+rng.uniform(0,.13,len(vals)),vals,s=5,color=col,alpha=.40,lw=0)
        q=np.quantile(vals,[.25,.5,.75]);ax.plot([i,i],q[[0,2]],color=INK,lw=2.5)
        ax.scatter(i,q[1],s=17,facecolors='white',edgecolors=INK,lw=.7,zorder=4)
    ax.set_xlim(-.48,2.45);ax.set_xticks(range(3),[f'{LABELS[g]}\nn = {sum(base.DX_STD.eq(g))}' for g in ['CN','MCI','AD']]);ax.tick_params(axis='x',length=0)
    ax.set_ylabel('AD-like JSD score (ADNI SD units)',labelpad=5)
    ax=fig.add_axes([.675,.615,.265,.285]);rows=[]
    for g in ['CN','MCI','AD']:
        r=effects[effects.CONTRAST.eq(g)].iloc[0].copy();r['LABEL']=LABELS[g];rows.append(r)
    forest(ax,rows,'Annual score change\n(SD units/year)',[CN,MCI,AD],True)
    diff=effects[effects.CONTRAST.eq('AD-CN')].iloc[0]
    fig.text(.675,.495,f'AD–CN difference: {diff.ESTIMATE:.3f}\n95% CI {diff.CI_LOW:.3f} to {diff.CI_HIGH:.3f}',fontsize=8,color=GRAY,linespacing=1.5,va='top')
    panel(fig,.025,.46,'C','Observed trajectories and adjusted group estimates')
    ylo=min(tr.JSD_SCORE_Z.min(),curves.CI_LOW.min());yhi=max(tr.JSD_SCORE_Z.max(),curves.CI_HIGH.max());pad=.05*(yhi-ylo)
    for k,g in enumerate(['CN','MCI','AD']):
        ax=style(fig.add_axes([.115+k*.305,.105,.245,.265]),True);q=tr[tr.BASELINE_DX.eq(g)];col=COLORS[g]
        for _,sub in q.groupby('PTID'):
            sub=sub.sort_values('TIME');ax.plot(sub.TIME,sub.JSD_SCORE_Z,color=col,alpha=.16,lw=.35)
            ax.scatter(sub.TIME,sub.JSD_SCORE_Z,color=col,alpha=.18,s=1,lw=0)
        c=curves[curves.GROUP.eq(g)]
        ax.fill_between(c.TIME.to_numpy(),c.CI_LOW.to_numpy(),c.CI_HIGH.to_numpy(),color=col,alpha=.17,lw=0,zorder=4)
        ax.plot(c.TIME,c.ESTIMATE,color='white',lw=3.5,zorder=5);ax.plot(c.TIME,c.ESTIMATE,color=col,lw=1.7,zorder=6)
        ax.set_xlim(0,q.TIME.max()*1.015);ax.set_ylim(ylo-pad,yhi+pad);ax.xaxis.set_major_locator(MaxNLocator(4,integer=True))
        ax.set_title(f'{LABELS[g]}\n{q.PTID.nunique()} people; {len(q):,} MRI',fontsize=8.5,pad=8,linespacing=1.5)
        ax.set_xlabel('Years from first MRI',fontsize=8.5)
        if k==0:ax.set_ylabel('AD-like JSD score',labelpad=6)
        else:ax.tick_params(labelleft=False)
    fig.text(.115,.028,'Group-specific observed follow-up ranges; shaded bands show 95% CIs for adjusted fixed effects.',fontsize=8,color=GRAY)
    save(fig,1)
if __name__=='__main__':main()
