import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'src'))
from plot_common import *

def main():
    inc=read('adni_cognition_included.csv');effects=read('adni_cognition_effects.csv');fig=figure(180)
    for k,(ep,term,letter,left,bottom) in enumerate([
       ('ADAS13','SCORE_BETWEEN','A',.115,.64),('ADAS13','SCORE_WITHIN','B',.625,.64),
       ('MMSE','SCORE_BETWEEN','C',.115,.17),('MMSE','SCORE_WITHIN','D',.625,.17)]):
        r=effects[effects.OUTCOME.eq(ep)&effects.TERM.eq(term)].iloc[0];s=inc[inc.ANALYSIS.eq(r.ANALYSIS)].copy()
        between=term=='SCORE_BETWEEN'
        if between:
            v=s.groupby('PTID',sort=False).agg(x=('SCORE_BETWEEN','first'),y=('COGNITION','mean'))
            xl='Participant mean AD-like JSD score';yl=f'Participant mean {ep}';col=CN
        else:
            v=pd.DataFrame({'x':s.SCORE_WITHIN,'y':s.COGNITION-s.groupby('PTID').COGNITION.transform('mean')})
            xl='Score deviation from\nparticipant mean';yl=f'{ep} deviation from\nparticipant mean';col=MCI
        panel(fig,left-.09,bottom+.31,letter,ep+': '+('between-person' if between else 'within-person'))
        ax=style(fig.add_axes([left,bottom,.33,.265]))
        if not between:zero(ax,True,True)
        ax.scatter(v.x,v.y,s=6 if between else 3.5,color=col,alpha=.32 if between else .25,lw=0)
        ax.margins(x=.07,y=.07);ax.set_xlabel(xl,fontsize=8.5,labelpad=5);ax.set_ylabel(yl,fontsize=8.5,labelpad=5)
                                                                              
        fig.text(left,bottom-.095,f'β = {r.BETA:.3f} (95% CI {r.CI_LOW:.3f} to {r.CI_HIGH:.3f})',fontsize=8,fontweight='bold')
        unit='one point per person' if between else f'{len(v):,} paired visits'
        fig.text(left,bottom-.122,f'{int(r.N_PATIENTS):,} people; {int(r.N_MRI):,} paired visits',fontsize=8,color=GRAY)
    save(fig,2)
if __name__=='__main__':main()
