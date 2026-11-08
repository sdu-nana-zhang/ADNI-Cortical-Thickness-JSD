import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'src'))
from plot_common import *
from matplotlib.collections import PolyCollection
from matplotlib.lines import Line2D
import nibabel as nib

def brain_image(hem,view):
    cache=ROOT/'data/anatomy'/f'{hem}_{view}_degree.png'
    if cache.exists():return Image.open(cache)
    v,f=nib.freesurfer.read_geometry(DATA/f'anatomy/{hem}.pial')
    lab,ctab,names=nib.freesurfer.read_annot(DATA/f'anatomy/{hem}.aparc.DKTatlas.annot')
    names=[n.decode() for n in names];roi=read('regional_involvement.csv')
    ri=np.zeros(len(v),int)
    for r in roi[roi.HEMISPHERE.eq(hem)].itertuples():
        assert r.REGION in names,r.REGION
        ri[lab==names.index(r.REGION)]=r.INDEX
    tri=v[f];normal=np.cross(tri[:,1]-tri[:,0],tri[:,2]-tri[:,0])
    normal/=np.maximum(np.linalg.norm(normal,axis=1)[:,None],1e-20)
    sign=-1 if hem=='lh' else 1
    if view=='medial':sign=-sign
    camera=np.array([sign,0.,0.]);u=np.array([0.,sign,0.]);up=np.array([0.,0.,1.])
    xy=np.c_[v@u,v@up];depth=v@camera
    visible=np.flatnonzero(normal@camera>0)
    visible=visible[np.argsort(depth[f[visible]].mean(axis=1))]
    light=camera+.45*up-.35*u;light/=np.linalg.norm(light)
                                                                                             
    shade=.57+.43*np.maximum(normal@light,0)
    fr=ri[f]
                                                                                                 
    label=np.where(fr[:,0]==fr[:,1],fr[:,0],np.where(fr[:,1]==fr[:,2],fr[:,1],fr[:,0]))
    values=np.r_[0,roi.sort_values('INDEX').DEGREE.to_numpy()][label]
    rgba=DEGREE(Normalize(0,roi.DEGREE.max())(values));rgba[label==0,:3]=[.82,.84,.85]
    rgba[:,:3]*=shade[:,None]
    fig=plt.figure(figsize=(3.1,2.25),dpi=450);ax=fig.add_axes([0,0,1,1])
    ax.add_collection(PolyCollection(xy[f[visible]],facecolors=rgba[visible],edgecolors='none',antialiaseds=False))
    ax.set_xlim(xy[:,0].min()-2,xy[:,0].max()+2);ax.set_ylim(xy[:,1].min()-2,xy[:,1].max()+2)
    ax.set_aspect('equal');ax.axis('off');fig.savefig(cache,dpi=450,transparent=True);plt.close(fig)
    print('Rendered actual DKT surface',hem,view,flush=True)
    return Image.open(cache)

def main():
    sig=read('significant277_edges.csv');roi=read('regional_involvement.csv');hom=read('homotopic4_edges.csv')
    fig=figure(215)
    panel(fig,.025,.987,'A','Regional involvement in differential connections')
    for i,(hem,view) in enumerate([('lh','lateral'),('lh','medial'),('rh','medial'),('rh','lateral')]):
        ax=fig.add_axes([.025+i*.241,.806,.235,.158]);ax.imshow(brain_image(hem,view));ax.axis('off')
        fig.text(.1425+i*.241,.803,('Left' if hem=='lh' else 'Right')+' '+view,ha='center',fontsize=8)
    cax=fig.add_axes([.34,.772,.34,.010]);cb=fig.colorbar(matplotlib.cm.ScalarMappable(norm=Normalize(0,roi.DEGREE.max()),cmap=DEGREE),cax=cax,orientation='horizontal')
    cb.set_label('Differential connections per region',fontsize=8,labelpad=3);cb.ax.tick_params(labelsize=8,length=2,pad=2);cb.outline.set_linewidth(.4)
    panel(fig,.025,.713,'B','AD–CN slope differences')
    panel(fig,.545,.713,'C','Group-specific JSD slopes')
    ax=fig.add_axes([.115,.385,.35,.287]);mat=np.full((62,62),np.nan)
    for r in sig.itertuples():
        i,j=int(r.ROI1_INDEX)-1,int(r.ROI2_INDEX)-1
        mat[i,j]=mat[j,i]=r.ESTIMATE*1000
    lim=np.nanmax(abs(mat));cmap=DIVERGE.copy();cmap.set_bad('#F0F2F3')
    ax.set_facecolor('#F0F2F3')
                                                                                             
    mesh=ax.pcolormesh(np.arange(63),np.arange(63),np.ma.masked_invalid(mat),cmap=cmap,vmin=-lim,vmax=lim,shading='flat',edgecolors='none')
    ax.set_xlim(0,62);ax.set_ylim(62,0);ax.set_aspect('equal');ax.axvline(31,color='white',lw=.8);ax.axhline(31,color='white',lw=.8)
    ax.set_xticks([15.5,46.5],['Left','Right']);ax.set_yticks([15.5,46.5],['Left','Right']);ax.tick_params(length=0)
    for s in ax.spines.values():s.set_visible(False)
    ax.set_ylabel('Cortical regions',labelpad=5)
    cax=fig.add_axes([.115,.352,.35,.012]);cb=fig.colorbar(mesh,cax=cax,orientation='horizontal');cb.outline.set_linewidth(.4)
    cb.set_label('AD–CN difference (10$^{-3}$ JSD/year)',fontsize=8,labelpad=3);cb.ax.tick_params(labelsize=8,length=2,pad=2)
    cb.locator=MaxNLocator(5);cb.update_ticks()
    ax=style(fig.add_axes([.64,.385,.32,.287]));zero(ax,True,True)
    categories=['Both decreasing','Both increasing','CN decreasing; AD increasing','CN increasing; AD decreasing']
    cols=[CN,AD,PURPLE,MCI];marks=['o','s','^','D'];handles=[]
    for name,col,marker in zip(categories,cols,marks):
        q=sig[sig.QUADRANT.eq(name)]
        ax.scatter(q.CN_SLOPE*1000,q.AD_SLOPE*1000,s=10,c=col,marker=marker,alpha=.75,linewidths=.15,edgecolors='white')
        legend_name=name.replace('; ', ';\n')
        handles.append(Line2D([],[],ls='',marker=marker,color=col,markersize=4,label=f'{legend_name} ({len(q)})'))
    vals=np.r_[sig.CN_SLOPE,sig.AD_SLOPE]*1000;lo=min(vals.min(),0);hi=max(vals.max(),0);pad=.06*(hi-lo)
    ax.set_xlim(lo-pad,hi+pad);ax.set_ylim(lo-pad,hi+pad)
    ax.plot([lo-pad,hi+pad],[lo-pad,hi+pad],color='#C3CDD4',ls='--',lw=.65,zorder=0)
    ax.set_xlabel('CN slope (10$^{-3}$ JSD/year)',labelpad=4);ax.set_ylabel('Clinical AD slope\n(10$^{-3}$ JSD/year)',labelpad=3)
    ax.legend(handles=handles,loc='upper center',bbox_to_anchor=(.5,-.25),ncol=2,frameon=False,fontsize=7.7,handletextpad=.4,columnspacing=.6,borderaxespad=0)
    panel(fig,.025,.252,'D','Four homotopic connections')
    names={'caudalmiddlefrontal':'Caudal middle frontal','inferiorparietal':'Inferior parietal','precentral':'Precentral','superiortemporal':'Superior temporal'}
    order=list(names);hom['name']=hom.ROI1.str.replace('lh_','',regex=False);hom=hom.set_index('name').loc[order].reset_index()
    rows=[]
    for _,r in hom.iterrows():
        r=r.copy();r['LABEL']=names[r['name']]
        for col in ['ESTIMATE','CI_LOW','CI_HIGH']:r[col]*=1000
        rows.append(r)
    ax=fig.add_axes([.335,.078,.47,.14]);forest(ax,rows,'AD–CN slope difference (10$^{-3}$ JSD/year)',colors=[PURPLE,GREEN,MCI,AD])
    q_labels=[]
    for i,r in enumerate(rows):q_labels.append(ax.text(1.10,i,f'q = {r.Q_BH1891:.3f}',transform=ax.get_yaxis_transform(),va='center',fontsize=8))
    fig.canvas.draw()
    boxes=[t.get_window_extent(fig.canvas.get_renderer()) for t in q_labels]
    q_center=fig.transFigure.inverted().transform(((min(b.x0 for b in boxes)+max(b.x1 for b in boxes))/2,0))[0]
    fig.text(q_center,.231,'BH q',fontsize=8,color=GRAY,ha='center')
    save(fig,3)
if __name__=='__main__':main()
