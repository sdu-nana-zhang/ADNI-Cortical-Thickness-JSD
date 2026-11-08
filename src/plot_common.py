from pathlib import Path
import json,os
from project import stage
os.environ.setdefault('MPLCONFIGDIR',str(stage('figures')/'.mplconfig'))
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from matplotlib.colors import LinearSegmentedColormap, Normalize, to_rgb
from PIL import Image

ROOT=stage('figures')
DATA=ROOT/'data'
INK='#263746';GRAY='#738390';GRID='#E9EEF1'
CN='#236E98';MCI='#D29A31';AD='#C9586E';GREEN='#527D70';PURPLE='#8574AF'
COLORS={'CN':CN,'MCI':MCI,'AD':AD}
LABELS={'CN':'CN','MCI':'MCI','AD':'Clinical AD'}
DEGREE=LinearSegmentedColormap.from_list('degree',['#EDF2F5','#BCD2E0','#719DB9','#285D7B'])
DIVERGE=LinearSegmentedColormap.from_list('slope_difference',['#246F99','#8CB4CD','#FFFFFF','#E4B39D','#C36E4D'])
plt.rcParams.update({'font.family':'Arial','font.size':8.5,'axes.labelsize':9,'axes.titlesize':9.5,
 'axes.titleweight':'bold','axes.labelcolor':INK,'text.color':INK,'xtick.color':INK,'ytick.color':INK,
 'axes.edgecolor':'#93A1AA','axes.linewidth':.65,'xtick.labelsize':8,'ytick.labelsize':8,
 'xtick.major.width':.6,'ytick.major.width':.6,'xtick.major.size':3,'ytick.major.size':3,
 'pdf.fonttype':42,'ps.fonttype':42,'svg.fonttype':'none','figure.facecolor':'white','savefig.facecolor':'white',
 'lines.linewidth':1.2,'axes.unicode_minus':True})
def read(name):return pd.read_csv(DATA/name)
def figure(height_mm):return plt.figure(figsize=(180/25.4,height_mm/25.4))
def panel(fig,x,y,letter,title):
    fig.text(x,y,letter,fontsize=12,fontweight='bold',va='top')
    fig.text(x+.036,y-.001,title,fontsize=9.5,fontweight='bold',va='top')
def style(ax,grid=False):
    for s in ['top','right']:ax.spines[s].set_visible(False)
    ax.tick_params(pad=3)
    ax.set_axisbelow(True)
    if grid:ax.grid(axis='y',color=GRID,lw=.5)
    ax.xaxis.set_major_locator(MaxNLocator(4))
    ax.yaxis.set_major_locator(MaxNLocator(5))
    return ax
def zero(ax,x=True,y=False):
    if x:ax.axvline(0,color='#ABB9C4',lw=.7,ls=(0,(3,2)),zorder=0)
    if y:ax.axhline(0,color='#ABB9C4',lw=.7,ls=(0,(3,2)),zorder=0)
def interval(ax,y,est,lo,hi,color=CN,marker='o',size=4.7):
    ax.errorbar(est,y,xerr=np.array([[est-lo],[hi-est]]),fmt=marker,ms=size,color=color,
       mfc=color,mec='white',mew=.5,lw=1.5,capsize=0,zorder=3)
def forest(ax,rows,xlabel,colors=None,show_values=False):
    style(ax);zero(ax)
    low=min(0,min(r.CI_LOW for r in rows));high=max(0,max(r.CI_HIGH for r in rows));span=high-low
    ax.set_xlim(low-.10*span,high+(.27 if show_values else .12)*span)
    for i,r in enumerate(rows):
        col=colors[i] if colors else GREEN
        interval(ax,i,r.get('ESTIMATE',r.get('BETA')),r.CI_LOW,r.CI_HIGH,col)
        if show_values:ax.text(r.CI_HIGH+.04*span,i,f'{r.get("ESTIMATE",r.get("BETA")):.3f}',va='center',fontsize=8)
    ax.set_ylim(len(rows)-.5,-.5);ax.set_yticks(range(len(rows)),[r.LABEL for r in rows]);ax.set_xlabel(xlabel,labelpad=6)
    ax.tick_params(axis='y',length=0,pad=5);ax.spines['left'].set_visible(False)
def save(fig,num):
    out=ROOT/f'figure{num}';out.mkdir(exist_ok=True)
                                                                     
    for ext in ['pdf','svg','png']:
        fig.savefig(out/f'Figure{num}.{ext}',dpi=600)
    with Image.open(out/f'Figure{num}.png') as im:
        im.save(out/f'Figure{num}.tiff',dpi=(600,600),compression='tiff_lzw')
        im.thumbnail((1400,1600));im.save(out/f'Figure{num}_preview.png')
    fig.canvas.draw()
    renderer=fig.canvas.get_renderer();w,h=fig.canvas.get_width_height();overflow=[]
                                                                                  
                                                                                
    offrange=set()
    for ax in fig.axes:
        for axis,limits in [(ax.xaxis,ax.get_xlim()),(ax.yaxis,ax.get_ylim())]:
            for tick in axis.get_major_ticks()+axis.get_minor_ticks():
                if not min(limits)-1e-9<=tick.get_loc()<=max(limits)+1e-9:
                    offrange.update([id(tick.label1),id(tick.label2)])
    for text in fig.findobj(matplotlib.text.Text):
        if id(text) in offrange or not text.get_visible() or not text.get_text():continue
        b=text.get_window_extent(renderer)
        if b.x0 < -1 or b.y0 < -1 or b.x1>w+1 or b.y1>h+1:overflow.append(text.get_text())
    (out/'layout_check.json').write_text(json.dumps({'text_outside_canvas':overflow,'width_mm':180,'height_mm':round(fig.get_figheight()*25.4,2)},indent=2),encoding='utf-8')
    plt.close(fig)
    print(f'Figure {num}: exported PDF, SVG, PNG and TIFF; canvas overflow: {overflow}',flush=True)
