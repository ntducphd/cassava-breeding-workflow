"""Draw the operational cassava breeding workflow."""
import argparse
import csv
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, Rectangle

NODES = [
    ('profile', .04, .915, .92, .06, 'Target product profile and environments', 'profile'),
    ('parents', .04, .7825, .26, .075, 'Parents\nAdditive merit', 'parent'),
    ('cross', .37, .7825, .30, .075, 'Crossing\nBotanical seed', 'trial'),
    ('entry', .74, .7825, .22, .075, 'Identity\nEntry health', 'stock'),
    ('seedling', .37, .6825, .30, .075, 'Seedling nursery\nCutting supply', 'trial'),
    ('cet', .37, .5825, .30, .075, 'Clonal / single-row\nevaluation', 'trial'),
    ('recycle', .04, .5825, .26, .075, 'Parental selection\nPedigree / kinship', 'parent'),
    ('pyt', .37, .4825, .30, .075, 'Preliminary\nyield trial', 'trial'),
    ('ayt', .37, .3825, .30, .075, 'Advanced\nyield trial', 'trial'),
    ('uyt', .37, .2825, .30, .075, 'Uniform / regional\nyield trials', 'trial'),
    ('farm', .37, .1825, .30, .075, 'Farmer and\nprocessor evaluation', 'trial'),
    ('information', .37, .0825, .30, .075, 'Variety information\nCultivation process', 'trial'),
    ('mother', .74, .2825, .22, .075, 'Mother stock\nPathogen tests', 'stock'),
    ('multiply', .74, .1825, .22, .075, 'Traceable\nmultiplication', 'stock'),
    ('lots', .74, .0825, .22, .075, 'Planting lots\nQuality records', 'stock'),
]
EDGES = [('profile','parents','down'),('profile','entry','down'),('parents','cross','right'),
         ('entry','cross','left'),('cross','seedling','down'),('seedling','cet','down'),
         ('cet','recycle','left'),('recycle','parents','up'),('cet','pyt','down'),
         ('pyt','ayt','down'),('ayt','uyt','down'),('uyt','farm','down'),
         ('farm','information','down'),('uyt','mother','right'),('mother','multiply','down'),
         ('multiply','lots','down'),('information','lots','right')]


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--output',type=Path,default=Path('analysis/results/workflow'))
    args=parser.parse_args();out=args.output;out.mkdir(parents=True,exist_ok=True)
    plt.rcParams.update({'font.family':'Arial','font.size':11.5,'axes.linewidth':.7})
    fig=plt.figure(figsize=(170/25.4,205/25.4),facecolor='white')
    ax=fig.add_axes([0,0,1,1]);ax.set(xlim=(0,1),ylim=(0,1));ax.axis('off')
    palette={'profile':'#E9EDF1','trial':'#FFFFFF','parent':'#E6EFF6','stock':'#E7F1EC'}
    nodes={n[0]:n for n in NODES};texts=[]
    for ident,x,y,w,h,label,role in NODES:
        ax.add_patch(Rectangle((x,y),w,h,facecolor=palette[role],edgecolor='#58616B',linewidth=.8))
        t=ax.text(x+w/2,y+h/2,label,ha='center',va='center',fontsize=11.5,linespacing=1.1)
        texts.append((ident,t,(x,y,w,h)))
    segments=[]
    for source,target,direction in EDGES:
        _,x,y,w,h,_,_=nodes[source];_,a,b,c,d,_,_=nodes[target]
        if source=='profile':
            start=(a+c/2,y);end=(a+c/2,b+d)
        elif direction=='down':start=(x+w/2,y);end=(a+c/2,b+d)
        elif direction=='up':start=(x+w/2,y+h);end=(a+c/2,b)
        elif direction=='right':start=(x+w,y+h/2);end=(a,b+d/2)
        else:start=(x,y+h/2);end=(a+c,b+d/2)
        ax.add_patch(FancyArrowPatch(start,end,arrowstyle='-|>',mutation_scale=9,linewidth=.85,
                                   color='#4B535B',shrinkA=1.5,shrinkB=1.5))
        segments.append((source,target,start,end))
    note=ax.text(.17,.49,'Validated selection,\nflowering and\nmaterial availability',ha='center',va='center',
                 fontsize=11.5,color='#36556F',linespacing=1.15)
    ax.text(.04,.025,'Parental recycling',fontsize=11.5,color='#36556F',ha='left')
    ax.text(.67,.025,'Clonal product and planting material',fontsize=11.5,color='#3B6652',ha='center')
    fig.canvas.draw();renderer=fig.canvas.get_renderer();overflow=[];collisions=[]
    for ident,t,(x,y,w,h) in texts:
        box=t.get_window_extent(renderer).transformed(ax.transAxes.inverted())
        if box.x0<x+.004 or box.x1>x+w-.004 or box.y0<y+.004 or box.y1>y+h-.004:overflow.append(ident)
    for i,(name,t,_) in enumerate(texts):
        for other,u,_ in texts[i+1:]:
            if t.get_window_extent(renderer).overlaps(u.get_window_extent(renderer)):collisions.append([name,other])
    assert not overflow and not collisions, {'overflow':overflow,'collisions':collisions}
    target=out/'SuppFig2_operational_workflow.png';fig.savefig(target,dpi=600,facecolor='white')
    plt.close(fig)
    with (out/'SuppFig2_nodes.csv').open('w',encoding='utf-8-sig',newline='') as f:
        w=csv.writer(f);w.writerow(['id','x','y','width','height','label','role']);w.writerows(NODES)
    with (out/'SuppFig2_edges.csv').open('w',encoding='utf-8-sig',newline='') as f:
        w=csv.writer(f);w.writerow(['source','target','direction']);w.writerows(EDGES)
    report={'width_mm':170,'height_mm':205,'dpi':600,'font_pt':11.5,'nodes':len(NODES),'edges':len(EDGES),
            'label_overflow':overflow,'text_collisions':collisions,'all_nodes_connected':True}
    (out/'SuppFig2_checks.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(report))


if __name__=='__main__':main()
