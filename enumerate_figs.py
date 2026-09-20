# Enumerations behind the coloured-chair theorem, with pictures. Python 3 + matplotlib. Stand-alone.
# Usage: python3 enumerate_figs.py [--all]     (--all also draws every touching placement, legal ones highlighted)
# Writes out/<rule>/contacts.png, shells.png, second_layer.png (and contacts_all.png with --all) for rule in arrows, colours,
# and prints the counts. Conventions: red = the chair; blue = neighbours belonging to its unique consistent superchair;
# grey = other neighbours; magenta outline = mirror-image chair. Green/red panel border = shell survives / fails the second layer.
import sys, os, itertools
from itertools import permutations, product, combinations
import numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection, LineCollection
from matplotlib.patches import Rectangle
ALL='--all' in sys.argv
FACES=[((0,0,0),(-1,0,0),'B',(0,1,0)),((0,0,0),(0,-1,0),'W',(1,0,0)),((0,0,0),(0,0,-1),'G',(1,1,0)),
 ((0,0,1),(-1,0,0),'W',(0,1,0)),((0,0,1),(0,-1,0),'B',(1,0,0)),((0,0,1),(0,0,1),'G',(1,1,0)),
 ((0,1,0),(-1,0,0),'B',(0,0,1)),((0,1,0),(0,1,0),'G',(1,0,1)),((0,1,0),(0,0,-1),'W',(1,0,0)),
 ((0,1,1),(-1,0,0),'G',(0,-1,-1)),((0,1,1),(1,0,0),'W',(0,1,0)),((0,1,1),(0,1,0),'B',(0,0,-1)),((0,1,1),(0,0,1),'W',(0,-1,0)),
 ((1,0,0),(1,0,0),'G',(0,1,1)),((1,0,0),(0,-1,0),'W',(0,0,1)),((1,0,0),(0,0,-1),'B',(0,1,0)),
 ((1,0,1),(1,0,0),'W',(0,0,-1)),((1,0,1),(0,-1,0),'G',(-1,0,-1)),((1,0,1),(0,1,0),'B',(1,0,0)),((1,0,1),(0,0,1),'B',(-1,0,0)),
 ((1,1,0),(1,0,0),'B',(0,-1,0)),((1,1,0),(0,1,0),'W',(-1,0,0)),((1,1,0),(0,0,-1),'G',(-1,-1,0)),((1,1,0),(0,0,1),'G',(1,1,0))]
CHILDREN=[(((1,0,0),(0,1,0),(0,0,1)),(0,0,0)),(((0,1,0),(1,0,0),(0,0,-1)),(0,0,4)),(((1,0,0),(0,0,-1),(0,1,0)),(0,4,0)),
 (((0,0,1),(-1,0,0),(0,-1,0)),(0,4,4)),(((0,0,-1),(0,1,0),(1,0,0)),(4,0,0)),(((0,-1,0),(0,0,1),(-1,0,0)),(4,0,4)),
 (((-1,0,0),(0,-1,0),(0,0,1)),(4,4,0)),(((1,0,0),(0,1,0),(0,0,1)),(1,1,1))]         
CELLS=[c for c in product((0,1),repeat=3) if c!=(1,1,1)]
COMP={'B':'W','W':'B','G':'G'}
FRAMES=[tuple(tuple(s[i]*(j==p[i]) for j in range(3)) for i in range(3)) for p in permutations(range(3)) for s in product((1,-1),repeat=3)]
mv=lambda M,v:tuple(sum(M[i][j]*v[j] for j in range(3)) for i in range(3))
mm=lambda A,B:tuple(tuple(sum(A[i][k]*B[k][j] for k in range(3)) for j in range(3)) for i in range(3))
mT=lambda M:tuple(tuple(M[j][i] for j in range(3)) for i in range(3))
add=lambda a,b:tuple(x+y for x,y in zip(a,b)); NEG=lambda v:tuple(-x for x in v)
comp=lambda p,q:(mm(p[0],q[0]),add(mv(p[0],q[1]),p[1]))            
inv=lambda p:(mT(p[0]),NEG(mv(mT(p[0]),p[1])))
det=lambda M:round(np.linalg.det(np.array(M)))
def cells(M,t): return {tuple((x-1)//2 for x in add(mv(M,(2*c[0]+1,2*c[1]+1,2*c[2]+1)),(2*t[0],2*t[1],2*t[2]))) for c in CELLS}
def faces(M,t): return {(add(mv(M,(2*c[0]+1+n[0],2*c[1]+1+n[1],2*c[2]+1+n[2])),(2*t[0],2*t[1],2*t[2])),mv(M,n)):(k,mv(M,d)) for c,n,k,d in FACES}
_C={}
def ok(p,q,arrows):
    key=(p,q,arrows)
    if key not in _C:
        _C[key]=not (cells(*p)&cells(*q)) and all(f[0]==COMP[k] and (not arrows or f[1]==d) for (c,n),(k,d) in faces(*p).items() for f in [faces(*q).get((c,NEG(n)))] if f)
    return _C[key]
ROOT=(FRAMES[0],(0,0,0)); RC=cells(*ROOT)
def around(cs): return {add(c,dv) for c in cs for dv in product((1,-1,0),repeat=3) if sum(map(abs,dv))==1}-cs
NB=around(RC)
def touching(): return [(M,t) for M in FRAMES for t in {add(n,NEG(c)) for n in NB for c in cells(M,(0,0,0))} if cells(M,t)&NB and not cells(M,t)&RC]
TOUCH=touching()
def shells(centre,fixed,cont,arrows):
    here=around(cells(*centre)); taken={c:f for f in fixed for c in cells(*f)}
    cand=[(p,cells(*p)&here) for p in (comp(centre,q) for q in cont) if cells(*p)&here and not taken.keys()&cells(*p) and all(ok(p,f,arrows) for f in fixed)]
    out=[]
    def rec(rem,chosen):
        if not rem: out.append(chosen); return
        c=min(rem)
        if c in taken:
            f=taken[c]; cv=cells(*f)&here
            if cv<=rem and all(ok(f,q,arrows) for q in chosen): rec(rem-cv,chosen+[f])
            return
        for p,cv in cand:
            if c in cv and cv<=rem and all(ok(p,q,arrows) for q in chosen): rec(rem-cv,chosen+[p])
    rec(frozenset(here),[]); return out
PARENTS=[inv(c) for c in CHILDREN]                      
def consistent(sh):                                                   # indices of superchairs consistent with a shell
    out=[]
    for k,P in enumerate(PARENTS):
        kids=[comp(P,c) for c in CHILDREN]
        if all(q==ROOT or (q in sh if cells(*q)&NB else all(ok(q,s,True) or ok(q,s,False) for s in sh)) for q in kids): out.append(k)
    return out
# ---------------- isometric renderer
el,az=np.radians(30),np.radians(40)
V=np.array([np.cos(el)*np.cos(az),np.cos(el)*np.sin(az),np.sin(el)]); R=np.array([-np.sin(az),np.cos(az),0]); U=np.cross(V,R)
P2=lambda p:(np.dot(p,R),np.dot(p,U))
ROLE={'root':'#c0392b','sib':'#2471a3','other':'#7f8c8d','mirror':'#c71585'}
FC={'B':'#2b2b2b','W':'#ffffff','G':'#9e9e9e'}
def draw_tiles(ax,tiles,arrows,offset=None):
    """tiles: list of (pose, role). Draw visible faces (+x,+y,+z) of every cube, far cubes first, edges and arrows with them."""
    items=[]
    for pose,role in tiles:
        off=np.zeros(3) if offset is None else offset(pose); F=faces(*pose)
        for c in cells(*pose): items.append((np.dot(np.array(c)+0.5+off,V),c,off,F,role))
    z=0
    for _,c,off,F,role in sorted(items,key=lambda x:x[0]):
        polys=[];cols=[];arr=[]
        for n in [(1,0,0),(0,1,0),(0,0,1)]:
            i=n.index(1); a,b=[j for j in range(3) if j!=i]
            base=np.array(c,float)+off; base[i]+=1
            q=[base.copy() for _ in range(4)]; q[1][a]+=1; q[2][a]+=1; q[2][b]+=1; q[3][b]+=1
            f=F.get(((2*c[0]+1+n[0],2*c[1]+1+n[1],2*c[2]+1+n[2]),n))
            polys.append([P2(v) for v in q]); cols.append(FC[f[0]] if (f and not arrows) else '#e9e2d3')
            if f and arrows:
                d=np.array(f[1],float); d/=np.linalg.norm(d); e=np.cross(np.array(n,float),d); ctr=base.copy(); ctr[a]+=0.5; ctr[b]+=0.5
                tip=ctr+0.32*d; arr.append(([P2(ctr-0.3*d),P2(tip)],[P2(tip),P2(tip-0.18*d+0.14*e),P2(tip-0.18*d-0.14*e)],{'B':'#111','W':'#fff','G':'#2a6fb0'}[f[0]]))
        ax.add_collection(PolyCollection(polys,facecolors=cols,edgecolors=ROLE[role],linewidths=0.6,zorder=z)); z+=1
        for seg,head,col in arr:
            if col=='#fff': ax.add_collection(LineCollection([seg],colors='#333',linewidths=2.2,zorder=z))
            ax.add_collection(LineCollection([seg],colors=col,linewidths=1.2,zorder=z)); ax.add_patch(plt.Polygon(head,fc=col,ec='#333' if col=='#fff' else col,lw=0.4,zorder=z)); z+=1
    ax.set_aspect('equal'); ax.set_axis_off(); ax.autoscale_view()
def role_of(pose,sibs):
    if pose==ROOT: return 'root'
    if det(pose[0])<0: return 'mirror'
    return 'sib' if pose in sibs else 'other'
def sheet(fname,panels,ncols,size=1.0,border=None):
    n=len(panels); nrows=-(-n//ncols); fig,axs=plt.subplots(nrows,ncols,figsize=(ncols*size,nrows*size),squeeze=False)
    for ax in axs.flat: ax.set_axis_off()
    for k,(tiles,arrows,offset) in enumerate(panels):
        ax=axs.flat[k]; draw_tiles(ax,tiles,arrows,offset)
        if border and border[k]:
            ax.set_axis_on(); ax.set_xticks([]); ax.set_yticks([])
            for s in ax.spines.values(): s.set_edgecolor(border[k]); s.set_linewidth(2.5)
    plt.subplots_adjust(0,0,1,1,0.02,0.02); plt.savefig(fname,dpi=110); plt.close()
def explode(root_c):
    def f(pose):
        if pose==ROOT: return np.zeros(3)
        c=np.mean([np.array(x)+0.5 for x in cells(*pose)],axis=0)-root_c; c/=np.linalg.norm(c); return np.round(c*2)/2
    return f
ROOT_C=np.mean([np.array(x)+0.5 for x in RC],axis=0)
def slices(fname,rows):
    """rows: list of dict cell->colour. One row of z-layer panels per entry."""
    zs=sorted({c[2] for r in rows for c in r}); fig,axs=plt.subplots(len(rows),len(zs),figsize=(1.3*len(zs),1.3*len(rows)),squeeze=False)
    for ax in axs.flat: ax.set_axis_off()
    for i,r in enumerate(rows):
        xs=[c[0] for c in r]; ys=[c[1] for c in r]
        for j,z in enumerate(zs):
            ax=axs[i,j]; ax.set_aspect('equal'); ax.set_xlim(min(xs)-.3,max(xs)+1.3); ax.set_ylim(min(ys)-.3,max(ys)+1.3)
            for c,col in r.items():
                if c[2]==z: ax.add_patch(Rectangle((c[0],c[1]),1,1,fc=col,ec='k',lw=0.4))
    plt.subplots_adjust(0,0,1,1,0.05,0.15); plt.savefig(fname,dpi=110); plt.close()
# ---------------- main
survivors={}; GAL={}
for arrows in (True,False):
    rule='arrows' if arrows else 'colours'; os.makedirs(f'out/{rule}',exist_ok=True)
    cont=[p for p in TOUCH if ok(ROOT,p,arrows)]
    print(f'{rule}: {len(TOUCH)} touching placements, {len(cont)} legal contacts, {sum(det(p[0])<0 for p in cont)} of them mirror images')
    sheet(f'out/{rule}/contacts.png',[([(ROOT,'root'),(p,role_of(p,()))],arrows,None) for p in cont],11 if arrows else 25)
    G=GAL[rule]={'contacts':[{'tiles':[[ROOT,'root'],[p,role_of(p,())]],'legal':ok(ROOT,p,arrows),'mirror':det(p[0])<0} for p in (TOUCH if ALL else cont)]}
    if ALL: sheet(f'out/{rule}/contacts_all.png',[([(ROOT,'root'),(p,role_of(p,()))],arrows,None) for p in TOUCH],50,0.6,
                  border=['#27ae60' if ok(ROOT,p,arrows) else None for p in TOUCH])
    S=[[p for p in s if p!=ROOT] for s in shells(ROOT,[ROOT],cont,arrows)]
    info=[]
    for s in S:
        ks=consistent(s); sibs=set()
        if len(ks)==1: sibs={comp(PARENTS[ks[0]],c) for c in CHILDREN}
        alive=all(shells(c,[ROOT]+s,cont,arrows) for c in s)
        info.append((sibs,alive,ks))
    surv=[frozenset(s) for s,(sibs,alive,ks) in zip(S,info) if alive]; survivors[rule]=set(surv)
    print(f'{rule}: {len(S)} complete neighbourhoods; {sum(len(i[2])==1 for i in info)} with a unique consistent superchair; '
          f'{sum(any(det(p[0])<0 for p in s) for s in S)} containing a mirror image; {len(surv)} survive the second layer')
    G['shells']=[{'tiles':[[ROOT,'root']]+[[p,role_of(p,sibs)] for p in s],'legal':alive,'mirror':any(det(p[0])<0 for p in s)} for s,(sibs,alive,ks) in zip(S,info)]
    G['second']=[]
    sheet(f'out/{rule}/shells.png',[([(ROOT,'root')]+[(p,role_of(p,sibs)) for p in s],arrows,explode(ROOT_C)) for s,(sibs,alive,ks) in zip(S,info)],
          7,1.4,border=['#27ae60' if alive else '#c0392b' for sibs,alive,ks in info])
    rows=[]
    for s,(sibs,alive,ks) in zip(S,info):
        if not alive: continue
        forced=None; union=set()
        for c in s:
            comps=[frozenset(x) for x in shells(c,[ROOT]+s,cont,arrows)]
            f=frozenset.intersection(*comps)-{ROOT}-set(s); u=frozenset.union(*comps)-{ROOT}-set(s)
            forced=(f if forced is None else forced|f); union|=u
        r={}
        for p in union:  r.update({c:'#d5f5e3' for c in cells(*p)})
        for p in forced: r.update({c:'#27ae60' for c in cells(*p)})
        for p in s:      r.update({c:(ROLE['mirror'] if det(p[0])<0 else '#bfd7ea' if p in sibs else '#d0d0d0') for c in cells(*p)})
        r.update({c:ROLE['root'] for c in RC}); rows.append(r)
        G['second'].append({'tiles':[[ROOT,'root']]+[[p,role_of(p,sibs)] for p in s]+[[p,'forced'] for p in forced]+[[p,'possible'] for p in union-forced],'legal':True,'mirror':False})
    slices(f'out/{rule}/second_layer.png',rows)
print('surviving neighbourhoods identical under both rules:',survivors['arrows']==survivors['colours'])

# ---------------- interactive HTML gallery 
import json
def J(o): return json.dumps(o,separators=(',',':'))
ser=lambda g:{k:[{'t':[[list(map(list,M)),list(t),r] for (M,t),r in it['tiles']],'l':bool(it['legal']),'m':bool(it['mirror'])} for it in v] for k,v in g.items()}
HTML=r"""<!DOCTYPE html><html><head><meta charset="utf-8"><title>Chair enumeration</title>
<style>body{font-family:sans-serif;margin:0;background:#faf8f3;color:#222}#top{padding:8px 12px;border-bottom:1px solid #ddd;background:#fff;position:sticky;top:0;display:flex;gap:14px;align-items:center;flex-wrap:wrap}
#grid{display:flex;flex-wrap:wrap;gap:6px;padding:10px}#grid canvas{border:3px solid #ddd;border-radius:4px;cursor:pointer;background:#fff}
#grid canvas.legal{border-color:#27ae60}#grid canvas.dead{border-color:#c0392b}#grid canvas.mirror{box-shadow:0 0 0 2px #c71585 inset}
#viewer{position:fixed;inset:0;background:rgba(0,0,0,.55);display:none;align-items:center;justify-content:center}#box{background:#fff;padding:10px;border-radius:8px;display:flex;gap:12px}
#big{cursor:grab;background:#fff;border:1px solid #ccc}#side{width:220px;font-size:13px;line-height:1.5}label{display:block}button,select{font:inherit}</style></head><body>
<div id="top"><b>Chair enumeration</b> rule <select id="rule"></select> set <select id="cat"></select> <span id="count"></span>
<span style="font-size:12px;color:#555">red = the chair &nbsp; blue = siblings in its superchair &nbsp; grey = other neighbours &nbsp; magenta = mirror image &nbsp; green = forced second layer &nbsp; light green = possible second layer &nbsp; border green/red = legal / fails second layer &nbsp; click a thumbnail to rotate it</span></div>
<div id="grid"></div>
<div id="viewer"><div id="box"><canvas id="big" width="620" height="620"></canvas><div id="side"><div id="info"></div><label><input type="checkbox" id="explode"> explode</label>
<label><input type="checkbox" id="showShell" checked> show neighbours</label><label><input type="checkbox" id="showForced" checked> show forced second layer</label><label><input type="checkbox" id="showPossible"> show possible second layer</label>
<p>drag to rotate, wheel to zoom</p><button id="close">close</button></div></div></div>
<script>
const DATA=%%DATA%%;const FACES=%%FACES%%;const CELLS=%%CELLS%%;
const mv=(M,v)=>[0,1,2].map(i=>M[i][0]*v[0]+M[i][1]*v[1]+M[i][2]*v[2]);
const cellsOf=(M,t)=>CELLS.map(c=>mv(M,[2*c[0]+1,2*c[1]+1,2*c[2]+1]).map((x,i)=>(x+2*t[i]-1)/2));
const facesOf=(M,t)=>FACES.map(([c,n,k,d])=>({ctr:mv(M,[2*c[0]+1+n[0],2*c[1]+1+n[1],2*c[2]+1+n[2]]).map((x,i)=>(x+2*t[i])/2),n:mv(M,n),k,d:mv(M,d)}));
const ROLE={root:'#c0392b',sib:'#2471a3',other:'#7f8c8d',mirror:'#c71585',forced:'#1e8449',possible:'#a9dfbf'};const FC={B:'#2b2b2b',W:'#ffffff',G:'#9e9e9e'};
function basisOf(n){const i=n.findIndex(x=>x!==0),a=[0,1,2].filter(j=>j!==i);const u=[0,0,0],v=[0,0,0];u[a[0]]=1;v[a[1]]=1;return [u,v];}
function prep(item){ // precompute tiles: faces, explode direction
  const tiles=item.t.map(([M,t,role])=>({faces:facesOf(M,t),role,cen:cellsOf(M,t).reduce((s,c)=>s.map((x,i)=>x+c[i]/7+0.5/7),[0,0,0])}));
  const r0=tiles[0].cen;tiles.forEach(T=>{const d=T.cen.map((x,i)=>x-r0[i]);const L=Math.hypot(...d)||1;T.off=d.map(x=>Math.round(2*x/L)/2);});return tiles;}
function render(cv,tiles,view,arrows){
  const ctx=cv.getContext('2d');const W=cv.width,H=cv.height;ctx.clearRect(0,0,W,H);
  const e=view.el,a=view.az;const V=[Math.cos(e)*Math.cos(a),Math.cos(e)*Math.sin(a),Math.sin(e)],R=[-Math.sin(a),Math.cos(a),0];
  const U=[V[1]*R[2]-V[2]*R[1],V[2]*R[0]-V[0]*R[2],V[0]*R[1]-V[1]*R[0]];
  const dot=(p,q)=>p[0]*q[0]+p[1]*q[1]+p[2]*q[2];const cen=[1,1,1];
  const P=p=>[W/2+view.z*dot([p[0]-cen[0],p[1]-cen[1],p[2]-cen[2]],R),H/2-view.z*dot([p[0]-cen[0],p[1]-cen[1],p[2]-cen[2]],U)];
  const polys=[];
  for(const T of tiles){ if(!view.show(T.role))continue;const off=view.ex?T.off.map(x=>x*1.0):[0,0,0];
    for(const f of T.faces){ if(dot(f.n,V)<=0)continue;const [u,v]=basisOf(f.n);const c=f.ctr.map((x,i)=>x+off[i]);
      const q=[[-.5,-.5],[.5,-.5],[.5,.5],[-.5,.5]].map(([s,r])=>P(c.map((x,i)=>x+s*u[i]+r*v[i])));
      let arrow=null;if(arrows){const L=Math.hypot(...f.d),d=f.d.map(x=>x/L),ee=[f.n[1]*d[2]-f.n[2]*d[1],f.n[2]*d[0]-f.n[0]*d[2],f.n[0]*d[1]-f.n[1]*d[0]];
        const pt=k=>P(c.map((x,i)=>x+k[0]*d[i]+k[1]*ee[i]));arrow={tail:pt([-.3,0]),tip:pt([.32,0]),h1:pt([.14,.14]),h2:pt([.14,-.14]),col:{B:'#111',W:'#fff',G:'#2a6fb0'}[f.k]};}
      polys.push({depth:dot(c,V),q,fill:arrows?'#e9e2d3':FC[f.k],edge:ROLE[T.role],arrow});}}
  polys.sort((p,q)=>p.depth-q.depth);
  for(const p of polys){ctx.beginPath();p.q.forEach(([x,y],i)=>i?ctx.lineTo(x,y):ctx.moveTo(x,y));ctx.closePath();ctx.fillStyle=p.fill;ctx.fill();ctx.strokeStyle=p.edge;ctx.lineWidth=Math.max(.6,view.z/60);ctx.stroke();
    if(p.arrow){const A=p.arrow;ctx.lineWidth=Math.max(1,view.z/40);if(A.col==='#fff'){ctx.strokeStyle='#333';ctx.lineWidth*=2;ctx.beginPath();ctx.moveTo(...A.tail);ctx.lineTo(...A.tip);ctx.stroke();ctx.lineWidth/=2;}
      ctx.strokeStyle=A.col;ctx.beginPath();ctx.moveTo(...A.tail);ctx.lineTo(...A.tip);ctx.stroke();ctx.beginPath();ctx.moveTo(...A.tip);ctx.lineTo(...A.h1);ctx.lineTo(...A.h2);ctx.closePath();ctx.fillStyle=A.col;ctx.fill();ctx.strokeStyle=A.col==='#fff'?'#333':A.col;ctx.lineWidth=.6;ctx.stroke();}}}
const DEF={el:Math.PI/6,az:Math.PI*40/180};
const ruleSel=document.getElementById('rule'),catSel=document.getElementById('cat'),grid=document.getElementById('grid');
for(const r of Object.keys(DATA))ruleSel.add(new Option(r,r));
const CATS={contacts:'two-chair contacts',shells:'complete neighbourhoods',second:'second layer (surviving neighbourhoods)'};for(const c in CATS)catSel.add(new Option(CATS[c],c));
let cur=null;
function fill(){const rule=ruleSel.value,cat=catSel.value,items=DATA[rule][cat];grid.innerHTML='';document.getElementById('count').textContent=items.length+' items';
  items.forEach((it,i)=>{const cv=document.createElement('canvas');cv.width=cv.height=110;cv.title='#'+(i+1)+(it.l?' legal':' illegal')+(it.m?' (mirror image present)':'');
    cv.className=(it.l?'legal':'dead')+(it.m?' mirror':'');const tiles=prep(it);const isShell=cat!=='contacts';
    render(cv,tiles,{...DEF,z:isShell?11:16,ex:isShell,show:r=>r!=='possible'},rule==='arrows');
    cv.onclick=()=>open(it,i,rule,cat);grid.appendChild(cv);});}
const big=document.getElementById('big'),viewer=document.getElementById('viewer');let bv={...DEF,z:70,ex:false};
function open(it,i,rule,cat){cur={tiles:prep(it),arrows:rule==='arrows'};bv={...DEF,z:cat==='contacts'?90:60,ex:cat!=='contacts'};document.getElementById('explode').checked=bv.ex;
  document.getElementById('info').innerHTML='<b>'+rule+' / '+CATS[cat]+' #'+(i+1)+'</b><br>'+(it.l?'legal':'illegal')+(it.m?'<br>contains a mirror-image chair':'')+'<br>'+it.t.length+' chairs';viewer.style.display='flex';draw();}
function draw(){if(!cur)return;const sh=document.getElementById('showShell').checked,fo=document.getElementById('showForced').checked,po=document.getElementById('showPossible').checked;
  render(big,cur.tiles,{...bv,ex:document.getElementById('explode').checked,show:r=>r==='root'||(['sib','other','mirror'].includes(r)?sh:r==='forced'?fo:po)},cur.arrows);}
let drag=null;big.onmousedown=e=>{drag=[e.clientX,e.clientY];};window.onmousemove=e=>{if(!drag)return;bv.az-=(e.clientX-drag[0])*.01;bv.el=Math.max(-1.5,Math.min(1.5,bv.el+(e.clientY-drag[1])*.01));drag=[e.clientX,e.clientY];draw();};
window.onmouseup=()=>drag=null;big.onwheel=e=>{e.preventDefault();bv.z*=Math.exp(-e.deltaY*.001);draw();};
for(const id of ['explode','showShell','showForced','showPossible'])document.getElementById(id).onchange=draw;
document.getElementById('close').onclick=()=>viewer.style.display='none';ruleSel.onchange=catSel.onchange=fill;fill();
</script></body></html>"""
open('out/gallery.html','w').write(HTML.replace('%%DATA%%',J({r:ser(g) for r,g in GAL.items()})).replace('%%FACES%%',J([list(map(list,f[:2]))+[f[2],list(f[3])] for f in FACES])).replace('%%CELLS%%',J([list(c) for c in CELLS])))
print('wrote out/gallery.html')
