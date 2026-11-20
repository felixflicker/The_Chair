# growth.py: three growth experiments on the three-colour Chair (Chair611), each drawn as an interactive page.
#
# Requires Python 3.8 or later and numpy. Run as
#
#   python3 growth.py MODE [key=value ...]
#
# MODE:
#   base      Forced growth from the base of a level-L superchair. The chairs resting on one of its three square
#             faces are set on the plane z=0. A chair is then added wherever it is the only legal chair covering
#             an empty cube next to the cluster, staying at z>=0, until no cube is forced. Legal means that it
#             overlaps nothing and that black faces meet white and blue faces meet blue.
#   skeleton  Growth from a single chair in 3D by the R=1 rule: a candidate must be legal and must leave every
#             empty cube it touches coverable by some other legal chair. A cube with exactly one candidate is
#             forced; if no cube is forced, a random candidate is placed (a guess).
#   anneal    The chair as a lattice model. Each cube of a level-L superchair holds one of 168 states: one of
#             the 7 cubes of a chair in one of its 24 orientations. Two neighbouring cubes cost energy 1 unless
#             they continue the same chair or show each other legal colours. The cubes just outside the
#             superchair are frozen as in the perfect tiling, the inside starts random, and it is cooled
#             through nT temperatures from Thi to Tlo by heat-bath updates of single cubes.
#
# Options (defaults)
#   L=2         level of the superchair (base, anneal)
#   seed=1      random seed
#   face=z      base: the square face the superchair rests on (x, y or z)
#   faces=1     base: least number of square faces two touching chairs must share (1 is the colour rule alone)
#   n=600       skeleton: number of chairs to place
#   Thi=0.5     anneal: first temperature
#   Tlo=0.4     anneal: last temperature
#   nT=21       anneal: number of temperatures, equally spaced
#   sweeps=1000 anneal: sweeps at each temperature (a sweep updates every cube once, in two interleaved halves)
#   quench=300  anneal: final sweeps at T=0.1, which strip the thermal noise before the result is drawn
#
# Output, in the working directory; progress is shown on one line of the terminal while a run is under way
#   NAME.txt   the line printed at the end. For anneal it also holds the table E(T): temperature, energy
#              per cube averaged over the second half of the sweeps, and fraction of cubes in the correct state.
#   NAME.dat   final state. base, skeleton: one chair per line, orientation 1-48 (FRAMES of enumerate_figs.py)
#              then x y z. anneal: one cube per line, x y z, orientation, and which cube of the chair (1-7).
#   NAME.html  the result (needs a connection, to load three.js). Drag to rotate, wheel to zoom.
#              burst moves the pieces apart at their own size, level by level; the menu beside it chooses the
#              pieces: chairs, supertiles of one level or, for anneal, single Wang cubes. Zoom out to follow it.
#              colour by level shows the largest complete supertile each chair belongs to.
#              base, skeleton: play or step through the growth; forced cubes are marked green.
#              skeleton: "difference from the perfect tiling" shows the chairs of the tiling that are absent
#              and enclosed by the structure, coloured by their place in the hierarchy.
#              anneal: the same toggle shows the cubes in the wrong state in red. The cubes are coloured as
#              in the Wang-cube figure: black, white and blue outside; inside, one violet-to-rose hue for each
#              orientation, strong where the chair's corner cube joins its three neighbours, pale elsewhere.
#              Keys: space play, arrows step, C difference, H hide the panel.
#
# Run with no mode, or with a mistyped mode or option, and this text is printed.
import sys,json,random
from collections import Counter
from itertools import permutations,product
PS=[(p,s) for p in permutations(range(3)) for s in product((1,-1),repeat=3)]
R=[a for a,(p,s) in enumerate(PS) if s[0]*s[1]*s[2]*(3-2*((p[1]-p[0])%3))>0]
MUL=[[PS.index((tuple(q[j] for j in p),tuple(s[i]*u[p[i]] for i in range(3)))) for q,u in PS] for p,s in PS]
INV=[m.index(0) for m in MUL]
CELLS=[c for c in product((0,1),repeat=3) if sum(c)<3]
DIRS=[(1,0,0),(-1,0,0),(0,1,0),(0,-1,0),(0,0,1),(0,0,-1)]
COLS=['010203','020130','013002','231020','300201','201310','102033']
KID=list(zip([0,17,10,35,44,29,6,0],[tuple(4*x for x in c) for c in CELLS]+[(1,1,1)]))
add=lambda u,v:tuple(x+y for x,y in zip(u,v))
mv=lambda a,v:tuple(s*v[j] for j,s in zip(*PS[a]))
comp=lambda P,Q:(MUL[P[0]][Q[0]],add(mv(P[0],Q[1]),P[1]))
inv=lambda P:(INV[P[0]],mv(INV[P[0]],tuple(-x for x in P[1])))
kid=lambda i,k:(KID[i][0],tuple(2**(k-1)*x for x in KID[i][1]))
FT=[[(tuple((x-1)//2 for x in mv(a,tuple(2*x+1 for x in c))),[int(k[DIRS.index(mv(INV[a],d))]) for d in DIRS]) for c,k in zip(CELLS,COLS)] for a in range(48)]
K=1024
E=lambda v:v[0]+K*v[1]+K*K*v[2]
dec=lambda n:tuple((n+E((512,512,512)))//K**i%K-512 for i in range(3))
O=[[E(o) for o,_ in f] for f in FT]
FA=[[(E(add(o,DIRS[d]))*6+(d^1),3-k[d] or 3,E(o)*6+d,k[d]) for o,k in f for d in range(6) if k[d]] for f in FT]
NB=[sorted({u//6 for u,*_ in f}) for f in FA]
B4,B6,B9=[[E(d) for d in product(range(a,b),repeat=3)] for a,b in ((-4,4),(-6,6),(-4,5))]
PL=[(a,o) for a in R for o in O[a]]
BAD={}
def tick(*a): print(*a,' '*30,end='\r',file=sys.stderr,flush=True)
def sup(k,P,p=''): return [x for i in range(8) for x in sup(k-1,comp(P,kid(i,k)),p+str(i))] if k else [(P,p)]
def grow(seed,n,rng,look=0,m=1,floor=-K**3):
    lo=-K**3; occ,col,cov,adm,frc,memo,kp,front,ch,S={},{},{},{},{},{},{},set(),[],[]
    legal=lambda a,t:not any(t+o in occ or t+o<lo for o in O[a]) and all(col.get(t*6+u,k)==k for u,k,w,c in FA[a]) and (m<2 or all(v>=m for v in Counter(occ[t+u//6] for u,*_ in FA[a] if t+u//6 in occ).values()))
    def keep(a,d):
        if (a,d) not in kp: kp[a,d]=sum(1<<j for j,(b,o) in enumerate(PL) if (b,d-o) not in bad[a])
        return kp[a,d]
    def cv(q):
        if q not in cov: cov[q]=sum(1<<j for j,(a,o) in enumerate(PL) if legal(a,q-o))
        return cov[q]
    def put(a,t,f):
        for o in O[a]: occ[t+o]=len(ch); front.discard(t+o); cov.pop(t+o,0)
        for u,k,w,c in FA[a]: col[t*6+w]=c
        ch.append((a,t,f)); front.update(q for q in (t+x for x in NB[a]) if q not in occ and q>=lo)
        for d in B4:
            if t+d in cov: cov[t+d]&=keep(a,d)
        return [q for q in (t+d for d in B6) if q in front or q in adm]
    def ok(p):
        if p not in memo: a,t=p; memo[p]=all(t+q in occ or cv(t+q)&keep(a,q) for q in NB[a])
        return memo[p]
    bad=BAD.setdefault(m,{a:() for a in R})
    for a in R:
        if not bad[a]:
            put(a,0,0); bad[a]={(b,d) for b in R for d in B9 if not legal(b,d)}
            for x in occ,col,cov,front,ch: x.clear()
    lo=floor; dirty=[q for a,t in seed for q in put(a,E(t),0)]
    while len(ch)-len(seed)<n:
        for c in dirty:
            if c in front: adm[c]=[p for p in ((a,c-o) for j,(a,o) in enumerate(PL) if cv(c)>>j&1) if not look or ok(p)]
            else: adm.pop(c,0)
            if c in front and len(adm[c])==1: frc.setdefault(c,len(ch)-len(seed))
            elif c in frc: S.append((c,frc.pop(c),len(ch)-len(seed)))
        if (len(ch)-len(seed))%100==0: tick(len(ch)-len(seed),'placed')
        if frc: p,f=adm[rng.choice(list(frc))][0],1
        elif look and (pool:=sorted({p for v in adm.values() for p in v})): p,f=rng.choice(pool),2
        else: break
        memo.clear(); dirty=put(*p,f)
    return [(a,dec(t),f) for a,t,f in ch],[(*dec(c),s,e) for c,s,e in S+[(c,s,n+1) for c,s in frc.items()]]
def inside(k,P,g):
    l=[(x-1)//2 for x in mv(INV[P[0]],tuple(2*g[i]+1-2*P[1][i] for i in range(3)))]
    return all(0<=x<2<<k for x in l) and min(l)<1<<k
def at(k,P,g,p=''):
    while k: i,P=next((i,Q) for i in range(8) for Q in [comp(P,kid(i,k))] if inside(k-1,Q,g)); p+=str(i); k-=1
    return P,p
def nest(ch):
    def fit(k,Q,I,bad):
        out=[]; B=[(min(x,x+y)-1,max(x,x+y)+1) for x,y in zip(Q[1],mv(Q[0],(2<<k,)*3))]
        for i in I:
            a,t,_=ch[i]
            if any(x<l or x>h for x,(l,h) in zip(t,B)): out.append(i); continue
            c=[add(t,o) for o,_ in FT[a]]; n=sum(inside(k,Q,g) for g in c)
            if n==0: out.append(i)
            elif n<7 or at(k,Q,c[0])[0]!=(a,t): bad=bad+[i]
        return len(bad),len(out),k,Q,bad,out
    live=[fit(0,ch[0][:2],range(len(ch)),[])]
    while live[0][1] and live[0][2]<16:
        tick('fit level',live[0][2]+1)
        new=sorted(fit(k+1,comp(Q,inv(kid(r,k+1))),out,bad) for *_,k,Q,bad,out in live for r in range(8))
        live=[x for x in new if x[0]==new[0][0]][:8]
    return live[0][2:5]
def voids(ch,k,Q):
    occ={add(t,o) for a,t,_ in ch for o,_ in FT[a]}; S,V,done={},{},set()
    for g in occ:
        for i in range(3): S.setdefault((i,g[i-1],g[i-2]),[]).append(g[i])
    I=lambda g:g not in occ and all(min(v)<g[i]<max(v) for i in range(3) for v in [S.get((i,g[i-1],g[i-2]),[g[i]])])
    for g in product(*(range(min(v),max(v)) for v in zip(*occ))):
        if g not in done and I(g) and inside(k,Q,g):
            c,p=at(k,Q,g); cs=[add(c[1],o) for o,_ in FT[c[0]]]; done.update(cs)
            if all(map(I,cs)): V[c]=p
    return V
TT=[[[*o,*k] for o,k in f] for f in FT]+[[[0,0,0]+[k[d] or 4+2*(R.index(a) if a in R else 0)+(i>0 and any(add(CELLS[i],mv(INV[a],DIRS[d])))) for d in range(6)]] for a in range(48) for i,(_,k) in enumerate(FT[a])]
HTML=r'''<!DOCTYPE html><html><head><meta charset="utf-8"><title>__T__</title><script src="https://cdnjs.cloudflare.com/ajax/libs/three.js/r128/three.min.js"></script>
<style>body{margin:0;background:#000;color:#ddd;font:14px sans-serif;overflow:hidden}#ui{position:absolute;left:10px;top:10px;background:rgba(0,0,0,.6);padding:8px 12px;border-radius:6px}input,select{vertical-align:middle}</style></head><body><div id="ui"><b>__T__</b><br>
<span id="pt"><button id="play">play</button> <button id="bk">&lt;</button> <button id="fw">&gt;</button> speed <input id="sp" type="range" min="1" max="40" value="8" style="width:80px"> <input id="sl" type="range" min="0" value="0" style="width:300px"> <span id="lb"></span> &nbsp; <label><input id="hl" type="checkbox" checked> forced cubes (green), last chair (white; magenta if a guess)</label><br></span>
<span id="pv"><label><input id="df" type="checkbox"> difference from the perfect tiling:</label> <span id="pk">absent <label><input class="k" type="checkbox" checked> <span style="color:#ff8c33">central chairs</span></label> <label><input class="k" type="checkbox" checked> <span style="color:#33aaff">chairs opposite the notch</span></label> <label><input class="k" type="checkbox" checked> <span style="color:#55dd55">central superchairs</span></label> <label><input class="k" type="checkbox" checked> <span style="color:#dddd33">other</span></label>, </span><span style="color:#ff3333">wrong pieces</span> &nbsp; faint <input id="fa" type="range" min="0" max="100" value="12" style="width:80px"><br></span>
colour by <select id="cm"><option>faces<option>level</select> &nbsp; burst <input id="ex" type="range" min="0" max="100" value="0" style="width:120px"> into <select id="kp"></select></div>
<script>
const D=__D__,$=i=>document.getElementById(i),N=D.N,Q=D.q||0,DIR=[[1,0,0],[-1,0,0],[0,1,0],[0,-1,0],[0,0,1],[0,0,-1]],SH=[.85,.65,.75,.55,1,.45],FC=Q?[0,0x141414,0xf8f8f8,0x378add]:[0,0x2b2b2b,0xffffff,0x2a6fb0],KC=[0xff8c33,0x33aaff,0x55dd55,0xdddd33],LC=[0x777777,0x44bb55,0x3377ee,0xee8833,0xcc3399,0x33bbbb,0xdddd33,0xff3333],
hsl=(h,s,l)=>new THREE.Color().setHSL(h,s,l).getHex(),MT={},mat=o=>MT[o]=MT[o]||new THREE.ShaderMaterial({uniforms:{o:{value:o}},side:2,transparent:o<1,depthWrite:o==1,extensions:{derivatives:true},vertexShader:'attribute vec3 kc;attribute vec2 ke;varying vec3 c;varying vec2 v;void main(){c=kc;v=ke;gl_Position=projectionMatrix*modelViewMatrix*vec4(position,1.);}',fragmentShader:'uniform float o;varying vec3 c;varying vec2 v;void main(){float d=max(abs(v.x),abs(v.y)),t=smoothstep(.93-fwidth(d),.93,d);gl_FragColor=vec4(mix(c,max(c.r,max(c.g,c.b))<.3?c+.2:vec3(0.),t),o);}'}),ck=q=>q[0]+512+1024*(q[1]+512+1024*(q[2]+512)),
mk=(c,i)=>({i,a:c[0],s:c[4],f:c[5],p:c[6],w:c[6]==null||c[5]==3,q:D.T[c[0]].map(r=>[r[0]+c[1],r[1]+c[2],r[2]+c[3]])}),C=D.C.map(mk),V=D.V.map((c,i)=>mk([...c.slice(0,4),0,0,c[4]],i)),
key=(c,k)=>c.p==null?'#'+c.i:c.p.slice(0,Math.max(0,c.p.length-k)),cl=p=>p.slice(-1)=='7'?0:p.slice(-1)=='0'?1:p.slice(-2,-1)=='7'?2:3,K=Math.max(1,...D.C.map(c=>c[6]==null?0:c[6].length)),U={},cen=[0,0,0],
sc=new THREE.Scene(),cam=new THREE.PerspectiveCamera(40,1,.1,1e4),ren=new THREE.WebGLRenderer({antialias:true}),fit=()=>{cam.aspect=innerWidth/innerHeight;cam.updateProjectionMatrix();ren.setSize(innerWidth,innerHeight)};
let cur=N,play=0,acc=0,dirty=1,ms=[],th=.7,ph=1.05,R=0,dr=0;
for(const c of C.concat(V))c.m=[0,1,2].map(j=>c.q.reduce((s,q)=>s+q[j]+.5,0)/c.q.length);
for(const c of C){c.m.forEach((v,j)=>cen[j]+=v/C.length);for(let k=0;k<=K;k++){const u=U[k+key(c,k)]=U[k+key(c,k)]||{n:0,t:0,s:0,m:[0,0,0]};u.n+=!c.w;u.t++;u.s=Math.max(u.s,c.s);c.m.forEach((v,j)=>u.m[j]+=v)}}
for(const c of C)for(const q of c.q)R=Math.max(R,...q.map((v,j)=>3.2*Math.abs(v+.5-cen[j])+6));
$('sl').max=N;for(let k=-Q;k<K-Q;k++)$('kp').add(new Option(k<0?'Wang cubes':k?'level-'+k+' supertiles':'chairs',k));$('kp').value=0;$('df').checked=!!D.df;$('pt').hidden=!N;$('pk').hidden=!D.V.length;$('pv').hidden=!D.V.length&&C.every(c=>!c.w);
ren.setPixelRatio(devicePixelRatio);document.body.appendChild(ren.domElement);fit();onresize=fit;
function draw(){for(const m of ms){sc.remove(m);m.geometry.dispose()}ms=[];
const df=$('df').checked,e=$('ex').value/100,kp=Math.min(K,Math.max(0,+$('kp').value+Q)),cm=$('cm').value,on=[...document.querySelectorAll('.k')].map(b=>b.checked);
if(df)cur=N;
const vis=C.filter(c=>c.s<=cur),hl=$('hl').checked&&!df&&N,
uc=(c,k)=>{const u=U[k+key(c,k)];return u?u.m.map(v=>v/u.t):c.m},
off=c=>{const d=[0,0,0];if(e)for(let k=kp;k<=K;k++){const a=uc(c,k),b=k<K?uc(c,k+1):cen;for(let j=0;j<3;j++)d[j]+=e*(1+.3*Math.max(0,k-Q))*(a[j]-b[j])}return d},
lev=c=>{let l=0;for(let k=1;k<=K;k++){const u=U[k+key(c,k)];if(u.n==(Q?7*8**(k-1):8**k)&&u.s<=cur)l=k}return Math.min(7,Math.max(0,l-Q))},
mesh=(L,col,op=1,h=.5)=>{const P=[],A=[],B=[],own=new Set();
if(!e&&h==.5)for(const c of L)for(const q of c.q)own.add(ck(q));
for(const c of L){const o=off(c),k0=col(c);c.q.forEach((q,i)=>{const r=c.r||D.T[c.a][i];
for(let d=0;d<6;d++){const n=DIR[d],x=r[3+d],a=d>>1;if(!x||own.has(ck(q.map((v,j)=>v+n[j]))))continue;
const k=new THREE.Color(k0<0?x<4?FC[x]:hsl((285+(x-4>>1)/24*80)%360/360,x&1?.7:.8,x&1?.8:.58):k0).multiplyScalar(SH[d]);
for(const[s,t] of [[1,1],[-1,1],[-1,-1],[1,1],[-1,-1],[1,-1]]){const p=q.map((v,j)=>v+.5+o[j]+n[j]*h-cen[j]);p[(a+1)%3]+=s*h;p[(a+2)%3]+=t*h;P.push(p[0],p[2],p[1]);A.push(k.r,k.g,k.b);B.push(s,t)}}})}
const g=new THREE.BufferGeometry(),m=new THREE.Mesh(g,mat(op));
g.setAttribute('position',new THREE.Float32BufferAttribute(P,3));g.setAttribute('kc',new THREE.Float32BufferAttribute(A,3));g.setAttribute('ke',new THREE.Float32BufferAttribute(B,2));m.frustumCulled=false;sc.add(m);ms.push(m)};
if(df){mesh(V.filter(v=>on[cl(v.p)]),v=>KC[cl(v.p)]);mesh(vis.filter(c=>c.w),()=>0xff3333);mesh(vis.filter(c=>!c.w),()=>0x999999,$('fa').value/100)}
else mesh(vis,c=>hl&&cur&&c.s==cur?c.f==2?0xff44ff:0xffffff:cm=='faces'?-1:LC[lev(c)]);
if(hl&&!e)mesh(D.S.filter(s=>s[3]<=cur&&cur<s[4]).map(s=>({q:[s],r:[0,0,0,1,1,1,1,1,1],m:cen})),()=>0x33dd55,.85,.25);
$('sl').value=cur;$('lb').textContent=cur+' / '+N}
addEventListener('input',e=>{if(e.target.id=='sl')cur=+e.target.value;dirty=1});
$('play').onclick=()=>{play=!play;if(play&&cur==N)cur=0;dirty=1;$('play').textContent=play?'pause':'play'};$('bk').onclick=()=>{cur=Math.max(0,cur-1);dirty=1};$('fw').onclick=()=>{cur=Math.min(N,cur+1);dirty=1};
onkeydown=e=>{const k=e.key;if(/INPUT|SELECT/.test(e.target.tagName))return;if(k==' '){e.preventDefault();$('play').click()}if(k=='h')$('ui').hidden=!$('ui').hidden;if(k=='c'){$('df').checked=!$('df').checked;dirty=1}if(k=='ArrowRight')$('fw').click();if(k=='ArrowLeft')$('bk').click()};
onmousedown=e=>{dr=e.target.tagName=='CANVAS'&&[e.clientX,e.clientY]};onmouseup=()=>dr=0;
onmousemove=e=>{if(dr){th+=(e.clientX-dr[0])*.01;ph=Math.min(3,Math.max(.1,ph-(e.clientY-dr[1])*.01));dr=[e.clientX,e.clientY]}};
onwheel=e=>R*=Math.exp(e.deltaY*.001);
(function loop(){requestAnimationFrame(loop);if(play){acc+=$('sp').value**2/64*Math.max(1,N/2000);const n=Math.floor(acc);acc-=n;cur=Math.min(N,cur+n);dirty=dirty||n;if(cur==N)$('play').click()}
if(dirty){dirty=0;draw()}cam.position.set(R*Math.sin(ph)*Math.cos(th),R*Math.cos(ph),R*Math.sin(ph)*Math.sin(th));cam.lookAt(0,0,0);ren.render(sc,cam)})();
</script></body></html>'''
def page(name,title,C,V=[],S=[],**kw): open(name+'.html','w').write(HTML.replace('__T__',title).replace('__D__',json.dumps(dict(T=TT,C=C,V=V,S=S,N=max(c[4] for c in C),**kw),separators=(',',':'))))
def say(name,msg,m='w'): print(msg,' '*12,flush=True); open(name+'.txt',m).write(msg+'\n')
def show(name,title,label,note,ch,S):
    k,Q,bad=nest(ch); tick('absent chairs'); V=voids(ch,k,Q); tick('page'); n0=sum(not c[2] for c in ch)
    cl=Counter(0 if p[-1:]=='7' else 1 if p[-1:]=='0' else 2 if p[-2:-1]=='7' else 3 for p in V.values())
    open(name+'.dat','w').write(''.join(f'{a+1} {t[0]} {t[1]} {t[2]}\n' for a,t,_ in ch))
    page(name,title,[[a,*t,max(0,i-n0+1),f,None if i in bad else at(k,Q,add(t,FT[a][0][0]))[1]] for i,(a,t,f) in enumerate(ch)],[[a,*t,p] for (a,t),p in V.items()],S)
    say(name,f'{label}: {len(ch)-n0} placed on {n0}, {note}, {len(bad)} wrong'+f', absent {cl[0]} central, {cl[1]} notch, {cl[2]} central-superchair, {cl[3]} other'*bool(V)+f' -> {name}.html')
def anneal(name,L,Thi,Tlo,nT,sweeps,quench,seed):
    import numpy as np
    g=np.random.default_rng(seed); h=2**L; ix={a:7*r for r,a in enumerate(R)}
    cells=[(add(t,o),ix[a]+i,p+str(i)) for (a,t),p in sup(L,(0,(h,h,h))) for i,(o,_) in enumerate(FT[a])]
    idx={c[0]:n for n,c in enumerate(cells)}; N=len(cells); wall={}
    nb=np.array([[idx[q] if q in idx else wall.setdefault(q,N+len(wall)) for q in (add(c[0],d) for d in DIRS)] for c in cells])
    def frozen(q): a,t=at(L+2,(0,(0,0,0)),q)[0]; return ix[a]+[add(t,o) for o,_ in FT[a]].index(q)
    def cost(s,t,d):
        a,i,b,j=R[s//7],s%7,R[t//7],t%7; x,y=FT[a][i][1][d],FT[b][j][1][d^1]
        return not(x==0 and a==b and add(FT[a][i][0],DIRS[d])==FT[a][j][0] or x*y and (x+y)%3==0)
    ref=np.array([c[1] for c in cells]); spin=np.concatenate([g.integers(168,size=N),[frozen(q) for q in wall]])
    C=np.array([[[cost(s,t,d) for s in range(168)] for t in range(168)] for d in range(6)],np.int8)
    P=[(C[2*k][:,None,:]+C[2*k+1][None,:,:]).reshape(-1,168) for k in range(3)]
    half=[np.flatnonzero(np.array([sum(c[0])%2 for c in cells])==p) for p in (0,1)]
    I,D=np.nonzero(nb>np.arange(N)[:,None]); J=nb[I,D]
    energy=lambda:int(C[D,spin[J],spin[I]].sum())
    def sweep(T):
        w=np.exp(-np.arange(7)/T)
        for p in half:
            s=spin[nb[p]]; c=np.cumsum(w[P[0][s[:,0]*168+s[:,1]]+P[1][s[:,2]*168+s[:,3]]+P[2][s[:,4]*168+s[:,5]]],1)
            spin[p]=(c<g.random((len(p),1))*c[:,-1:]).sum(1)
    say(name,'# T, energy per cube, fraction correct')
    for T in np.linspace(Thi,Tlo,nT):
        e=[]
        for s in range(sweeps):
            sweep(T)
            if s%100==0: tick(f'T={T:.4f} sweep {s}')
            if s>=sweeps//2 and s%10==0: e.append(energy())
        say(name,f'{T:.4f} {sum(e)/len(e)/N:.5f} {(spin[:N]==ref).mean():.4f}','a')
    for s in range(quench): sweep(.1)
    ok=spin[:N]==ref; title=f'Level-{L} superchair cooled from random cubes: energy {energy()}, {N-int(ok.sum())} of {N} cubes wrong'
    open(name+'.dat','w').write(''.join(f'{c[0][0]} {c[0][1]} {c[0][2]} {R[s//7]+1} {s%7+1}\n' for c,s in zip(cells,spin.tolist())))
    page(name,title,[[48+7*R[s//7]+s%7,*c[0],0,3*(s!=c[1]),c[2]] for c,s in zip(cells,spin.tolist())],q=1,df=int(not ok.all()))
    say(name,f'# anneal L={L}: energy {energy()}, {N-int(ok.sum())}/{N} wrong -> {name}.html','a')
if __name__=='__main__':
    A=dict(L=2,seed=1,face='z',faces=1,n=600,Thi=.5,Tlo=.4,nT=21,sweeps=1000,quench=300)
    try:
        mode=sys.argv[1]; A.update((k,type(A[k])(v)) for k,v in (x.split('=') for x in sys.argv[2:])); assert mode in ('base','skeleton','anneal') and A['face'] in ('x','y','z')
    except Exception: sys.exit(''.join(l[2:] or '\n' for l in open(__file__) if l[0]=='#'))
    L,sd,rng=A['L'],A['seed'],random.Random(A['seed'])
    if mode=='base':
        a=next(a for a in R if mv(INV[a],(0,0,-1))==tuple(-(i=='xyz'.index(A['face'])) for i in range(3)))
        unit=[comp((a,(0,0,0)),c) for c,_ in sup(L,(0,(0,0,0)))]; z=min(add(t,o)[2] for b,t in unit for o,_ in FT[b]); unit=[(b,(t[0],t[1],t[2]-z)) for b,t in unit]
        ch,S=grow([u for u in unit if min(add(u[1],o)[2] for o,_ in FT[u[0]])==0],10**9,rng,0,A['faces'],E((-512,-512,0)))
        show(f'base_L{L}{A["face"]}_f{A["faces"]}_s{sd}',f'Forced growth from the base of a level-{L} superchair',f'base L={L} {A["face"]}',f'superchair {"complete" if set(unit)<={c[:2] for c in ch} else "incomplete"}',ch,S)
    if mode=='skeleton': ch,S=grow([(0,(0,0,0))],A['n'],rng,1); show(f'skeleton_n{A["n"]}_s{sd}','Growth from a single chair by the R=1 rule',f'skeleton n={A["n"]}',f'{sum(c[2]==2 for c in ch)} guesses',ch,S)
    if mode=='anneal': anneal(f'anneal_L{L}_s{sd}',L,A['Thi'],A['Tlo'],A['nT'],A['sweeps'],A['quench'],sd)
