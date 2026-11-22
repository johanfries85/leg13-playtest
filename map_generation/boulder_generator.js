// DOM-free structured boulder generator. Also usable as a classic offline script.
(() => {
const defaults=Object.freeze({
  "deflection": 0.9,
  "density": 8000,
  "max-ratio": 2,
  "scale": 0.91,
  "short": 152,
  "contrast": 2.4,
  "isolation": 2,
  "separation": 10,
  "radius": 44,
  "removal": 19,
  "background": 60
});
const tau=2*Math.PI;
function createGenerator(radius,mapNumber,parameters={}){
 if(typeof radius!=='number'||!Number.isFinite(radius)||radius<=0)throw new RangeError('Radius must be a positive finite number of meters');
 if(!Number.isInteger(mapNumber)||mapNumber<1||mapNumber>4294967295)throw new RangeError('Map number must be an integer from 1 to 4294967295');
 const seed=mapNumber,settings={...defaults,...parameters};
 const setting=key=>settings[key];
 let deflected;
function rng(salt) {
  let state = (seed ^ salt) >>> 0;
  return () => {
    state = (state + 0x6d2b79f5) | 0;
    let bits = Math.imul(state ^ (state >>> 15),state | 1);
    bits ^= bits + Math.imul(bits ^ (bits >>> 7),bits | 61);
    return ((bits ^ (bits >>> 14)) >>> 0)/4294967296;
  };
}
function wavelengthRatio(maximalRatio, stretch) {
  return 1 + (maximalRatio - 1) * stretch ** 2;
}

function waves(ratio) {
  const random = rng(0x1e337a), scale = setting('scale');
  function wave(angle, wavelength, displacement) {
    const f = tau/wavelength;
    return {x:f*Math.cos(angle),y:f*Math.sin(angle),
      amplitude:displacement === undefined ? 1/f**2 : displacement/f,phase:tau*random()};
  }
  const long = Array.from({length:7},(_,i) =>
    wave(i*Math.PI/7+.36*random()-.18,350+230*random(),13+20*random()));
  const short = Array.from({length:16},() => wave(Math.PI*random(),175+115*random()));
  const power = list => list.reduce((sum,w) => sum+(w.amplitude*(w.x*w.x+w.y*w.y))**2,0);
  const gain = setting('short')/100*Math.sqrt(power(long)/power(short));
  return [...long,...(gain ? short.map(w=>({...w,amplitude:w.amplitude*gain})) : [])]
    .map(w=>({x:w.x/(scale*ratio),y:w.y/scale,amplitude:w.amplitude*scale,phase:w.phase}));
}
function inCircle(p) { return p.x*p.x+p.y*p.y <= radius*radius; }

function makeCandidates(field) {
  const step=8, side=2*Math.ceil(radius/step)+6, edge=side*step/2;
  const k=setting('deflection')*setting('scale'), gamma=setting('contrast');
  const pad=k*field.reduce((sum,w)=>sum+w.amplitude*Math.hypot(w.x,w.y),0);
  function land(x,y) {
    let dx=0,dy=0;
    for(const w of field) {
      const slope=w.amplitude*Math.sin(w.x*x+w.y*y+w.phase);
      dx+=slope*w.x;
      dy+=slope*w.y;
    }
    return {x:x+k*dx,y:y+k*dy};
  }
  const source=rng(0xc31729), extent=edge+pad;
  const n=Math.round(4*extent**2*setting('density')/1e6), reference=[];
  for(let i=0;i<n;i++) {
    const x=(2*source()-1)*extent,y=(2*source()-1)*extent,p=land(x,y);
    if(!inCircle(p)) continue;
    const col=Math.max(0,Math.min(2,Math.floor(3*(x+radius)/(2*radius))));
    deflected.push({...p,colorIndex:(y>0?0:3)+col});
    reference.push(p);
  }
  if(k===0||gamma===1) return reference;

  const bins=new Float64Array(side*side), probe=rng(0x21601d);
  const index=p=>Math.floor((p.y+edge)/step)*side+Math.floor((p.x+edge)/step);
  for(let y=-extent;y<extent;y+=5) for(let x=-extent;x<extent;x+=5) {
    const p=land(x+5*probe(),y+5*probe());
    if(Math.abs(p.x)<edge&&Math.abs(p.y)<edge) bins[index(p)]++;
  }
  const weights=[1,Math.exp(-.5),Math.exp(-2),Math.exp(-4.5)];
  function smooth(data,horizontal) {
    const stride=horizontal?1:side;
    return data.map((_,i)=>{
      const at=horizontal?i%side:Math.floor(i/side);
      let sum=0;
      for(let d=-3;d<=3;d++) if(at+d>=0&&at+d<side) sum+=weights[Math.abs(d)]*data[i+d*stride];
      return sum;
    });
  }
  const density=smooth(smooth(bins,true),false);
  for(let i=0;i<density.length;i++) {
    if(!inCircle({x:(i%side+.5)*step-edge,y:(Math.floor(i/side)+.5)*step-edge})) density[i]=0;
  }
  const weight=density.map(d=>d**gamma);
  const baseline=weight.reduce((a,b)=>a+b,0)/density.reduce((a,b)=>a+b,0);
  let total=0;
  const cumulative=weight.map((v,i)=>total+=Math.max(0,v-density[i]*baseline));
  const sample=rng(0x42a17b);
  return reference.map(p=>{
    const keep=sample(),pick=sample(),x=sample(),y=sample(),i=index(p);
    if(total===0||keep*density[i]*baseline<weight[i]) return p;
    let low=0,high=cumulative.length-1;
    while(low<high) {
      const mid=(low+high)>>>1;
      if(cumulative[mid]<=pick*total) low=mid+1;
      else high=mid;
    }
    return {x:(low%side+x)*step-edge,y:(Math.floor(low/side)+y)*step-edge};
  }).filter(inCircle);
}

function proximity(spacing) {
  const buckets=new Map();
  const key=(x,y)=>x+','+y;
  function coordinates(p) { return [Math.floor(p.x/spacing),Math.floor(p.y/spacing)]; }
  return {
    add(p) {
      const id=key(...coordinates(p));
      if(!buckets.has(id)) buckets.set(id,[]);
      buckets.get(id).push(p);
    },
    *neighbors(p) {
      const [x,y]=coordinates(p);
      for(let j=-1;j<=1;j++) for(let i=-1;i<=1;i++) yield* buckets.get(key(x+i,y+j)) ?? [];
    },
    blocks(p) {
      const limit=spacing**2;
      for(const q of this.neighbors(p)) if(distance2(p,q)<limit) return true;
      return false;
    }
  };
}
function distance2(p,q) { return (p.x-q.x)**2+(p.y-q.y)**2; }

function filter(prepared) {
  const {candidates,nearest,order}=prepared;
  const isolation=setting('isolation')**2, separation=setting('separation');
  const spacing=separation ? proximity(separation) : null;
  const retained=[];
  for(const i of order) {
    const p=candidates[i];
    if((isolation&&nearest[i]>isolation)||spacing?.blocks(p)) continue;
    spacing?.add(p);
    retained.push(p);
  }
  const remove=Math.round(retained.length*setting('removal')/100);
  const removed=new Set();
  if(remove) {
    const r=setting('radius'), grid=proximity(r), random=rng(0x2e89b1), priorities=[];
    retained.forEach(p=>grid.add(p));
    for(const p of retained) {
      let count=0;
      for(const q of grid.neighbors(p)) if(p!==q&&distance2(p,q)<=r*r) count++;
      const draw=random();
      if(count) priorities.push([p,-Math.log1p(-draw)/count**2]);
    }
    priorities.sort((a,b)=>a[1]-b[1]);
    for(const [p] of priorities.slice(0,remove)) removed.add(p);
  }
  const angle=tau*rng(0x3dd95b)(), c=Math.cos(angle), s=Math.sin(angle);
  const caustic=retained.filter(p=>!removed.has(p)).map(p=>({x:c*p.x-s*p.y,y:s*p.x+c*p.y}));
  return caustic;
}

function background(points) {
  const separation=setting('separation');
  const grid=separation ? proximity(separation) : null;
  if(grid) points.forEach(p=>grid.add(p));
  const expected=setting('background')*Math.PI*radius**2/1e6;
  const random=rng(0x74c2ed), accepted=[];
  for(let t=-Math.log1p(-random());t<expected;t-=Math.log1p(-random())) {
    const distance=radius*Math.sqrt(random()), angle=tau*random();
    const p={x:distance*Math.cos(angle),y:distance*Math.sin(angle)};
    if(grid?.blocks(p)) continue;
    grid?.add(p);
    accepted.push(p);
  }
  return accepted;
}

function prepare() {
  deflected=[];
  const stretch = rng(0x5c20f3)();
  const ratio = wavelengthRatio(setting('max-ratio'), stretch);
  const field=waves(ratio);
  const candidates=makeCandidates(field);
  const grid=proximity(10);
  candidates.forEach(p=>grid.add(p));
  const nearest=candidates.map(p=>{
    let best=Infinity;
    for(const q of grid.neighbors(p)) if(q!==p) best=Math.min(best,distance2(p,q));
    return best;
  });
  const order=candidates.map((_,i)=>i);
  const random=rng(0x57c3a1);
  for(let i=order.length-1;i>0;i--) {
    const j=Math.floor(random()*(i+1));
    [order[i],order[j]]=[order[j],order[i]];
  }
  return {stretch,ratio,field,candidates,deflected,nearest,order};
}

 return {prepare,filter,background};
}
/** Return all boulders as {x,y} in meters relative to the circle center.
 * Uses fixed default parameters and a local seed; never consumes Math.random(). */
function get_boulder_positions(radius,mapNumber){
 const generator=createGenerator(radius,mapNumber);
 const prepared=generator.prepare(),caustic=generator.filter(prepared);
 return [...caustic,...generator.background(caustic)];
}
globalThis.BoulderGenerator=Object.freeze({get_boulder_positions,createGenerator,defaults});
globalThis.get_boulder_positions=get_boulder_positions;
})();
