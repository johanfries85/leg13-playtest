import webfft from "webfft";
"use strict";

// Shared grid primitives. Coordinates are normalized by LEG_M = 1000 m.
const LEG_M=1000;
// Race-time cost uses requested distance plus a fixed map-reading cost per click.
const RUNNING_TIME_SECONDS_PER_METER=0.285;
const MAP_READING_TIME_SECONDS_PER_CLICK=3;
const MAX_RACE_TIME_SECONDS=60*60;
let raceTimeLimitSeconds=MAX_RACE_TIME_SECONDS;

export function normalize(a){let z=0;for(let i=0;i<a.length;i++)z+=a[i];if(!(z>0))return false;let iz=1/z;for(let i=0;i<a.length;i++)a[i]*=iz;return true}

export function makeGrid(N,hM,cx,cy,data=null){
 return {N,hM,h:hM/LEG_M,widthM:N*hM,width:N*hM/LEG_M,c:{x:cx,y:cy},a:data||new Float32Array(N*N)};
}
export function xAt(g,ix){return g.c.x-g.width/2+(ix+.5)*g.h}
export function yAt(g,iy){return g.c.y-g.width/2+(iy+.5)*g.h}
export function idx(g,ix,iy){return iy*g.N+ix}
export function deltaGrid(N, hM, x, y) {
    const g = makeGrid(N, hM, x, y);
    const ix = Math.floor(N / 2);
    const iy = Math.floor(N / 2);

    g.a[idx(g, ix, iy)] = 1;

    return g;
}
export function moments(g) {
    let mx = 0, my = 0;

    for (let iy = 0; iy < g.N; iy++) {
        const y = yAt(g, iy);

        for (let ix = 0; ix < g.N; ix++) {
            const q = g.a[idx(g, ix, iy)];
            mx += xAt(g, ix) * q;
            my += y * q;
        }
    }

    let vx = 0, vy = 0;

    for (let iy = 0; iy < g.N; iy++) {
        const y = yAt(g, iy);

        for (let ix = 0; ix < g.N; ix++) {
            const q = g.a[idx(g, ix, iy)];
            vx += (xAt(g, ix) - mx) ** 2 * q;
            vy += (y - my) ** 2 * q;
        }
    }

    return {
        x: mx,
        y: my,
        sx: Math.sqrt(Math.max(0, vx)),
        sy: Math.sqrt(Math.max(0, vy))
    };
}
export function supportEdges(g, P) {
 const N=g.N,row=new Float64Array(N),col=new Float64Array(N);
 for(let iy=0;iy<N;iy++)for(let ix=0;ix<N;ix++){const q=g.a[idx(g,ix,iy)];row[iy]+=q;col[ix]+=q}
 let cum=0,left=0;for(left=0;left<N;left++){if(cum+col[left]>=P)break;cum+=col[left]}
 cum=0;let right=N-1;for(right=N-1;right>=0;right--){if(cum+col[right]>=P)break;cum+=col[right]}
 cum=0;let bottom=0;for(bottom=0;bottom<N;bottom++){if(cum+row[bottom]>=P)break;cum+=row[bottom]}
 cum=0;let top=N-1;for(top=N-1;top>=0;top--){if(cum+row[top]>=P)break;cum+=row[top]}
 return {xmin:xAt(g,left)-g.h/2,xmax:xAt(g,right)+g.h/2,ymin:yAt(g,bottom)-g.h/2,ymax:yAt(g,top)+g.h/2};
}
export function sample(g, x, y) {
    const fx = (x - (g.c.x - g.width / 2)) / g.h - 0.5;
    const fy = (y - (g.c.y - g.width / 2)) / g.h - 0.5;

    const x0 = Math.floor(fx);
    const y0 = Math.floor(fy);

    if (x0 < 0 || y0 < 0 || x0 >= g.N - 1 || y0 >= g.N - 1) {
        return 0;
    }

    const tx = fx - x0;
    const ty = fy - y0;
    const k = y0 * g.N + x0;

    const a = g.a;
    const N = g.N;

    return (1 - ty) * (
        (1 - tx) * a[k] + tx * a[k + 1]
    ) + ty * (
        (1 - tx) * a[k + N] + tx * a[k + N + 1]
    );
}

// Conservative cell-overlap regridding.
// Treat each source cell as having uniform density within its area and transfer
// probability by geometric overlap. This preserves probability mass under
// h changes instead of sampling only at the new cell centres.
export function resample(src,N,hM,cx,cy){
 const out=makeGrid(N,hM,cx,cy);
 const sh=src.h, oh=out.h;
 const sx0=src.c.x-src.width/2, sy0=src.c.y-src.width/2;
 const ox0=out.c.x-out.width/2, oy0=out.c.y-out.width/2;

 for(let sy=0;sy<src.N;sy++){
   const y0=sy0+sy*sh, y1=y0+sh;
   let j0=Math.floor((y0-oy0)/oh), j1=Math.floor((y1-oy0-1e-14)/oh);
   if(j1<0||j0>=N) continue;
   j0=Math.max(0,j0); j1=Math.min(N-1,j1);

   for(let sx=0;sx<src.N;sx++){
     const mass=src.a[sy*src.N+sx];
     if(!(mass>0)) continue;
     const x0=sx0+sx*sh, x1=x0+sh;
     let i0=Math.floor((x0-ox0)/oh), i1=Math.floor((x1-ox0-1e-14)/oh);
     if(i1<0||i0>=N) continue;
     i0=Math.max(0,i0); i1=Math.min(N-1,i1);

     for(let j=j0;j<=j1;j++){
       const yy0=oy0+j*oh, yy1=yy0+oh;
       const fy=Math.max(0,Math.min(y1,yy1)-Math.max(y0,yy0))/sh;
       if(!(fy>0)) continue;
       for(let i=i0;i<=i1;i++){
         const xx0=ox0+i*oh, xx1=xx0+oh;
         const fx=Math.max(0,Math.min(x1,xx1)-Math.max(x0,xx0))/sh;
         if(fx>0) out.a[j*N+i]+=mass*fx*fy;
       }
     }
   }
 }
 let z=0; for(let i=0;i<out.a.length;i++) z+=out.a[i];
 if(!(z>0)) throw new Error("regridding window does not overlap probability support");
 // Tiny tail loss at a deliberately finite window is renormalized.
 const iz=1/z; for(let i=0;i<out.a.length;i++) out.a[i]*=iz;
 return out;
}




const H_SET_M=[.25,.5,1,2,4,8,16];
const LOCALIZATION_REFINEMENT_POWER=1;
const PPI_PADDING_SIGMAS=6;
const TAIL_P=1e-6;
export const Z_TAIL=4.753424308822899; // Phi^-1(1-1e-6)
const MOVEMENT_ALONG_TRACK_SD_RATIO=.10;
const MOVEMENT_CROSS_TRACK_SD_RATIO=.05;
const SIGMA_PPI_M=3;
const CANDIDATE_CLOUD_PEAK_ALPHA=105;
const CONTROL_RADIUS_M=15,SIGMA_FIND_M=3;
const DETECT_PEAK=.9,DETECT_SIGMA_M=25;
const DETECT_CUTOFF_M=75;
const DISPLAY_N_X=160,DISPLAY_N_Y=340;
const MAX_LOCALIZATION_N=2048;
const MAX_PPI_N=2048;
const MEMORY_LIMIT_MB=450;
// Course references are resolved only when the game is initialized.
// The localization lab imports this module without loading map_course.js.

let fixedN=256,bouldersOn=true,actual={x:0,y:0},finished=false;
let disqualified=false;
let actualTrajectory=[];
let moveCount=0,raceTimeSec=0,lastDetections=[],hoverS={x:0,y:0},movementHovered=false;
const KEYBOARD_DISTANCE_STEP_M=10;
const KEYBOARD_TURN_STEP_RAD=Math.PI/180,KEYBOARD_REPEAT_TURN_STEP_RAD=5*Math.PI/180;
let movementAngle=Math.PI/2,keyboardPreview=false;
let mapDrag=null,mapPreviewMode=null;
let state=null,candidate=null,numericalMessage="";
let fftLocalization=null;
let lastTimes={candidate:0,localization:0,ppi:0,control:0,commit:0};
const fftCache=new Map();

function activeBoulders(){return MapCourse.activeBoulders(bouldersOn)}
function normal(){let u=0,v=0;while(!u)u=Math.random();while(!v)v=Math.random();return Math.sqrt(-2*Math.log(u))*Math.cos(2*Math.PI*v)}
function gaussian2(dx,dy,s){if(s<=0)return 0;return Math.exp(-.5*(dx*dx+dy*dy)/(s*s))/(2*Math.PI*s*s)}
function nextPow2(x){let n=1;while(n<x)n<<=1;return n}
function movementCov(sx,sy){
 const d=Math.hypot(sx,sy);if(d<1e-15)return {d:0,Sxx:0,Syy:0,Sxy:0,sxSD:0,sySD:0};
 const ex=sx/d,ey=sy/d,tx=-ey,ty=ex;
 const sr=MOVEMENT_ALONG_TRACK_SD_RATIO*d,st=MOVEMENT_CROSS_TRACK_SD_RATIO*d;
 const Sxx=sr*sr*ex*ex+st*st*tx*tx,Syy=sr*sr*ey*ey+st*st*ty*ty,Sxy=sr*sr*ex*ey+st*st*tx*ty;
 return {d,Sxx,Syy,Sxy,sxSD:Math.sqrt(Sxx),sySD:Math.sqrt(Syy),sr,st,ex,ey,tx,ty};
}

function chooseCandidateH(sx, sy) {
  if (!document.getElementById("autoH").checked) {
    return Number(document.getElementById("manualH").value);
  }

  const e = supportEdges(state, TAIL_P);
  const m = movementCov(sx, sy);

  const xmin = e.xmin - Z_TAIL * m.sxSD;
  const xmax = e.xmax + Z_TAIL * m.sxSD;
  const ymin = e.ymin - Z_TAIL * m.sySD;
  const ymax = e.ymax + Z_TAIL * m.sySD;

  const needM = LEG_M * Math.max(xmax - xmin, ymax - ymin);

  return selectGridSpacing(needM, fixedN);
}

export function selectGridSpacing(needM, fixedN) {
  for (const hM of H_SET_M) {
    if (fixedN * hM >= needM * 1.03) {
      return hM;
    }
  }
  return H_SET_M[H_SET_M.length - 1];
}

function getFFT(N,hM){
 const key=N+"|"+hM;
 if(fftCache.has(key))return fftCache.get(key);
 const eng=new webfft(N),k=new Float64Array(N),k2=new Float64Array(N);
 const width=(N*hM)/LEG_M,f=2*Math.PI/width;
 for(let i=0;i<N;i++){const v=f*(i<=N/2?i:i-N);k[i]=v;k2[i]=v*v}
 const o={eng,k,k2};fftCache.set(key,o);return o;
}
function realRows(g){
 const rows=new Array(g.N);
 for(let y=0;y<g.N;y++){const r=new Float32Array(2*g.N),off=y*g.N;for(let x=0;x<g.N;x++)r[2*x]=g.a[off+x];rows[y]=r}
 return rows;
}
function inverseRows(rows,eng,N){
 for(let y=0;y<N;y++)for(let x=0;x<N;x++)rows[y][2*x+1]*=-1;
 const R=eng.fft2d(rows),out=new Float32Array(N*N),inv=1/(N*N);
 for(let y=0;y<N;y++)for(let x=0;x<N;x++)out[y*N+x]=Math.max(0,R[y][2*x]*inv);
 normalize(out);return out;
}
function convolveMovement(src,sx,sy,targetH){
 const t0=performance.now(),base=resample(src,fixedN,targetH,src.c.x,src.c.y);
 const F0=getFFT(fixedN,targetH),F=F0.eng.fft2d(realRows(base)),m=movementCov(sx,sy);
 for(let y=0;y<fixedN;y++){const ky=F0.k[y],yy=m.Syy*F0.k2[y],cross=2*m.Sxy*ky,row=F[y];
   for(let x=0;x<fixedN;x++){const amp=Math.exp(-.5*(m.Sxx*F0.k2[x]+cross*F0.k[x]+yy)),q=2*x;row[q]*=amp;row[q+1]*=amp}}
 const a=inverseRows(F,F0.eng,fixedN);
 const out=makeGrid(fixedN,targetH,src.c.x+sx,src.c.y+sy,a);
 lastTimes.candidate=performance.now()-t0;return out;
}
function gaussianSmooth(g,sigmaM){
 const t0=performance.now(),F0=getFFT(g.N,g.hM),F=F0.eng.fft2d(realRows(g)),s=sigmaM/LEG_M;
 for(let y=0;y<g.N;y++){const gy=Math.exp(-.5*s*s*F0.k2[y]),row=F[y];for(let x=0;x<g.N;x++){const a=gy*Math.exp(-.5*s*s*F0.k2[x]),q=2*x;row[q]*=a;row[q+1]*=a}}
 const a=inverseRows(F,F0.eng,g.N),out=makeGrid(g.N,g.hM,g.c.x,g.c.y,a);
 lastTimes.ppi=performance.now()-t0;return out;
}
// Keep localization and PPI grids separate. FFT smoothing is periodic, so
// zero-padding prevents appreciable wrap-around into the posterior support.
function ppiGrid(raw,sigmaM){
 const pad=Math.ceil(PPI_PADDING_SIGMAS*sigmaM/raw.hM);
 const N=nextPow2(raw.N+2*pad);
 if(N>MAX_PPI_N)throw new Error(`PPI grid ${N} exceeds ${MAX_PPI_N} cells per axis`);
 const out=makeGrid(N,raw.hM,raw.c.x,raw.c.y);
 const offset=(N-raw.N)/2;
 for(let iy=0;iy<raw.N;iy++){
   out.a.set(raw.a.subarray(iy*raw.N,(iy+1)*raw.N),(iy+offset)*N+offset);
 }
 return out;
}
function smoothLocalizationPosterior(raw,sigmaM){
 return gaussianSmooth(ppiGrid(raw,sigmaM),sigmaM);
}
function displaySmoothed(g){return gaussianSmooth(g,SIGMA_PPI_M)}

function sampleDetections(position, bs, random=Math.random) {
  const D = [];

  for (const b of bs) {
    const rx = b.x - position.x;
    const ry = b.y - position.y;
    const dM = LEG_M * Math.hypot(rx, ry);

    const u = random(); // Always consume one random number per boulder

    if (dM > DETECT_CUTOFF_M) continue;

    const pd = DETECT_PEAK *
      Math.exp(-0.5 * (dM / DETECT_SIGMA_M) ** 2);

    if (u < pd) {
      D.push({
        rx,
        ry,
        d: Math.hypot(rx, ry)
      });
    }
  }

  return D;
}


/**
 * Reference implementation of Localization B.
 *
 * Directly evaluates the localization likelihood at each candidate
 * runner position. Retained as an independent reference for verifying
 * optimized implementations of the same mathematical procedure.
 *
 * Not intended for gameplay due to its computational cost.
 */
function localizeBReference(g,D,bs){
 const t0=performance.now();if(!D.length){lastTimes.localization=performance.now()-t0;return g}
 const logL = new Float32Array(g.a.length);
 let max = -Infinity;
 for(let iy=0;iy<g.N;iy++){const y=yAt(g,iy);for(let ix=0;ix<g.N;ix++){const x=xAt(g,ix);let ll=0;
   for(const d of D){const sig=Math.max(d.d/5,g.h/2),cut=5*sig;let Li=0;
     for(const b of bs){const ux=b.x-x,uy=b.y-y;if(Math.abs(ux-d.rx)>cut||Math.abs(uy-d.ry)>cut)continue;Li+=gaussian2(ux-d.rx,uy-d.ry,sig)}
     if(!(Li>0)){ll=-Infinity;break}ll+=Math.log(Li)}
   const k=idx(g,ix,iy);logL[k]=ll;if(ll>max)max=ll}}
 const out=makeGrid(g.N,g.hM,g.c.x,g.c.y);
 if(!Number.isFinite(max)){out.a.set(g.a);return out}
 for(let k=0;k<out.a.length;k++)out.a[k]=g.a[k]*Math.exp(logL[k]-max);
 if(!normalize(out.a))out.a.set(g.a);
 lastTimes.localization=performance.now()-t0;return out;
}

/**
 * Optimized implementation of Localization B.
 *
 * Accumulates Gaussian likelihood contributions on a grid using
 * separable Gaussian factors for efficient evaluation.
 *
 * Mathematically equivalent to localizeBReference(), within
 * floating-point precision. Used for gameplay.
 */
export function localizeB(g, D, bs) {
  const t0 = performance.now();

  if (!D.length) {
    lastTimes.localization = performance.now() - t0;
    return g;
  }

  const N = g.N;
  const nCells = g.a.length;

  // Accumulate log likelihoods across detections.
  const logL = new Float64Array(nCells);
  const field = new Float64Array(nCells);

  // Reusable arrays for separable Gaussian factors.
  const gaussianX = new Float64Array(N);
  const gaussianY = new Float64Array(N);

  // Coordinates of the first grid point.
  const x0 = xAt(g, 0);
  const y0 = yAt(g, 0);
  const h = g.h;

  for (const d of D) {
    field.fill(0);

    const sig = Math.max(d.d / 5, h / 2);
    const cut = 5 * sig;
    const invTwoSig2 = -0.5 / (sig * sig);
    const gaussianNorm = 1 / (2 * Math.PI * sig * sig);

    for (const b of bs) {
      // Gaussian center in runner-position coordinates.
      const cx = b.x - d.rx;
      const cy = b.y - d.ry;

      // Candidate grid cells inside the cutoff square.
      // Slightly expand the index bounds; the exact cutoff
      // is checked below.
      const ixMin = Math.max(
        0, Math.ceil((cx - cut - x0) / h) - 1
      );
      const ixMax = Math.min(
        N - 1, Math.floor((cx + cut - x0) / h) + 1
      );
      const iyMin = Math.max(
        0, Math.ceil((cy - cut - y0) / h) - 1
      );
      const iyMax = Math.min(
        N - 1, Math.floor((cy + cut - y0) / h) + 1
      );

      // Precompute horizontal Gaussian factors for this boulder.
      for (let ix = ixMin; ix <= ixMax; ix++) {
        const x = xAt(g, ix);
        const dx = x - cx;

        gaussianX[ix] = Math.exp(dx * dx * invTwoSig2);
      }

      // Precompute vertical Gaussian factors for this boulder.
      for (let iy = iyMin; iy <= iyMax; iy++) {
        const y = yAt(g, iy);
        const dy = y - cy;

        gaussianY[iy] = Math.exp(dy * dy * invTwoSig2);
      }

      for (let iy = iyMin; iy <= iyMax; iy++) {
        const y = yAt(g, iy);

        for (let ix = ixMin; ix <= ixMax; ix++) {
          const x = xAt(g, ix);

          const ux = b.x - x;
          const uy = b.y - y;

          // Identical cutoff to direct Localization B.
          if (
            Math.abs(ux - d.rx) > cut ||
            Math.abs(uy - d.ry) > cut
          ) {
            continue;
          }

          field[idx(g, ix, iy)] +=
            gaussianX[ix] * gaussianY[iy] * gaussianNorm;
        }
      }
    }

    // Combine this detection with previous detections.
    for (let k = 0; k < nCells; k++) {
      if (field[k] > 0 && Number.isFinite(logL[k])) {
        logL[k] += Math.log(field[k]);
      } else {
        logL[k] = -Infinity;
      }
    }
  }

// Match the Float32 storage used by direct Localization B.
const roundedLogL = new Float32Array(logL);

let max = -Infinity;

for (let k = 0; k < nCells; k++) {
  if (roundedLogL[k] > max) max = roundedLogL[k];
}

  const out = makeGrid(g.N, g.hM, g.c.x, g.c.y);

  if (!Number.isFinite(max)) {
    out.a.set(g.a);
    lastTimes.localization = performance.now() - t0;
    return out;
  }

  for (let k = 0; k < nCells; k++) {
   out.a[k] = g.a[k] * Math.exp(roundedLogL[k] - max);
  }

  if (!normalize(out.a)) {
    out.a.set(g.a);
  }

  lastTimes.localization = performance.now() - t0;
  return out;
}

function localizationGrid(
  prior,
  refinementPower = LOCALIZATION_REFINEMENT_POWER
) {
  if (!Number.isInteger(refinementPower) || refinementPower < 0) {
    throw new Error(
      "Localization refinement power must be a nonnegative integer"
    );
  }

  const factor = 2 ** refinementPower;
  const N = prior.N * factor;

  if (!Number.isSafeInteger(N) || N > MAX_LOCALIZATION_N) {
    throw new Error(
      `Localization grid ${N} exceeds ${MAX_LOCALIZATION_N} cells per axis`
    );
  }
 const out=makeGrid(N,prior.hM/factor,prior.c.x,prior.c.y);
 // Exact conservative subdivision: each source-cell probability mass is
 // distributed equally among its factor × factor child cells.
 const fraction=1/(factor*factor);
 for(let iy=0;iy<prior.N;iy++)for(let ix=0;ix<prior.N;ix++){
   const mass=prior.a[iy*prior.N+ix]*fraction;
   for(let dy=0;dy<factor;dy++)for(let dx=0;dx<factor;dx++){
     out.a[(iy*factor+dy)*N+ix*factor+dx]=mass;
   }
 }
 return out;
}

// q(d) = probability that a 2-D N(0,sigma_find^2 I) perturbation of a point at
// distance d lies inside the R_find disk. Deterministic polar quadrature table.
let qTable=null;const QSTEP=.1,QMAX=35;
function buildQ(){
 const n=Math.ceil(QMAX/QSTEP)+1,A=96,R=100,s=SIGMA_FIND_M,rmax=6*s,dr=rmax/R;qTable=new Float32Array(n);
 for(let ii=0;ii<n;ii++){const d=ii*QSTEP;let sum=0;
  for(let ir=0;ir<R;ir++){const r=(ir+.5)*dr,wr=(r/(s*s))*Math.exp(-.5*(r/s)*(r/s))*dr;let hit=0;
   for(let ia=0;ia<A;ia++){const th=2*Math.PI*(ia+.5)/A,rr=Math.sqrt(d*d+r*r+2*d*r*Math.cos(th));if(rr<=CONTROL_RADIUS_M)hit++}
   sum+=wr*hit/A}qTable[ii]=Math.max(0,Math.min(1,sum))}
}
function qFind(dM){if(!qTable)buildQ();if(dM>=QMAX)return 0;const f=dM/QSTEP,i=Math.floor(f),t=f-i;return qTable[i]*(1-t)+qTable[i+1]*t}
function findsControl(position, control) {
  const dM = LEG_M * Math.hypot(
    position.x - control.x,
    position.y - control.y
  );
  return Math.random() < qFind(dM);
}
function conditionControl(g, control, found) {
  const t0 = performance.now();
  const out = makeGrid(
    g.N, g.hM, g.c.x, g.c.y,
    new Float32Array(g.a)
  );

  for (let iy = 0; iy < g.N; iy++) {
    const y = yAt(g, iy);

    for (let ix = 0; ix < g.N; ix++) {
      const x = xAt(g, ix);
      const dM = LEG_M * Math.hypot(
        x - control.x,
        y - control.y
      );

      const q = qFind(dM);
      const likelihood = found ? q : 1 - q;

      out.a[idx(g, ix, iy)] *= likelihood;
    }
  }

  normalize(out.a);
  lastTimes.control = performance.now() - t0;
  return out;
}
// Diagnostics only: inspect the posterior immediately before fixed-grid resampling.
let lastGridDiagnostic = "No movement committed yet.";
function diagnoseGridSelection(g, m, e, needM, selectedH) {
  let total=0, outside128=0, outside256=0, outside512=0;
  let positive=0, minPositive=Infinity, maxMass=0;
  for(let iy=0;iy<g.N;iy++){
    const y=yAt(g,iy), dy=Math.abs(y-m.y)*LEG_M;
    for(let ix=0;ix<g.N;ix++){
      const q=g.a[idx(g,ix,iy)];
      if(!(q>0))continue;
      const dx=Math.abs(xAt(g,ix)-m.x)*LEG_M;
      total+=q;positive++;minPositive=Math.min(minPositive,q);maxMass=Math.max(maxMass,q);
      if(dx>128||dy>128)outside128+=q;
      if(dx>256||dy>256)outside256+=q;
      if(dx>512||dy>512)outside512+=q;
    }
  }
  const fmt=v=>Number(v).toExponential(3);
  return `
GRID DIAGNOSTICS (before return to fixed grid):
 input grid: N=${g.N}, h=${g.hM} m, width=${g.widthM} m
 support x: ${(e.xmax-e.xmin)*LEG_M|0} m; y: ${(e.ymax-e.ymin)*LEG_M|0} m (TAIL_P=${TAIL_P})
 support bounds x: ${(e.xmin*LEG_M).toFixed(1)} to ${(e.xmax*LEG_M).toFixed(1)} m
 support bounds y: ${(e.ymin*LEG_M).toFixed(1)} to ${(e.ymax*LEG_M).toFixed(1)} m
 required width (+3%): ${needM.toFixed(1)} m; selected h: ${selectedH} m
 input total mass: ${fmt(total)}; positive cells: ${positive}
 min positive cell mass: ${fmt(minPositive)}; max cell mass: ${fmt(maxMass)}
 mass outside mean +/-128 m square: ${fmt(outside128)}
 mass outside mean +/-256 m square: ${fmt(outside256)}
 mass outside mean +/-512 m square: ${fmt(outside512)}`;
}
// PPI visibility cutoff: alpha=2 in the position-cloud heatmap (peak alpha=145).
const PPI_CUTOFF_ALPHA=2, POSITION_CLOUD_ALPHA=145;
const PPI_RELATIVE_DENSITY_CUTOFF=(PPI_CUTOFF_ALPHA/POSITION_CLOUD_ALPHA)**2;
let lastPpiDiscarded=0, largestPpiDiscards=[];
function truncatePpiDensity(g, cutoffAlpha=PPI_CUTOFF_ALPHA){
  // Uniform-area cells: comparing cell mass to peak mass is equivalent to
  // comparing probability densities. Mutate only the final conditioned grid.
  const a=g.a;
  let peak=0, total=0;
  for(let i=0;i<a.length;i++){
    const q=a[i];
    if(q>peak)peak=q;
    total+=q;
  }
  if(!(peak>0)||!(total>0))throw new Error("PPI truncation: empty posterior");
  const cutoff=peak*(cutoffAlpha/POSITION_CLOUD_ALPHA)**2;
  let removed=0;
  for(let i=0;i<a.length;i++){
    const q=a[i];
    if(q>0&&q<cutoff){removed+=q;a[i]=0;}
  }
  if(!normalize(a))throw new Error("PPI truncation removed all probability");
  const fraction=removed/total;
  return fraction;
}
function recordPpiDiscard(fraction,move){
  lastPpiDiscarded=fraction;
  largestPpiDiscards.push({move,fraction});
  largestPpiDiscards.sort((a,b)=>b.fraction-a.fraction||a.move-b.move);
  if(largestPpiDiscards.length>5)largestPpiDiscards.length=5;
}
function recenterAndReturnToFixed(g){
 const m=moments(g);
 // After PPI/control, choose finest h that contains the Player-facing support.
 const e=supportEdges(g,TAIL_P),needM=LEG_M*Math.max(e.xmax-e.xmin,e.ymax-e.ymin)*1.03;
 let hM=H_SET_M[H_SET_M.length-1];for(const h of H_SET_M)if(fixedN*h>=needM){hM=h;break}
 lastGridDiagnostic=diagnoseGridSelection(g,m,e,needM,hM);
 return resample(g,fixedN,hM,m.x,m.y);
}

function consumeRaceTime(fullTimeSeconds){
  const remaining=Math.max(0,raceTimeLimitSeconds-raceTimeSec);
  const timedOut=fullTimeSeconds>remaining;
  // Spread the full click cost across its realized displacement. This includes
  // map-reading overhead for requests, and running time only for punching.
  const fraction=timedOut?remaining/fullTimeSeconds:1;
  raceTimeSec=timedOut?raceTimeLimitSeconds:raceTimeSec+fullTimeSeconds;
  return {fraction,timedOut};
}

function commit(sx,sy,previewAfterCommit){
  if(finished)return;
  const t0=performance.now();

  const prior=candidate||convolveMovement(
    state,sx,sy,chooseCandidateH(sx,sy)
  );

  const d=Math.hypot(sx,sy);
  const fullTimeSeconds=RUNNING_TIME_SECONDS_PER_METER*d*LEG_M+MAP_READING_TIME_SECONDS_PER_CLICK;
  const timing=consumeRaceTime(fullTimeSeconds);

  if(d>0){
    const m=movementCov(sx,sy);
    const eR=normal()*m.sr;
    const eT=normal()*m.st;
    actual.x+=timing.fraction*(sx+eR*m.ex+eT*m.tx);
    actual.y+=timing.fraction*(sy+eR*m.ey+eT*m.ty);
  }

  actualTrajectory.push({...actual});
  moveCount++;
  disqualified=timing.timedOut;
  if(disqualified){
    // An interrupted request never reaches its observation/localization step.
    // Keep the last posterior rather than applying the uncompleted full prior.
    finished=true;candidate=null;lastTimes.commit=performance.now()-t0;
    draw();return;
  }

  // Observe both controls at R's actual position.
  const foundC12=findsControl(actual,MapCourse.C12);
  const foundC13=findsControl(actual,MapCourse.C13);

  // Boulder observations and localization.
  lastDetections = sampleDetections(actual, activeBoulders());
  let fine=localizationGrid(prior);
  const localize=document.getElementById("fftLocalization")?.checked&&fftLocalization?fftLocalization:localizeB;
  fine=localize(fine,lastDetections,activeBoulders());
  fine=smoothLocalizationPosterior(fine,SIGMA_PPI_M);

  // Condition on both control observations.
  fine=conditionControl(fine,MapCourse.C12,foundC12);
  fine=conditionControl(fine,MapCourse.C13,foundC13);

  // Truncate after all control observations, before support-based resampling.
  const discarded=truncatePpiDensity(fine);
  state=recenterAndReturnToFixed(fine);
  recordPpiDiscard(discarded,moveCount);

  // Once the flag is found, R runs directly to C13 to punch automatically.
  // This is not another requested movement: no map-reading cost, random draw,
  // localization, or move-count increment. Record its distance in the route.
  if(foundC13&&!disqualified){
    const punchingDistanceM=LEG_M*Math.hypot(MapCourse.C13.x-actual.x,MapCourse.C13.y-actual.y);
    const punchingTimeSec=punchingDistanceM*RUNNING_TIME_SECONDS_PER_METER;
    const punchingTiming=consumeRaceTime(punchingTimeSec);
    actual=punchingTiming.timedOut?{
      x:actual.x+punchingTiming.fraction*(MapCourse.C13.x-actual.x),
      y:actual.y+punchingTiming.fraction*(MapCourse.C13.y-actual.y)
    }:{...MapCourse.C13};
    if(punchingDistanceM>0&&punchingTiming.fraction>0)actualTrajectory.push({...actual});
    disqualified=punchingTiming.timedOut;
  }
  // Finish only after posterior processing and the automatic punching move.
  finished=disqualified||foundC13;

  candidate=null;
  lastTimes.commit=performance.now()-t0;
  if(previewAfterCommit==="clear"){
    hoverS={x:0,y:0};
    movementAngle=Math.PI/2;
    keyboardPreview=false;
    mapPreviewMode=null;
  }else if(previewAfterCommit==="refresh"&&!finished)updateMovementPreview(sx,sy);
  draw();
}
function startConfiguration(index){
 try{
  MapCourse.setConfiguration(index);
  reset();return true;
 }catch(error){
  numericalMessage="CONFIGURATION ERROR: "+error.message;
  document.getElementById("configurationError").textContent=error.message;
  draw();return false;
 }
}
function nextConfiguration(){if(MapCourse.configuration.index<Number.MAX_SAFE_INTEGER)return startConfiguration(MapCourse.configuration.index+1)}
let choosingReplay=false;
function reset(){
 choosingReplay=false;
 disqualified=false;
 document.getElementById("configurationForm").hidden=true;
 document.getElementById("configurationError").textContent="";
 document.getElementById("playByNumber").setAttribute("aria-expanded","false");
 movementHovered=false;
 mapDrag=null;
 mapPreviewMode=null;
 movementAngle=Math.PI/2;
 keyboardPreview=false;
 actual={x:MapCourse.C12.x,y:MapCourse.C12.y};finished=false;moveCount=0;raceTimeSec=0;lastDetections=[];hoverS={x:0,y:0};candidate=null;numericalMessage="";lastGridDiagnostic="No movement committed yet.";lastPpiDiscarded=0;largestPpiDiscards=[];
 actualTrajectory=[{...actual}];
 const hM=Number(document.getElementById("manualH").value)||.5;state=deltaGrid(fixedN,hM,MapCourse.C12.x,MapCourse.C12.y);draw();
}

let map,mapVector,move,moveCloud,mc,moveCtx;
function fit(c){const r=c.getBoundingClientRect(),q=devicePixelRatio||1,w=Math.round(r.width*q),h=Math.round(r.height*q);if(c.width!==w||c.height!==h){c.width=w;c.height=h;c.getContext("2d").setTransform(q,0,0,q,0,0)}}
// Optional rendering target/viewport lets the lab use the same heatmap algorithm.
// Existing game calls (without options) retain their original behavior.
function drawHeatmap(g,w,h,r,gc,b,alpha,alreadyPPI=false,options={}){
 const shown=alreadyPPI ? g : displaySmoothed(g);
 let maxD=0;
 for(let i=0;i<shown.a.length;i++)maxD=Math.max(maxD,shown.a[i]);
 maxD/=shown.h*shown.h;
 if(!(maxD>0))return;
 const ctx=options.ctx||mc;
 const bounds=options.bounds||{
   xmin:MapCourse.DISPLAY_X_MIN,xmax:MapCourse.DISPLAY_X_MAX,
   ymin:MapCourse.DISPLAY_Y_MIN,ymax:MapCourse.DISPLAY_Y_MAX
 };
 const nx=options.nx||DISPLAY_N_X,ny=options.ny||DISPLAY_N_Y;
 const off=document.createElement("canvas");off.width=nx;off.height=ny;
 const oc=off.getContext("2d"),im=oc.createImageData(nx,ny);
 for(let py=0;py<ny;py++)for(let px=0;px<nx;px++){
   const x=bounds.xmin+(px+.5)/nx*(bounds.xmax-bounds.xmin);
   const y=bounds.ymax-(py+.5)/ny*(bounds.ymax-bounds.ymin);
   const v=(sample(shown,x,y)/(shown.h*shown.h))/maxD;
   const k=4*(py*nx+px);
   im.data[k]=r;im.data[k+1]=gc;im.data[k+2]=b;
   im.data[k+3]=Math.round(alpha*Math.sqrt(Math.max(0,v)));
 }
 oc.putImageData(im,0,0);ctx.drawImage(off,0,0,w,h);
}
function rebuildMap(){const svgEl=MapCourse.svgEl;MapCourse.render(mapVector,{bouldersOn,beforeBoulders:document.getElementById("debug").checked?(target)=>{const q=svgEl("g");target.appendChild(q);q.appendChild(svgEl("circle",{cx:MapCourse.sxm(actual.x),cy:MapCourse.sym(actual.y),r:3,fill:"#c00"}))}:null})}
function drawMap(){fit(map);const r=map.getBoundingClientRect();mc.clearRect(0,0,r.width,r.height);
 // The initial delta needs display smoothing to remain visible on the heatmap.
 // After a move, state already includes substantive PPI smoothing.
 drawHeatmap(state,r.width,r.height,255,0,0,145,moveCount>0);
 // candidate is an internal unsmoothed prior, so smooth only its displayed copy.
 if(candidate&&!finished)drawHeatmap(candidate,r.width,r.height,255,0,255,CANDIDATE_CLOUD_PEAK_ALPHA,false);
 rebuildMap();
 // Never construct route geometry in the DOM while the leg is active.
 if(finished){
  const points=actualTrajectory.filter((p,i)=>i===0||p.x!==actualTrajectory[i-1].x||p.y!==actualTrajectory[i-1].y);
  if(points.length>1)mapVector.appendChild(MapCourse.svgEl("polyline",{
   points:points.map(p=>`${MapCourse.sxm(p.x)},${MapCourse.sym(p.y)}`).join(" "),
   fill:"none",stroke:"#ff2428","stroke-width":4,"stroke-linecap":"round","stroke-linejoin":"round",
   "pointer-events":"none","data-actual-trajectory":"","aria-label":"Actual runner trajectory"
  }));
 }
 const overlayOn=document.getElementById("worstCaseLeg").checked;
 const slider=document.getElementById("overlayShift"),shift=Number(slider.value);
 slider.disabled=!overlayOn;
 document.getElementById("overlayShiftValue").textContent=`${shift>0?"+":""}${shift} m`;
 if(overlayOn){
  const overlay=MapCourse.drawWorstCaseOverlay(mapVector,shift);
  drawPanelFrames(overlay);
 }
}
function drawPanelFrames(overlay){
 const bounds=map.getBoundingClientRect();
 const form=document.getElementById("configurationForm"),formHidden=form.hidden;
 form.hidden=true; // Compare the compact choice panel, before number entry.
 try {
  for(const [id,color,dash] of [["finishPanel","#ff2428","none"],["replayPanel","#000","12 8"]]){
   const panel=document.getElementById(id),hidden=panel.hidden,visibility=panel.style.visibility;
   let rect,texts=[],buttons=[];
   try {
    panel.style.visibility="hidden";panel.hidden=false;rect=panel.getBoundingClientRect();
    if(id==="finishPanel")texts=Array.from(panel.querySelectorAll(".finishTitle, #finishTime, th, td, button"),el=>{
     const box=el.getBoundingClientRect(),style=getComputedStyle(el);
     if(el.tagName.toLowerCase()==="button")buttons.push({rect:box,id:el.id});
     const align=style.textAlign;
     const x=align==="right"?box.left+box.width-parseFloat(style.paddingRight||0):
      align==="center"?box.left+box.width/2:box.left+parseFloat(style.paddingLeft||0);
     return {text:el.textContent,x,anchor:align==="right"?"end":align==="center"?"middle":"start",
      y:box.top+box.height/2,fontSize:parseFloat(style.fontSize),fontWeight:style.fontWeight};
    });
   }
   finally {panel.hidden=hidden;panel.style.visibility=visibility;}
   overlay.appendChild(MapCourse.svgEl("rect",{
    x:(rect.left-bounds.left)*800/bounds.width,y:(rect.top-bounds.top)*1700/bounds.height,
    width:rect.width*800/bounds.width,height:rect.height*1700/bounds.height,
    rx:5*800/bounds.width,fill:"none",stroke:color,"stroke-width":2,
    "stroke-dasharray":dash,"data-panel-frame":id
   }));
   for(const t of texts)overlay.appendChild(MapCourse.svgEl("text",{
    x:(t.x-bounds.left)*800/bounds.width,y:(t.y-bounds.top)*1700/bounds.height,
    fill:"#ff2428","font-family":"Arial, sans-serif","font-size":t.fontSize*800/bounds.width,
    "font-weight":t.fontWeight,"text-anchor":t.anchor,"dominant-baseline":"central","data-panel-text":"finishPanel"
   },t.text));
   for(const {rect:buttonRect,id:buttonId} of buttons)overlay.appendChild(MapCourse.svgEl("rect",{
    x:(buttonRect.left-bounds.left)*800/bounds.width,y:(buttonRect.top-bounds.top)*1700/bounds.height,
    width:buttonRect.width*800/bounds.width,height:buttonRect.height*1700/bounds.height,
    rx:2*800/bounds.width,fill:"none",stroke:"#ff2428","stroke-width":2,
    "data-panel-button":buttonId
   }));
  }
 } finally {form.hidden=formHidden;}
}
function drawMovementCloud(){
 const c=moveCtx,w=moveCloud.width,h=moveCloud.height;c.clearRect(0,0,w,h);const sx=hoverS.x,sy=hoverS.y,d=Math.hypot(sx,sy);if(!d)return;
 const m=movementCov(sx,sy),cx=w/2,cy=h/2,scale=w/2,img=c.createImageData(w,h),a=img.data;
 // Apply the same display smoothing and square-root opacity as the candidate cloud.
 const displaySigma=SIGMA_PPI_M/LEG_M;
 const alongSigma=Math.hypot(m.sr,displaySigma),crossSigma=Math.hypot(m.st,displaySigma);
 for(let j=0;j<h;j++){
  const dy=(cy-j-.5)/scale-sy;
  for(let i=0;i<w;i++){
   const dx=(i+.5-cx)/scale-sx;
   const u=(dx*m.ex+dy*m.ey)/alongSigma,v=(dx*m.tx+dy*m.ty)/crossSigma;
   const alpha=Math.round(CANDIDATE_CLOUD_PEAK_ALPHA*Math.exp(-.25*(u*u+v*v)));
   if(!alpha)continue;
   const k=4*(j*w+i);a[k]=255;a[k+2]=255;a[k+3]=alpha;
  }
 }
 c.putImageData(img,0,0);
}
export function movementDirectionDegrees(x,y,northDegrees){
 if(x===0&&y===0)return null;
 const degrees=Math.round(Math.atan2(x,y)*180/Math.PI-northDegrees);
 return ((degrees+180)%360+360)%360-180;
}
export function formatLegPace(timeSeconds,lengthM){
 const seconds=Math.round(timeSeconds*1000/lengthM);
 return `${Math.floor(seconds/60)}:${String(seconds%60).padStart(2,"0")}/km`;
}
function formatRaceTime(timeSeconds){
 return `${Math.floor(timeSeconds/60)}:${String(Math.floor(timeSeconds%60)).padStart(2,"0")}`;
}
export function routeLengthMeters(points){
 let distance=0;
 for(let i=1;i<points.length;i++)distance+=Math.hypot(points[i].x-points[i-1].x,points[i].y-points[i-1].y);
 return distance*LEG_M;
}
function drawMove(){
 const svgEl=MapCourse.svgEl;
 drawMovementCloud();move.replaceChildren();const S=2*LEG_M,cx=LEG_M,cy=LEG_M,g=svgEl("g");move.appendChild(g);
 g.appendChild(svgEl("line",{x1:cx,y1:0,x2:cx,y2:S,stroke:"#aaa"}));g.appendChild(svgEl("line",{x1:0,y1:cy,x2:S,y2:cy,stroke:"#aaa"}));
 for(let m=100;m<=LEG_M;m+=100){g.appendChild(svgEl("circle",{cx,cy,r:m,fill:"none",stroke:"#777","stroke-dasharray":"5 5"}))}
 for(let m=200;m<=LEG_M;m+=200){
  const offset=m/Math.SQRT2;
  g.appendChild(svgEl("text",{x:cx+offset,y:cy-offset,fill:"#ff2428","font-family":"Arial, sans-serif","font-size":37,"text-anchor":"middle","dominant-baseline":"central"},`${m} m`));
 }
 // Both SVGs display 100 units per CSS cm, so match the map arrow dimensions.
 const arrowWidth=56,arrowHeight=103;
 const north=svgEl("g",{transform:`rotate(${MapCourse.MAP_ORIENTATION_DEG} ${cx} ${cy})`});
 north.appendChild(MapCourse.createNorthArrow(cx-arrowWidth/2,cy-arrowHeight/2,arrowWidth,arrowHeight));
 g.appendChild(north);
 const px=cx+hoverS.x*LEG_M,py=cy-hoverS.y*LEG_M;g.appendChild(svgEl("line",{x1:cx,y1:cy,x2:px,y2:py,stroke:"#111","stroke-width":2}));g.appendChild(svgEl("circle",{cx:px,cy:py,r:5,fill:"#111"}));move.appendChild(svgEl("circle",{cx,cy,r:LEG_M-.5,fill:"none",stroke:"#777"}))
 const readout=document.getElementById("movementReadout");
 readout.hidden=!candidate||(!movementHovered&&!keyboardPreview&&!mapDrag&&!mapPreviewMode)||finished;
 if(!readout.hidden){
  const direction=movementDirectionDegrees(hoverS.x,hoverS.y,MapCourse.MAP_ORIENTATION_DEG);
  const distance=Math.hypot(hoverS.x,hoverS.y)*1000;
  const lines=[
   `Direction: ${direction===null?"—":direction+"°"}`,
   `Distance: ${Math.round(distance)} m`,
   `Time cost: ${Math.round(RUNNING_TIME_SECONDS_PER_METER*distance+MAP_READING_TIME_SECONDS_PER_CLICK)} s`
  ];
  document.getElementById("movementReadoutText").replaceChildren(
   ...lines.map((line,index)=>svgEl("tspan",{x:0,y:29+index*36},line))
  );
 }
}
function status(){
 const m=moments(state),d=Math.hypot(hoverS.x,hoverS.y),mins=Math.floor(raceTimeSec/60),secs=Math.floor(raceTimeSec%60);
 let cand=candidate?`${candidate.hM} m (${candidate.widthM} m window)`:"—";
 return `${finished?(disqualified?"MAX TIME EXCEEDED — disqualified\n":"C13 FOUND — run finished\n"):""}${numericalMessage?numericalMessage+"\n":""}race time: ${mins}:${String(secs).padStart(2,"0")}
moves: ${moveCount}
positive detections last move: ${lastDetections.length}
fixed N: ${fixedN}
stored position grid: h=${state.hM} m, width=${state.widthM} m
candidate grid: ${cand}
automatic h: ${document.getElementById("autoH").checked?"on":"off"}
σ_PPI: ${SIGMA_PPI_M} m
localization grid: N=${fixedN*2**LOCALIZATION_REFINEMENT_POWER}, k=${LOCALIZATION_REFINEMENT_POWER}, h=prior h/${2**LOCALIZATION_REFINEMENT_POWER}
PPI padding: at least ${PPI_PADDING_SIGMAS}σ per side (FFT grid cap ${MAX_PPI_N})
cloud mean: (${(m.x*1000).toFixed(0)}, ${(m.y*1000).toFixed(0)}) m
cloud SD x/y: ${(m.sx*1000).toFixed(1)} / ${(m.sy*1000).toFixed(1)} m
requested movement: (${(hoverS.x*1000).toFixed(0)}, ${(hoverS.y*1000).toFixed(0)}) m
requested length: ${(d*1000).toFixed(0)} m

last timings:
 candidate prior: ${lastTimes.candidate.toFixed(1)} ms
 Localization B: ${lastTimes.localization.toFixed(1)} ms
 PPI smoothing: ${lastTimes.ppi.toFixed(1)} ms
 smooth control miss: ${lastTimes.control.toFixed(1)} ms
 full commit: ${lastTimes.commit.toFixed(1)} ms

PPI density truncation (after controls, before resampling):
 cutoff alpha: ${PPI_CUTOFF_ALPHA} / 255 (peak ${POSITION_CLOUD_ALPHA} / 255)
 relative density cutoff: ${PPI_RELATIVE_DENSITY_CUTOFF.toFixed(7)}
 discarded last move: ${(lastPpiDiscarded*100).toFixed(5)}%
 five largest discarded masses this game:
${largestPpiDiscards.length?largestPpiDiscards.map(x=>`  move ${x.move}: ${(x.fraction*100).toFixed(5)}%`).join("\n"):"  —"}${lastGridDiagnostic}`;
}
function draw(){
 drawMove();document.getElementById("status").textContent=status();
 document.getElementById("configurationNumber").textContent=`Configuration: ${MapCourse.configuration.index}`;
 document.getElementById("legLength").textContent=`${Math.round(MapCourse.configuration.legLengthM)} m`;
 document.getElementById("legPace").textContent=disqualified?"-":formatLegPace(raceTimeSec,MapCourse.configuration.legLengthM);
 document.getElementById("finishTitle").textContent=disqualified?"Disqualified.":"Leg 13 finished!";
 const routeM=finished?routeLengthMeters(actualTrajectory):0;
 document.getElementById("routeLength").textContent=finished?`${Math.round(routeM)} m`:"—";
 document.getElementById("routePace").textContent=finished&&routeM>0?formatLegPace(raceTimeSec,routeM):"—";
 document.getElementById("next").disabled=MapCourse.configuration.index===Number.MAX_SAFE_INTEGER;
 document.getElementById("finishPanel").hidden=!finished||choosingReplay;
 document.getElementById("replayPanel").hidden=!choosingReplay;
 document.getElementById("playNew").disabled=MapCourse.configuration.index===Number.MAX_SAFE_INTEGER;
 document.getElementById("playNext").disabled=MapCourse.configuration.index===Number.MAX_SAFE_INTEGER;
 const time=formatRaceTime(raceTimeSec);
 document.getElementById("elapsedPanel").hidden=finished;
 document.getElementById("elapsedTime").textContent=`Elapsed time: ${time}`;
 const finishTime=document.getElementById("finishTime");
 finishTime.className=disqualified?"finishReason":"";
 finishTime.textContent=disqualified?`Max time exceeded (limit: ${formatRaceTime(raceTimeLimitSeconds)}).`:`Time: ${time}`;
 drawMap();
}

function updateMovementPreview(x,y){
 hoverS={x,y};
 try{candidate=convolveMovement(state,x,y,chooseCandidateH(x,y));numericalMessage=""}
 catch(err){candidate=null;numericalMessage="HOVER ERROR: "+err.message}
}

function commitPreview(previewAfterCommit=mapPreviewMode||"refresh"){
 try{commit(hoverS.x,hoverS.y,previewAfterCommit)}
 catch(err){numericalMessage="COMMIT ERROR: "+err.message;candidate=null;draw()}
}

function focusMovementControls(){
 const target=document.getElementById("movementAxes").hidden?document.getElementById("movementPanel"):move;
 target.focus();
}

function handleMovementKey(event){
 if(finished||mapDrag||!document.getElementById("keyboardMovement")?.checked||
    event.defaultPrevented||event.altKey||event.ctrlKey||event.metaKey||event.shiftKey||
    event.target?.closest?.('input, textarea, select, button, summary, a, [contenteditable]'))return;
 const {key}=event;
 if(key===" "){
  event.preventDefault();
  if(!event.repeat&&candidate){
   keyboardPreview=true;
   commitPreview();
  }
  return;
 }
 if(!["ArrowUp","ArrowDown","ArrowLeft","ArrowRight"].includes(key))return;
 event.preventDefault();
 let distance=Math.hypot(hoverS.x,hoverS.y);
 const step=KEYBOARD_DISTANCE_STEP_M/LEG_M;
 if(key==="ArrowUp")distance=Math.min(1,distance+step);
 if(key==="ArrowDown")distance=Math.max(0,distance-step);
 const turnStep=event.repeat?KEYBOARD_REPEAT_TURN_STEP_RAD:KEYBOARD_TURN_STEP_RAD;
 if(key==="ArrowLeft")movementAngle+=turnStep;
 if(key==="ArrowRight")movementAngle-=turnStep;
 keyboardPreview=true;
 mapPreviewMode=null;
 updateMovementPreview(distance*Math.cos(movementAngle),distance*Math.sin(movementAngle));
 draw();
}

function cancelMapDrag(){
 if(!mapDrag)return;
 ({hoverS,candidate,movementAngle,keyboardPreview,numericalMessage}=mapDrag.previous);
 mapDrag=null;
 draw();
}

function updateMapDrag(event){
 const {width,height}=map.getBoundingClientRect();
 const x=mapDrag.baseX+(event.clientX-mapDrag.startX)/width*(MapCourse.DISPLAY_X_MAX-MapCourse.DISPLAY_X_MIN);
 const y=mapDrag.baseY+(mapDrag.startY-event.clientY)/height*(MapCourse.DISPLAY_Y_MAX-MapCourse.DISPLAY_Y_MIN);
 const distance=Math.hypot(x,y);
 const limit=Math.max(1,distance);
 if(distance>0)movementAngle=Math.atan2(y,x);
 keyboardPreview=false;
 updateMovementPreview(x/limit,y/limit);
 mapDrag.lastX=event.clientX;
 mapDrag.lastY=event.clientY;
}

export function initializeGame({fftLocalizer=null,maxRaceTimeSeconds=MAX_RACE_TIME_SECONDS}={}){
 if(!Number.isFinite(maxRaceTimeSeconds)||maxRaceTimeSeconds<=0)throw new RangeError('Race time limit must be positive and finite');
 raceTimeLimitSeconds=maxRaceTimeSeconds;
 fftLocalization=fftLocalizer;
 map=document.getElementById("map");
 mapVector=document.getElementById("mapVector");
 move=document.getElementById("move");
 moveCloud=document.getElementById("moveCloud");
 mc=map.getContext("2d");
 moveCtx=moveCloud.getContext("2d");
const movementAxes=document.getElementById("movementAxes");
const movementAxesEnabled=document.getElementById("movementAxesEnabled");
const keyboardMovement=document.getElementById("keyboardMovement");
const dragOnMap=document.getElementById("dragOnMap");
const mapDragSettings=document.getElementById("mapDragSettings");
const mapDragAdjust=document.getElementById("mapDragAdjust");
const mapDragRelease=document.getElementById("mapDragRelease");
const mapDragClick=document.getElementById("mapDragClick");
const mapDragKeepPreview=document.getElementById("mapDragKeepPreview");
function restoreMovementControlDefaults(){
 movementAxesEnabled.checked=movementAxesEnabled.defaultChecked;
 movementAxes.hidden=!movementAxesEnabled.checked;
 if(keyboardMovement)keyboardMovement.checked=keyboardMovement.defaultChecked;
 dragOnMap.checked=dragOnMap.defaultChecked;
 for(const input of mapDragSettings.querySelectorAll("input"))input.checked=input.defaultChecked;
 mapDragSettings.disabled=!dragOnMap.checked;
}
restoreMovementControlDefaults();
movementAxesEnabled.onchange=()=>{
 movementAxes.hidden=!movementAxesEnabled.checked;
 movementHovered=false;
 drawMove();
 if(keyboardMovement?.checked)focusMovementControls();
};
move.addEventListener("mousemove",e=>{
 if(movementAxes.hidden||mapDrag)return;
 movementHovered=true;
 if(finished)return;const r=move.getBoundingClientRect();let x=2*(e.clientX-r.left)/r.width-1,y=1-2*(e.clientY-r.top)/r.height,d=Math.hypot(x,y);if(d>1){x/=d;y/=d}
 keyboardPreview=false;
 mapPreviewMode=null;
 if(d>0)movementAngle=Math.atan2(y,x);
 updateMovementPreview(x,y);
 draw();
});
move.addEventListener("mouseleave",()=>{movementHovered=false;drawMove()});
move.addEventListener("click",()=>{if(!movementAxes.hidden&&!mapDrag)commitPreview()});
dragOnMap.onchange=()=>{
 mapDragSettings.disabled=!dragOnMap.checked;
 if(!dragOnMap.checked){
  cancelMapDrag();
  mapPreviewMode=null;
  drawMove();
 }
 if(keyboardMovement?.checked)focusMovementControls();
};
mapDragSettings.onchange=()=>{
 cancelMapDrag();
 if(mapPreviewMode)mapPreviewMode=mapDragKeepPreview.checked?"refresh":"clear";
};
map.addEventListener("pointerdown",event=>{
 if(!dragOnMap.checked||finished||mapDrag||event.pointerType!=="mouse"||event.button!==0)return;
 event.preventDefault();
 map.setPointerCapture(event.pointerId);
 const adjust=mapDragAdjust.checked&&candidate;
 mapDrag={
  pointerId:event.pointerId,startX:event.clientX,startY:event.clientY,
  lastX:event.clientX,lastY:event.clientY,hasDragged:false,
  baseX:adjust?hoverS.x:0,
  baseY:adjust?hoverS.y:0,
  previous:{hoverS,candidate,movementAngle,keyboardPreview,numericalMessage}
 };
});
map.addEventListener("pointermove",event=>{
 if(mapDrag?.pointerId!==event.pointerId)return;
 if(event.clientX===mapDrag.lastX&&event.clientY===mapDrag.lastY)return;
 mapDrag.hasDragged=true;
 updateMapDrag(event);
 draw();
});
map.addEventListener("pointerup",event=>{
 if(mapDrag?.pointerId!==event.pointerId||event.button!==0)return;
 const {hasDragged,startX,startY}=mapDrag;
 const release=mapDragRelease.checked;
 const previewAfterCommit=mapDragKeepPreview.checked?"refresh":"clear";
 const movedFromStart=event.clientX!==startX||event.clientY!==startY;
 if(!hasDragged&&!movedFromStart){
  mapDrag=null;
  if(mapDragClick.checked&&candidate)commitPreview(previewAfterCommit);
  return;
 }
 if(release&&!movedFromStart){cancelMapDrag();return}
 if(event.clientX!==mapDrag.lastX||event.clientY!==mapDrag.lastY)updateMapDrag(event);
 mapDrag=null;
 mapPreviewMode=candidate?previewAfterCommit:null;
 if(release&&candidate)commitPreview(previewAfterCommit);
 else draw();
});
for(const eventName of ["pointercancel","lostpointercapture"]){
 map.addEventListener(eventName,event=>{if(mapDrag?.pointerId===event.pointerId)cancelMapDrag()});
}
if(keyboardMovement)keyboardMovement.onchange=()=>{
 if(keyboardMovement.checked)focusMovementControls();
 else{keyboardPreview=false;drawMove()}
};
window.addEventListener("keydown",handleMovementKey);
document.getElementById("reset").onclick=()=>{
 restoreMovementControlDefaults();
 document.getElementById("instructions").open=true;
 document.getElementById("movementOptions").open=false;
 if(startConfiguration(1)&&keyboardMovement?.checked)focusMovementControls();
};
document.getElementById("next").onclick=nextConfiguration;
document.getElementById("playAgain").onclick=()=>{
 if(!finished)return;
 choosingReplay=true;draw();document.getElementById("playSame").focus();
};
document.getElementById("playSame").onclick=()=>{reset();focusMovementControls()};
document.getElementById("replayBack").onclick=()=>{
 choosingReplay=false;
 document.getElementById("configurationForm").hidden=true;
 document.getElementById("configurationError").textContent="";
 document.getElementById("playByNumber").setAttribute("aria-expanded","false");
 draw();document.getElementById("playAgain").focus();
};
document.getElementById("playNew").onclick=()=>{nextConfiguration();focusMovementControls()};
document.getElementById("playNext").onclick=()=>{
 if(!finished)return;
 nextConfiguration();focusMovementControls();
};
document.getElementById("playByNumber").onclick=()=>{
 document.getElementById("configurationForm").hidden=false;
 document.getElementById("playByNumber").setAttribute("aria-expanded","true");
 const input=document.getElementById("configurationInput");
 input.value=String(MapCourse.configuration.index);input.focus();input.select();
};
document.getElementById("configurationForm").onsubmit=e=>{
 e.preventDefault();
 const input=document.getElementById("configurationInput"),value=input.value.trim(),index=Number(value);
 if(!/^\d+$/.test(value)||!Number.isSafeInteger(index)||index<1){
  document.getElementById("configurationError").textContent=`Enter a whole number from 1 to ${Number.MAX_SAFE_INTEGER}.`;
  input.focus();return;
 }
 if(startConfiguration(index))focusMovementControls();
};
document.getElementById("gridN").onchange=e=>{fixedN=Number(e.target.value);reset()};
document.getElementById("manualH").onchange=()=>{if(!document.getElementById("autoH").checked)reset()};
document.getElementById("autoH").onchange=()=>{candidate=null;draw()};
document.getElementById("boulders").onchange=e=>{bouldersOn=e.target.checked;reset()};
document.getElementById("debug").onchange=draw;
document.getElementById("worstCaseLeg").onchange=draw;
document.getElementById("overlayShift").oninput=draw;
window.addEventListener("resize",draw);
window.addEventListener("pageshow",draw);

fixedN=256;buildQ();startConfiguration(1);

}

export { movementCov, convolveMovement, gaussianSmooth, ppiGrid, smoothLocalizationPosterior, sampleDetections, localizeBReference, localizationGrid, qFind, findsControl, conditionControl, truncatePpiDensity, recenterAndReturnToFixed, drawHeatmap };
