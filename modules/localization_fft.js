// Experimental numerical approximation to Localization B. Not used by gameplay.
import webfft from "webfft";
import { makeGrid, normalize, xAt, yAt } from "./pdf_engine.js";
import { BufferedFFT2D, complexRows } from "./buffered_fft.js";

const engines = new Map();
const bufferedEngines=new Map();
function bufferedFor(N){
  if(!bufferedEngines.has(N))bufferedEngines.set(N,{
    fft:new BufferedFFT2D(N),work:complexRows(N),spectrum:complexRows(N),inverse:complexRows(N),
    raster:new Float32Array(N*N),kernelInput:new Float32Array(2*N),
    kernels:[new Float32Array(2*N),new Float32Array(2*N)]
  });
  return bufferedEngines.get(N);
}
const MAX_FFT_N = 2048;
function engineFor(N) {
  if (!engines.has(N)) engines.set(N, new webfft(N));
  return engines.get(N);
}
function nextPow2(n) { let N=1; while(N<n) N*=2; return N; }

// Samples lie on the world lattice (i*h,j*h). Each boulder has unit total weight.
export function rasterizeBoulders(bs, N, h, ix0, iy0, storage=null) {
  const a=storage||new Float32Array(N*N);if(storage)a.fill(0);
  for(const b of bs) {
    const x=b.x/h-ix0,y=b.y/h-iy0,ix=Math.floor(x),iy=Math.floor(y);
    if(ix < -1 || iy < -1 || ix >= N || iy >= N) continue;
    const tx=x-ix,ty=y-iy;
    for(let dy=0;dy<2;dy++)for(let dx=0;dx<2;dx++) {
      const i=ix+dx,j=iy+dy;
      if(i>=0&&j>=0&&i<N&&j<N)a[j*N+i]+=(dx?tx:1-tx)*(dy?ty:1-ty);
    }
  }
  return a;
}

function directField(g,d,bs,sigma) {
  const out=new Float64Array(g.N*g.N),cut=5*sigma,inv=-.5/(sigma*sigma);
  const gx=new Float64Array(g.N),gy=new Float64Array(g.N);
  const x0=xAt(g,0),y0=yAt(g,0),norm=1/(2*Math.PI*sigma*sigma);
  for(const b of bs) {
    const cx=b.x-d.rx,cy=b.y-d.ry;
    const left=Math.max(0,Math.ceil((cx-cut-x0)/g.h)),right=Math.min(g.N-1,Math.floor((cx+cut-x0)/g.h));
    const bottom=Math.max(0,Math.ceil((cy-cut-y0)/g.h)),top=Math.min(g.N-1,Math.floor((cy+cut-y0)/g.h));
    for(let x=left;x<=right;x++)gx[x]=Math.exp((xAt(g,x)-cx)**2*inv);
    for(let y=bottom;y<=top;y++)gy[y]=Math.exp((yAt(g,y)-cy)**2*inv);
    for(let y=bottom;y<=top;y++)for(let x=left;x<=right;x++)out[y*g.N+x]+=gx[x]*gy[y]*norm;
  }
  return out;
}

// An exact support mask prevents FFT roundoff from inventing likelihood where
// no boulder lies within Localization B's square 5-sigma cutoff.
function supportMask(g,d,bs,cut) {
  const stride=g.N+1,diff=new Int32Array(stride*stride),x0=xAt(g,0),y0=yAt(g,0);
  for(const b of bs) {
    const left=Math.max(0,Math.ceil((b.x-d.rx-cut-x0)/g.h));
    const right=Math.min(g.N-1,Math.floor((b.x-d.rx+cut-x0)/g.h));
    const bottom=Math.max(0,Math.ceil((b.y-d.ry-cut-y0)/g.h));
    const top=Math.min(g.N-1,Math.floor((b.y-d.ry+cut-y0)/g.h));
    if(left>right||bottom>top)continue;
    diff[bottom*stride+left]++;diff[bottom*stride+right+1]--;
    diff[(top+1)*stride+left]--;diff[(top+1)*stride+right+1]++;
  }
  for(let y=0;y<g.N;y++)for(let x=0;x<g.N;x++) {
    const k=y*stride+x;
    diff[k]+=(x?diff[k-1]:0)+(y?diff[k-stride]:0)-(x&&y?diff[k-stride-1]:0);
  }
  return {a:diff,stride};
}

function boulderSpectrum(bs,N,h,ix0,iy0,buffered=null) {
  const a=rasterizeBoulders(bs,N,h,ix0,iy0,buffered?.raster);
  if(buffered){
    for(let y=0;y<N;y++)for(let x=0;x<N;x++){
      buffered.work[y][2*x]=a[y*N+x];buffered.work[y][2*x+1]=0;
    }
    return buffered.fft.transform(buffered.work,buffered.spectrum);
  }
  const rows=Array.from({length:N},(_,y)=>{
    const row=new Float32Array(2*N);
    for(let x=0;x<N;x++)row[2*x]=a[y*N+x];
    return row;
  });
  return engineFor(N).fft2d(rows);
}

function fftField(g,d,bs,sigma,N,pad,shared=null) {
  // Requested output is the localization rectangle translated by the detection.
  // Include a full kernel halo plus deposition/interpolation guard cells.
  const x0=xAt(g,0)+d.rx,y0=yAt(g,0)+d.ry;
  const ix0=shared?shared.ix0:Math.floor(x0/g.h)-pad;
  const iy0=shared?shared.iy0:Math.floor(y0/g.h)-pad;
  const eng=engineFor(N);
  // Filtering mutates F, so preserve the common forward transform across kernels.
  const F=shared?shared.spectrum.map(row=>new Float32Array(row)):boulderSpectrum(bs,N,g.h,ix0,iy0);
  // Separable sampled Gaussian, truncated at 5 sigma, circularly embedded.
  // Only the interior output is read: the halo makes that convolution linear.
  const K=kernelSpectrum(N,g.h,sigma);
  for(let y=0;y<N;y++)for(let x=0;x<N;x++) {
    const kr=K[2*x]*K[2*y]-K[2*x+1]*K[2*y+1];
    const ki=K[2*x]*K[2*y+1]+K[2*x+1]*K[2*y];
    const re=F[y][2*x],im=F[y][2*x+1];
    F[y][2*x]=re*kr-im*ki;
    F[y][2*x+1]=-(re*ki+im*kr); // conjugate for inverse via forward FFT
  }
  return extractField(eng.fft2d(F),g,d,bs,sigma,N,ix0,iy0);
}

function kernelSpectrum(N,h,sigma,buffered=null,slot=0){
  const kernel=buffered?buffered.kernelInput:new Float32Array(2*N),radius=Math.floor(5*sigma/h);
  if(buffered)kernel.fill(0);
  for(let i=-radius;i<=radius;i++)kernel[2*((i+N)%N)]=Math.exp(-.5*(i*h/sigma)**2)/(Math.sqrt(2*Math.PI)*sigma);
  return buffered?buffered.fft.transform1d(kernel,buffered.kernels[slot]):engineFor(N).fft(kernel);
}

function extractField(inverse,g,d,bs,sigma,N,ix0,iy0,component=0,sign=1,fieldH=g.h){
  const invN2=sign/(N*N),x0=xAt(g,0)+d.rx,y0=yAt(g,0)+d.ry;
  const fx=x0/fieldH-ix0,fy=y0/fieldH-iy0,ox=Math.floor(fx),oy=Math.floor(fy),tx=fx-ox,ty=fy-oy;
  let xs,ys,xt,yt;
  if(fieldH!==g.h){
    xs=new Int32Array(g.N);ys=new Int32Array(g.N);xt=new Float64Array(g.N);yt=new Float64Array(g.N);
    for(let k=0;k<g.N;k++){
      const x=fx+k*g.h/fieldH,y=fy+k*g.h/fieldH;
      xs[k]=Math.floor(x);ys[k]=Math.floor(y);xt[k]=x-xs[k];yt[k]=y-ys[k];
    }
  }
  const support=supportMask(g,d,bs,5*sigma),out=new Float64Array(g.N*g.N);
  for(let y=0;y<g.N;y++)for(let x=0;x<g.N;x++) {
    if(!support.a[y*support.stride+x])continue;
    const j=ys?ys[y]:y+oy,i=2*(xs?xs[x]:x+ox)+component;
    const u=xt?xt[x]:tx,v=yt?yt[y]:ty;
    out[y*g.N+x]=Math.max(0,invN2*((1-v)*((1-u)*inverse[j][i]+u*inverse[j][i+2])+
      v*((1-u)*inverse[j+1][i]+u*inverse[j+1][i+2])));
  }
  return out;
}

function pairedFields(g,first,second,bs,N,shared){
  // K spectra may use backend-owned buffers; preserve K1 before computing K2.
  const buffered=shared.buffered;
  const fieldH=shared.fieldH??g.h;
  const K1=buffered?kernelSpectrum(N,fieldH,first.sigma,buffered,0):new Float32Array(kernelSpectrum(N,fieldH,first.sigma));
  const K2=second?kernelSpectrum(N,fieldH,second.sigma,buffered,1):null;
  const packed=buffered?buffered.work:complexRows(N);
  for(let y=0;y<N;y++)for(let x=0;x<N;x++){
    const re=shared.spectrum[y][2*x],im=shared.spectrum[y][2*x+1];
    const kr1=K1[2*x]*K1[2*y]-K1[2*x+1]*K1[2*y+1];
    const ki1=K1[2*x]*K1[2*y+1]+K1[2*x+1]*K1[2*y];
    let ar=re*kr1-im*ki1,ai=re*ki1+im*kr1;
    if(K2){
      const kr2=K2[2*x]*K2[2*y]-K2[2*x+1]*K2[2*y+1];
      const ki2=K2[2*x]*K2[2*y+1]+K2[2*x+1]*K2[2*y];
      const br=re*kr2-im*ki2,bi=re*ki2+im*kr2;
      ar-=bi;ai+=br; // H = F*K1 + i*F*K2
    }
    packed[y][2*x]=ar;packed[y][2*x+1]=-ai;
  }
  // forward(conj(H))/N^2 = conj(inverse(H)): negate the second output channel.
  const inverse=buffered?buffered.fft.transform(packed,buffered.inverse):engineFor(N).fft2d(packed);
  return [extractField(inverse,g,first.d,bs,first.sigma,N,shared.ix0,shared.iy0,0,1,fieldH),
    second?extractField(inverse,g,second.d,bs,second.sigma,N,shared.ix0,shared.iy0,1,-1,fieldH):null];
}

/** Approximate FFT evaluation of Localization B, preserving its definitions of
 * sigma, observations, and likelihood combination. Rasterization/interpolation
 * introduce numerical differences; existing exact implementations are retained.
 * diagnostics is optional and never contributes to the returned probability grid.
 */
export function localizeBFFT(g,D,bs,{diagnostics=null}={}) {
  if(diagnostics)Object.assign(diagnostics,{fftDetections:0,directDetections:0,fftSizes:[],forwardTransforms:0,inverseTransforms:0});
  if(!D.length)return g;
  const logL=new Float64Array(g.a.length);
  for(const d of D) {
    const sigma=Math.max(d.d/5,g.h/2),pad=Math.ceil(5*sigma/g.h)+2;
    const N=nextPow2(g.N+2*pad+2);
    const direct=sigma/g.h<2||N>MAX_FFT_N;
    const field=direct?directField(g,d,bs,sigma):fftField(g,d,bs,sigma,N,pad);
    if(diagnostics) {
      diagnostics[direct?'directDetections':'fftDetections']++;
      if(!direct){diagnostics.fftSizes.push(N);diagnostics.forwardTransforms++;diagnostics.inverseTransforms++;}
    }
    for(let k=0;k<logL.length;k++)logL[k]=field[k]>0?logL[k]+Math.log(field[k]):-Infinity;
  }
  return posteriorFromLogLikelihood(g,logL);
}

function posteriorFromLogLikelihood(g,logL) {
  const rounded=new Float32Array(logL),out=makeGrid(g.N,g.hM,g.c.x,g.c.y);
  let max=-Infinity;for(const v of rounded)max=Math.max(max,v);
  if(!Number.isFinite(max)){out.a.set(g.a);return out;}
  for(let k=0;k<out.a.length;k++)out.a[k]=g.a[k]*Math.exp(rounded[k]-max);
  if(!normalize(out.a))out.a.set(g.a);
  return out;
}

/** Shared-forward experiment. Reuse one world-aligned boulder spectrum, with
 * one inverse transform for each eligible detection. Retain the baseline when
 * the union would exceed the existing memory cap or no reuse is possible.
 */
export function localizeBFFTShared(g,D,bs,{diagnostics=null}={}) {
  return localizeShared(g,D,bs,diagnostics,false);
}

export function localizeBFFTPaired(g,D,bs,{diagnostics=null}={}) {
  return localizeShared(g,D,bs,diagnostics,true);
}

export function localizeBFFTBuffered(g,D,bs,{diagnostics=null}={}) {
  return localizeShared(g,D,bs,diagnostics,true,true);
}

/** Additional approximation: sample each broad likelihood on a dyadic field
 * lattice no coarser than sigma/samplesPerSigma, then interpolate to the original
 * position cells. Does not change, crop, or threshold the prior grid.
 */
export function localizeBFFTCoarse(g,D,bs,{diagnostics=null,samplesPerSigma=4}={}){
  if(!Number.isFinite(samplesPerSigma)||samplesPerSigma<2)throw new RangeError('samplesPerSigma must be >= 2');
  if(diagnostics)Object.assign(diagnostics,{samplesPerSigma,groups:[],forwardTransforms:0,inverseTransforms:0,fftDetections:0,directDetections:0});
  if(!D.length)return g;
  const groups=new Map(),logL=new Float64Array(g.a.length);
  const accumulate=field=>{
    for(let k=0;k<logL.length;k++)logL[k]=field[k]>0?logL[k]+Math.log(field[k]):-Infinity;
  };
  for(const d of D){
    const sigma=Math.max(d.d/5,g.h/2);
    if(sigma/g.h<2){accumulate(directField(g,d,bs,sigma));if(diagnostics)diagnostics.directDetections++;continue;}
    const power=Math.max(0,Math.floor(Math.log2(sigma/(samplesPerSigma*g.h))));
    const fieldH=g.h*2**power;
    if(!groups.has(fieldH))groups.set(fieldH,[]);
    groups.get(fieldH).push({d,sigma,pad:Math.ceil(5*sigma/fieldH)+2});
  }
  for(const [fieldH,specs] of groups){
    const ix0=Math.min(...specs.map(s=>Math.floor((xAt(g,0)+s.d.rx)/fieldH)-s.pad));
    const iy0=Math.min(...specs.map(s=>Math.floor((yAt(g,0)+s.d.ry)/fieldH)-s.pad));
    const ixEnd=Math.max(...specs.map(s=>Math.ceil((xAt(g,g.N-1)+s.d.rx)/fieldH)+s.pad+2));
    const iyEnd=Math.max(...specs.map(s=>Math.ceil((yAt(g,g.N-1)+s.d.ry)/fieldH)+s.pad+2));
    const N=nextPow2(Math.max(4,ixEnd-ix0,iyEnd-iy0));
    if(N>MAX_FFT_N){
      for(const s of specs)accumulate(directField(g,s.d,bs,s.sigma));
      if(diagnostics)diagnostics.directDetections+=specs.length;
      continue;
    }
    const buffered=bufferedFor(N),shared={ix0,iy0,fieldH,buffered,
      spectrum:boulderSpectrum(bs,N,fieldH,ix0,iy0,buffered)};
    if(diagnostics){
      diagnostics.groups.push({fieldSpacingM:fieldH*1000,N,detections:specs.length});
      diagnostics.forwardTransforms++;diagnostics.inverseTransforms+=Math.ceil(specs.length/2);
      diagnostics.fftDetections+=specs.length;
    }
    for(let i=0;i<specs.length;i+=2){
      const [a,b]=pairedFields(g,specs[i],specs[i+1],bs,N,shared);
      accumulate(a);if(b)accumulate(b);
    }
  }
  return posteriorFromLogLikelihood(g,logL);
}

function localizeShared(g,D,bs,diagnostics,paired,buffered=false){
  const specs=D.map(d=>{
    const sigma=Math.max(d.d/5,g.h/2),pad=Math.ceil(5*sigma/g.h)+2;
    const N=nextPow2(g.N+2*pad+2);
    return {d,sigma,pad,N,direct:sigma/g.h<2||N>MAX_FFT_N};
  });
  const eligible=specs.filter(s=>!s.direct);
  const ix0=Math.min(...eligible.map(s=>Math.floor((xAt(g,0)+s.d.rx)/g.h)-s.pad));
  const iy0=Math.min(...eligible.map(s=>Math.floor((yAt(g,0)+s.d.ry)/g.h)-s.pad));
  const ixEnd=Math.max(...eligible.map(s=>Math.floor((xAt(g,0)+s.d.rx)/g.h)+g.N+s.pad+2));
  const iyEnd=Math.max(...eligible.map(s=>Math.floor((yAt(g,0)+s.d.ry)/g.h)+g.N+s.pad+2));
  const N=eligible.length?nextPow2(Math.max(ixEnd-ix0,iyEnd-iy0)):0;
  if(eligible.length<2||N>MAX_FFT_N){
    const result=localizeBFFT(g,D,bs,{diagnostics});
    if(diagnostics)Object.assign(diagnostics,{usedShared:false,pairedInverse:false,buffered:false,sharedFFTSize:N,
      fallbackReason:eligible.length<2?'fewer than two FFT detections':'shared patch exceeds size limit'});
    return result;
  }
  if(diagnostics)Object.assign(diagnostics,{usedShared:true,pairedInverse:paired,buffered,sharedFFTSize:N,
    fftDetections:eligible.length,directDetections:specs.length-eligible.length,
    fftSizes:eligible.map(()=>N),forwardTransforms:1,inverseTransforms:paired?Math.ceil(eligible.length/2):eligible.length});
  const workspace=buffered?bufferedFor(N):null;
  const shared={ix0,iy0,buffered:workspace,spectrum:boulderSpectrum(bs,N,g.h,ix0,iy0,workspace)};
  const logL=new Float64Array(g.a.length);
  const accumulate=field=>{
    for(let k=0;k<logL.length;k++)logL[k]=field[k]>0?logL[k]+Math.log(field[k]):-Infinity;
  };
  if(paired){
    for(const s of specs)if(s.direct)accumulate(directField(g,s.d,bs,s.sigma));
    for(let i=0;i<eligible.length;i+=2){
      const [a,b]=pairedFields(g,eligible[i],eligible[i+1],bs,N,shared);
      accumulate(a);if(b)accumulate(b);
    }
  }else for(const s of specs){
    accumulate(s.direct?directField(g,s.d,bs,s.sigma):fftField(g,s.d,bs,s.sigma,N,s.pad,shared));
  }
  return posteriorFromLogLikelihood(g,logL);
}

