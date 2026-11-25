"use strict";
// Shared map/course generation and rendering. Probability fields and game/debug state stay outside.
const MapCourse=(()=>{
 // Rectangle-corner coordinates from get_leg13_configuration, in 1000 m/unit.
 const GAME_PARAMETERS=Object.freeze({mapWidth:800,mapHeight:1700,minimalLegLength:500,maximalLegLength:1000,mapCenter:150});
 const MAP_REGION=Object.freeze({xmin:0,xmax:.8,ymin:0,ymax:1.7});
 const DISPLAY_X_MIN=MAP_REGION.xmin, DISPLAY_X_MAX=MAP_REGION.xmax;
 const DISPLAY_Y_MIN=MAP_REGION.ymin, DISPLAY_Y_MAX=MAP_REGION.ymax;
 const LEG_CENTER_Y=(GAME_PARAMETERS.mapHeight/2-GAME_PARAMETERS.mapCenter)/1000;
 const C12={x:0,y:0}, C13={x:0,y:1};
 const COURSE_PURPLE='rgb(166,38,255)';
 function rng(seed){let s=seed>>>0;return()=>{s=(1664525*s+1013904223)>>>0;return s/4294967296}}
 function northDirectionToOrientationDeg(northDirectionRadians){
  // get_leg13_configuration uses mathematical radians: zero along +x, increasing
  // counterclockwise. The game uses degrees clockwise from up:
  // orientationDeg = -180 * (northDirectionRadians - pi/2) / pi.
  const degrees=-180*(northDirectionRadians-Math.PI/2)/Math.PI;
  return ((degrees%360)+360)%360;
 }
 function generateExternalLegs(index){
  if(!Number.isSafeInteger(index)||index<1)throw new RangeError('Configuration must be a positive safe integer');
  // External C11/C14 lines have their own stream, independent of terrain/north.
  let seed=2166136261;
  const sourceIndex=Leg13Configuration.sourceConfigurationNumber(index);
  for(const ch of `Leg13:external:${sourceIndex}`)seed=Math.imul(seed^ch.charCodeAt(0),16777619)>>>0;
  const random=rng(seed);
  return {angle12Deg:90+60*random(),side12:random()<.5?-1:1,
   angle13Deg:90+60*random(),side13:random()<.5?-1:1};
 }
 function generateConfiguration(index){
  if(!Number.isSafeInteger(index)||index<1)throw new RangeError('Configuration must be a positive safe integer');
  const p=GAME_PARAMETERS;
  const generated=get_leg13_configuration(index,p.mapWidth,p.mapHeight,p.minimalLegLength,p.maximalLegLength,p.mapCenter);
  const scaled=point=>Object.freeze({x:point.x/1000,y:point.y/1000});
  const c12=scaled(generated.C12),c13=scaled(generated.C13);
  // The public terrain already contains the two control boulders. Keep those
  // out of this list because activeBoulders adds the shared control objects.
  const atControl=b=>Math.hypot(b.x-c12.x,b.y-c12.y)<1e-12||Math.hypot(b.x-c13.x,b.y-c13.y)<1e-12;
  const boulders=Object.freeze(generated.boulders.map(scaled).filter(b=>!atControl(b)));
  return Object.freeze({index,C12:c12,C13:c13,boulders,
   northDirectionRadians:generated.northDirectionRadians,
   orientationDeg:northDirectionToOrientationDeg(generated.northDirectionRadians),
   legLengthM:1000*Math.hypot(c13.x-c12.x,c13.y-c12.y),...generateExternalLegs(index)});
 }
 let configuration;
 function setConfiguration(index){
  if(configuration?.index===index)return;
  const next=generateConfiguration(index); // A failure leaves the current game intact.
  configuration=next;
  Object.assign(C12,next.C12);Object.assign(C13,next.C13);
  rulerPlacement=null;
 }
 function activeBoulders(bouldersOn=true){const a=[C12,C13];if(bouldersOn)a.push(...configuration.boulders);return a}
 function rot(v,a){return {x:v.x*Math.cos(a)-v.y*Math.sin(a),y:v.x*Math.sin(a)+v.y*Math.cos(a)}}
 const RAY12_OUT={x:0,y:1};
 const RAY13_IN={x:0,y:-1};
 function rayToBoundary(c,v){let ts=[];if(v.x>0)ts.push((DISPLAY_X_MAX-c.x)/v.x);if(v.x<0)ts.push((DISPLAY_X_MIN-c.x)/v.x);if(v.y>0)ts.push((DISPLAY_Y_MAX-c.y)/v.y);if(v.y<0)ts.push((DISPLAY_Y_MIN-c.y)/v.y);return Math.min(...ts.filter(t=>t>0))}
 function unitOppositeBisector(a,b){let x=a.x+b.x,y=a.y+b.y,n=Math.hypot(x,y);return {x:-x/n,y:-y/n}}
 function svgEl(name,attrs={},text=''){const e=document.createElementNS('http://www.w3.org/2000/svg',name);for(const [k,v] of Object.entries(attrs))e.setAttribute(k,v);if(text)e.textContent=text;return e}
 function createNorthArrow(x,y,width,height){
  const arrow=svgEl('svg',{x,y,width,height,viewBox:'0 0 161 295',fill:'#ff2428',role:'img','aria-label':'North arrow'});
  // Original paths from the supplied north-arrow.svg, shared by both axes.
  arrow.appendChild(svgEl('path',{d:'M78 20 145 100H94V285H66V100H11Z'}));
  arrow.appendChild(svgEl('path',{d:'M11 109 121 203V111H147V257L38 166V255H11Z'}));
  return arrow;
 }
function drawBoulders(parent, boulders, transform) {
  const group = svgEl('g', { fill: '#000' });
  parent.appendChild(group);

  for (const b of boulders) {
    const p = transform(b);

    group.appendChild(svgEl('circle', {
      cx: p.x,
      cy: p.y,
      r: 3
    }));
  }

  return group;
}
 function sxm(x){return (x-DISPLAY_X_MIN)/(DISPLAY_X_MAX-DISPLAY_X_MIN)*800}
 function sym(y){return (DISPLAY_Y_MAX-y)/(DISPLAY_Y_MAX-DISPLAY_Y_MIN)*1700}
 function northUpGroup(cx,cy){return svgEl('g',{transform:`rotate(${configuration.orientationDeg} ${cx} ${cy})`})}
 // Cache placement for this generated course, independently of terrain and redraws.
 let rulerPlacement=null;
 setConfiguration(1);
 function bottomRightRulerPlacement(ruler){
  const box=ruler.getBBox(),angle=configuration.orientationDeg*Math.PI/180;
  const corners=[];
  for(const x of [box.x,box.x+box.width])for(const y of [box.y,box.y+box.height]){
   corners.push(rot({x,y},angle));
  }
  const right=Math.max(...corners.map(point=>point.x));
  const bottom=Math.max(...corners.map(point=>point.y));
  // Reserve the rightmost 48 units for the configuration number, plus a 12-unit gap.
  return {x:800-48-12-right,y:1700-24-bottom};
 }
 function render(mapVector,{bouldersOn=true,beforeBoulders=null}={}){
  const MAP_ORIENTATION_DEG=configuration.orientationDeg;
  const RAY12_IN=rot(RAY12_OUT,configuration.side12*configuration.angle12Deg*Math.PI/180);
  const RAY13_OUT=rot(RAY13_IN,configuration.side13*configuration.angle13Deg*Math.PI/180);
  mapVector.replaceChildren();
  const gC=svgEl('g',{'fill':'none','stroke':COURSE_PURPLE,'stroke-width':'5.2','stroke-linecap':'butt'});mapVector.appendChild(gC);
  const stop=47.5;
  const p12={x:sxm(C12.x),y:sym(C12.y)},p13={x:sxm(C13.x),y:sym(C13.y)};
  gC.appendChild(svgEl('line',{x1:p12.x,y1:p12.y-stop,x2:p13.x,y2:p13.y+stop}));
  for(const [c,ray,pix] of [[C12,RAY12_IN,p12],[C13,RAY13_OUT,p13]]){const t=rayToBoundary(c,ray),end={x:c.x+t*ray.x,y:c.y+t*ray.y};const ex=sxm(end.x),ey=sym(end.y),dx=ex-pix.x,dy=ey-pix.y,n=Math.hypot(dx,dy);gC.appendChild(svgEl('line',{x1:pix.x+dx/n*stop,y1:pix.y+dy/n*stop,x2:ex,y2:ey}))}
  for(const c of [C12,C13])gC.appendChild(svgEl('circle',{cx:sxm(c.x),cy:sym(c.y),r:37.5}));
  const b12=unitOppositeBisector(RAY12_IN,RAY12_OUT),b13=unitOppositeBisector(RAY13_IN,RAY13_OUT),numDist=108;
  const gN=svgEl('g',{'fill':COURSE_PURPLE,'font-family':'Arial, sans-serif','font-size':'85.725','text-anchor':'middle','dominant-baseline':'middle'});mapVector.appendChild(gN);
  const n12x=p12.x+b12.x*numDist,n12y=p12.y-b12.y*numDist,n13x=p13.x+b13.x*numDist,n13y=p13.y-b13.y*numDist;
  gN.appendChild(svgEl('text',{x:n12x,y:n12y,transform:`rotate(${MAP_ORIENTATION_DEG} ${n12x} ${n12y})`},'12'));
  gN.appendChild(svgEl('text',{x:n13x,y:n13y,transform:`rotate(${MAP_ORIENTATION_DEG} ${n13x} ${n13y})`},'13'));
  const rulerCx=0,rulerCy=0,rw=100,rh=24,labelGap=7;
  const rulerAnchor=svgEl('g');mapVector.appendChild(rulerAnchor);
  const rg=northUpGroup(rulerCx,rulerCy);rulerAnchor.appendChild(rg);
  const rulerLines=svgEl('g',{'stroke':'#ff2428','stroke-width':'2.7','fill':'none'});rg.appendChild(rulerLines);
  rulerLines.appendChild(svgEl('line',{x1:rulerCx-rw/2,y1:rulerCy-rh/2,x2:rulerCx-rw/2,y2:rulerCy+rh/2}));
  rulerLines.appendChild(svgEl('line',{x1:rulerCx+rw/2,y1:rulerCy-rh/2,x2:rulerCx+rw/2,y2:rulerCy+rh/2}));
  rulerLines.appendChild(svgEl('line',{x1:rulerCx-rw/2,y1:rulerCy,x2:rulerCx+rw/2,y2:rulerCy}));
  rg.appendChild(svgEl('text',{x:rulerCx,y:rulerCy+rh/2+labelGap,'dominant-baseline':'hanging','font-family':'Arial, sans-serif','font-size':'29.1','text-anchor':'middle','fill':'#ff2428','stroke':'none'},'100 m'));
  rg.appendChild(createNorthArrow(-28,-111.74,56,103));
  if(!rulerPlacement)rulerPlacement=bottomRightRulerPlacement(rg);
  rulerAnchor.setAttribute('transform',`translate(${rulerPlacement.x} ${rulerPlacement.y})`);
  if(beforeBoulders)beforeBoulders(mapVector);
const visibleBoulders = activeBoulders(bouldersOn).filter(b =>
  b.x >= DISPLAY_X_MIN &&
  b.x <= DISPLAY_X_MAX &&
  b.y >= DISPLAY_Y_MIN &&
  b.y <= DISPLAY_Y_MAX
);

drawBoulders(mapVector, visibleBoulders, b => ({
  x: sxm(b.x),
  y: sym(b.y)
}));
  // Screen-aligned, reading bottom-to-top; longer numbers expand without truncation.
  mapVector.appendChild(svgEl('text',{
   x:780,y:1676,transform:'rotate(-90 780 1676)',
   'font-family':'Arial, sans-serif','font-size':'24',
   style:'font-variant-numeric:tabular-nums',fill:'#666',
   stroke:'#fff','stroke-width':'5','stroke-linejoin':'round','paint-order':'stroke',
   'data-configuration-number':''
  },String(configuration.index).padStart(6,'0')));
 }
 // Development overlay only: independent of the active course and hidden runner.
 function drawWorstCaseOverlay(target,shiftM=0){
  const red='#ff2428',strokeWidth=5.2,stop=47.5;
  const overlay=svgEl('g',{'data-worst-case-leg':'','pointer-events':'none'});
  target.appendChild(overlay);
  // Positive slider values move only this diagnostic leg downward.
  const centerY=LEG_CENTER_Y-shiftM/1000;
  const centerX=GAME_PARAMETERS.mapWidth/2000;
  const controls=[{x:centerX,y:centerY-.5},{x:centerX,y:centerY+.5}];
  // At C13, 90 degrees puts the complementary bisector furthest north:
  // 108*cos(theta/2) is maximal at the smallest allowed theta.
  const rays=[rot(RAY12_OUT,-Math.PI/2),rot(RAY13_IN,Math.PI/2)];
  const lines=svgEl('g',{fill:'none',stroke:red,'stroke-width':strokeWidth});overlay.appendChild(lines);
  lines.appendChild(svgEl('line',{x1:sxm(centerX),y1:sym(controls[0].y)-stop,x2:sxm(centerX),y2:sym(controls[1].y)+stop}));
  const numbers=[];
  for(let i=0;i<2;i++){
   const c=controls[i],ray=rays[i],cx=sxm(c.x),cy=sym(c.y);
   const t=rayToBoundary(c,ray),ex=sxm(c.x+t*ray.x),ey=sym(c.y+t*ray.y);
   const length=Math.hypot(ex-cx,ey-cy);
   lines.appendChild(svgEl('line',{x1:cx+(ex-cx)*stop/length,y1:cy+(ey-cy)*stop/length,x2:ex,y2:ey}));
   lines.appendChild(svgEl('circle',{cx,cy,r:37.5}));
   const bisector=unitOppositeBisector(ray,i===0?RAY12_OUT:RAY13_IN);
   const x=cx+108*bisector.x,y=cy-108*bisector.y;
   const text=svgEl('text',{x,y,fill:red,'font-family':'Arial, sans-serif','font-size':85.725,
    'text-anchor':'middle','dominant-baseline':'middle'},String(12+i));
   overlay.appendChild(text);numbers.push({text,x,y});
  }
  // Rotate the shared number orientation so C13's most distant text-box corner
  // points straight up. This is a conservative bound using browser font metrics.
  const top=numbers[1],box=top.text.getBBox();
  let radius=-1,angle=0;
  for(const x of [box.x,box.x+box.width])for(const y of [box.y,box.y+box.height]){
   const dx=x-top.x,dy=y-top.y,r=Math.hypot(dx,dy);
   if(r>radius){radius=r;angle=-90-Math.atan2(dy,dx)*180/Math.PI;}
  }
  for(const n of numbers){
   n.text.setAttribute('transform',`rotate(${angle} ${n.x} ${n.y})`);
  }
  return overlay;
 }
return {
  get MAP_ORIENTATION_DEG(){return configuration.orientationDeg},
  get configuration(){return configuration},
  generateConfiguration,
  generateExternalLegs,
  northDirectionToOrientationDeg,
  setConfiguration,
  MAP_REGION,
  GAME_PARAMETERS,
  DISPLAY_X_MIN, DISPLAY_X_MAX,
  DISPLAY_Y_MIN, DISPLAY_Y_MAX,
  C12, C13,
  activeBoulders,
  render,
  svgEl,
  createNorthArrow,
  drawBoulders,
  drawWorstCaseOverlay,
  sxm, sym
};
})();
// Preserve the existing game's constant names without duplicating their definitions.
const {DISPLAY_X_MIN,DISPLAY_X_MAX,DISPLAY_Y_MIN,DISPLAY_Y_MAX,C12,C13}=MapCourse;
globalThis.MapCourse = MapCourse;
