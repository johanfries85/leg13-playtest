// A separate seeded stream for course generation, independent of boulders.
(() => {
  function createSeededRandom(legNumber,purpose){
    if(!Number.isSafeInteger(legNumber)||legNumber<1)throw new RangeError('Configuration number must be a positive safe integer');
    let state=2166136261;
    for(const ch of `Leg13:${purpose}:${legNumber}`)state=Math.imul(state^ch.charCodeAt(0),16777619)>>>0;
    return ()=>{
      state=(state+0x6d2b79f5)|0;
      let bits=Math.imul(state^(state>>>15),state|1);
      bits^=bits+Math.imul(bits^(bits>>>7),bits|61);
      return ((bits^(bits>>>14))>>>0)/4294967296;
    };
  }
  function createCourseRandom(legNumber){return createSeededRandom(legNumber,'course');}
  function northAngleForConfiguration(number){return -Math.PI+2*Math.PI*createSeededRandom(number,'north')();}
  function mapNumberForLeg(legNumber){
    // Uniform draw among valid map numbers excluding the leg number itself.
    return drawMapNumber(legNumber,createSeededRandom(legNumber,'map'));
  }
  function drawMapNumber(legNumber,random){
    if(legNumber>4294967295)return 1+Math.floor(random()*4294967295);
    const pick=1+Math.floor(random()*4294967294);
    return pick>=legNumber?pick+1:pick;
  }
  function generateCourse(legNumber,minimumLength,maximumLength){
    return sampleCourse(createCourseRandom(legNumber),minimumLength,maximumLength);
  }
  function sampleCourse(random,minimumLength,maximumLength){
    if(!Number.isFinite(minimumLength)||minimumLength<100||minimumLength>1000||!Number.isFinite(maximumLength)||maximumLength<minimumLength||maximumLength>2000)throw new RangeError('Invalid leg length range');
    const lengthM=minimumLength+(maximumLength-minimumLength)*random();
    const angleDeg=-180+360*random(),alpha=angleDeg*Math.PI/180;
    const x=lengthM/2*Math.cos(alpha),y=lengthM/2*Math.sin(alpha);
    return {lengthM,angleDeg,C12:{x:-x,y:-y},C13:{x,y}};
  }
  function nearestBoulder(points,target){
    let nearest=null,distance=Infinity;
    for(const point of points){
      const d=(point.x-target.x)**2+(point.y-target.y)**2;
      if(d<distance){distance=d;nearest=point;}
    }
    return nearest;
  }
  function generateLeg(legNumber,radiusM,minimumLength,maximumLength,getBoulders){
    const courseRandom=createCourseRandom(legNumber),mapRandom=createSeededRandom(legNumber,'map');
    if(!Number.isFinite(radiusM)||radiusM<=0)throw new RangeError('Map radius must be positive and finite');
    let rejection='';
    for(let trial=1;trial<=10;trial++){
      const mapNumber=drawMapNumber(legNumber,mapRandom);
      const course=sampleCourse(courseRandom,minimumLength,maximumLength);
      const points=getBoulders(radiusM,mapNumber);
      const C12=nearestBoulder(points,course.C12),C13=nearestBoulder(points,course.C13);
      if(!C12||!C13){rejection='No boulders available';continue;}
      if(C12===C13){rejection='Both controls selected the same boulder';continue;}
      const lengthM=Math.hypot(C13.x-C12.x,C13.y-C12.y);
      if(lengthM<minimumLength||lengthM>maximumLength){rejection=`Snapped length ${lengthM.toFixed(1)} m is outside ${minimumLength}–${maximumLength} m`;continue;}
      return {number:legNumber,mapNumber,radiusM,points,course,controls:{C12,C13},lengthM,trials:trial};
    }
    throw new Error(`Could not generate a valid leg within 10 trials. ${rejection}. Try more forgiving bounds for the leg length.`);
  }
  function initialMapRadius(widthM,heightM,offsetM){
    if(!Number.isFinite(widthM)||widthM<=0||!Number.isFinite(heightM)||heightM<=0||!Number.isFinite(offsetM)||offsetM<0)throw new RangeError('Invalid map dimensions or center offset');
    // Original 100 m snapping allowance plus the extra 100 m initial margin.
    return Math.hypot(widthM/2,heightM/2+offsetM)+200;
  }
  function rectangleForLeg(leg,widthM,heightM,offsetM){
    const {C12,C13}=leg.controls;
    const legAngle=Math.atan2(C13.y-C12.y,C13.x-C12.x);
    const vx=Math.cos(legAngle),vy=Math.sin(legAngle);
    const point={x:(C12.x+C13.x)/2+offsetM*vx,y:(C12.y+C13.y)/2+offsetM*vy};
    const corners=[[-1,-1],[1,-1],[1,1],[-1,1]].map(([w,h])=>({
      x:point.x+w*widthM/2*vy+h*heightM/2*vx,
      y:point.y-w*widthM/2*vx+h*heightM/2*vy
    }));
    return {point,legAngle,corners,furthestDistanceM:Math.max(...corners.map(p=>Math.hypot(p.x,p.y)))};
  }
  function generateFramedLeg(legNumber,minimumLength,maximumLength,widthM,heightM,offsetM,getBoulders){
    const initialRadiusM=initialMapRadius(widthM,heightM,offsetM);
    const northAngleRad=northAngleForConfiguration(legNumber);
    let radiusM=initialRadiusM,violationM=0;
    for(let radiusTrial=1;radiusTrial<=10;radiusTrial++){
      // generateLeg recreates both streams from the same configuration each time.
      const leg=generateLeg(legNumber,radiusM,minimumLength,maximumLength,getBoulders);
      const frame=rectangleForLeg(leg,widthM,heightM,offsetM);
      violationM=frame.furthestDistanceM-radiusM;
      if(violationM<=0)return {...leg,frame,initialRadiusM,northAngleRad,radiusTrials:radiusTrial};
      radiusM+=2*violationM;
    }
    throw new Error('Could not fit the map rectangle within 10 radius trials. Try smaller map dimensions or a smaller map-center offset.');
  }
  function bouldersInsideRectangle(points,corners){
    const origin=corners[0];
    const wx=corners[1].x-origin.x,wy=corners[1].y-origin.y;
    const hx=corners[3].x-origin.x,hy=corners[3].y-origin.y;
    const width=Math.hypot(wx,wy),height=Math.hypot(hx,hy),epsilon=1e-9;
    return points.filter(p=>{
      const dx=p.x-origin.x,dy=p.y-origin.y;
      const x=(dx*wx+dy*wy)/width,y=(dx*hx+dy*hy)/height;
      return x>=-epsilon&&x<=width+epsilon&&y>=-epsilon&&y<=height+epsilon;
    });
  }
  globalThis.LegCourseGeneration=Object.freeze({createCourseRandom,northAngleForConfiguration,mapNumberForLeg,generateCourse,generateLeg,initialMapRadius,rectangleForLeg,generateFramedLeg,bouldersInsideRectangle});
})();
