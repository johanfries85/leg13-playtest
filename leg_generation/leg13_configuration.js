// Reusable classic-script API. Load boulder_generator.js and course_random.js first.
(() => {
  // Diagnostic world coordinates stay outside the public configuration result.
  // Weak keys let discarded experiment configurations release their source maps.
  const experimentDetails=new WeakMap();
  const CONFIGURATION_OFFSET=74n;
  function sourceConfigurationNumber(configurationNumber){
    if(!Number.isSafeInteger(configurationNumber)||configurationNumber<1)throw new RangeError('Configuration number must be a positive safe integer');
    // Renumber all configurations: new 1 is former 75. BigInt avoids rounding,
    // and wrapping preserves the existing full positive-safe-integer range.
    return Number((BigInt(configurationNumber)-1n+CONFIGURATION_OFFSET)%BigInt(Number.MAX_SAFE_INTEGER)+1n);
  }

  function toRectangleCoordinates(point,frame){
    const origin=frame.corners[0],angle=frame.legAngle;
    const dx=point.x-origin.x,dy=point.y-origin.y;
    return {x:dx*Math.sin(angle)-dy*Math.cos(angle),y:dx*Math.cos(angle)+dy*Math.sin(angle)};
  }

  function get_leg13_configuration(configurationNumber,mapWidth,mapHeight,minimalLegLength,maximalLegLength,mapCenter){
    const generated=LegCourseGeneration.generateFramedLeg(sourceConfigurationNumber(configurationNumber),minimalLegLength,maximalLegLength,mapWidth,mapHeight,mapCenter,get_boulder_positions);
    const leg={...generated,number:configurationNumber};
    const local=point=>toRectangleCoordinates(point,leg.frame);
    const northAngle=leg.northAngleRad+Math.PI/2-leg.frame.legAngle;
    const configuration={
      boulders:LegCourseGeneration.bouldersInsideRectangle(leg.points,leg.frame.corners).map(local),
      C12:local(leg.controls.C12),
      C13:local(leg.controls.C13),
      northDirectionRadians:Math.atan2(Math.sin(northAngle),Math.cos(northAngle))
    };
    experimentDetails.set(configuration,leg);
    return configuration;
  }

  globalThis.Leg13Configuration=Object.freeze({
    get_leg13_configuration,
    sourceConfigurationNumber,
    toRectangleCoordinates,
    getExperimentDetails:configuration=>experimentDetails.get(configuration)
  });
  globalThis.get_leg13_configuration=get_leg13_configuration;
})();
