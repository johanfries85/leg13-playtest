export function readEnvironmentCapabilities(browser=window){
 const finePointer=browser.matchMedia?.("(any-pointer: fine)");
 return {
  supportsTouch:browser.navigator?.maxTouchPoints>0,
  hasPrimaryFinePointer:Boolean(browser.matchMedia?.("(pointer: fine)").matches),
  get hasFinePointer(){return finePointer?.matches},
  onFinePointerChange(listener){finePointer?.addEventListener?.("change",listener)}
 };
}
