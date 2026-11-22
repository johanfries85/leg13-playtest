// Experimental adapter for the Indutny backend bundled with WebFFT 1.0.3.
// Isolate this internal dependency here; no node_modules files are modified.
import FFT from "../node_modules/webfft/lib/indutny/fft.js";

export function complexRows(N){return Array.from({length:N},()=>new Float32Array(2*N));}

export class BufferedFFT2D {
  constructor(N){
    if(!Number.isInteger(N)||N<4||(N&(N-1)))throw new RangeError('FFT size must be a power of two >= 4');
    this.N=N;this.plan=new FFT(N);
    this.rows=complexRows(N);
    this.columnIn=new Float32Array(2*N);this.columnOut=new Float32Array(2*N);
  }
  transform1d(input,output){this.plan.transform(output,input);return output;}
  // Caller owns input/output; neither the output nor an immutable spectrum is
  // silently overwritten by a later transform. Scratch storage is reused.
  transform(input,output){
    const N=this.N;
    if(input.length!==N||output.length!==N||input.some(r=>r.length!==2*N)||output.some(r=>r.length!==2*N))
      throw new RangeError('2D transform requires N rows of 2*N interleaved values');
    for(let y=0;y<N;y++)this.plan.transform(this.rows[y],input[y]);
    for(let x=0;x<N;x++){
      for(let y=0;y<N;y++){
        this.columnIn[2*y]=this.rows[y][2*x];this.columnIn[2*y+1]=this.rows[y][2*x+1];
      }
      this.plan.transform(this.columnOut,this.columnIn);
      for(let y=0;y<N;y++){
        output[y][2*x]=this.columnOut[2*y];output[y][2*x+1]=this.columnOut[2*y+1];
      }
    }
    return output;
  }
}
