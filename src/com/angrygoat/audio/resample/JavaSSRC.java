/*******************************************************************************
 * Copyright (c) 2013 Wayne Tam.
 * All rights reserved. This program (except for SplitRadixFft and I0Bessel)
 * is made available under the terms of the GNU Lesser Public License v2.1
 * which accompanies this distribution, and is available at
 * http://www.gnu.org/licenses/old-licenses/gpl-2.0.html
 * 
 * Please the source files of SplitRadixFft and I0Bessel for their respective
 * licenses.
 * 
 * Contributors:
 *     Wayne Tam - initial API and implementation
 ******************************************************************************/
package com.angrygoat.audio.resample;

import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.DataInputStream;
import java.io.EOFException;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.nio.DoubleBuffer;
import java.util.Arrays;
import java.util.Random;

import vavi.util.I0Bessel;
import vavi.util.SplitRadixFft;

public class JavaSSRC {
	
	public interface ProgressListener{
		public void onChanged(double progress);
		public void onShowMessage(String message);
	}
//	private static final boolean DEBUG = false; 

	public static final String VERSION = "1.30";
	private static final int RANDBUFLEN = 65536;

	private static final int[] scoeflen =  new int[]{1,16,20,16,16,15,16,15};
	private static final int[] samp =  new int[]{8,18,27,8,8,8,10,9};
	private static final int[] scoeffreq = new int[]{0,48000,44100,37800,32000,22050,48000,44100};

	private static final double[][] shapercoefs = new double[][]{
		{-1}, /* triangular dither */

		{-2.8720729351043701172,   5.0413231849670410156,  -6.2442994117736816406,   5.8483986854553222656,
			-3.7067542076110839844,   1.0495119094848632812,   1.1830236911773681641,  -2.1126792430877685547,
			1.9094531536102294922,  -0.99913084506988525391,  0.17090806365013122559,  0.32615602016448974609,
			-0.39127644896507263184,  0.26876461505889892578, -0.097676105797290802002, 0.023473845794796943665,
		}, /* 48k, N=16, amp=18 */

		{-2.6773197650909423828,   4.8308925628662109375,  -6.570110321044921875,    7.4572014808654785156,
			-6.7263274192810058594,   4.8481650352478027344,  -2.0412089824676513672,  -0.7006359100341796875,
			2.9537565708160400391,  -4.0800385475158691406,   4.1845216751098632812,  -3.3311812877655029297,
			2.1179926395416259766,  -0.879302978515625,       0.031759146600961685181, 0.42382788658142089844,
			-0.47882103919982910156,  0.35490813851356506348, -0.17496839165687561035,  0.060908168554306030273,
		}, /* 44.1k, N=20, amp=27 */

		{-1.6335992813110351562,   2.2615492343902587891,  -2.4077029228210449219,   2.6341717243194580078,
			-2.1440362930297851562,   1.8153258562088012695,  -1.0816224813461303711,   0.70302653312683105469,
			-0.15991993248462677002, -0.041549518704414367676, 0.29416576027870178223, -0.2518316805362701416,
			0.27766478061676025391, -0.15785403549671173096,  0.10165894031524658203, -0.016833892092108726501,
		}, /* 37.8k, N=16 */

		{-0.82901298999786376953,  0.98922657966613769531, -0.59825712442398071289,  1.0028809309005737305,
			-0.59938216209411621094,  0.79502451419830322266, -0.42723315954208374023,  0.54492527246475219727,
			-0.30792605876922607422,  0.36871799826622009277, -0.18792048096656799316,  0.2261127084493637085,
			-0.10573341697454452515,  0.11435490846633911133, -0.038800679147243499756, 0.040842197835445404053,
		}, /* 32k, N=16 */

		{-0.065229974687099456787, 0.54981261491775512695,  0.40278548002243041992,  0.31783768534660339355,
			0.28201797604560852051,  0.16985194385051727295,  0.15433363616466522217,  0.12507140636444091797,
			0.08903945237398147583,  0.064410120248794555664, 0.047146003693342208862, 0.032805237919092178345,
			0.028495194390416145325, 0.011695005930960178375, 0.011831838637590408325,
		}, /* 22.05k, N=15 */

		{-2.3925774097442626953,   3.4350297451019287109,  -3.1853709220886230469,   1.8117271661758422852,
			0.20124770700931549072, -1.4759907722473144531,   1.7210904359817504883,  -0.97746700048446655273,
			0.13790138065814971924,  0.38185903429985046387, -0.27421241998672485352, -0.066584214568138122559,
			0.35223302245140075684, -0.37672343850135803223,  0.23964276909828186035, -0.068674825131893157959,
		}, /* 48k, N=16, amp=10 */

		{-2.0833916664123535156,   3.0418450832366943359,  -3.2047898769378662109,   2.7571926116943359375,
			-1.4978630542755126953,   0.3427594602108001709,   0.71733748912811279297, -1.0737057924270629883,
			1.0225815773010253906,  -0.56649994850158691406,  0.20968692004680633545,  0.065378531813621520996,
			-0.10322438180446624756,  0.067442022264003753662, 0.00495197344571352005,
		}, /* 44.1k, N=15, amp=9 */
	};

	private static int RINT(double x){
		return ((x) >= 0 ? ((int)((x) + 0.5)) : ((int)((x) - 0.5)));
	}

	private static double[] noiseAmpPresets = new double[]{0.7,0.9,0.18};

	private static final int POOLSIZE = 97;

	private final SplitRadixFft FFT = new SplitRadixFft();

	private ProgressListener listener = null;
	private int srcChannels = 2;
	private int dstChannels = 2;
	private int monoChannel = -1;
	private ByteOrder srcByteOrder = ByteOrder.LITTLE_ENDIAN;
	private int srcBPS = 16;
	private int dstBPS = 16;
	private int srcSamplingRate = 44100;
	private int dstSamplingRate = 44100;
	private double gain = 1f;
	private int ditherType = 0;
	private int pdfType = 0;
	private double noiseAmplitude = 0.18f;
	private boolean twoPass = false;
	private boolean normalize = false;
	private boolean fast = false;
	private String tempFilename = null;
	
	private ResampleContext rCtx = null;
	
	private static class ResampleContext{
		protected int rnch = 2;
		protected int mono = -1;
		protected int nch = 2;
		protected int dnch = 2;
		protected int bps = 16;
		protected int dbps = 16;
		protected int sfrq = 44100;
		protected int dfrq = 44100;
		protected double gain = 1f;
		protected int dither = 0;
		protected int pdf = 0;
		protected double noiseamp = 0.18f;
		protected boolean twopass = false;
		protected boolean normalize = false;
		protected double AA=170; /* stop band attenuation(dB) */
		protected double DF=100;
		protected int FFTFIRLEN=65536;
		protected int ditherSample = 0;
		protected ByteOrder srcByteOrder = ByteOrder.LITTLE_ENDIAN;
		protected String tmpFn = null;
		protected SplitRadixFft FFT; 
		protected double[][] shapebuf;
		protected int shaper_type,shaper_len,shaper_clipmin,shaper_clipmax;
		protected double[] randbuf;
		protected int randptr = 0;
		protected byte[] outBytes = null;

		
		int osf;
		int fs1,fs2;
		int n1,n2;
		int nx,ny,nb,nb2;
		int[] fOrder,fInc;
		int[] fft_ip;
		double[] fft_w, stageA;
		double[][] stageB;
		byte[] rawinbuf,rawoutbuf;
		ByteBuffer inBuffer, outBuffer;
		double[] inbuf,outbuf;
		double[][] buf1,buf2;
		int frqgcd;
		int ip;
		int inbuflen = 0;
		int sp = 0;
		int rps = 0;
		int ds  = 0;
		int rp = 0;
		int delay = 0;
		int osc = 0;
		double peak = 0;
		long sumread = 0;
		long sumwrite = 0;
		boolean init = true;
	}
	
	private static void reset(ResampleContext rCtx)
	{
		rCtx.init = true;
		rCtx.sumread = rCtx.sumwrite = 0;
		rCtx.sp = rCtx.rp = rCtx.rps = rCtx.ds = rCtx.osc = 0;
		rCtx.peak = 0;

		rCtx.fft_ip[0] = 0;
		rCtx.FFT.rdft(rCtx.nb,1,rCtx.stageA,rCtx.fft_ip,rCtx.fft_w);
	    
	    for(int i=0;i<rCtx.nch;i++){
	    	Arrays.fill(rCtx.buf1[i],0);
	    	Arrays.fill(rCtx.buf2[i],0);
	    }
	    rCtx.inBuffer.clear();
	    rCtx.outBuffer.clear();
	    
		if (rCtx.sfrq < rCtx.dfrq){
			rCtx.inbuflen = rCtx.n1/2/(rCtx.fs1/rCtx.sfrq)+1;
			rCtx.delay = (int)((double)rCtx.n2/2/(rCtx.fs2/rCtx.dfrq));
		}else if (rCtx.sfrq > rCtx.dfrq){
			rCtx.inbuflen = 0;
			rCtx.delay = (int) ((double)rCtx.n1/2/((double)rCtx.fs1/rCtx.dfrq)+(double)rCtx.n2/2/((double)rCtx.fs2/rCtx.dfrq));
		}
	}
	
	public byte[] getOutBytes()
	{
		return rCtx.outBytes;
	}
	
	public double getPeak(){
		if(rCtx == null)
			return 0;
		return rCtx.peak;
	}
	
	public void setOnProgressListener(ProgressListener listener){
		this.listener = listener;
	}
	
	public void setFastProfile(boolean fast){
		this.fast = fast;
	}
	
	public void setSrcChannels(int numChannels) {
		this.srcChannels = numChannels;
	}

	public void setDstChannels(int numChannels) {
		this.dstChannels = numChannels;
	}
	
	public void setMonoChannel(int channel){
		this.monoChannel = channel;
	}

	public void setSrcByteOrder(ByteOrder bo) {
		this.srcByteOrder = bo;
	}
	public void setSrcBPS(int srcBPS) {
		if (srcBPS != 8 && srcBPS != 16 && srcBPS != 24 && srcBPS != 32)
			throw new IllegalArgumentException("Src BPS type must be 8, 16, 24, or 32 bits (input: "+srcBPS+")");
		this.srcBPS = srcBPS;
	}

	public void setDstBPS(int dstBPS) {
		if (dstBPS != 8 && dstBPS != 16 && dstBPS != 24)
			throw new IllegalArgumentException("Dst BPS type must be 8, 16, or 24 bits (input: "+dstBPS+")");
		this.dstBPS = dstBPS;
	}

	public void setSrcSamplingRate(int srcSamplingRate) {
		this.srcSamplingRate = srcSamplingRate;
	}

	public void setDstSamplingRate(int dstSamplingRate) {
		this.dstSamplingRate = dstSamplingRate;
	}

	public void setAttenuation(double attenuation) {
		this.gain = dBToGain(-attenuation);
	}

	public void setGain(double gain) {
		this.gain = gain;
	}

	public void setDitherType(int ditherType) {
		if (ditherType < 0 || ditherType > 4)
			throw new IllegalArgumentException("Dither type must be 0, 1, 2, 3, or 4");
		this.ditherType = ditherType;
	}

	public void setPdfType(int pdfType) {
		if (pdfType < 0 || pdfType > 2)
			throw new IllegalArgumentException("PDF type must be 0, 1, or 2");
		this.pdfType = pdfType;
		this.noiseAmplitude = noiseAmpPresets[pdfType];
	}

	public void setNoiseAmplitude(double noiseAmplitude) {
		this.noiseAmplitude = noiseAmplitude;
	}

	public void setTwoPass(boolean twoPass) {
		this.twoPass = twoPass;
	}

	public void setNormalize(boolean normalize) {
		this.normalize = normalize;
	}

	public void setTempFilename(String tempFilename) {
		this.tempFilename = tempFilename;
	}

	public void initialize(){
		rCtx = new ResampleContext();
		rCtx.FFT = FFT;
		rCtx.rnch = srcChannels;
		rCtx.mono = monoChannel < srcChannels?monoChannel:srcChannels-1;
		if(rCtx.rnch > 1 && rCtx.mono > -1)
			rCtx.nch = 1;
		else
			rCtx.nch = rCtx.rnch;
		rCtx.dnch = dstChannels;
		rCtx.srcByteOrder = srcByteOrder;
		rCtx.bps = srcBPS/8;
		rCtx.dbps = dstBPS/8;
		rCtx.sfrq = srcSamplingRate;
		rCtx.dfrq = dstSamplingRate;
		rCtx.gain = gain;
		if(rCtx.bps == rCtx.dbps)
			rCtx.dither = 0;
		else
			rCtx.dither = ditherType;
		rCtx.pdf = pdfType;
		rCtx.noiseamp = noiseAmplitude;
		rCtx.normalize = normalize;
		if(rCtx.sfrq == rCtx.dfrq && rCtx.dither == 0 && rCtx.gain == 1)
			rCtx.twopass = false;
		else{
			if(rCtx.normalize)
				rCtx.twopass = true;
			else
				rCtx.twopass = twoPass;
		}
		rCtx.tmpFn = tempFilename; 
		if(fast){
			rCtx.AA = 96;
			rCtx.DF = 8000;
			rCtx.FFTFIRLEN = 1024;
		}else{
			rCtx.AA=170;
			rCtx.DF=100;
			rCtx.FFTFIRLEN=65536;
		}

		if (rCtx.dither > 0)
			init_shaper();
		
		if (rCtx.sfrq < rCtx.dfrq){
			initUpSample(rCtx);
		}else if (rCtx.sfrq > rCtx.dfrq){
			initDownSample(rCtx);
		}else{
			rCtx.rawinbuf = new byte[16384*rCtx.dbps*rCtx.rnch];
			if(rCtx.twopass)
				rCtx.rawoutbuf = new byte[16384*8*rCtx.nch];
			else
				rCtx.rawoutbuf = new byte[(int)(16384*rCtx.dbps*rCtx.dnch*((double)rCtx.dbps/rCtx.bps))];
			rCtx.inBuffer = ByteBuffer.wrap(rCtx.rawinbuf).order(rCtx.srcByteOrder);
			rCtx.outBuffer = ByteBuffer.wrap(rCtx.rawoutbuf).order(ByteOrder.LITTLE_ENDIAN);
			rCtx.outBytes = new byte[rCtx.rawoutbuf.length];
		}
	}
	
	private void init_shaper()
	{
		int i;
		int[] pool = new int[POOLSIZE];

		for(i=1;i<6;i++) if (rCtx.dfrq == scoeffreq[i]) break;
		if ((rCtx.dither == 3 || rCtx.dither == 4) && i == 6) {
			showMessage(String.format("Warning: ATH based noise shaping for destination frequency %dHz is not available, using triangular dither",rCtx.dfrq));
		}
		if (rCtx.dither == 2 || i == 6) i = 0;
		if (rCtx.dither == 4 && (i == 1 || i == 2)) i += 5;

		rCtx.shaper_type = i;

		rCtx.shaper_len = scoeflen[rCtx.shaper_type];
		rCtx.shapebuf = new double[rCtx.nch][rCtx.shaper_len];

		if (rCtx.dbps == 1) {rCtx.shaper_clipmin = -0x80; rCtx.shaper_clipmax = 0x7f;}
		if (rCtx.dbps == 2) {rCtx.shaper_clipmin = -0x8000; rCtx.shaper_clipmax = 0x7fff;}
		if (rCtx.dbps == 3) {rCtx.shaper_clipmin = -0x800000; rCtx.shaper_clipmax = 0x7fffff;}
		if (rCtx.dbps == 4) {rCtx.shaper_clipmin = -0x80000000; rCtx.shaper_clipmax = 0x7fffffff;}

		rCtx.randbuf = new double[RANDBUFLEN];
		Random random = new Random(System.nanoTime());

		for(i=0;i<POOLSIZE;i++) pool[i] = random.nextInt(Integer.MAX_VALUE); 

		int r1,r2,p;
		double r;

		switch(rCtx.pdf){
		case 0: // rectangular
			for(i=0;i<RANDBUFLEN;i++){
				p = random.nextInt(Integer.MAX_VALUE) % POOLSIZE;
				r1 = pool[p]; pool[p] = random.nextInt(Integer.MAX_VALUE);
				rCtx.randbuf[i] = rCtx.noiseamp * (((double)r1)/Integer.MAX_VALUE-0.5);
			}
			break;
		case 1: // triangular
			for(i=0;i<RANDBUFLEN;i++){
				p = random.nextInt(Integer.MAX_VALUE) % POOLSIZE;
				r1 = pool[p]; pool[p] = random.nextInt(Integer.MAX_VALUE);
				p = random.nextInt(Integer.MAX_VALUE) % POOLSIZE;
				r2 = pool[p]; pool[p] = random.nextInt(Integer.MAX_VALUE);
				rCtx.randbuf[i] = rCtx.noiseamp * ((((double)r1)/Integer.MAX_VALUE)-(((double)r2)/Integer.MAX_VALUE));
			}
			break;
		case 2: // gaussian
			int sw = 0;
			double t = 0, u = 0;

			for(i=0;i<RANDBUFLEN;i++)
			{
				if (sw == 0) {
					sw = 1;

					p = random.nextInt(Integer.MAX_VALUE) % POOLSIZE;
					r = ((double)pool[p])/Integer.MAX_VALUE; pool[p] = random.nextInt(Integer.MAX_VALUE);
					if (r == 1.0) r = 0.0;

					t = Math.sqrt(-2 * Math.log(1-r));

					p = random.nextInt(Integer.MAX_VALUE) % POOLSIZE;
					r = ((double)pool[p])/Integer.MAX_VALUE; pool[p] = random.nextInt(Integer.MAX_VALUE);

					u = 2 * Math.PI * r;

					rCtx.randbuf[i] = rCtx.noiseamp * t * Math.cos(u);
				} else {
					sw = 0;
					rCtx.randbuf[i] = rCtx.noiseamp * t * Math.sin(u);
				}
			}
			break;
		}

		rCtx.randptr = 0;

		if (rCtx.dither == 0 || rCtx.dither == 1)
			rCtx.ditherSample = 1;
		else
			rCtx.ditherSample = samp[rCtx.shaper_type];
	}

	private void showProgress(double p)
	{
		if(listener != null){
			listener.onChanged(p);
		}
	}

	private void showMessage(String msg)
	{
		if(listener != null)
			listener.onShowMessage(msg);
	}

	private static void initUpSample(ResampleContext rCtx){
		int filter2len = rCtx.FFTFIRLEN; /* stage 2 filter length */

		/* Make stage 1 filter */

		double lpf,d,df,alp,iza;
		double guard = 2;
		int i;
		
		rCtx.frqgcd = gcd(rCtx.sfrq,rCtx.dfrq);
		
		rCtx.fs1 = rCtx.sfrq / rCtx.frqgcd * rCtx.dfrq;

		if (rCtx.fs1/rCtx.dfrq == 1) rCtx.osf = 1;
		else if (rCtx.fs1/rCtx.dfrq % 2 == 0) rCtx.osf = 2;
		else if (rCtx.fs1/rCtx.dfrq % 3 == 0) rCtx.osf = 3;
		else {
			throw new UnsupportedOperationException(String.format("Resampling from %dHz to %dHz is not supported.\n" +
					"%d/gcd(%d,%d)=%d must be divisible by 2 or 3.",
					rCtx.sfrq,rCtx.dfrq,rCtx.sfrq,rCtx.sfrq,rCtx.dfrq,rCtx.fs1/rCtx.dfrq));
		}

		df = (rCtx.dfrq*rCtx.osf/2 - rCtx.sfrq/2) * 2 / guard;
		lpf = rCtx.sfrq/2 + (rCtx.dfrq*rCtx.osf/2 - rCtx.sfrq/2)/guard;

		if (rCtx.AA <= 21) d = 0.9222; else d = (rCtx.AA-7.95)/14.36;

		rCtx.n1 = (int)(rCtx.fs1/df*d+1);
		if (rCtx.n1 % 2 == 0) rCtx.n1++;

		alp = alpha(rCtx);
		iza = I0Bessel.value(alp);
		//printf("iza = %g\n",iza);

		rCtx.ny = rCtx.fs1/rCtx.sfrq;
		rCtx.nx = rCtx.n1/rCtx.ny+1;

		rCtx.fOrder = new int[rCtx.ny*rCtx.osf];
		for(i=0;i<rCtx.ny*rCtx.osf;i++) {
			rCtx.fOrder[i] = rCtx.fs1/rCtx.sfrq-(i*(rCtx.fs1/(rCtx.dfrq*rCtx.osf)))%(rCtx.fs1/rCtx.sfrq);
			if (rCtx.fOrder[i] == rCtx.fs1/rCtx.sfrq) rCtx.fOrder[i] = 0;
		}

		rCtx.fInc = new int[rCtx.ny*rCtx.osf];
		for(i=0;i<rCtx.ny*rCtx.osf;i++) {
			rCtx.fInc[i] = rCtx.fOrder[i] < rCtx.fs1/(rCtx.dfrq*rCtx.osf) ? rCtx.nch : 0;
			if (rCtx.fOrder[i] == rCtx.fs1/rCtx.sfrq) rCtx.fOrder[i] = 0;
		}

		rCtx.stageB = new double[rCtx.ny][rCtx.nx];

		for(i=-(rCtx.n1/2);i<=rCtx.n1/2;i++)
		{
			rCtx.stageB[(i+rCtx.n1/2)%rCtx.ny][(i+rCtx.n1/2)/rCtx.ny] = win(i,rCtx.n1,alp,iza)*hn_lpf(i,lpf,rCtx.fs1)*rCtx.fs1/rCtx.sfrq;
		}

		/* Make stage 2 filter */

		int ipsize,wsize;

		if (rCtx.AA <= 21) d = 0.9222; else d = (rCtx.AA-7.95)/14.36;

		rCtx.fs2 = rCtx.dfrq*rCtx.osf;

		for(i=1;;i = i * 2)
		{
			rCtx.n2 = filter2len * i;
			if (rCtx.n2 % 2 == 0) rCtx.n2--;
			df = (rCtx.fs2*d)/(rCtx.n2-1);
			lpf = rCtx.sfrq/2;
			if (df < rCtx.DF) break;
		}

		alp = alpha(rCtx);

		iza = I0Bessel.value(alp);

		for(rCtx.nb=1;rCtx.nb<rCtx.n2;rCtx.nb*=2);
		rCtx.nb *= 2;

		rCtx.stageA = new double[rCtx.nb];

		for(i=-(rCtx.n2/2);i<=rCtx.n2/2;i++) {
			rCtx.stageA[i+rCtx.n2/2] = win(i,rCtx.n2,alp,iza)*hn_lpf(i,lpf,rCtx.fs2)/rCtx.nb*2;
		}

		ipsize    = (int)(2+Math.sqrt(rCtx.nb));
		rCtx.fft_ip    = new int[ipsize];
		wsize     = rCtx.nb/2;
		rCtx.fft_w     = new double[wsize];

		rCtx.FFT.rdft(rCtx.nb,1,rCtx.stageA,rCtx.fft_ip,rCtx.fft_w);

		/* Apply filters */

		rCtx.nb2 = rCtx.nb/2;

		rCtx.buf1 = new double[rCtx.nch][rCtx.nb2/rCtx.osf+1];
		rCtx.buf2 = new double[rCtx.nch][rCtx.nb];

		rCtx.rawinbuf  = new byte[rCtx.rnch*(rCtx.nb2+rCtx.nx)*rCtx.bps];

		if(rCtx.twopass)
			rCtx.rawoutbuf = new byte[rCtx.nch*(rCtx.nb2/rCtx.osf+1)*8];
		else
			rCtx.rawoutbuf = new byte[rCtx.dnch*(rCtx.nb2/rCtx.osf+1)*rCtx.dbps];

		rCtx.inbuf  = new double[rCtx.nch*(rCtx.nb2+rCtx.nx)];
		rCtx.outbuf = new double[rCtx.nch*(rCtx.nb2/rCtx.osf+1)];

		rCtx.inbuflen = rCtx.n1/2/(rCtx.fs1/rCtx.sfrq)+1;
		rCtx.delay = (int)((double)rCtx.n2/2/(rCtx.fs2/rCtx.dfrq));

		rCtx.inBuffer = ByteBuffer.wrap(rCtx.rawinbuf).order(rCtx.srcByteOrder);
		rCtx.outBuffer = ByteBuffer.wrap(rCtx.rawoutbuf).order(ByteOrder.LITTLE_ENDIAN);
	}
	
	private static void initDownSample(ResampleContext rCtx){
		int filter1len = rCtx.FFTFIRLEN; // stage 1 filter length 

		// Make stage 1 filter 

		double lpf,d,df,alp,iza;
		int ipsize,wsize,i;

		rCtx.frqgcd = gcd(rCtx.sfrq,rCtx.dfrq);

		if (rCtx.dfrq/rCtx.frqgcd == 1) rCtx.osf = 1;
		else if (rCtx.dfrq/rCtx.frqgcd % 2 == 0) rCtx.osf = 2;
		else if (rCtx.dfrq/rCtx.frqgcd % 3 == 0) rCtx.osf = 3;
		else {
			throw new UnsupportedOperationException(String.format("Resampling from %dHz to %dHz is not supported.\n" +
					"%d/gcd(%d,%d)=%d must be divisible by 2 or 3.",
					rCtx.sfrq,rCtx.dfrq,rCtx.dfrq,rCtx.sfrq,rCtx.dfrq,rCtx.dfrq/rCtx.frqgcd));
		}

		rCtx.fs1 = rCtx.sfrq*rCtx.osf;

		if (rCtx.AA <= 21) d = 0.9222; else d = (rCtx.AA-7.95)/14.36;

		rCtx.n1 = filter1len;
		for(i=1;;i = i * 2)
		{
			rCtx.n1 = filter1len * i;
			if (rCtx.n1 % 2 == 0) rCtx.n1--;
			df = (rCtx.fs1*d)/(rCtx.n1-1);
			lpf = (rCtx.dfrq-df)/2;
			if (df < rCtx.DF) break;
		}

		alp = alpha(rCtx);

		iza = I0Bessel.value(alp);

		for(rCtx.nb=1;rCtx.nb<rCtx.n1;rCtx.nb*=2){}
		rCtx.nb *= 2;

		rCtx.stageA = new double[rCtx.nb];

		for(i=-(rCtx.n1/2);i<=rCtx.n1/2;i++) {
			rCtx.stageA[i+rCtx.n1/2] = win(i,rCtx.n1,alp,iza)*hn_lpf(i,lpf,rCtx.fs1)*rCtx.fs1/rCtx.sfrq/rCtx.nb*2;
		}

		ipsize    = (int)(2+Math.sqrt(rCtx.nb));
		rCtx.fft_ip    = new int[ipsize];
		rCtx.fft_ip[0] = 0;
		wsize     = rCtx.nb/2;
		rCtx.fft_w     = new double[wsize];

		rCtx.FFT.rdft(rCtx.nb,1,rCtx.stageA,rCtx.fft_ip,rCtx.fft_w);

		// Make stage 2 filter 

		if (rCtx.osf == 1) {
			rCtx.fs2 = rCtx.sfrq/rCtx.frqgcd*rCtx.dfrq;
			rCtx.n2 = 1;
			rCtx.ny = rCtx.nx = 1;
			rCtx.fOrder = new int[rCtx.ny];
			rCtx.fInc = new int[rCtx.ny];
			rCtx.fInc[0] = rCtx.sfrq/rCtx.dfrq;
			rCtx.stageB = new double[rCtx.ny][rCtx.nx];
			rCtx.stageB[0][0] = 1;
		} else {
			double guard = 2;

			rCtx.fs2 = rCtx.sfrq / rCtx.frqgcd * rCtx.dfrq ;

			df = (rCtx.fs1/2 - rCtx.sfrq/2) * 2 / guard;
			lpf = rCtx.sfrq/2 + (rCtx.fs1/2 - rCtx.sfrq/2)/guard;

			if (rCtx.AA <= 21) d = 0.9222; else d = (rCtx.AA-7.95)/14.36;

			rCtx.n2 = (int) (rCtx.fs2/df*d+1);
			if (rCtx.n2 % 2 == 0) rCtx.n2++;

			alp = alpha(rCtx);
			iza = I0Bessel.value(alp);

			rCtx.ny = rCtx.fs2/rCtx.fs1; // 0でないサンプルがfs2で何サンプルおきにあるか？
			rCtx.nx = rCtx.n2/rCtx.ny+1;

			rCtx.fOrder = new int[rCtx.ny];
			for(i=0;i<rCtx.ny;i++) {
				rCtx.fOrder[i] = rCtx.fs2/rCtx.fs1-(i*(rCtx.fs2/rCtx.dfrq))%(rCtx.fs2/rCtx.fs1);
				if (rCtx.fOrder[i] == rCtx.fs2/rCtx.fs1) rCtx.fOrder[i] = 0;
			}

			rCtx.fInc = new int[rCtx.ny];
			for(i=0;i<rCtx.ny;i++) {
				rCtx.fInc[i] = (rCtx.fs2/rCtx.dfrq-rCtx.fOrder[i])/(rCtx.fs2/rCtx.fs1)+1;
				if (rCtx.fOrder[i+1==rCtx.ny ? 0 : i+1] == 0) rCtx.fInc[i]--;
			}

			rCtx.stageB = new double[rCtx.ny][rCtx.nx];
			for(i=-(rCtx.n2/2);i<=rCtx.n2/2;i++){
				rCtx.stageB[(i+rCtx.n2/2)%rCtx.ny][(i+rCtx.n2/2)/rCtx.ny] = win(i,rCtx.n2,alp,iza)*hn_lpf(i,lpf,rCtx.fs2)*rCtx.fs2/rCtx.fs1;
			}
		}

		rCtx.nb2 = rCtx.nb/2;

		rCtx.buf1 = new double[rCtx.nch][rCtx.nb];
		rCtx.buf2 = new double[rCtx.nch][rCtx.nx+1+rCtx.nb2];

		rCtx.rawinbuf = new byte[(rCtx.nb2/rCtx.osf+rCtx.osf+1)*rCtx.rnch*rCtx.bps];
		
		if(rCtx.twopass)
			rCtx.rawoutbuf = new byte[(int)(((double)rCtx.nb2*rCtx.sfrq/rCtx.dfrq+1)*8*rCtx.nch)];
		else
			rCtx.rawoutbuf = new byte[(int)(((double)rCtx.nb2*rCtx.sfrq/rCtx.dfrq+1)*rCtx.dbps*rCtx.dnch)];
		
		rCtx.inbuf = new double[rCtx.nch*(rCtx.nb2/rCtx.osf+rCtx.osf+1)];
		rCtx.outbuf = new double[(int)(rCtx.nch*((double)rCtx.nb2*rCtx.sfrq/rCtx.dfrq+1))];

		rCtx.inBuffer = ByteBuffer.wrap(rCtx.rawinbuf).order(rCtx.srcByteOrder);
		rCtx.outBuffer = ByteBuffer.wrap(rCtx.rawoutbuf).order(ByteOrder.LITTLE_ENDIAN);
	
		rCtx.delay = (int) ((double)rCtx.n1/2/((double)rCtx.fs1/rCtx.dfrq)+(double)rCtx.n2/2/((double)rCtx.fs2/rCtx.dfrq));
	}

	private static int do_shaping(ResampleContext rCtx,double s,int ch)
	{
		double u,h;
		int i;

		if (rCtx.dither == 1) {
			s += rCtx.randbuf[rCtx.randptr++ & (RANDBUFLEN-1)];

			if (s < rCtx.shaper_clipmin) {
				double d = (double)s/rCtx.shaper_clipmin;
				rCtx.peak = rCtx.peak < d ? d : rCtx.peak;
				s = rCtx.shaper_clipmin;
			}
			if (s > rCtx.shaper_clipmax) {
				double d = (double)s/rCtx.shaper_clipmax;
				rCtx.peak = rCtx.peak < d ? d : rCtx.peak;
				s = rCtx.shaper_clipmax;
			}

			return RINT(s);
		}

		h = 0;
		for(i=0;i<rCtx.shaper_len;i++) h += shapercoefs[rCtx.shaper_type][i]*rCtx.shapebuf[ch][i];
		s += h;
		u = s;
		s += rCtx.randbuf[rCtx.randptr++ & (RANDBUFLEN-1)];

		for(i=rCtx.shaper_len-2;i>=0;i--) rCtx.shapebuf[ch][i+1] = rCtx.shapebuf[ch][i];

		if (s < rCtx.shaper_clipmin) {
			double d = (double)s/rCtx.shaper_clipmin;
			rCtx.peak = rCtx.peak < d ? d : rCtx.peak;
			s = rCtx.shaper_clipmin;
			rCtx.shapebuf[ch][0] = s-u;

			if (rCtx.shapebuf[ch][0] >  1) rCtx.shapebuf[ch][0] =  1;
			if (rCtx.shapebuf[ch][0] < -1) rCtx.shapebuf[ch][0] = -1;
		} else if (s > rCtx.shaper_clipmax) {
			double d = (double)s/rCtx.shaper_clipmax;
			rCtx.peak = rCtx.peak < d ? d : rCtx.peak;
			s = rCtx.shaper_clipmax;
			rCtx.shapebuf[ch][0] = s-u;

			if (rCtx.shapebuf[ch][0] >  1) rCtx.shapebuf[ch][0] =  1;
			if (rCtx.shapebuf[ch][0] < -1) rCtx.shapebuf[ch][0] = -1;
		} else {
			s = RINT(s);
			rCtx.shapebuf[ch][0] = s-u;
		}

		return (int)s;
	}

	private static double alpha(ResampleContext rCtx)
	{
		if (rCtx.AA <= 21) return 0;
		if (rCtx.AA <= 50) return 0.5842*Math.pow(rCtx.AA-21,0.4)+0.07886*(rCtx.AA-21);
		return 0.1102*(rCtx.AA-8.7);
	}

	private static double win(double n,int len,double alp,double iza)
	{
		return I0Bessel.value(alp*Math.sqrt(1-4*n*n/(((double)len-1)*((double)len-1))))/iza;
	}

	private static double sinc(double x)
	{
		return x == 0 ? 1 : Math.sin(x)/x;
	}

	private static double hn_lpf(int n,double lpf,double fs)
	{
		double t = 1/fs;
		double omega = 2*Math.PI*lpf;
		return 2*lpf*t*sinc(n*omega*t);
	}

	private static int gcd(int x, int y)
	{
		int t;

		while (y != 0) {
			t = x % y;  x = y;  y = t;
		}
		return x;
	}

	private static int fillInBuf(ResampleContext rCtx, int[][] samples, int offset, int length){
		int ibOffset = rCtx.nch * rCtx.inbuflen;
		int len = length * rCtx.nch;
		switch(rCtx.bps){
		case 1:
			for(int i=0;i<len;i++) {
				rCtx.inbuf[ibOffset+i] = (1/(double)0x7f) * ((samples[i%rCtx.nch][offset+i/rCtx.nch] & 0xff)-128);
			}
			break;
		case 2:
			for(int i=0;i<len;i++) {
				rCtx.inbuf[ibOffset+i] = (1/(double)0x7fff) * samples[i%rCtx.nch][offset+i/rCtx.nch];
			}
			break;
		case 3:
			for(int i=0;i<len;i++) {
				rCtx.inbuf[ibOffset+i] = (1/(double)0x7fffff) * samples[i%rCtx.nch][offset+i/rCtx.nch];
			}
			break;
		case 4:
			for(int i=0;i<len;i++) {
				rCtx.inbuf[ibOffset+i] = (1/(double)0x7fffffff) * samples[i%rCtx.nch][offset+i/rCtx.nch];
			}
			break;
		}
		return length;
	}
	
	private static void fillInBuf(ResampleContext rCtx, int length){
		int offset = rCtx.nch * rCtx.inbuflen;
		int len = length * rCtx.nch;
		int j = 0,jCnt = rCtx.bps;
		if(rCtx.nch == 1 && rCtx.rnch != rCtx.nch){
			j = rCtx.mono*rCtx.bps;
			jCnt = rCtx.rnch*rCtx.bps;
		}
		switch(rCtx.bps){
		case 1:
			for(int i = 0; i < len; i++,j+=jCnt){
				rCtx.inbuf[offset + i] = (1 / (double)0x7f) * 
					(double)(((short)rCtx.inBuffer.get(j) & 0xff) - 128);
			}
			break;
		case 2:
			for(int i=0;i<len;i++,j+=jCnt) {
				rCtx.inbuf[offset+i] = (1/(double)0x7fff) * (double)rCtx.inBuffer.getShort(j);
			}
			break;
		case 3:
			if(rCtx.srcByteOrder == ByteOrder.LITTLE_ENDIAN){
				for(int i=0;i<len;i++,j+=jCnt) {
					rCtx.inbuf[offset+i] = (1/(double)0x7fffff)*
							(double)(((int)rCtx.inBuffer.getShort(j) & 0xffff) | ((int)rCtx.inBuffer.get(j+2) << 24) >> 8);
				}
			}else{
				for(int i=0;i<len;i++) {
					rCtx.inbuf[offset+i] = (1/(double)0x7fffff)*
							(double)((((int)rCtx.inBuffer.get(j) << 24) >> 8) | ((int)rCtx.inBuffer.getShort(j+1) & 0xffff));
				}
			}
			break;
		case 4:
			for(int i=0;i<len;i++,j+=jCnt) {
				rCtx.inbuf[offset+i] = (1/(double)0x7fffffff) * (double)rCtx.inBuffer.getInt(j);
			}
			break;
		}
	}

	private static void fillOutBuf(ResampleContext rCtx, int dbps, double gain, int nsmplwrt){
		int i,j,s;
		double gain2,d;
		if (rCtx.twopass) {
			for(i=0;i<nsmplwrt*rCtx.nch;i++)
			{
				d = rCtx.outbuf[i];
				rCtx.outBuffer.putDouble(d);
				d = Math.abs(d);
				rCtx.peak = rCtx.peak < d ? d : rCtx.peak;
			}
		} else {
			switch(dbps){
			case 1:
				gain2 = gain * (double)0x7f;
				if (rCtx.dither > 0) {
					for(i=0;i<nsmplwrt;i++)
					{
						for(j=0;j<rCtx.dnch;j++){
							s = do_shaping(rCtx,rCtx.outbuf[i*rCtx.nch+(j%rCtx.nch)]*gain2,(j%rCtx.nch));
							rCtx.outBuffer.put((byte)(s + 0x80)); //	((unsigned char *)rawoutbuf)[i] = s + 0x80;
						}
					}
				} else {
					for(i=0;i<nsmplwrt;i++)
					{
						for(j=0;j<rCtx.dnch;j++){
							s = RINT(rCtx.outbuf[i*rCtx.nch+(j%rCtx.nch)]*gain2);

							if (s < -0x80) {
								d = (double)s/-0x80;
								rCtx.peak = rCtx.peak < d ? d : rCtx.peak;
								s = -0x80;
							}
							if (0x7f <  s) {
								d = (double)s/ 0x7f;
								rCtx.peak = rCtx.peak < d ? d : rCtx.peak;
								s =  0x7f;
							}
							rCtx.outBuffer.put((byte)(s + 0x80)); //	((unsigned char *)rawoutbuf)[i] = s + 0x80;
						}
					}
				}
				break;
			case 2:
				gain2 = gain * (double)0x7fff;
				if (rCtx.dither > 0) {
					for(i=0;i<nsmplwrt;i++)	{
						for(j=0;j<rCtx.dnch;j++){
							s = do_shaping(rCtx,rCtx.outbuf[i*rCtx.nch+(j%rCtx.nch)]*gain2,(j%rCtx.nch));
							rCtx.outBuffer.putShort((short)s);
						}
					}
				} else {
					for(i=0;i<nsmplwrt;i++)	{
						for(j=0;j<rCtx.dnch;j++){
							s = RINT(rCtx.outbuf[i*rCtx.nch+(j%rCtx.nch)]*gain2);
		
							if (s < -0x8000) {
								d = (double)s/-0x8000;
								rCtx.peak = rCtx.peak < d ? d : rCtx.peak;
								s = -0x8000;
							}
							if (0x7fff <  s) {
								d = (double)s/ 0x7fff;
								rCtx.peak = rCtx.peak < d ? d : rCtx.peak;
								s = 0x7fff;
							}
							rCtx.outBuffer.putShort((short)s);
						}
					}
				}
				break;
			case 3:
				gain2 = gain * (double)0x7fffff;
				if (rCtx.dither > 0) {
					for(i=0;i<nsmplwrt;i++)
					{
						for(j=0;j<rCtx.dnch;j++){
							s = do_shaping(rCtx,rCtx.outbuf[i*rCtx.nch+(j%rCtx.nch)]*gain2,(j%rCtx.nch));
							rCtx.outBuffer.putShort((short)s);
							s >>= 16;
							rCtx.outBuffer.put((byte)s);
						}
					}
				} else {
					for(i=0;i<nsmplwrt;i++)
					{
						for(j=0;j<rCtx.dnch;j++){
							s = RINT(rCtx.outbuf[i*rCtx.nch+(j%rCtx.nch)]*gain2);
		
							if (s < -0x800000) {
								d = (double)s/-0x800000;
								rCtx.peak = rCtx.peak < d ? d : rCtx.peak;
								s = -0x800000;
							}
							if (0x7fffff <  s) {
								d = (double)s/ 0x7fffff;
								rCtx.peak = rCtx.peak < d ? d : rCtx.peak;
								s =  0x7fffff;
							}
							rCtx.outBuffer.putShort((short)s);
							s >>= 16;
							rCtx.outBuffer.put((byte)s);
						}
					}
				}
				break;
			}
		}
	}

	private static int writeOutBytes(ResampleContext rCtx, int nsmplwrt, int dbps, int offset, boolean isLast){
		int writeLen = 0,writeOffset = 0;
		boolean isDone = false;
		int nch = rCtx.twopass?rCtx.nch:rCtx.dnch;

		if (!rCtx.init) {
			if (isLast) {
				if ((double)rCtx.sumread*rCtx.dfrq/rCtx.sfrq+2 > rCtx.sumwrite+nsmplwrt) {
					writeLen = dbps*nch*nsmplwrt;
				} else {
					writeLen = dbps*nch*(int)(Math.floor((double)rCtx.sumread*rCtx.dfrq/rCtx.sfrq)+2-rCtx.sumwrite);
					reset(rCtx);
					isDone = true;
				}
			} else {
				writeLen = dbps*nch*nsmplwrt;
			}
		} else {
			if (nsmplwrt < rCtx.delay) {
				rCtx.delay -= nsmplwrt;
			} else {
				writeOffset = dbps*nch*rCtx.delay;
				if (isLast) {
					if ((double)rCtx.sumread*rCtx.dfrq/rCtx.sfrq+2 > rCtx.sumwrite+nsmplwrt-rCtx.delay) {
						writeLen = dbps*nch*(nsmplwrt-rCtx.delay);
					} else {
						writeLen = (int)(dbps*nch*(Math.floor((double)rCtx.sumread*rCtx.dfrq/rCtx.sfrq)+2-rCtx.sumwrite-rCtx.delay));
						reset(rCtx);
						isDone = true;
					}
				} else {
					writeLen = dbps*nch*(nsmplwrt-rCtx.delay);
					rCtx.init = false;
				}
			}
		}

		if(writeLen > 0){
			if(rCtx.outBytes == null)
				rCtx.outBytes = new byte[rCtx.outBuffer.limit()];
			if(rCtx.outBytes.length - offset < rCtx.outBuffer.limit()){
				byte[] tmpBytes = new byte[offset + rCtx.outBuffer.limit()];
				System.arraycopy(rCtx.outBytes, 0, tmpBytes, 0, offset);
				rCtx.outBytes = tmpBytes;
			}
			rCtx.outBuffer.position(writeOffset);
			rCtx.outBuffer.get(rCtx.outBytes, offset, writeLen);
		}
			
		return isDone?-1:writeLen;
	}

	private static boolean writeOutStream(ResampleContext rCtx, OutputStream fpo, int nsmplwrt, int dbps, boolean isLast) throws IOException{
		int nch = rCtx.twopass?rCtx.nch:rCtx.dnch;

		if (!rCtx.init) {
			if (isLast) {
				if ((double)rCtx.sumread*rCtx.dfrq/rCtx.sfrq+2 > rCtx.sumwrite+nsmplwrt) {
					fpo.write(rCtx.rawoutbuf,0,dbps*nch*nsmplwrt);
					rCtx.sumwrite += nsmplwrt;
				} else {
                    int limitData = (int) (dbps * nch * (Math.floor((double) rCtx.sumread * rCtx.dfrq / rCtx.sfrq) + 2 - rCtx.sumwrite));
                    if (limitData > 0)
                    	fpo.write(rCtx.rawoutbuf,0,limitData);
					reset(rCtx);
					return true;
				}
			} else {
				fpo.write(rCtx.rawoutbuf,0,dbps*nch*nsmplwrt);
				rCtx.sumwrite += nsmplwrt;
			}
		} else {
			if (nsmplwrt < rCtx.delay) {
				rCtx.delay -= nsmplwrt;
			} else {
				if (isLast) {
					if ((double)rCtx.sumread*rCtx.dfrq/rCtx.sfrq+2 > rCtx.sumwrite+nsmplwrt-rCtx.delay) {
						fpo.write(rCtx.rawoutbuf,dbps*nch*rCtx.delay,dbps*nch*(nsmplwrt-rCtx.delay));
						rCtx.sumwrite += nsmplwrt-rCtx.delay;
					} else {
						fpo.write(rCtx.rawoutbuf,dbps*nch*rCtx.delay,(int)(dbps*nch*(Math.floor((double)rCtx.sumread*rCtx.dfrq/rCtx.sfrq)+2-rCtx.sumwrite-rCtx.delay)));
						reset(rCtx);
						return true;
					}
				} else {
					fpo.write(rCtx.rawoutbuf,dbps*nch*rCtx.delay,dbps*nch*(nsmplwrt-rCtx.delay));
					rCtx.sumwrite += nsmplwrt-rCtx.delay;
					rCtx.init = false;
				}
			}
		}
		return false;
	}
	
	private static int upSample(ResampleContext rCtx, int nsmplwrt1)
	{
		// apply stage 1 filter
		int ch,p,i,j,s1o,no,ip2,nsmplwrt2 = 0;
		int s1p_backup = rCtx.sp;
		int ip_backup  = rCtx.ip;
		int osc_backup = rCtx.osc;
		double re,im,d,tmp;

		for(ch=0;ch<rCtx.nch;ch++)
		{
			no = rCtx.ny*rCtx.osf;

			rCtx.sp = s1p_backup;
			rCtx.ip = ip_backup+ch;

			switch(rCtx.nx)
			{
			case 7:
				for(p=0;p<nsmplwrt1;p++)
				{
					s1o = rCtx.fOrder[rCtx.sp];

					rCtx.buf2[ch][p] =
							rCtx.stageB[s1o][0] * rCtx.inbuf[rCtx.ip+0*rCtx.nch]+
							rCtx.stageB[s1o][1] * rCtx.inbuf[rCtx.ip+1*rCtx.nch]+
							rCtx.stageB[s1o][2] * rCtx.inbuf[rCtx.ip+2*rCtx.nch]+
							rCtx.stageB[s1o][3] * rCtx.inbuf[rCtx.ip+3*rCtx.nch]+
							rCtx.stageB[s1o][4] * rCtx.inbuf[rCtx.ip+4*rCtx.nch]+
							rCtx.stageB[s1o][5] * rCtx.inbuf[rCtx.ip+5*rCtx.nch]+
							rCtx.stageB[s1o][6] * rCtx.inbuf[rCtx.ip+6*rCtx.nch];

					rCtx.ip += rCtx.fInc[rCtx.sp];

					rCtx.sp++;
					if (rCtx.sp == no) rCtx.sp = 0;
				}
				break;
			case 9:
				for(p=0;p<nsmplwrt1;p++)
				{
					s1o = rCtx.fOrder[rCtx.sp];

					rCtx.buf2[ch][p] =
							rCtx.stageB[s1o][0] * rCtx.inbuf[rCtx.ip+0*rCtx.nch]+
							rCtx.stageB[s1o][1] * rCtx.inbuf[rCtx.ip+1*rCtx.nch]+
							rCtx.stageB[s1o][2] * rCtx.inbuf[rCtx.ip+2*rCtx.nch]+
							rCtx.stageB[s1o][3] * rCtx.inbuf[rCtx.ip+3*rCtx.nch]+
							rCtx.stageB[s1o][4] * rCtx.inbuf[rCtx.ip+4*rCtx.nch]+
							rCtx.stageB[s1o][5] * rCtx.inbuf[rCtx.ip+5*rCtx.nch]+
							rCtx.stageB[s1o][6] * rCtx.inbuf[rCtx.ip+6*rCtx.nch]+
							rCtx.stageB[s1o][7] * rCtx.inbuf[rCtx.ip+7*rCtx.nch]+
							rCtx.stageB[s1o][8] * rCtx.inbuf[rCtx.ip+8*rCtx.nch];

					rCtx.ip += rCtx.fInc[rCtx.sp];

					rCtx.sp++;
					if (rCtx.sp == no) rCtx.sp = 0;
				}
				break;
			default:
				for(p=0;p<nsmplwrt1;p++)
				{
					tmp = 0;
					ip2=rCtx.ip;

					s1o = rCtx.fOrder[rCtx.sp];

					for(i=0;i<rCtx.nx;i++)
					{
						tmp += rCtx.stageB[s1o][i] * rCtx.inbuf[ip2];
						ip2 += rCtx.nch;
					}
					rCtx.buf2[ch][p] = tmp;

					rCtx.ip += rCtx.fInc[rCtx.sp];

					rCtx.sp++;
					if (rCtx.sp == no) rCtx.sp = 0;
				}
				break;
			}

			rCtx.osc = osc_backup;

			// apply stage 2 filter

			Arrays.fill(rCtx.buf2[ch], nsmplwrt1, rCtx.nb, 0);

			rCtx.FFT.rdft(rCtx.nb,1,rCtx.buf2[ch],rCtx.fft_ip,rCtx.fft_w);

			rCtx.buf2[ch][0] = rCtx.stageA[0]*rCtx.buf2[ch][0];
			rCtx.buf2[ch][1] = rCtx.stageA[1]*rCtx.buf2[ch][1]; 

			for(i=1;i<rCtx.nb/2;i++)
			{

				re = rCtx.stageA[i*2  ]*rCtx.buf2[ch][i*2] - rCtx.stageA[i*2+1]*rCtx.buf2[ch][i*2+1];
				im = rCtx.stageA[i*2+1]*rCtx.buf2[ch][i*2] + rCtx.stageA[i*2  ]*rCtx.buf2[ch][i*2+1];

				//System.out.println(String.format("%d : %g %g %g %g %g %g\n",i,rCtx.stageA[i*2],rCtx.stageA[i*2+1],rCtx.buf2[ch][i*2],rCtx.buf2[ch][i*2+1],re,im));

				rCtx.buf2[ch][i*2  ] = re;
				rCtx.buf2[ch][i*2+1] = im;
			}

			rCtx.FFT.rdft(rCtx.nb,-1,rCtx.buf2[ch],rCtx.fft_ip,rCtx.fft_w);

			for(i=rCtx.osc,j=0;i<rCtx.nb2;i+=rCtx.osf,j++)
			{
				d = (rCtx.buf1[ch][j] + rCtx.buf2[ch][i]);
				rCtx.outbuf[ch + j*rCtx.nch] = d;
			}

			nsmplwrt2 = j;

			rCtx.osc = i - rCtx.nb2;

			for(j=0;i<rCtx.nb;i+=rCtx.osf,j++)
				rCtx.buf1[ch][j] = rCtx.buf2[ch][i];
		}
		
		return nsmplwrt2;
	}
	
	private static int downSample(ResampleContext rCtx){
		int rps_backup = rCtx.rps;
		int s2p_backup = rCtx.sp;
		int i,j,k,ch,t1,bp,nsmplwrt2 = 0;
		double re,im,tmp;
		int bp2;
		int s2o;


		for( ch=0;ch<rCtx.nch;ch++)
		{
			rCtx.rps = rps_backup;

			
			for(k=0;k<rCtx.rps;k++) rCtx.buf1[ch][k] = 0;

			for(i=rCtx.rps,j=0;i<rCtx.nb2;i+=rCtx.osf,j++)
			{
//				assert(j < ((rCtx.nb2-rCtx.rps-1)/rCtx.osf+1));
				rCtx.buf1[ch][i] = rCtx.inbuf[j*rCtx.nch+ch];
				for(k=i+1;k<i+rCtx.osf;k++) rCtx.buf1[ch][k] = 0;
			}

//			assert(j == ((rCtx.nb2-rCtx.rps-1)/rCtx.osf+1));

			for(k=rCtx.nb2;k<rCtx.nb;k++) rCtx.buf1[ch][k] = 0;
			
			rCtx.rps = i - rCtx.nb2;
			rCtx.FFT.rdft(rCtx.nb,1,rCtx.buf1[ch],rCtx.fft_ip,rCtx.fft_w);
			rCtx.buf1[ch][0] = rCtx.stageA[0]*rCtx.buf1[ch][0];
			rCtx.buf1[ch][1] = rCtx.stageA[1]*rCtx.buf1[ch][1]; 
			for(i=1;i<rCtx.nb2;i++)
			{
				re = rCtx.stageA[i*2  ]*rCtx.buf1[ch][i*2] - rCtx.stageA[i*2+1]*rCtx.buf1[ch][i*2+1];
				im = rCtx.stageA[i*2+1]*rCtx.buf1[ch][i*2] + rCtx.stageA[i*2  ]*rCtx.buf1[ch][i*2+1];

				rCtx.buf1[ch][i*2  ] = re;
				rCtx.buf1[ch][i*2+1] = im;
			}

			rCtx.FFT.rdft(rCtx.nb,-1,rCtx.buf1[ch],rCtx.fft_ip,rCtx.fft_w);
			for(i=0;i<rCtx.nb2;i++) {
				rCtx.buf2[ch][rCtx.nx+1+i] += rCtx.buf1[ch][i];
			}

			t1 = rCtx.rp/(rCtx.fs2/rCtx.fs1);
			if (rCtx.rp%(rCtx.fs2/rCtx.fs1) != 0) t1++;

			bp = t1; // bp = &(buf2[ch][t1]);
			rCtx.sp = s2p_backup;

			for(j=0;bp<rCtx.nb2+1;j++)
			{
				tmp = 0;
				bp2 = bp;
				s2o = rCtx.fOrder[rCtx.sp];
				bp += rCtx.fInc[rCtx.sp];
				rCtx.sp++;

				if (rCtx.sp == rCtx.ny) rCtx.sp = 0;

//				assert(bp2*(rCtx.fs2/rCtx.fs1)-(rCtx.rp+j*(rCtx.fs2/rCtx.dfrq)) == s2o);

				for(i=0;i<rCtx.nx;i++)
					tmp += rCtx.stageB[s2o][i] * rCtx.buf2[ch][bp2++];

				rCtx.outbuf[j*rCtx.nch+ch] = tmp;
			}

			nsmplwrt2 = j;
		}
		return nsmplwrt2;
	}

	private static int upsample(ResampleContext rCtx, int[][] samples, int length, double gain, boolean isLast){
		int nsmplwrt1,nsmplwrt2;
		int writeLen = 0;
		int dbps = rCtx.twopass?8:rCtx.dbps;
		int toberead = (int)Math.floor((double)rCtx.nb2*rCtx.sfrq/(rCtx.dfrq*rCtx.osf))+1+rCtx.nx-rCtx.inbuflen;

		if(length < toberead && !isLast){
			rCtx.inbuflen += fillInBuf(rCtx,samples,0,length);
			return 0;
		}

		if(length == 0 && rCtx.inbuflen > 0 && isLast){
			rCtx.sumread += rCtx.inbuflen;
			nsmplwrt1 = rCtx.nb2;
			rCtx.ip = ((rCtx.sfrq*(rCtx.rp-1)+rCtx.fs1)/rCtx.fs1)*rCtx.nch;
			nsmplwrt2 = upSample(rCtx,nsmplwrt1);
			rCtx.rp += nsmplwrt1 * (rCtx.sfrq / rCtx.frqgcd) / rCtx.osf;
			rCtx.outBuffer.clear();
			fillOutBuf(rCtx, dbps, gain, nsmplwrt2);
			rCtx.outBuffer.flip();
			return writeOutBytes(rCtx,nsmplwrt2,dbps,0,isLast);
		}

		int outBytesWritten = 0;
		int lenUsed = 0;
		
		while(lenUsed < length){
			toberead = (int)Math.floor((double)rCtx.nb2*rCtx.sfrq/(rCtx.dfrq*rCtx.osf))+1+rCtx.nx-rCtx.inbuflen;
			
			if(length - lenUsed < toberead && !isLast){
				rCtx.inbuflen += fillInBuf(rCtx,samples,lenUsed,length - lenUsed);
				return outBytesWritten;
			}
			rCtx.inbuflen += fillInBuf(rCtx,samples,lenUsed,toberead);
			lenUsed += toberead;
			rCtx.sumread += toberead;

			nsmplwrt1 = rCtx.nb2;

			rCtx.ip = ((rCtx.sfrq*(rCtx.rp-1)+rCtx.fs1)/rCtx.fs1)*rCtx.nch;

			nsmplwrt2 = upSample(rCtx,nsmplwrt1);
			
			rCtx.rp += nsmplwrt1 * (rCtx.sfrq / rCtx.frqgcd) / rCtx.osf;

			rCtx.outBuffer.clear();
			fillOutBuf(rCtx, dbps, gain, nsmplwrt2);
			rCtx.outBuffer.flip();

			writeLen = writeOutBytes(rCtx,nsmplwrt2,dbps,outBytesWritten,isLast);
			if(writeLen < 0)
				break;
			outBytesWritten += writeLen;

			rCtx.sumwrite += writeLen/(rCtx.bps*(rCtx.twopass?rCtx.nch:rCtx.dnch));

			rCtx.ds = (rCtx.rp-1)/(rCtx.fs1/rCtx.sfrq);

			System.arraycopy(rCtx.inbuf, rCtx.nch*rCtx.ds, rCtx.inbuf, 0, rCtx.nch*(rCtx.inbuflen-rCtx.ds));

			rCtx.inbuflen -= rCtx.ds;
			rCtx.rp -= rCtx.ds*(rCtx.fs1/rCtx.sfrq);
		}
		return outBytesWritten;
	}
	
	private static int upsample(ResampleContext rCtx, byte[] samples, int offset, int length, double gain, boolean isLast){
		int nsmplread,nsmplwrt1,nsmplwrt2;
		int writeLen = 0;
		int dbps = rCtx.twopass?8:rCtx.dbps;
		int toberead = (int)Math.floor((double)rCtx.nb2*rCtx.sfrq/(rCtx.dfrq*rCtx.osf))+1+rCtx.nx-rCtx.inbuflen;
		
		if(rCtx.inBuffer.position() + length < toberead*rCtx.bps*rCtx.rnch && !isLast){
			rCtx.inBuffer.put(samples, offset, length);
			return 0;
		}

		if(length == 0 && rCtx.inBuffer.hasRemaining() && isLast){
			nsmplread = rCtx.inBuffer.position()/(rCtx.bps*rCtx.rnch);
			rCtx.inBuffer.flip();
			fillInBuf(rCtx,nsmplread);
			Arrays.fill(rCtx.inbuf, rCtx.nch*rCtx.inbuflen+nsmplread*rCtx.nch, rCtx.nch*rCtx.inbuflen+rCtx.nch*toberead, 0);
			rCtx.inBuffer.clear();
			rCtx.inbuflen += toberead;
			rCtx.sumread += nsmplread;
			nsmplwrt1 = rCtx.nb2;
			rCtx.ip = ((rCtx.sfrq*(rCtx.rp-1)+rCtx.fs1)/rCtx.fs1)*rCtx.nch;
			nsmplwrt2 = upSample(rCtx,nsmplwrt1);
			rCtx.rp += nsmplwrt1 * (rCtx.sfrq / rCtx.frqgcd) / rCtx.osf;
			rCtx.outBuffer.clear();
			fillOutBuf(rCtx, dbps, gain, nsmplwrt2);
			rCtx.outBuffer.flip();
			return writeOutBytes(rCtx,nsmplwrt2,dbps,0,isLast);
		}

		int outBytesWritten = 0;
		int lenUsed = 0;
		
		while(lenUsed < length){
			toberead = (int)Math.floor((double)rCtx.nb2*rCtx.sfrq/(rCtx.dfrq*rCtx.osf))+1+rCtx.nx-rCtx.inbuflen;
			nsmplread = toberead*rCtx.bps*rCtx.rnch - rCtx.inBuffer.position();
			if(nsmplread > length - lenUsed){
				nsmplread = length - lenUsed;
			}

			rCtx.inBuffer.put(samples, offset+lenUsed, nsmplread);
			lenUsed += nsmplread;
			
			if(rCtx.inBuffer.position() < toberead*rCtx.bps*rCtx.rnch && !isLast)
				return outBytesWritten;

			rCtx.inBuffer.flip();
			fillInBuf(rCtx,toberead);
			rCtx.inBuffer.clear();

			rCtx.inbuflen += toberead;
			
			rCtx.sumread += toberead;

			nsmplwrt1 = rCtx.nb2;

			rCtx.ip = ((rCtx.sfrq*(rCtx.rp-1)+rCtx.fs1)/rCtx.fs1)*rCtx.nch;

			nsmplwrt2 = upSample(rCtx,nsmplwrt1);
			
			rCtx.rp += nsmplwrt1 * (rCtx.sfrq / rCtx.frqgcd) / rCtx.osf;

			rCtx.outBuffer.clear();
			fillOutBuf(rCtx, dbps, gain, nsmplwrt2);
			rCtx.outBuffer.flip();

			writeLen = writeOutBytes(rCtx,nsmplwrt2,dbps,outBytesWritten,isLast);
			if(writeLen < 0)
				break;
			outBytesWritten += writeLen;

			rCtx.sumwrite += writeLen/(rCtx.bps*(rCtx.twopass?rCtx.nch:rCtx.dnch));

			rCtx.ds = (rCtx.rp-1)/(rCtx.fs1/rCtx.sfrq);

			System.arraycopy(rCtx.inbuf, rCtx.nch*rCtx.ds, rCtx.inbuf, 0, rCtx.nch*(rCtx.inbuflen-rCtx.ds));

			rCtx.inbuflen -= rCtx.ds;
			rCtx.rp -= rCtx.ds*(rCtx.fs1/rCtx.sfrq);
		}		
		return outBytesWritten;
	}
	
	private void upsample(InputStream fpi,OutputStream fpo,double gain,long length) throws IOException
	{
		int spcount = 0;
		boolean ending = false;
		int nsmplread,toberead,toberead2,tmpLen,readLen,nsmplwrt1;
		boolean EOF = false;
		int nsmplwrt2 = 0;
		int dbps = rCtx.twopass?8:rCtx.dbps;
		long chanklen = length/rCtx.bps/rCtx.rnch;

		rCtx.sumread = rCtx.sumwrite = 0;

		for(;;)
		{
			toberead2 = toberead = (int)Math.floor((double)rCtx.nb2*rCtx.sfrq/(rCtx.dfrq*rCtx.osf))+1+rCtx.nx-rCtx.inbuflen;
			if (toberead+rCtx.sumread > chanklen) {
				toberead = (int)(chanklen-rCtx.sumread);
			}

			rCtx.inBuffer.clear();
			nsmplread = 0;
			readLen = rCtx.bps*rCtx.rnch*toberead;
			try{
				while(nsmplread < readLen && !EOF){
					tmpLen = fpi.read(rCtx.rawinbuf,nsmplread,readLen - nsmplread);
					if(tmpLen < 0)
						EOF = true;
					else
						nsmplread += tmpLen;
				}
			}catch(Exception e){
				EOF = true;
			}
			rCtx.inBuffer.limit(nsmplread);
			nsmplread /= rCtx.bps*rCtx.rnch;
			fillInBuf(rCtx,nsmplread);
			Arrays.fill(rCtx.inbuf, rCtx.nch*rCtx.inbuflen+nsmplread*rCtx.nch, rCtx.nch*rCtx.inbuflen+rCtx.nch*toberead2, 0);

			rCtx.inbuflen += toberead2;

			rCtx.sumread += nsmplread;
			ending = EOF || rCtx.sumread >= chanklen;

			//nsmplwrt1 = ((rp-1)*srcSamplingRate/fs1+inbuflen-n1x)*dstSamplingRate*osf/srcSamplingRate;
			//if (nsmplwrt1 > n2b2) nsmplwrt1 = n2b2;
			nsmplwrt1 = rCtx.nb2;

			rCtx.ip = ((rCtx.sfrq*(rCtx.rp-1)+rCtx.fs1)/rCtx.fs1)*rCtx.nch;
			nsmplwrt2 = upSample(rCtx,nsmplwrt1);
			rCtx.rp += nsmplwrt1 * (rCtx.sfrq / rCtx.frqgcd) / rCtx.osf;

			rCtx.outBuffer.clear();
			fillOutBuf(rCtx, dbps, gain, nsmplwrt2);
			rCtx.outBuffer.flip();

			if(writeOutStream(rCtx, fpo, nsmplwrt2, dbps, ending))
				break;

			rCtx.ds = (rCtx.rp-1)/(rCtx.fs1/rCtx.sfrq);
			assert(rCtx.inbuflen >= rCtx.ds);
			System.arraycopy(rCtx.inbuf, rCtx.nch*rCtx.ds, rCtx.inbuf, 0, rCtx.nch*(rCtx.inbuflen-rCtx.ds));
			rCtx.inbuflen -= rCtx.ds;
			rCtx.rp -= rCtx.ds*(rCtx.fs1/rCtx.sfrq);
			if ((spcount++ & 7) == 7) showProgress((double)rCtx.sumread / chanklen);
		}
		showProgress(1);
	}

	private static int downsample(ResampleContext rCtx, int[][] samples, int length, double gain, boolean isLast)
	{
		int nsmplwrt;
		int writeLen = 0;
		int dbps = rCtx.twopass?8:rCtx.dbps;
		int toberead = ((rCtx.nb2-rCtx.rps-1)/rCtx.osf+1);
		
		
		if(rCtx.inbuflen + length < toberead && !isLast){
			rCtx.inbuflen += fillInBuf(rCtx,samples,0,length);
			return 0;
		}

		if(length == 0 && rCtx.inbuflen > 0 && isLast){
			Arrays.fill(rCtx.inbuf, rCtx.inbuflen*rCtx.nch, toberead*rCtx.nch, 0);
			rCtx.sumread += rCtx.inbuflen;
			nsmplwrt = downSample(rCtx);
			rCtx.inbuflen = 0;
			rCtx.rp += nsmplwrt * (rCtx.fs2 / rCtx.dfrq);
			rCtx.outBuffer.clear();
			fillOutBuf(rCtx, dbps, gain, nsmplwrt);
			rCtx.outBuffer.flip();
			return writeOutBytes(rCtx,nsmplwrt,dbps,0,isLast);
		}
		
		int outBytesWritten = 0;
		int lenUsed = 0;
		
		while(lenUsed < length){
			toberead = ((rCtx.nb2-rCtx.rps-1)/rCtx.osf+1) - rCtx.inbuflen;

			if(length - lenUsed < toberead && !isLast){
				rCtx.inbuflen += fillInBuf(rCtx,samples,lenUsed,length - lenUsed);
				return outBytesWritten;
			}
			rCtx.inbuflen += fillInBuf(rCtx,samples,lenUsed,toberead);
			lenUsed += toberead;

			rCtx.sumread += toberead;

			nsmplwrt = downSample(rCtx);
			rCtx.inbuflen = 0;
			rCtx.rp += nsmplwrt * (rCtx.fs2 / rCtx.dfrq);
			
			rCtx.outBuffer.clear();
			fillOutBuf(rCtx, dbps, gain, nsmplwrt);
			rCtx.outBuffer.flip();

			writeLen = writeOutBytes(rCtx,nsmplwrt,dbps,outBytesWritten,isLast);
			if(writeLen < 0)
				break;
			outBytesWritten += writeLen;

			rCtx.sumwrite += writeLen/(rCtx.bps*(rCtx.twopass?rCtx.nch:rCtx.dnch));

			rCtx.ds = (rCtx.rp-1)/(rCtx.fs2/rCtx.fs1);

			if (rCtx.ds > rCtx.nb2) 
				rCtx.ds = rCtx.nb2;

			int ch;
			for(ch=0;ch<rCtx.nch;ch++)
				System.arraycopy(rCtx.buf2[ch], rCtx.ds, rCtx.buf2[ch], 0, rCtx.nx+1+rCtx.nb2-rCtx.ds);

			rCtx.rp -= rCtx.ds*(rCtx.fs2/rCtx.fs1);

			for(ch=0;ch<rCtx.nch;ch++)
				System.arraycopy(rCtx.buf1[ch], rCtx.nb2, rCtx.buf2[ch], rCtx.nx+1, rCtx.nb2);
			
		}
		return outBytesWritten;
	}
	private static int downsample(ResampleContext rCtx, byte[] samples, int offset, int length, double gain, boolean isLast)
	{
		int nsmplread,nsmplwrt;
		int writeLen = 0;
		int dbps = rCtx.twopass?8:rCtx.dbps;
		int toberead = ((rCtx.nb2-rCtx.rps-1)/rCtx.osf+1);
		
		if(rCtx.inBuffer.position() + length < toberead*rCtx.bps*rCtx.rnch && !isLast){
			rCtx.inBuffer.put(samples, offset, length);
			return 0;
		}

		if(length == 0 && rCtx.inBuffer.hasRemaining() && isLast){
			nsmplread = rCtx.inBuffer.position()/(rCtx.bps*rCtx.rnch);
			rCtx.inBuffer.flip();
			fillInBuf(rCtx,nsmplread);
			Arrays.fill(rCtx.inbuf, nsmplread*rCtx.nch, toberead*rCtx.nch, 0);
			rCtx.inBuffer.clear();
			rCtx.sumread += nsmplread;
			nsmplwrt = downSample(rCtx);
			rCtx.rp += nsmplwrt * (rCtx.fs2 / rCtx.dfrq);
			rCtx.outBuffer.clear();
			fillOutBuf(rCtx, dbps, gain, nsmplwrt);
			rCtx.outBuffer.flip();
			return writeOutBytes(rCtx,nsmplwrt,dbps,0,isLast);
		}
		
		int outBytesWritten = 0;
		int lenUsed = 0;
		
		while(lenUsed < length){
			toberead = ((rCtx.nb2-rCtx.rps-1)/rCtx.osf+1);
			nsmplread = toberead * rCtx.bps * rCtx.rnch - rCtx.inBuffer.position();
			if(nsmplread > length - lenUsed){
				nsmplread = length - lenUsed;
			}
			rCtx.inBuffer.put(samples, offset+lenUsed, nsmplread);
			lenUsed += nsmplread;
			
			if(rCtx.inBuffer.position() < toberead*rCtx.bps*rCtx.rnch && !isLast)
				return outBytesWritten;
			
			rCtx.inBuffer.flip();
			fillInBuf(rCtx,toberead);
			rCtx.inBuffer.clear();
			
			rCtx.sumread += toberead;

			nsmplwrt = downSample(rCtx);
			rCtx.rp += nsmplwrt * (rCtx.fs2 / rCtx.dfrq);
			
			rCtx.outBuffer.clear();
			fillOutBuf(rCtx, dbps, gain, nsmplwrt);
			rCtx.outBuffer.flip();

			writeLen = writeOutBytes(rCtx,nsmplwrt,dbps,outBytesWritten,isLast);
			if(writeLen < 0)
				break;
			outBytesWritten += writeLen;

			rCtx.sumwrite += writeLen/(rCtx.bps*(rCtx.twopass?rCtx.nch:rCtx.dnch));

			rCtx.ds = (rCtx.rp-1)/(rCtx.fs2/rCtx.fs1);

			if (rCtx.ds > rCtx.nb2) 
				rCtx.ds = rCtx.nb2;

			int ch;
			for(ch=0;ch<rCtx.nch;ch++)
				System.arraycopy(rCtx.buf2[ch], rCtx.ds, rCtx.buf2[ch], 0, rCtx.nx+1+rCtx.nb2-rCtx.ds);

			rCtx.rp -= rCtx.ds*(rCtx.fs2/rCtx.fs1);

			for(ch=0;ch<rCtx.nch;ch++)
				System.arraycopy(rCtx.buf1[ch], rCtx.nb2, rCtx.buf2[ch], rCtx.nx+1, rCtx.nb2);
		}
		return outBytesWritten;
	}
	
	private void downsample(InputStream fpi,OutputStream fpo, double gain, long length) throws IOException
	{

		int spcount = 0;
		int nsmplwrt2 = 0; 
		boolean ending;
		int ch;
		int dbps = rCtx.twopass?8:rCtx.dbps;
		long chanklen = length/rCtx.bps/rCtx.rnch;

		ending = false;

		rCtx.sumread = rCtx.sumwrite = 0;

		int nsmplread,toberead,readLen,tmpLen;
		boolean EOF = false;

		for(;;)
		{
			toberead = (rCtx.nb2-rCtx.rps-1)/rCtx.osf+1;
			if (toberead+rCtx.sumread > chanklen) {
				toberead = (int) (chanklen-rCtx.sumread);
			}

			rCtx.inBuffer.clear();
			nsmplread = 0;
			readLen = rCtx.bps*rCtx.rnch*toberead;
			try{
				while(nsmplread < readLen && !EOF){
					tmpLen = fpi.read(rCtx.rawinbuf,nsmplread,readLen - nsmplread);
					if(tmpLen < 0)
						EOF = true;
					else
						nsmplread += tmpLen;
				}
			}catch(Exception e){
				EOF = true;
			}
			rCtx.inBuffer.limit(nsmplread);
			nsmplread /= rCtx.bps*rCtx.rnch;
			fillInBuf(rCtx,nsmplread);
			Arrays.fill(rCtx.inbuf, nsmplread*rCtx.nch, rCtx.nch*toberead, 0);

			rCtx.sumread += nsmplread;
			ending = EOF || rCtx.sumread >= chanklen;

			nsmplwrt2 = downSample(rCtx);
            
			rCtx.rp += nsmplwrt2 * (rCtx.fs2 / rCtx.dfrq);
			
			rCtx.outBuffer.clear();
			fillOutBuf(rCtx, dbps, gain, nsmplwrt2);
			rCtx.outBuffer.flip();
			if(writeOutStream(rCtx, fpo, nsmplwrt2, dbps, ending))
				break;
			
			rCtx.ds = (rCtx.rp-1)/(rCtx.fs2/rCtx.fs1);
			if (rCtx.ds > rCtx.nb2) 
				rCtx.ds = rCtx.nb2;
			for(ch=0;ch<rCtx.nch;ch++)
				System.arraycopy(rCtx.buf2[ch], rCtx.ds, rCtx.buf2[ch], 0, rCtx.nx+1+rCtx.nb2-rCtx.ds);
			rCtx.rp -= rCtx.ds*(rCtx.fs2/rCtx.fs1);
			for(ch=0;ch<rCtx.nch;ch++)
				System.arraycopy(rCtx.buf1[ch], rCtx.nb2, rCtx.buf2[ch], rCtx.nx+1, rCtx.nb2);
			if ((spcount++ & 7) == 7)
				showProgress((double)rCtx.sumread / chanklen);
		}

		showProgress(1);
	}

	private static double readFromInBuffer(ResampleContext rCtx, int i){
		switch(rCtx.bps) {
		case 1:
			return (1 / (double)0x7f) * (double)(((short)rCtx.inBuffer.get(i) & 0xff) - 128);
		case 2:
			return (1 / (double)0x7fff) * rCtx.inBuffer.getShort(i);
		case 3:
			if(rCtx.srcByteOrder == ByteOrder.LITTLE_ENDIAN){
				return (1 / (double)0x7fffff) * (double)(((int)rCtx.inBuffer.getShort(i) & 0xffff) | ((int)rCtx.inBuffer.get(i+2) << 24) >> 8);
			}else{
				return (1 / (double)0x7fffff) * (double)((((int)rCtx.inBuffer.get(i) << 24) >> 8) | ((int)rCtx.inBuffer.getShort(i+1) & 0xffff));
			}
		case 4:
			return (1 / (double)0x7fffffff) * (double)rCtx.inBuffer.getInt(i); 
		};
		return 0;
	}
	
	private static double intSampleToDouble(ResampleContext rCtx, int sample){
		switch(rCtx.bps) {
		case 1:
			return (1 / (double)0x7f) * (double)((sample & 0xff) - 128);
		case 2:
			return (1 / (double)0x7fff) * sample;
		case 3:
			return (1 / (double)0x7fffff) * sample;
		case 4:
			return (1 / (double)0x7fffffff) * sample; 
		};
		return 0;
	}
	private static void writeToOutBuffer(ResampleContext rCtx, double f, int ch)
	{
		int s;
		switch(rCtx.dbps) {
		case 1:
			f *= 0x7f;
			s = rCtx.dither > 0? do_shaping(rCtx, f,ch) : RINT(f);
			rCtx.outBuffer.put((byte) (s + 128));
			break;
		case 2:
			f *= 0x7fff;
			s = rCtx.dither > 0? do_shaping(rCtx, f,ch) : RINT(f);
			rCtx.outBuffer.putShort((short)s);
			break;
		case 3:
			f *= 0x7fffff;
			s = rCtx.dither > 0? do_shaping(rCtx, f,ch) : RINT(f);
			rCtx.outBuffer.putShort((short)s);
			s >>= 16;
			rCtx.outBuffer.put((byte)s);
			break;
		}
	}

	private static void writeIntToBuffer(ResampleContext rCtx, int sample, double gain)
	{
		int s = (int)(sample * gain);
		switch(rCtx.dbps) {
		case 1:
			rCtx.outBuffer.put((byte)s);
			break;
		case 2:
			rCtx.outBuffer.putShort((short)s);
			break;
		case 3:
			rCtx.outBuffer.putShort((short)s);
			s >>= 16;
			rCtx.outBuffer.put((byte)s);
			break;
		}
	}
	
	private static int no_src(ResampleContext rCtx, int[][] samples, int length, double gain)
	{
		int i,ch;
		double f = 0,p;
		int len;
		
		int outBytesWritten = 0;
		int lenUsed = 0;

		while(lenUsed < length){
			len = length - lenUsed;
			rCtx.outBuffer.clear();
			if(len > rCtx.outBuffer.limit()/rCtx.nch)
				len = rCtx.outBuffer.limit()/rCtx.nch;
			lenUsed += len;
			
			if (rCtx.twopass) {
				for(i=0;i<len;i++){
					for(ch=0;ch<rCtx.nch;ch++){
						f = intSampleToDouble(rCtx, samples[ch%rCtx.nch][i]) * gain;
						p = f > 0 ? f : -f;
						rCtx.peak = rCtx.peak < p ? p : rCtx.peak;
						rCtx.outBuffer.putDouble(f);
					}
				}
			} else if(rCtx.dbps == rCtx.bps){
				for(i=0;i<len;i++){
					for(ch=0;ch<rCtx.dnch;ch++){
						writeIntToBuffer(rCtx, samples[ch%rCtx.nch][i],gain);
					}
				}
			} else {
				for(i=0;i<len;i++){
					for(ch=0;ch<rCtx.dnch;ch++){
						f = intSampleToDouble(rCtx, samples[ch%rCtx.nch][i]) * gain;
						writeToOutBuffer(rCtx, f,(ch%rCtx.nch));
					}
				}
			}
			rCtx.outBuffer.flip();
			if(rCtx.outBytes.length - outBytesWritten < rCtx.outBuffer.limit()){
				byte[] tmpBytes = new byte[outBytesWritten + rCtx.outBuffer.limit()];
				System.arraycopy(rCtx.outBytes, 0, tmpBytes, 0, outBytesWritten);
				rCtx.outBytes = tmpBytes;
			}
			rCtx.outBuffer.get(rCtx.outBytes, outBytesWritten, rCtx.outBuffer.limit());
			outBytesWritten += rCtx.outBuffer.limit();
		}
		return outBytesWritten;
	}
	
	private static int no_src(ResampleContext rCtx, byte[] samples, int offset, int length, double gain)
	{
		int i,ch;
		double f = 0,p;
		int len = length;

		if(len >= rCtx.inBuffer.remaining())
			len = rCtx.inBuffer.remaining();
		
		if(rCtx.inBuffer.position() + len < rCtx.bps*rCtx.nch){
			rCtx.inBuffer.put(samples, offset, len);
			return 0;
		}
		
		int outBytesWritten = 0;
		int lenUsed = 0;

		int j = 0;
		if(rCtx.nch == 1 && rCtx.rnch != rCtx.nch)
			j = rCtx.mono;
		
		while(lenUsed < length){
			len = rCtx.inBuffer.remaining();
			if(len > length - lenUsed)
				len = length - lenUsed;
			
			rCtx.inBuffer.put(samples, offset+lenUsed, len);

			if(rCtx.inBuffer.position() < rCtx.bps*rCtx.nch)
				break;

			rCtx.inBuffer.flip();
			rCtx.outBuffer.clear();
			
			lenUsed += len;

			if (rCtx.twopass) {
				for(i=0;i<rCtx.inBuffer.limit()-rCtx.bps*rCtx.rnch;i+=rCtx.bps*rCtx.rnch){
						for(ch=0;ch<rCtx.nch;ch++){
							f = readFromInBuffer(rCtx, i+(ch+j)*rCtx.bps) * gain;
							p = f > 0 ? f : -f;
							rCtx.peak = rCtx.peak < p ? p : rCtx.peak;
							rCtx.outBuffer.putDouble(f);
						}
				}
			} else {
				for(i=0;i<rCtx.inBuffer.limit()-rCtx.bps*rCtx.rnch;i+=rCtx.bps*rCtx.rnch){
					for(ch=0;ch<rCtx.dnch;ch++){
						f = readFromInBuffer(rCtx, i+((ch%rCtx.nch)+j)*rCtx.bps) * gain;
						writeToOutBuffer(rCtx, f,(ch%rCtx.nch));
					}
				}
			}
			rCtx.inBuffer.position(i);
			rCtx.inBuffer.compact();
			rCtx.outBuffer.flip();
			if(rCtx.outBytes.length - outBytesWritten < rCtx.outBuffer.limit()){
				byte[] tmpBytes = new byte[outBytesWritten + rCtx.outBuffer.limit()];
				System.arraycopy(rCtx.outBytes, 0, tmpBytes, 0, outBytesWritten);
				rCtx.outBytes = tmpBytes;
			}
			rCtx.outBuffer.get(rCtx.outBytes, outBytesWritten, rCtx.outBuffer.limit());
			outBytesWritten += rCtx.outBuffer.limit();
		}
		return outBytesWritten;
	}
	
	private void no_src(InputStream fpi, OutputStream fpo,double gain,long length) throws IOException
	{
		int ch=0,sumread=0,readLen,tmpLen;
		double f = 0,p;
		long chunklen = length/rCtx.bps/rCtx.rnch;
		int j = 0;
		if(rCtx.nch == 1 && rCtx.rnch != rCtx.nch)
			j = rCtx.mono;
		
		while(sumread < chunklen)
		{
			try{
				rCtx.inBuffer.clear();
				f = 0;
				readLen = 0;
				while(readLen < rCtx.bps*rCtx.rnch){
					tmpLen = fpi.read(rCtx.rawinbuf, readLen, rCtx.bps*rCtx.rnch - readLen);
					if(tmpLen < 0)
						throw new EOFException();
					readLen += tmpLen;
				}
			}catch(EOFException e){
			}
			rCtx.outBuffer.clear();

			if (rCtx.twopass) {
				for(ch=0;ch<rCtx.nch;ch++){
					f = readFromInBuffer(rCtx, (ch+j)*rCtx.bps) * gain;
					p = f > 0 ? f : -f;
					rCtx.peak = rCtx.peak < p ? p : rCtx.peak;
					rCtx.outBuffer.putDouble(f);
				}
				fpo.write(rCtx.rawoutbuf,0,8*rCtx.nch);
			} else {
				for(ch=0;ch<rCtx.dnch;ch++){
					f = readFromInBuffer(rCtx, (ch%rCtx.nch+j)*rCtx.bps) * gain;
					writeToOutBuffer(rCtx,f,(ch%rCtx.nch));
				}
				fpo.write(rCtx.rawoutbuf,0,rCtx.dbps*rCtx.dnch);
			}

			sumread++;

			if ((sumread & 0x3ffff) == 0) showProgress((double)sumread / chunklen);
		}

		fpo.flush();
		showProgress(1);
	}

	public void resetShaper()
	{
		rCtx.randptr = 0;
	}
	
	public void doubleToIntSample(double d, int ch, ByteBuffer ob)
	{
		int s;
		switch(rCtx.dbps) {
		case 1:
			s = rCtx.dither > 0? do_shaping(rCtx,d,ch) : RINT(d);
			ob.put((byte) (s + 128));
			break;
		case 2:
			s = rCtx.dither > 0? do_shaping(rCtx,d,ch) : RINT(d);
			ob.putShort((short) s);
			break;
		case 3:
			s = rCtx.dither > 0? do_shaping(rCtx,d,ch) : RINT(d);
			ob.putShort((short) s);
			s >>= 16;
			ob.put((byte)s);
			break;
		}
	}
	
	public double calcSecondPassGain()
	{
		if (!rCtx.normalize) {
			if (rCtx.peak < rCtx.gain)
				rCtx.peak = 1;
			else 
				rCtx.peak *= Math.pow(10,-Math.log10(rCtx.gain));
		} else 
			rCtx.peak *= Math.pow(10,-Math.log10(rCtx.gain));

		if (rCtx.dither > 0) {
			switch(rCtx.dbps)
			{
			case 1:
				return (rCtx.normalize || rCtx.peak >= (0x7f-rCtx.ditherSample)/(double)0x7f) ? 1/rCtx.peak*(0x7f-rCtx.ditherSample) : 1/rCtx.peak*0x7f;
			case 2:
				return (rCtx.normalize || rCtx.peak >= (0x7fff-rCtx.ditherSample)/(double)0x7fff) ? 1/rCtx.peak*(0x7fff-rCtx.ditherSample) : 1/rCtx.peak*0x7fff;
			case 3:
				return (rCtx.normalize || rCtx.peak >= (0x7fffff-rCtx.ditherSample)/(double)0x7fffff) ? 1/rCtx.peak*(0x7fffff-rCtx.ditherSample) : 1/rCtx.peak*0x7fffff;
			}
		} else {
			switch(rCtx.dbps)
			{
			case 1:
				return 1/rCtx.peak * 0x7f;
			case 2:
				return 1/rCtx.peak * 0x7fff;
			case 3:
				return 1/rCtx.peak * 0x7fffff;
			}
		}
		return 1;
	}

	public static double dBToGain(double att)
	{
		return Math.pow(10,att/20);
	}
	
	public int resample(int[][] samples, int length, boolean isLast){
		if(rCtx == null)
			throw new IllegalStateException("Resampler has not been initialized");

		if (rCtx.sfrq < rCtx.dfrq)
			return upsample(rCtx,samples,length,rCtx.gain, isLast);
		else if (rCtx.sfrq > rCtx.dfrq) 
			return downsample(rCtx,samples,length,rCtx.gain, isLast);
		else
			return no_src(rCtx,samples,length,rCtx.gain);
	}

	public int resample(byte[] samples, int offset, int length, boolean isLast){
		if(rCtx == null)
			throw new IllegalStateException("Resampler has not been initialized");

		if (rCtx.sfrq < rCtx.dfrq)
			return upsample(rCtx,samples,offset,length,rCtx.gain, isLast);
		else if (rCtx.sfrq > rCtx.dfrq) 
			return downsample(rCtx,samples,offset,length,rCtx.gain, isLast);
		else 
			return no_src(rCtx,samples,offset,length,rCtx.gain);
	}
	
	public double resample(InputStream fpi, OutputStream fpo, long length) throws IOException{
		if(rCtx == null)
			throw new IllegalStateException("Resampler has not been initialized");
		
		if (rCtx.twopass) {
			File tmpFile;
			FileOutputStream fpt;
			if (rCtx.tmpFn != null) {
				tmpFile = new File(rCtx.tmpFn);
			} else {
				tmpFile = File.createTempFile("JavaSSRC_", null);
				tmpFile.deleteOnExit();
			}
			fpt = new FileOutputStream(tmpFile);
			
			showMessage("Pass 1");
			try{
				if (rCtx.normalize) {
					if (rCtx.sfrq < rCtx.dfrq) 
						upsample(fpi,fpt,1,length);
					else if (rCtx.sfrq > rCtx.dfrq)
						downsample(fpi,fpt,1,length);
					else 
						no_src(new BufferedInputStream(fpi),new BufferedOutputStream(fpt),1,length);
				} else {
					if (rCtx.sfrq < rCtx.dfrq)
						upsample(fpi,fpt,rCtx.gain,length);
					else if (rCtx.sfrq > rCtx.dfrq)
						downsample(fpi,fpt,rCtx.gain,length);
					else 
						no_src(new BufferedInputStream(fpi),new BufferedOutputStream(fpt),rCtx.gain,length);
				}

				fpt.close();
			}catch(IOException e){
				throw new IOException("Error processing audio data",e);
			}

			showMessage(String.format("\npeak : %gdB",20*Math.log10(rCtx.peak)));

			showMessage("\nPass 2");
			
			double gain = calcSecondPassGain();

			resetShaper();

			long fptlen = tmpFile.length() / (8 * rCtx.nch);
			int sumread = 0;
			int ch = 0;
			DataInputStream inStrm = null;
			BufferedOutputStream outStrm = null;
			try {
				inStrm = new DataInputStream(new BufferedInputStream(new FileInputStream(tmpFile)));
				outStrm = new BufferedOutputStream(fpo);
				byte[] ibuf = new byte[8*rCtx.nch];
				ByteBuffer bb = ByteBuffer.wrap(ibuf).order(ByteOrder.LITTLE_ENDIAN);
				DoubleBuffer db = bb.asDoubleBuffer();
				double f;

				for(;;)
				{
					try{
						inStrm.readFully(ibuf); 
					}catch(EOFException e){
						break;
					}
					bb.clear();
					for(ch=0;ch<rCtx.dnch;ch++){
						f = db.get(ch%rCtx.nch) * gain;
						doubleToIntSample(f,ch%rCtx.nch,bb);
					}
					outStrm.write(bb.array(), 0, bb.position());
					sumread++;

					if ((sumread & 0x3ffff) == 0) listener.onChanged((double)sumread / fptlen);
				}
				showProgress(1);
			} catch (FileNotFoundException e1) {
				throw new FileNotFoundException("Error opening temp file");
			} catch (IOException e) {
				throw new IOException("Error processing temp file",e);
			}finally{
				try{
					if(outStrm != null)
						outStrm.flush();
				}catch(Exception e){}
				try{
					if(inStrm != null)
						inStrm.close();
				}catch(Exception e){}
				if (tmpFile != null) {
					try{
						tmpFile.delete();
					}catch(Exception e){
						showMessage(String.format("Failed to delete temp file %s",rCtx.tmpFn));
					}
				}
			}
		}else{
			BufferedOutputStream outStrm = new BufferedOutputStream(fpo);
			try {
				if (rCtx.sfrq < rCtx.dfrq)
					upsample(fpi,outStrm,rCtx.gain,length);
				else if (rCtx.sfrq > rCtx.dfrq) 
					downsample(fpi,outStrm,rCtx.gain,length);
				else 
					no_src(new BufferedInputStream(fpi),outStrm,rCtx.gain,length);
			} catch (IOException e) {
				throw new IOException("Error processing audio data",e);
			}finally{
				outStrm.flush();
			}
		}
		return rCtx.peak;
	}
}
