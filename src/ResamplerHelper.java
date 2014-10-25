/*******************************************************************************
 * Copyright (c) 2013, Wayne Tam
 * All rights reserved.
 * 
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 * 
 *     Redistributions of source code must retain the above copyright notice,
 *     this list of conditions and the following disclaimer.
 *     Redistributions in binary form must reproduce the above copyright notice,
 *     this list of conditions and the following disclaimer in the documentation
 *     and/or other materials provided with the distribution.
 * 
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 ******************************************************************************/
import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.DataInputStream;
import java.io.EOFException;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.nio.DoubleBuffer;
import java.nio.channels.FileChannel;

/*
import org.kc7bfi.jflac.FLACDecoder;
import org.kc7bfi.jflac.FrameListener;
import org.kc7bfi.jflac.frame.Frame;
import org.kc7bfi.jflac.metadata.Metadata;
import org.kc7bfi.jflac.metadata.StreamInfo;
*/

import com.angrygoat.audio.resample.JavaSSRC;


public class ResamplerHelper{
	public static boolean quiet = false;
	public static String srcFile,destFile,tmpFile = null;
	public static boolean twoPass,normalize,isFast;
	public static int dither,pdf;
	public static ByteOrder byteOrder = ByteOrder.LITTLE_ENDIAN;;
	public static int srcChannels,dstChannels,srcBPS;
	public static int srcSamplingRate,dstSamplingRate,dstBPS;
	public static double attentuation,noiseAmp;
	public static long length = 0;
	public static int monoChannel = -1;
	
	public static double processPCM(FileInputStream fpi, FileOutputStream fpo, JavaSSRC.ProgressListener listener) throws IOException
	{
		double peak = 0;
		
		JavaSSRC javaSSRC = new JavaSSRC();
		// Configure the resampler
		
		// Number of source & destination channels. Useful for going from mono to stereo. 
		javaSSRC.setSrcChannels(srcChannels);
		javaSSRC.setDstChannels(dstChannels);
		
		// Use a single channel from the source audio for every output channel.
		// 0 = 1st channel, 1 = 2nd channel, etc...
		// -1 to disable (default).
		javaSSRC.setMonoChannel(monoChannel);
		
		// Source byte order. WAV=little endian, AIFF=big endian
		// Output is always little endian.
		javaSSRC.setSrcByteOrder(byteOrder); 
		
		
		javaSSRC.setSrcBPS(srcBPS*8); // Source bits per sample. 8, 16, 24, or 32
		javaSSRC.setDstBPS(dstBPS*8); // Destination bits per sample. 8, 16, or 24
		
		// Source & Destination sampling rates. 22050, 44100, 48000, etc...
		// For upsampling:  src rate/gcd(src rate,dst rate) must be divisible by 2 or 3.
		// For downsampling:  dst rate/gcd(src rate,dst rate) must be divisible by 2 or 3.
		javaSSRC.setSrcSamplingRate(srcSamplingRate); 
		javaSSRC.setDstSamplingRate(dstSamplingRate);   


		// Dithering type. 0=none, 1=no noise shaping, 2=triangular spectral shape,
		// 3=ATH based noise shaping, 4=less dither amplitude than type 3
		javaSSRC.setDitherType(dither); 

		// p.d.f. of dither noise. 0=rectangular, 1=triangular, 2=Gaussian
		javaSSRC.setPdfType(pdf); 
		
		// Amplitude of dither noise. Choosing a pdf automatically sets a default amplitude.
		// To set a custom amplitude, make sure you call this after calling setPdfType().
		javaSSRC.setNoiseAmplitude(noiseAmp); 
		
		// Attenuation in dB or you can use setGain(...) gain = Math.pow(10,-attentuation/20)
		javaSSRC.setAttenuation(attentuation); 
		
		// Two pass mode to avoid clipping. Normalizing the audio automatically enables two pass.
		// For the byte array version of resample(...), you will need to manually do the second
		// pass. See below.
		javaSSRC.setTwoPass(twoPass); 
		javaSSRC.setNormalize(normalize); 
		
		// Name of temp file for two pass mode. Set to null for default.
		// Only used with the stream version of resample(...)
		javaSSRC.setTempFilename(tmpFile);
		
		
		javaSSRC.setFastProfile(isFast); // True for faster processing, but lower quality. 

		if(!quiet){
			// Set to null to disable listening for progress.
			// Progress is only reported by the stream version of resample(...)
			javaSSRC.setOnProgressListener(listener);
		}
		
		javaSSRC.initialize(); // Resampler must be initialized after changing settings or before first use.

		DataInputStream inStrm = null;
		try {
			/*
			 * The following section uses the byte array version of resample(...) specifically to
			 * demonstrate how to use it. 
			 * 
			 * For processing an InputStream from beginning to end, the stream version of resample(...)
			 * would be much simpler and more efficient.
			 */
			inStrm = new DataInputStream(new BufferedInputStream(fpi));
			byte[] buffer = new byte[4096];
			byte[] outBytes;
			int readLen,totalWritten = 0;
			double dLength = 0;

			BufferedOutputStream outStrm;
			FileOutputStream fpt = null;
			File file = null;
			if(javaSSRC.isTwoPass()){ // Two Pass mode
				if (tmpFile != null) {
					file = new File(tmpFile);
				} else {
					file = File.createTempFile("JavaSSRC_", null);
					file.deleteOnExit();
				}
				fpt = new FileOutputStream(file);
				outStrm = new BufferedOutputStream(fpt);
				
				listener.onShowMessage("Pass 1");
				
				// When two pass is enabled each sample is written as a double = 8 bytes
				if(ResamplerHelper.monoChannel > -1 && ResamplerHelper.srcChannels > 1){
					dLength = Math.floor(length*((double)dstSamplingRate/srcSamplingRate)*((double)8/srcBPS)/srcChannels);
				}else
					dLength = Math.floor(length*((double)dstSamplingRate/srcSamplingRate)*((double)8/srcBPS));
			}else{
				outStrm = new BufferedOutputStream(fpo);
				dLength = Math.floor(length*((double)dstSamplingRate/srcSamplingRate)*((double)dstBPS/srcBPS)*((double)dstChannels/srcChannels));
			}
			
			while(true){
				readLen = inStrm.read(buffer);
				if(readLen < 0){
					// Setting the last param of resample(...) to true tells the resampler to process any
					// remaining data in its internal buffers and then reset itself.
					// 
					// Does not affect the settings and intialize() does NOT need to be called.
					readLen = javaSSRC.resample(buffer, 0, 0, true); 
					outBytes = javaSSRC.getOutBytes();
					if(readLen > 0)
						outStrm.write(outBytes, 0, readLen);
					listener.onChanged(1.0);
					break;
				}
				readLen = javaSSRC.resample(buffer, 0, readLen, false);
				outBytes = javaSSRC.getOutBytes(); // If readLen = 0, getOutBytes() may return null.
				if(readLen > 0){
					outStrm.write(outBytes, 0, readLen);
					totalWritten += readLen;
				}
				listener.onChanged((double)totalWritten/dLength);
			}
			inStrm.close();
			peak = javaSSRC.getPeak();
			
			if(file != null){  // Second pass
				srcChannels = 1;
				fpt.close();
				listener.onShowMessage(String.format("\npeak : %gdB",20*Math.log10(peak)));

				listener.onShowMessage("\nPass 2");
				
				double gain = javaSSRC.calcSecondPassGain();

				javaSSRC.resetShaper();
				
				long fptlen = file.length() / (8 * srcChannels);
				int sumread = 0;
				int ch = 0;
				try {
					inStrm = new DataInputStream(new BufferedInputStream(new FileInputStream(file)));
					outStrm = new BufferedOutputStream(fpo);
					ByteBuffer bb = ByteBuffer.wrap(buffer).order(ByteOrder.LITTLE_ENDIAN);
					DoubleBuffer db = bb.asDoubleBuffer();
					double f;

					for(;;)
					{
						try{
							inStrm.readFully(buffer,0,8 * srcChannels); 
						}catch(EOFException e){
							break;
						}
						bb.clear();
						// First pass doesn't convert number of channels.
						// Convert to destination number of channels here.
						for(ch=0;ch<dstChannels;ch++){
							f = db.get(ch%srcChannels) * gain;
							javaSSRC.doubleToIntSample(f,ch%srcChannels,bb);
						}
						outStrm.write(bb.array(), 0, bb.position());
						sumread++;

						if ((sumread & 0x3ffff) == 0) listener.onChanged((double)sumread / fptlen);
					}
					listener.onChanged(1);
				} catch (FileNotFoundException e1) {
					throw new FileNotFoundException("Error opening temp file");
				} catch (IOException e) {
					throw new IOException("Error processing temp file",e);
				}finally{
					if (file != null) {
						try{
							file.delete();
						}catch(Exception e){
							listener.onShowMessage(String.format("Failed to delete temp file %s",tmpFile));
						}
					}
				}
			}
			/* end byte array resample(...) section */
			
			
			/*
			 *  Stream version of resample(...)
			 */
			//peak = javaSSRC.resample(fpi, fpo, length);
		}finally{
			try {
				if(inStrm != null)
					inStrm.close();
			} catch (IOException e) {}
		}
		return peak;
	}
	
/*	public static int totalWritten = 0;
	
	public static double processFlac(FileInputStream fpi, FileOutputStream fpo, final JavaSSRC.ProgressListener listener) throws IOException
	{
		final JavaSSRC javaSSRC = new JavaSSRC();
		final FLACDecoder flacDec = new FLACDecoder(fpi);
		double peak = 0;
		totalWritten = 0;
		
		try {
			javaSSRC.setSrcByteOrder(ByteOrder.LITTLE_ENDIAN); 
			javaSSRC.setSrcChannels(srcChannels);
			javaSSRC.setDstChannels(dstChannels);
			javaSSRC.setMonoChannel(monoChannel);
			javaSSRC.setSrcBPS(srcBPS*8);
			javaSSRC.setDstBPS(dstBPS*8);
			javaSSRC.setSrcSamplingRate(srcSamplingRate); 
			javaSSRC.setDstSamplingRate(dstSamplingRate);   
			javaSSRC.setDitherType(dither);  
			javaSSRC.setPdfType(pdf); 
			javaSSRC.setFastProfile(isFast);
			javaSSRC.setTwoPass(twoPass);  
			javaSSRC.setNormalize(normalize);  

			if(!quiet){
				javaSSRC.setOnProgressListener(listener);
			}
			javaSSRC.initialize();
			
			BufferedOutputStream outStrm;
			double dLen;
			File file = null;
			FileOutputStream fpt = null;
			if(twoPass){ // Two Pass mode
				if (tmpFile != null) {
					file = new File(tmpFile);
				} else {
					file = File.createTempFile("JavaSSRC_", null);
					file.deleteOnExit();
				}
				fpt = new FileOutputStream(file);
				outStrm = new BufferedOutputStream(fpt);
				
				listener.onShowMessage("Pass 1");
				
				// When two pass is enabled each sample is written as a double = 8 bytes
				if(ResamplerHelper.monoChannel > -1 && ResamplerHelper.srcChannels > 1){
					dLen = Math.floor(length*((double)dstSamplingRate/srcSamplingRate)*((double)8/srcBPS)/srcChannels);
				}else
					dLen = Math.floor(length*((double)dstSamplingRate/srcSamplingRate)*((double)8/srcBPS));
			}else{
				outStrm = new BufferedOutputStream(fpo);
				dLen = Math.floor(length*((double)dstSamplingRate/srcSamplingRate)*((double)dstBPS/srcBPS)*((double)dstChannels/srcChannels));
			}
			
			final BufferedOutputStream fOut = outStrm;
			final double dLength = dLen;
			
			flacDec.addFrameListener(new FrameListener(){

				@Override
				public void processMetadata(Metadata metadata) {
				}

				@Override
				public void processFrame(Frame frame) {
					int[][] samples;
					if(monoChannel > -1 && srcChannels > 1){
						samples = new int[1][];
						samples[0] = flacDec.getChannelData()[monoChannel].getOutput();
					}else{
						samples = new int[srcChannels][];
						for(int i=0;i<srcChannels;i++){
							samples[i] = flacDec.getChannelData()[i].getOutput();
						}
					}
						
					int readLen = javaSSRC.resample(samples, samples[0].length, false);
					if(readLen > 0)
						try {
							fOut.write(javaSSRC.getOutBytes(), 0, readLen);
							totalWritten += readLen;
							listener.onChanged((double)totalWritten/dLength);
						} catch (IOException e) {
								// TODO Auto-generated catch block
								e.printStackTrace(System.err);
						}
				}

				@Override
				public void processError(String msg) {
					System.err.println(msg);
				}
			});

			flacDec.decode();
			if(flacDec.isEOF()){
				int readLen = javaSSRC.resample(null, 0, true);
				if(readLen > 0)
					try {
						fOut.write(javaSSRC.getOutBytes(), 0, readLen);
					} catch (IOException e) {
						// TODO Auto-generated catch block
						e.printStackTrace(System.err);
					}
			}
			listener.onChanged(1);
			fpi.close();
			peak = javaSSRC.getPeak();

			if(file != null){  // Second pass
				srcChannels = 1;
				fpt.close();
				listener.onShowMessage(String.format("\npeak : %gdB",20*Math.log10(peak)));

				listener.onShowMessage("\nPass 2");
				
				double gain = javaSSRC.calcSecondPassGain();

				javaSSRC.resetShaper();
				
				long fptlen = file.length() / (8 * srcChannels);
				int sumread = 0;
				int ch = 0;
				DataInputStream inStrm = new DataInputStream(new BufferedInputStream(new FileInputStream(file)));
				try {
					outStrm = new BufferedOutputStream(fpo);
					byte[] buffer = new byte[8 * srcChannels];
					ByteBuffer bb = ByteBuffer.wrap(buffer).order(ByteOrder.LITTLE_ENDIAN);
					DoubleBuffer db = bb.asDoubleBuffer();
					double f;

					for(;;)
					{
						try{
							inStrm.readFully(buffer,0,8 * srcChannels); 
						}catch(EOFException e){
							break;
						}
						bb.clear();
						// First pass doesn't convert number of channels.
						// Convert to destination number of channels here.
						for(ch=0;ch<dstChannels;ch++){
							f = db.get(ch%srcChannels) * gain;
							javaSSRC.doubleToIntSample(f,ch%srcChannels,bb);
						}
						outStrm.write(bb.array(), 0, bb.position());
						sumread++;

						if ((sumread & 0x3ffff) == 0) listener.onChanged((double)sumread / fptlen);
					}
					listener.onChanged(1);
				} catch (FileNotFoundException e1) {
					throw new FileNotFoundException("Error opening temp file");
				} catch (IOException e) {
					throw new IOException("Error processing temp file",e);
				}finally{
					inStrm.close();
					if (file != null) {
						try{
							file.delete();
						}catch(Exception e){
							listener.onShowMessage(String.format("Failed to delete temp file %s",tmpFile));
						}
					}
				}
			}
		} catch (IOException e) {
			System.err.println("cannot processing flac file.");
			e.printStackTrace(System.err);
		}
		return peak;
	}//*/
	
	public static double extendedToDouble(byte[] buf) throws NumberFormatException
	{
		ByteBuffer b = ByteBuffer.wrap(buf,0,10).order(ByteOrder.BIG_ENDIAN);
		int e = b.getShort();
		double sign = ((e & 0x80) != 0)?-1:1;
		e &= 0x7fff;
		long m = b.getLong();
		
		double normal = (m & 0x8000000000000000L) != 0?1:0;
		m &= 0x7fffffffffffffffL;
		
		if(m == 0)
			return 0;
		else
		{
			if(e == 0)
				e = -16382;
			else if((e & 0x7fff) == 0x7fff)
				throw new NumberFormatException("NAN or Infinity");
			else
				e -= 16383;
		}
		// Shifting the divisor by 62 bits and then dividing the result by 2 deals
		// with the fact that java does not having unsign longs
		return sign * (normal + (double)m /(1L<<62)/2) * Math.pow(2.0, (float)e);
	}

	private static boolean readAiffHeader(DataInputStream inStrm, FileChannel fc){
		byte[] buf = new byte[10];
		boolean isAIFC = false;
		try {
			inStrm.readInt();

			inStrm.readFully(buf,0,4);
			if("AIFC".equals(new String(buf, 0,4)))
				isAIFC = true;
			else if(!"AIFF".equals(new String(buf, 0,4)))
				return false;
			
			String ckID;
			int ckSize;
			int offset;
			long dataPos = 0;
			int numSampleFrames = 0;
			long vTimestamp;
			boolean foundOne = false;
			double freq;
			while(true){
				inStrm.readFully(buf,0,4);
				ckSize = inStrm.readInt();

				ckID = new String(buf, 0,4);	
				if("COMM".equals(ckID)){
					srcChannels = inStrm.readShort();
					numSampleFrames = inStrm.readInt();
					srcBPS = inStrm.readShort();
					inStrm.readFully(buf,0,10);
					freq = extendedToDouble(buf);
					srcSamplingRate = (int)freq;
					if(isAIFC){
						ckID = new String(buf, 0,4);	
						if(!"NONE".equals(ckID)){
							System.err.println("Compressed AIFF files are not supported");
							return false;
						}
						if(ckSize - 22 > 0)
							inStrm.skip(ckSize -22);
					}
					if(foundOne)
						break;
					foundOne = true;
				}else if ("SSND".equals(ckID)){
					offset = inStrm.readInt();
					inStrm.readInt();
					dataPos = fc.position()+offset;
					if(foundOne)
						break;
					foundOne = true;
				}else if ("FVER".equals(ckID)){
					vTimestamp = inStrm.readInt();
					if(vTimestamp != 0xA2805140l)
					{
						System.err.println("AIFF format version not recognized");
						return false;
					}
				}
			}
			if (srcBPS != 8 && srcBPS != 16 && srcBPS != 24 && srcBPS != 32){
				System.err.println(String.format("Error : Only 8bit, 16bit, 24bit and 32bit PCM are supported."));
				return false;
			}
			srcBPS /= 8;
			length = numSampleFrames * srcChannels * srcBPS;
			if(length == 0){
				System.err.println(String.format("Couldn't find data chank"));
				return false;
			}
			fc.position(dataPos);
			
		} catch(EOFException eof){
			System.err.println(String.format("Couldn't find data chank"));
			return false;
		} catch (NumberFormatException e1) {
			System.err.println("Error reading file header");
			return false;
		} catch (IOException e1) {
			System.err.println("Error reading file header");
			return false;
		}
		return true;		
	}

	private static boolean readWavHeader(DataInputStream inStrm, FileChannel fc){
		byte[] buf = new byte[16];
		int size;
		ByteBuffer bb = ByteBuffer.wrap(buf).order(ByteOrder.LITTLE_ENDIAN);
		bb.clear();
		try {
			inStrm.readInt();

			inStrm.readFully(buf,0,8);
			if(!"WAVEfmt ".equals(new String(buf, 0,8)))
				return false;

			inStrm.readFully(buf,0,4);
			size = bb.getInt();
			if(size > buf.length){
				buf = new byte[size];
				bb = ByteBuffer.wrap(buf).order(ByteOrder.LITTLE_ENDIAN);
			}
			inStrm.readFully(buf,0,size);
			bb.clear();
			
			if (bb.getShort() != 1) {
				System.err.println(String.format("Error: Only PCM is supported."));
				return false;
			}
			srcChannels = bb.getShort();
			srcSamplingRate = bb.getInt();
			srcBPS = bb.getInt();
			if ((int)srcBPS % srcSamplingRate*srcChannels != 0) return false;

			srcBPS /= srcSamplingRate*srcChannels;

			byte[] cbuf = new byte[4];
			for(;;)
			{
				inStrm.readFully(cbuf,0,4);
				try{
					bb.clear();
					inStrm.readFully(buf,0,4);
					length = (long)bb.getInt();
				}catch(IOException ex){
					break;
				}
				if("data".equals(new String(cbuf, 0,4))) break;
				fc.position(fc.position() + length);
			}
			if (fc.position() > fc.size()) {
				System.err.println(String.format("Couldn't find data chank"));
				return false;
			}

			if (srcBPS != 1 && srcBPS != 2 && srcBPS != 3 && srcBPS != 4) {
				System.err.println(String.format("Error : Only 8bit, 16bit, 24bit and 32bit PCM are supported."));
				return false;
			}

		} catch (IOException e1) {
			System.err.println("Error reading file header");
			return false;
		}
		return true;		
	}
	
	public static int readFileHeader(FileInputStream fpi){
		DataInputStream inStrm = new DataInputStream(fpi);

		byte[] buf = new byte[4];
		try {
			inStrm.readFully(buf,0,4);
			String fType = new String(buf, 0,4);
			if("RIFF".equals(fType)){
				byteOrder = ByteOrder.LITTLE_ENDIAN;
				if(readWavHeader(inStrm, fpi.getChannel()))
					return 1;
			}else if("FORM".equals(fType)){
				byteOrder = ByteOrder.BIG_ENDIAN;
				if(readAiffHeader(inStrm, fpi.getChannel()))
					return 2;
			}else if("fLaC".equals(fType)){
/*				fpi.getChannel().position(0);
				FLACDecoder flacDec = new FLACDecoder(inStrm);
				StreamInfo strmInfo = flacDec.readStreamInfo();
				srcSamplingRate = strmInfo.getSampleRate();
				srcBPS = strmInfo.getBitsPerSample()/8;
				srcChannels = strmInfo.getChannels();
				length = strmInfo.getTotalSamples() * srcBPS * srcChannels;
				fpi.getChannel().position(0);
				return 3;*/
			}
		} catch (IOException e1) {
			System.err.println("Error reading file header");
		}
		System.err.println(String.format("Error: Unsupported file type."));
		return 0;
	}

	public static void writeWavHeader(FileOutputStream fpo) throws IOException{
		ByteBuffer bb = ByteBuffer.wrap(new byte[44]).order(ByteOrder.LITTLE_ENDIAN);
		bb.put("RIFF".getBytes());
		bb.putInt(0);

		bb.put("WAVEfmt ".getBytes());
		bb.putInt(16);
		bb.putShort((short)1); // PCM
		bb.putShort((short) dstChannels); // Channels
		bb.putInt(dstSamplingRate); // Sampling rate 
		bb.putInt(dstSamplingRate*dstChannels*dstBPS); // Bytes per sec
		bb.putShort((short) (dstBPS*dstChannels));// Block alignment
		bb.putShort((short) (dstBPS*8)); // Bits per sample

		bb.put("data".getBytes());
		bb.putInt(0);
		fpo.write(bb.array(), 0, bb.position());
	}

	public static boolean writeDataLengthsToHeader(FileOutputStream fpo){
		long len;
		ByteBuffer bb = ByteBuffer.allocate(4).order(ByteOrder.LITTLE_ENDIAN);
		FileChannel fc = fpo.getChannel();
		try {
			len = fc.size();

			fc.position(4);
			bb.putInt((int) (len-8));
			bb.flip();
			fc.write(bb);

			fc.position(40);
			bb.clear();
			bb.putInt((int) (len-44));
			bb.flip();
			fc.write(bb);
		} catch (IOException e) {
			System.err.println("error writing data lengths.");
			e.printStackTrace(System.err);
			return false;
		}
		return true;
	}
}
