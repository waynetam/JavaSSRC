
/**
 * *****************************************************************************
 * Copyright (c) 2013 Wayne Tam. All rights reserved. This program (except for
 * ResampleHelper) are made available under the terms of the GNU Lesser Public
 * License v2.1 which accompanies this distribution, and is available at
 * http://www.gnu.org/licenses/old-licenses/gpl-2.0.html
 *
 * Please see the source of ResampleHelper for its license.
 *
 * Contributors: Wayne Tam - initial API and implementation
 *****************************************************************************
 */
import java.io.BufferedOutputStream;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;

import com.angrygoat.audio.resample.JavaSSRC;

import java.io.FileNotFoundException;

public class Resample implements JavaSSRC.ProgressListener {

    private double[] presets = new double[]{0.7, 0.9, 0.18};
    private long starttime;
    private long lastshowed;
    private int lastshowed2;

    public BufferedOutputStream fOut;
    public int readLen, totalWritten = 0;
    public static double dLength;

    private void usage() {
        System.out.println("http://shibatch.sourceforge.net/\n");
        System.out.println("usage: java -jar JavaSSRC [<options>] <source wav file> <destination wav file>");
        System.out.println("usage: java Resample [<options>] <source wav file> <destination wav file>");
        System.out.println("options : --rate <sampling rate>     output sample rate");
        System.out.println("          --attentuation <attenuation(dB)>    attenuate signal");
        System.out.println("          --bits <number of bits>    output quantization bit length");
        System.out.println("          --channels <number of channels>    number output audio channels");
        System.out.println("          --monoChannel <channel>    use a single src channel");
        System.out.println("          --tmpfile <file name>      specify temporal file");
        System.out.println("          --twoPass                  two pass processing to avoid clipping");
        System.out.println("          --normalize                normalize the wave file");
        System.out.println("          --quiet                    nothing displayed except error");
        System.out.println("          --float                    floating point pcm output");
        System.out.println("          --dither [<type>]          dithering");
        System.out.println("                                       0 : no dither");
        System.out.println("                                       1 : no noise shaping");
        System.out.println("                                       2 : triangular spectral shape");
        System.out.println("                                       3 : ATH based noise shaping");
        System.out.println("                                       4 : less dither amplitude than type 3");
        System.out.println("          --pdf <type> [<amp>]       select p.d.f. of noise");
        System.out.println("                                       0 : rectangular");
        System.out.println("                                       1 : triangular");
        System.out.println("                                       2 : Gaussian");
    }

    private boolean parseArguments(String[] argv) {
        int i;

		// parse command line options
        ResamplerHelper.monoChannel = -1;
        ResamplerHelper.dstChannels = 2;
        ResamplerHelper.dstSamplingRate = -1;
        ResamplerHelper.attentuation = 0;
        ResamplerHelper.dstBPS = -1;
        ResamplerHelper.twoPass = false;
        ResamplerHelper.normalize = false;
        ResamplerHelper.isFast = false;
        ResamplerHelper.dither = 0;
        ResamplerHelper.pdf = 0;
        ResamplerHelper.noiseAmp = 0.18;

        for (i = 0; i < argv.length; i++) {
            if (argv[i].charAt(0) != '-') {
                break;
            }

            if ("--rate".equals(argv[i])) {
                ResamplerHelper.dstSamplingRate = Integer.parseInt(argv[++i]);
                continue;
            }

            if ("--attentuation".equals(argv[i])) {
                ResamplerHelper.attentuation = Double.parseDouble(argv[++i]);
                continue;
            }

            if ("--bits".equals(argv[i])) {
                ResamplerHelper.dstBPS = Integer.parseInt(argv[++i]);
                if (ResamplerHelper.dstBPS != 8 && ResamplerHelper.dstBPS != 16 && ResamplerHelper.dstBPS != 24) {
                    System.err.println("Error: Only 8bit, 16bit and 24bit PCM are supported.");
                    return false;
                }
                ResamplerHelper.dstBPS /= 8;
                continue;
            }

            if ("--monoChannel".equals(argv[i])) {
                ResamplerHelper.monoChannel = Integer.parseInt(argv[++i]);
                continue;
            }

            if ("--channels".equals(argv[i])) {
                ResamplerHelper.dstChannels = Integer.parseInt(argv[++i]);
                continue;
            }

            if ("--twoPass".equals(argv[i])) {
                ResamplerHelper.twoPass = true;
                continue;
            }

            if ("--float".equals(argv[i])) {
                ResamplerHelper.dstFloat = true;
                continue;
            }

            if ("--normalize".equals(argv[i])) {
                ResamplerHelper.twoPass = true;
                ResamplerHelper.normalize = true;
                continue;
            }

            if ("--dither".equals(argv[i])) {
                try {
                    ResamplerHelper.dither = Integer.parseInt(argv[i + 1], 10);
                    if (ResamplerHelper.dither < 0 || ResamplerHelper.dither > 4) {
                        System.err.println(String.format("unrecognized dither type : %s", argv[i + 1]));
                        return false;
                    }
                    i++;
                } catch (NumberFormatException e) {
                    ResamplerHelper.dither = -1;
                }
                continue;
            }

            if ("--pdf".equals(argv[i])) {
                try {
                    ResamplerHelper.pdf = Integer.parseInt(argv[i + 1], 10);
                    if (ResamplerHelper.pdf < 0 || ResamplerHelper.pdf > 2) {
                        System.err.println(String.format("unrecognized p.d.f. type : %s", argv[i + 1]));
                        return false;
                    }
                    i++;
                } catch (NumberFormatException e) {
                    System.err.println(String.format("unrecognized p.d.f. type : %s", argv[i + 1]));
                    return false;
                }

                try {
                    ResamplerHelper.noiseAmp = Double.parseDouble(argv[i + 1]);
                    i++;
                } catch (NumberFormatException e) {
                    ResamplerHelper.noiseAmp = presets[ResamplerHelper.pdf];
                }

                continue;
            }

            if ("--quiet".equals(argv[i])) {
                ResamplerHelper.quiet = true;
                continue;
            }

            if ("--tmpfile".equals(argv[i])) {
                ResamplerHelper.tmpFile = argv[++i];
                continue;
            }

            if ("--isFast".equals(argv[i])) {
                if ("isFast".equals(argv[i + 1])) {
                    ResamplerHelper.isFast = true;
                } else //noinspection StatementWithEmptyBody
                    if ("standard".equals(argv[i + 1])) {
                    /* nothing to do */
                } else {
                    System.err.println(String.format("unrecognized isFast : %s", argv[i + 1]));
                    return false;
                }
                i++;
                continue;
            }

            System.err.println(String.format("unrecognized option : %s\n", argv[i]));
            return false;
        }

        if (argv.length - i != 2) {
            return false;
        }

        ResamplerHelper.srcFile = argv[i];
        ResamplerHelper.destFile = argv[i + 1];

        return true;
    }

    private void showInfo() {
        String[] dtype = new String[]{
            "none", "no noise shaping", "triangular spectral shape", "ATH based noise shaping", "ATH based noise shaping(less amplitude)"
        };
        String[] ptype = new String[]{
            "rectangular", "triangular", "gaussian"
        };
        System.out.println(String.format("frequency : %d -> %d", ResamplerHelper.srcSamplingRate, ResamplerHelper.dstSamplingRate));
        System.out.println(String.format("attenuation : %gdB", ResamplerHelper.attentuation));
        if (ResamplerHelper.srcFloat && ResamplerHelper.srcFloat == ResamplerHelper.dstFloat) {
            System.out.println("bits per sample : float -> float");
        } else {
            if (ResamplerHelper.srcFloat) {
                System.out.println(String.format("bits per sample : float -> %d", ResamplerHelper.dstBPS * 8));
            } else if (ResamplerHelper.dstFloat) {
                System.out.println(String.format("bits per sample : %d -> float", ResamplerHelper.srcBPS * 8));
            } else {
                System.out.println(String.format("bits per sample : %d -> %d", ResamplerHelper.srcBPS * 8, ResamplerHelper.dstBPS * 8));
            }
        }
        System.out.println(String.format("channels : %d -> %d", ResamplerHelper.srcChannels, ResamplerHelper.dstChannels));
        if (ResamplerHelper.monoChannel > -1) {
            System.out.println(String.format("mono channel : %d", ResamplerHelper.monoChannel + 1));
        }
        System.out.println(String.format("length : %d bytes, %g secs", ResamplerHelper.length,
                (double) ResamplerHelper.length / ResamplerHelper.srcBPS / ResamplerHelper.srcChannels / ResamplerHelper.srcSamplingRate));
        if (ResamplerHelper.dither == 0) {
            System.out.println("dither type : none");
        } else {
            System.out.println(String.format("dither type : %s, %s p.d.f, amp = %g",
                    dtype[ResamplerHelper.dither], ptype[ResamplerHelper.pdf], ResamplerHelper.noiseAmp));
        }
        System.out.println();
    }

    private void setStartTime() {
        starttime = System.currentTimeMillis();
        lastshowed = 0;
        lastshowed2 = -1;
    }

    @Override
    public void onChanged(double progress) {
        if (!ResamplerHelper.quiet) {
            int eta, pc;
            long t;

            t = System.currentTimeMillis() - starttime;
            if (progress == 0) {
                eta = 0;
            } else {
                eta = (int) (t * (1 - progress) / progress);
            }

            pc = (int) (progress * 100);
            if (pc != lastshowed2 || t != lastshowed) {
                System.out.print(String.format(" %3d%% processed", pc));
                lastshowed2 = pc;
            }
            if (t != lastshowed) {
                System.out.print(String.format(", ETA = %7.3fsec  ", (float) eta / 1000));
                lastshowed = t;
            }
            System.out.print("\r");
            System.out.flush();

            if (progress == 1) {
                System.out.println(String.format("Total Time %7.3fsec                 ", (float) t / 1000));
                setStartTime();
            }
        }
    }

    @Override
    public void onShowMessage(String message) {
        System.out.println(message);
    }

    public static void main(String[] argv) {
        new Resample(argv);
    }

    public Resample(String[] argv) {
        FileInputStream fpi;
        FileOutputStream fpo = null;
        double peak;

        if (!parseArguments(argv)) {
            usage();
            return;
        }

        if (!ResamplerHelper.quiet) {
            System.err.println(String.format("Shibatch sampling rate converter version %s\n", JavaSSRC.VERSION));
        }

        try {
            fpi = new FileInputStream(ResamplerHelper.srcFile);
        } catch (FileNotFoundException e) {
            System.err.println("cannot open input file.");
            e.printStackTrace(System.err);
            return;
        }

        try {
            int fType = ResamplerHelper.readFileHeader(fpi);
            if (fType == 0) {
                return;
            }

            if (ResamplerHelper.dstBPS == -1) {
                if (ResamplerHelper.srcBPS != 1) {
                    ResamplerHelper.dstBPS = ResamplerHelper.srcBPS;
                } else {
                    ResamplerHelper.dstBPS = 2;
                }
                if (ResamplerHelper.dstBPS == 4) {
                    ResamplerHelper.dstBPS = 3;
                }
            }

            if (ResamplerHelper.dstSamplingRate == -1) {
                ResamplerHelper.dstSamplingRate = ResamplerHelper.srcSamplingRate;
            }

            if (ResamplerHelper.dither == -1) {
                if (ResamplerHelper.dstBPS < ResamplerHelper.srcBPS) {
                    if (ResamplerHelper.dstBPS == 1) {
                        ResamplerHelper.dither = 4;
                    } else {
                        ResamplerHelper.dither = 3;
                    }
                } else {
                    ResamplerHelper.dither = 1;
                }
            }

            if (!ResamplerHelper.quiet) {
                showInfo();
            }

            try {
                fpo = new FileOutputStream(ResamplerHelper.destFile);
                ResamplerHelper.writeWavHeader(fpo, ResamplerHelper.dstFloat);
            } catch (IOException e) {
                System.err.println("error writing output header.");
                e.printStackTrace(System.err);
                return;
            }

            setStartTime();

            try {
                /*				if(fType == 3)
                 peak = ResamplerHelper.processFlac(fpi, fpo, this);
                 else //*/
                peak = ResamplerHelper.processPCM(fpi, fpo, this);
            } catch (IOException e1) {
                e1.printStackTrace(System.err);
                return;
            }

            if (!ResamplerHelper.quiet) {
                System.out.println();
            }

            if (!ResamplerHelper.twoPass && peak > 1) {
                if (!ResamplerHelper.quiet) {
                    System.out.println(String.format("clipping detected : %gdB\n", 20 * Math.log10(peak)));
                }
            }

            ResamplerHelper.writeDataLengthsToHeader(fpo);
        } finally {
            try {
                fpi.close();
            } catch (IOException e) {
                e.printStackTrace(System.err);
            }
            if (fpo != null) {
                try {
                    fpo.close();
                } catch (IOException e) {
                    e.printStackTrace(System.err);
                }
            }
        }
        System.exit(0);
    }
}
