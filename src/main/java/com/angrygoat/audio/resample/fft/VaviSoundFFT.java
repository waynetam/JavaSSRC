package com.angrygoat.audio.resample.fft;

import vavi.util.SplitRadixFft;

/**
 * Created by Wayne on 5/19/2015.
 * VaviSoundFFT
 */
public class VaviSoundFFT implements FFT {

    private int[] ip;
    private double[] w;
    private int n;
    private SplitRadixFft fft;

    @Override
    public void init(int n) {
        this.n = n;
        this.ip = new int[(int) (2 + Math.sqrt(n))];
        this.w = new double[n / 2];
        this.fft = new SplitRadixFft();
        reset();
    }

    @Override
    public void reset() {
        this.ip[0] = 0;
    }

    @Override
    public void realDFT(double[] a) {
        fft.rdft(this.n,1,a,this.ip,this.w);
    }

    @Override
    public void realInverseDFT(double[] a) {
        fft.rdft(this.n,-1,a,this.ip,this.w);
    }
}
