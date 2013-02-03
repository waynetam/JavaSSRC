JavaSSRC
========

Resampling library in pure Java based on the SSRC library

This is based on the SSRC - High Quality Audio Sampling Rate Converter
by Naoki Shibata (http://shibatch.sourceforge.net/)

The original SSRC converter was written in C and was designed to work mainly with files.

This library started as a straight Java port of the C code, replacing the file I/O with
InputStream/OutputStream, but then reworked to add support for using byte and int arrays.

Look in Resample.java and ResampleHelper.java on how to use the library.


