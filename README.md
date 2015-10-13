# FFTW_example

Simple demonstration in C++ of how to use FFTW to calculate Fourier transforms.  Can be easily modified for specific applications.

First, you'll need to install [FFTW](http://www.fftw.org/download.html).  If you have Macports, type

    port install fftw-3

To compile the example code:

    g++ fftw_analysis.cpp -lfftw3


I've included an **example** folder with sample input files and the expected output after running the calculation.  The input text files give the height of a 2d surface as a function of x and y coordinates (\<x_coord\>, \<y_coord\>, \<height_surface\>).

