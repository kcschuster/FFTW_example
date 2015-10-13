/*
 *  Title: fftw_analysis.cpp
 *
 *  Description:  A simple demonstration of using FFTW to calculate
 *      two-dimensional Fourier transforms in C++.  Can be modified
 *      based on specific application.
 *
 *  Note:  Compile with "g++ fftw_analysis.cpp -lfftw3"
 *
 *  Author: Kelsey Schuster
 */


#include <iostream>
#include <cstdlib>
#include <stdio.h>
#include <math.h>
#include <sstream>
#include <fstream>
#include <map>
#include <string>
#include </usr/local/include/fftw3.h>


using namespace std;

//functions
void set_in(fftw_complex *, int, int, string, int);
void print_output(map<int, double>&, int, int, int);
void average_out(fftw_complex *, int, int, map<int, double>&);
void add_prefactor(map<int, double>&);
int get_index(int, int, int);

//data structures
map<int, double> out_avg;




//main function
int main()
{
    
    //=========== Input ============================================
    
    string inputFile = "example/interface";
    
    int N1 = 43;   //number of values input to fftw in x-dim
    int N2 = 43;   //number of values input to fftw in y-dim
    
    int nframes = 4;    //number of frames over which to average
    
    //==============================================================
    
    
    // define input and output matrices, make plan for 2d transform
    fftw_complex *in, *out;
    in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N1*N2);
    out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N1*N2);
    
    fftw_plan p;
    p = fftw_plan_dft_2d(N1, N2, in, out, FFTW_FORWARD, FFTW_MEASURE);
    
    
    // loop through all frames to get average
    for (int i=0; i<nframes; i++) {
        
        cout << "frame: " << i << endl;
        
        // fill input matrix with data from system
        set_in(in, N1, N2, inputFile, i);
    
        // do fourier transform of data
        fftw_execute(p);
        
        // add data to running average
        average_out(out, N1, N2, out_avg);
    }
    
    //if necessary, add prefactor
    add_prefactor(out_avg);
    
    // print output data
    print_output(out_avg, N1, N2, nframes);

    // free up space
    fftw_destroy_plan(p);
    fftw_free(in);
    fftw_free(out);
}


// adds prefactor to average and std dev
void add_prefactor(map<int, double>& avg)
{
    /*
     *  Add your own function here if you want to add prefactor to FFTW output.
     *  Dependent on your specific application.
     *
     */
}


// averages output from all frames (read in row-major order)
void average_out(fftw_complex *out, int N1, int N2, map<int, double>& avg)
{
    int index;
    double value;
    
    //add new data to cumulative average
    for (int x=0; x<N1; x++) {
        for (int y=0; y<N2; y++) {
            
            //get index and magnitude (real^2 + complex^2)
            index = get_index(x, y, N2);
            value = pow(out[index][0], 2.0) + pow(out[index][1], 2.0);
            
            avg[index] += value;
        }
    }
}


// prints output from FFTW
void print_output(map<int, double>& out_avg, int N1, int N2, int nframes) {
    
    ostringstream fn;
    fn << "FFTW_output.txt";
    ofstream file(fn.str().c_str());
    
    //loop through FFTW output
    for (int x=0; x<N1; x++) {
        for (int y=0; y<N2; y++) {
            
            //get row-major index
            int index = get_index(x, y, N2);
            
            //average over frames
            double value = out_avg[index]/((double)nframes);
        
            //print data to file
            file << x << "\t" << y << "\t" << value << "\n";
        }
    }
}


//returns appropriate row-major index for (x,y) coords
int get_index(int x, int y, int N2)
{
    return (y + N2*x);
}


// sets input for fftw
void set_in(fftw_complex *in, int N1, int N2, string file, int frame) {
    int x,y,h;
    
    //open input file (follows naming convention used in "example" folder)
    ifstream indata;
    ostringstream fn;
    fn << file << frame << ".txt";
    indata.open(fn.str().c_str());
    
    //exit program if file doens't exist
    if (!indata) {
        cerr << "Error: input file '" << fn.str() << "' could not be found" << endl;
        exit(1);
    }
    
    //read data from file and put into matrix (row-major order)
    for (;;) {
        indata >> x >> y >> h;
        if (indata.fail()) {
            break;
        }
        int index = get_index(x, y, N2);
        in[index][0] = h;   //real part
        in[index][1] = 0;   //imag part (if using real coords, just eq to 0)
    }
    
    /*
    //For testing: fill input matrix with sine wave
    int k = 0;
    for (int i=0; i<N1; i++) {
        for (int j=0; j<N2; j++) {
            in[k][0] = cos(2*M_PI*i/N1)*cos(2*M_PI*j/N1) - sin(2*M_PI*i/N1)*sin(2*M_PI*j/N1);
            in[k][1] = 0;
            k++;
        }
    }
    */
}



