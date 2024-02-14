/*
main.c
IrishSAT CLOVER SDR
February 2024
Rylan Paul

Notes:
- For memory leaks, every function that returns a pointer using malloc() or calloc(), mark with //! Returns a calloc ptr
            - For those using vscode, install the better comments extension, then you will see comments starting with --> ! <-- in bright red


Funtcions to add:
- np.sinc
- np.hamming
- np.sum
- np.convolve
- signal.resample_poly
- signal.upfirdn
- np.mean
- signal.fftconvolve
- demod.symbol_demod
- CRC.CRCcheck (double check this works)
- pulse_shaping.pulse_shaping

*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdbool.h>

// Our libraries --------------------------------
#include "standardArray.h" //#include <fftw3.h> //link flag -lfftw3, run: brew install fftw // to install library
#include "CRC.h"
#include "iq_imbalance.h"
//-----------------------------------------------

void encodepacket();

int main(){
    struct Array_Tuple real = defineArray((double[]){-0.1, 0.4, -0.5, -0.2, 0.1, -0.3, -6, -0.6, 0.1, -0.4}, 10);
    struct Array_Tuple imaj = defineArray((double[]){0, 0.6, -0.7, -0.3, 0.1, 0.9, 5.1, 1.1, -0.1, 0.2}, 10);
    struct Array_Tuple bits = defineArray((double[]){1,0,0,0,0,1,1,1,0,1,0,1,0}, 13);
    struct Complex_Array_Tuple complex = {real, imaj};
    struct Complex_Array_Tuple out = fftconvolve(complex, bits);
    printComplexArray("Convolve", out.real.array, out.imaginary.array, out.real.length);
    //printComplexArray("fft", out.real.array, out.imaginary.array, out.real.length);
    //encodepacket();
    return 0;
}

void encodepacket(){
    // Declare Variables
    const int data_length = 10;
    double fs = 2.4e9; // carrier frequency 
    double Ts = 1.0 / fs; // sampling period in seconds
    double f0 = 0.0; // homodyne (0 HZ IF)
    int M = 8; // oversampling factor
    double L = 8; // pulse shaping filter length (Is this the ideal length?)
    char pulse_shape[] = "rrc"; // type of pulse shaping, also "rect"
    char scheme[] = "BPSK"; // Modulation scheme 'OOK', 'QPSK', or 'QAM'
    double alpha = 0.5; // roll-off factor of the RRC pulse-shaping filter

    struct Array_Tuple preamble = defineArray((double[]){0,1,0,0,0,0,1,1,0,0,0,1,0,1,0,0,1,1,1,1,0,1,0,0,0,1,1,1,0,0,1,0,0,1,0,1,1,0,1,1,1,0,1,1,0,0,1,1,0,1,0,1,0,1,1,1,1,1,1,0}, 60);// optimal periodic binary code for N = 63 https://ntrs.nasa.gov/citations/19800017860
    struct Array_Tuple data = {randomArray(2,data_length), data_length}; // this points to the random array generated;  creates 256 random bits

    struct Array_Tuple CRC_key = defineArray((double[]){1,0,0,1,1,0,0,0,0,1,1,1}, 12);// Best CRC polynomials: https://users.ece.cmu.edu/~koopman/crc/
    struct Array_Tuple data_encoded = CRC_encodeData(data, CRC_key);
    freeArrayMemory(CRC_key); // if its not needed annymore
    struct Array_Tuple bits = append_array(preamble, data_encoded);

    int sps = M; // samples per symbol, equal to oversampling factor

    struct Array_Tuple matched_filter_coef = flip(preamble);

    struct Array_Tuple pulse_train = pulsetrain(bits, sps);
    
    //struct Array_Tuple testpacket = pulse_shape(pulse_train, sps, fs, pulse_shape, alpha, L);

    // Noise simulation ----------------------------
    //struct Array_Tuple testpacket_noise = generateNoise(testpacket);
    //printArray("T", testpacket_noise.array, testpacket_noise.length);
    // ---------------------------------------------
    return;
}
