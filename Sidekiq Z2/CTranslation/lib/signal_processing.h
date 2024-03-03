#ifndef SIGNAL_PROCESS_H
#define SIGNAL_PROCESS_H

#include "src/signal_processing.c"

// Using the fftw3 library, calculates discrete fft on an array np.fft.fft()
Complex_Array_Tuple fft(Complex_Array_Tuple array);
//np.hamming
Array_Tuple hamming(int M);
// Function to perform FFT shift along a single dimension
Array_Tuple fftshift(Array_Tuple data);
//fftshift for complex array
Complex_Array_Tuple complexfftshift(Complex_Array_Tuple input);
// Takes a bit array and converts it to a pulse train from -1 to 1 with sps-1 zeros between each bit
Array_Tuple pulsetrain(Array_Tuple bits, int sps);
Array_Tuple firwin(int M, double cutoff);
// same as scipy.resample_poly which upsamples an array of data, can also specify a down sample
Complex_Array_Tuple resample_poly(Complex_Array_Tuple a, int up, int down);
Complex_Array_Tuple generateComplexNoise(Complex_Array_Tuple testpacket, double std_dev, double phase_noise_strength, double noise_power);
#endif