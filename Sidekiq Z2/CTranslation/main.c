/*
main.c
IrishSAT CLOVER SDR
February 2024
Rylan Paul

Notes:
- For memory leaks, every function that returns a pointer using malloc() or calloc(), mark with //! Returns a calloc ptr
            - For those using vscode, install the better comments extension, then you will see comments starting with --> ! <-- in bright red
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdbool.h>

// Our libraries --------------------------------
#include "standardArray.h"
#include "CRC.h"
#include "iq_imbalance.h"
//-----------------------------------------------

void encodepacket();

struct Complex_Array_Tuple generateNoise(struct Complex_Array_Tuple testpacket);

int main(){
    encodepacket();
    return 0;
}

void stuff(){
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
    
    struct Array_Tuple testpacket = pulse_shape(pulse_train, sps, fs, pulse_shape, alpha, L);

    // Noise simulation ----------------------------
    struct Array_Tuple testpacket_noise = generateNoise(testpacket);
    printArray("T", testpacket_noise.array, testpacket_noise.length);
    // ---------------------------------------------
    return;
}

struct Complex_Array_Tuple generateNoise(struct Complex_Array_Tuple testpacket){
    // vars
    int max = 1;
    double mean = 0;
    double std_dev = 0.1;
    double noise_power = 0.2;
    int num_samples = testpacket.real.length;
    double phase_noise_strength = 0.1;

    // creating general normal distribution random arrays for real and complex
    double* awgn_comlpex_samples_real;
    double* awgn_comlpex_samples_imag;
    awgn_comlpex_samples_real = (double*)calloc(num_samples, sizeof(double));
    awgn_comlpex_samples_imag = (double*)calloc(num_samples, sizeof(double));
    for (int i = 0; i < num_samples; i++)
    {
        double phase_noise = rand_norm(mean, std_dev) * phase_noise_strength;
        awgn_comlpex_samples_real[i] = rand_norm(mean, std_dev) / sqrt(noise_power) * cos(phase_noise);
        awgn_comlpex_samples_imag[i] = rand_norm(mean, std_dev) / sqrt(noise_power) * sin(phase_noise);
    }

    struct Array_Tuple real = {awgn_comlpex_samples_real, num_samples};
    struct Array_Tuple imag = {awgn_comlpex_samples_imag, num_samples};
    struct Complex_Array_Tuple out = {real, imag};
    return out;
    //awgn_comlpex_samples_real = ;
    //awgn_comlpex_samples_imaginary = ;
}
