#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdbool.h>

#define MAX_ARRAY_LENGTH 1024

#include "standardArray.h"
#include "CRC.h"
#include "iq_imbalance.h"

void stuff();

struct Array_Tuple pulsetrain(struct Array_Tuple bits, int sps);
struct Complex_Array_Tuple generateNoise(struct Array_Tuple testpacket);
double rand_norm(double mu, double sigma);
double normal_distribution(double x, double mean, double stddev);
int main(){
    stuff();
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
    //freeArrayMemory(out.array); // when done with the array, need to free its memory
    //struct Array_Tuple testpacket = pulse_shape(pulse_train, sps, fs, pulse_shape, alpha, L);

    // Noise simulation ----------------------------
    struct Array_Tuple testpacket_noise = generateNoise(pulse_train);
    printArray("T", testpacket_noise.array, testpacket_noise.length);
    // ---------------------------------------------
    return;
}

struct Array_Tuple pulsetrain(struct Array_Tuple bits, int sps){
    double* array_ptr;
    int length = bits.length * sps;
    array_ptr = (double*)calloc(length, sizeof(double));
    int offset = 0;
    for (int i = 0; i < bits.length; i++)
    {
        for (int x = 0; x < sps; x++)
        {
            if (x == 0){
                array_ptr[x+offset] = bits.array[i]*2 - 1; // sets 0 to -1 and 1 to 1
            } else {
                array_ptr[x+offset] = 0; // puts zeros between the bits to make it easier to read
            }
        }
        offset += sps;
    }
    struct Array_Tuple t = {array_ptr, length};
    return t;
    
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
        double phase_noise = exp(rand_norm(mean, std_dev) * phase_noise_strength);
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
double rand_norm(double mu, double sigma){
    double U1, U2, W, mult;
    static double rand_norm_X1, rand_norm_X2;
    static int call_norm = 0;

    if (call_norm == 1)
    {
        call_norm = !call_norm;
        return (mu + sigma * (double)rand_norm_X2);
    }

    do
    {
        U1 = -1 + ((double)rand() / RAND_MAX) * 2;
        U2 = -1 + ((double)rand() / RAND_MAX) * 2;
        W = pow(U1, 2) + pow(U2, 2);
    } while (W >= 1 || W == 0);

    mult = sqrt((-2 * log(W)) / W);
    rand_norm_X1 = U1 * mult;
    rand_norm_X2 = U2 * mult;

    call_norm = !call_norm;

    return (mu + sigma * (double)rand_norm_X1);
}