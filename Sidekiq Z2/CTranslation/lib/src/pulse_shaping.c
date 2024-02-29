//pulse_shaping.c

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../standardArray.h" //! remove on build

// Function prototypes
Array_Tuple pulse_shaping(Array_Tuple a, int M, double fs, char pulse_shape[], double alpha, int L);
Array_Tuple rrcosfilter(int N, double alpha, double Ts, double Fs);

// save dynamic arrays into a static variable before freeing the array 
Array_Tuple pulse_shaping(Array_Tuple a, int M, double fs, char pulse_shape[], double alpha, int L) {
    if (strcmp(pulse_shape, "rrc") == 0) {
        int N = L * M;
        // Root Raised-Cosine span
        double T_symbol = 1.0 / (fs / M);

        // Generate RRC filter coefficients
        Array_Tuple h_rrc = rrcosfilter(N, alpha, T_symbol, fs);

        // Convolve the input with the RRC filter
        Array_Tuple a_imaj = zerosArray(a.length);
        Complex_Array_Tuple a_complex = {a, a_imaj};
        Complex_Array_Tuple convolve_out = convolve(a_complex, h_rrc);
        freeArrayMemory(a_imaj);
        freeArrayMemory(h_rrc);
        freeArrayMemory(convolve_out.imaginary);
        return convolve_out.real;
        //free(h_rrc);
    } 
    // else if (strcmp(pulse_shape, "rect") == 0) {//! This is not correct, goes off of scipy.upfirdn
    //     // Rectangular pulse
    //   for (int i = 0; i < a.length * M; i++) {
    //       baseband[i] = a.array[i];
    //   }
    //   Array_Tuple tuple = zerosArray(a.length * M);
    //   for (int i = 0; i < tuple.length; i++)
    //   {
    //     tuple.array[i] = baseband[i];
    //   }
    //   //free(baseband);
    //   return tuple;
    // }
    return a;
}

Array_Tuple rrcosfilter(int N, double alpha, double Ts, double Fs) {
    // Generate a root raised cosine (RRC) filter (FIR) impulse response
    double T_delta = 1.0 / Fs;

    // Allocate memory for time indices
    Array_Tuple time_idx_arange = arange(0, N, 1);
    Array_Tuple time_idx__minus = subtractDoubleFromArray(time_idx_arange, (double)N /2.0);
    Array_Tuple time_idx = multiplyDoubleFromArray(time_idx__minus, T_delta);

    Array_Tuple h_rrc = zerosArray(N);
    // Populate time indices
    for (int x = 0; x < N; x++) {
        time_idx.array[x] = (x - N / 2) * T_delta;
    }
    freeArrayMemory(time_idx_arange);
    freeArrayMemory(time_idx__minus);
    freeArrayMemory(time_idx);

    // Populate RRC filter coefficients
    for (int x = 0; x < N; x++) {
        double t = (x - (double)N / 2.0) * T_delta;
        if (t == 0.0) {
            h_rrc.array[x] = 1.0 - alpha + (4.0 * alpha / M_PI);
        } else if (alpha != 0 && t == Ts / (4 * alpha)) {
            h_rrc.array[x] = (alpha / sqrt(2)) * (((1 + 2.0 / M_PI) * (sin(M_PI / (4 * alpha)))) + ((1 - 2.0 / M_PI) * (cos(M_PI / (4 * alpha)))));
        } else if (alpha != 0 && t == -Ts / (4 * alpha)) {
            h_rrc.array[x] = (alpha / sqrt(2)) * (((1 + 2.0 / M_PI) * (sin(M_PI / (4 * alpha)))) + ((1 - 2.0 / M_PI) * (cos(M_PI / (4 * alpha)))));
        } else {
            h_rrc.array[x] = (sin(M_PI * t * (1 - alpha) / Ts) + 4 * alpha * (t / Ts) * cos(M_PI * t * (1 + alpha) / Ts)) / (M_PI * t * (1 - (4 * alpha * t / Ts) * (4 * alpha * t / Ts)) / Ts);
        }
    }

    return h_rrc;
}
