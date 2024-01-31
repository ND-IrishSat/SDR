#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define PI 3.14159265358979323846

struct array_baseband {
  int *baseband;
  int length;
}

// Function prototypes
void pulse_shaping(float *a, int M, float fs, char *pulse_shape, float alpha, int L, float *baseband);
void rrcosfilter(int N, float alpha, float Ts, float Fs, float *time_idx, float *h_rrc);

// save dynamic arrays into a static variable before freeing the array 
void pulse_shaping(float *a, int M, float fs, char *pulse_shape, float alpha, int L, float *baseband) {
    // Upsample by a factor of M
    // (Note: The convolution part has been simplified)

    if (strcmp(pulse_shape, "rrc") == 0) {
        // Root Raised-Cosine span
        int N = L * M;
        float T_symbol = 1.0 / (fs / M);

        // Generate RRC filter coefficients
        float *time_idx;
        float *h_rrc;
        rrcosfilter(N, alpha, T_symbol, fs, time_idx, h_rrc);

        // Convolve the input with the RRC filter
        // (Note: The convolution part has been simplified)
        for (int i = 0; i < N + (sizeof(a) / sizeof(float)) - 1; i++) {
          baseband[i] = 0.0;
          for (int j = 0; j < (sizeof(a) / sizeof(float)); j++) {
            if (i - j >= 0 && i - j < N) {
              baseband[i] += a[j] * h_rrc[i - j];
            }
          }
        }
      }
  
    if (strcmp(pulse_shape, "rect") == 0) {
        // Rectangular pulse
        Ts = 1 / fs;
      for (int i = 0; i < (sizeof(a) / sizeof(float)) * M; i++) {
          baseband[i] = a[i];
      }
    }
}

void rrcosfilter(int N, float alpha, float Ts, float Fs, float *time_idx, float *h_rrc) {
    // Generate a root raised cosine (RRC) filter (FIR) impulse response
    float T_delta = 1.0 / Fs;

    // Allocate memory for time indices
    time_idx = (float *)malloc(N * sizeof(float)); // DYNAMIC

    // Populate time indices
    for (int x = 0; x < N; x++) {
        time_idx[x] = (x - N / 2) * T_delta;
    }

    // Allocate memory for RRC filter coefficients
    h_rrc = (float *)malloc(N * sizeof(float)); // DYNAMIC

    // Populate RRC filter coefficients
    for (int x = 0; x < N; x++) {
        float t = (x - N / 2) * T_delta;
        if (t == 0.0) {
            h_rrc[x] = 1.0 - alpha + (4 * alpha / PI);
        } else if (alpha != 0 && t == Ts / (4 * alpha)) {
            h_rrc[x] = (alpha / sqrt(2)) * (((1 + 2 / PI) * (sin(PI / (4 * alpha)))) +
                                            ((1 - 2 / PI) * (cos(PI / (4 * alpha)))));
        } else if (alpha != 0 && t == -Ts / (4 * alpha)) {
            h_rrc[x] = (alpha / sqrt(2)) * (((1 + 2 / PI) * (sin(PI / (4 * alpha)))) +
                                            ((1 - 2 / PI) * (cos(PI / (4 * alpha)))));
        } else {
            h_rrc[x] = (sin(PI * t * (1 - alpha) / Ts) +
                        4 * alpha * (t / Ts) * cos(PI * t * (1 + alpha) / Ts)) /
                       (PI * t * (1 - (4 * alpha * t / Ts) * (4 * alpha * t / Ts)) / Ts);
        }
    }
    free(h_rrc); // not sure where to free this --> possibly in main() but we will put it here for now
    free(time_idx); //save to a static variable
  
}
