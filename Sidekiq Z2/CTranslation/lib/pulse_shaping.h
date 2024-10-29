#ifndef PULSE_SHAPING_H
#define PULSE_SHAPING_H

#include "src/pulse_shaping.c"


// Function prototypes
Array_Tuple pulse_shaping(Array_Tuple a, int M, double fs, char pulse_shape[], double alpha, int L);
Array_Tuple rrcosfilter(int N, double alpha, double Ts, double Fs);

#endif