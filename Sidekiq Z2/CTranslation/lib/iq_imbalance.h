#ifndef IQ_IMBALANCE_H // if CRC is not defined, then define it
#define IQ_IMBALANCE_H

/*
iq_imbalance.h
Rylan Paul
*/
#include "src/iq_imbalance.c"

void mean(Array_Tuple array, int period);
Complex_Array_Tuple IQImbalanceCorrect(Complex_Array_Tuple packet, int mean_period);

#endif
