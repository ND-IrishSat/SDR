/*
iq_imbalance.h
Rylan Paul
*/
#include "iq_imbalance.c"

void mean(struct Array_Tuple array, int period);
struct Complex_Array_Tuple IQImbalanceCorrect(struct Complex_Array_Tuple packet, int mean_period);