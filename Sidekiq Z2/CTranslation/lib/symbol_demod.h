#ifndef SYMBOL_DEMOD_H
#define SYMBOL_DEMOD_H

#include "src/symbol_demod_QAM.c"

Array_Tuple symbol_demod(Complex_Array_Tuple baseband_symbols_complex, char scheme[], double channel_gain, double preamble_len);

#endif