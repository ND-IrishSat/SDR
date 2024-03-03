#ifndef CRC_H // if CRC is not defined, then define it
#define CRC_H
// CRC.h
// Rylan Paul
#include "src/CRC.c"
// CRC

// CRC_xor 
Array_Tuple CRC_xor(Array_Tuple a_array, Array_Tuple b_array);
// CRC_mod2div
Array_Tuple CRC_mod2div(Array_Tuple dividend_array, Array_Tuple divisor_array);
// CRC_encodedData2
Array_Tuple CRC_encodeData(Array_Tuple data, Array_Tuple key);

int CRC_check(Array_Tuple codeword, Array_Tuple key);

#endif /* CRC.h end */