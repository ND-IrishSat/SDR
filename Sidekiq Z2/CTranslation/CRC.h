#include "CRC.c"
// CRC

// CRC_xor 
struct Array_Tuple CRC_xor(struct Array_Tuple a_array, struct Array_Tuple b_array);
// CRC_mod2div
struct Array_Tuple CRC_mod2div(struct Array_Tuple dividend_array, struct Array_Tuple divisor_array);
// CRC_encodedData2
struct Array_Tuple CRC_encodeData(struct Array_Tuple data, struct Array_Tuple key);

int CRC_check(struct Array_Tuple codeword, struct Array_Tuple key);
