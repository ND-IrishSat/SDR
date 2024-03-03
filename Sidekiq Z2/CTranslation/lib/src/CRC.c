// CRC.c
// Rylan Paul
#include "../standardArray.h"
// CRC_xor 
Array_Tuple CRC_xor(Array_Tuple a_array, Array_Tuple b_array){
    
    int length = b_array.length;
    Array_Tuple out_array = zerosArray(length-1);
    for (int i = 1; i < length; i++){
        double element1 = a_array.array[i];
        double element2 = b_array.array[i];
        if (element1 == element2){
            out_array.array[i-1] = 0;
        } else {
            out_array.array[i-1] = 1;
        }
    }
    return out_array;
}
// CRC_mod2div
Array_Tuple CRC_mod2div(Array_Tuple dividend_array, Array_Tuple divisor_array){
    int pick = divisor_array.length;
    Array_Tuple tmp = zerosArray(pick);
    for (int i=0; i < pick; i++){
        tmp.array[i] = dividend_array.array[i];
    }
    

    while (pick < dividend_array.length){
        if (tmp.array[0] == 1){
            Array_Tuple x = CRC_xor(divisor_array, tmp);
            Array_Tuple new_tmp = zerosArray(x.length + 1);
            int last_index;
            for (int i=0; i < x.length; i++){
                new_tmp.array[i] = x.array[i];
                last_index = i;
            }
            new_tmp.array[last_index + 1] = dividend_array.array[pick];
            for (int i=0; i < divisor_array.length; i++){
                tmp.array[i] = new_tmp.array[i];
            }
            freeArrayMemory(x);
            freeArrayMemory(new_tmp);
        } else {
            Array_Tuple ones_tuple = valueArray(pick, 1);
            Array_Tuple xor_out = CRC_xor(ones_tuple, tmp);
            Array_Tuple new_tmp_2 = zerosArray(xor_out.length + 1);
            int last_index;
            for (int i=0; i < xor_out.length; i++){
                new_tmp_2.array[i] = xor_out.array[i];
                last_index = i;
            }
            new_tmp_2.array[last_index + 1] = dividend_array.array[pick];
            for (int i=0; i < divisor_array.length; i++){
                tmp.array[i] = new_tmp_2.array[i];
            }
            freeArrayMemory(ones_tuple);
            freeArrayMemory(xor_out);
            freeArrayMemory(new_tmp_2);
        }

        pick += 1;
    }
    
    if (tmp.array[0] == 1){
        Array_Tuple xor_out_2 = CRC_xor(divisor_array, tmp);
        freeArrayMemory(tmp);
        return xor_out_2;
    } else {
        Array_Tuple zeros_tuple = zerosArray(pick);
        Array_Tuple xor_out_2 = CRC_xor(zeros_tuple, tmp);
        freeArrayMemory(zeros_tuple);
        freeArrayMemory(tmp);
        return xor_out_2;
    }
}
// CRC_encodedData2
Array_Tuple CRC_encodeData(Array_Tuple data, Array_Tuple key){
    Array_Tuple appended_data = zerosArray(data.length+key.length-1);
    for (int i=0; i<data.length; i++){
        appended_data.array[i] = data.array[i];
    }
    for (int i=data.length; i<data.length+key.length-1; i++){
        appended_data.array[i] = 0;
    }
    Array_Tuple mod2div_out = CRC_mod2div(appended_data, key);

    Array_Tuple codeword = zerosArray(data.length+mod2div_out.length);
    for (int i=0; i<data.length; i++){
        codeword.array[i] = data.array[i];
    }
    for (int i=data.length; i<data.length+mod2div_out.length; i++){
        codeword.array[i] = mod2div_out.array[i - data.length];
    }
    freeArrayMemory(appended_data);
    freeArrayMemory(mod2div_out);
    return codeword;
}

int CRC_check(Array_Tuple codeword, Array_Tuple key){
    Array_Tuple checkword = CRC_mod2div(codeword, key);
    int error = 0;
    for (int i = 0; i < checkword.length; i++) {
        if (checkword.array[i] == 0){
            error = 1;
            break;
        }
    }
    freeArrayMemory(checkword);
    return error;
}

