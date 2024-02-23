// CRC.c
// Rylan Paul

// CRC_xor 
#define MAX_ARRAY_LENGTH 1024
Array_Tuple CRC_xor(Array_Tuple a_array, Array_Tuple b_array){
    int length = b_array.length;
    static double out_array[MAX_ARRAY_LENGTH]; // Assuming maximum length of the array
    for (int i = 1; i < length; i++){
        double element1 = a_array.array[i];
        double element2 = b_array.array[i];
        if (element1 == element2){
            out_array[i-1] = 0;
        } else {
            out_array[i-1] = 1;
        }
    }
    Array_Tuple tuple = {out_array, length-1};
    return tuple;
}
// CRC_mod2div
Array_Tuple CRC_mod2div(Array_Tuple dividend_array, Array_Tuple divisor_array){
    int pick = divisor_array.length;
    static double temp[MAX_ARRAY_LENGTH]; // Assuming maximum length of the array
    
    for (int i=0; i < pick; i++){
        temp[i] = dividend_array.array[i];
    }
    Array_Tuple tmp = {temp, pick};

    while (pick < dividend_array.length){
        if (tmp.array[0] == 1){
            double new_tmp[MAX_ARRAY_LENGTH]; // Assuming maximum length of the array
            Array_Tuple x = CRC_xor(divisor_array, tmp);
            int last_index;
            for (int i=0; i < x.length; i++){
                new_tmp[i] = x.array[i];
                last_index = i;
            }
            new_tmp[last_index + 1] = dividend_array.array[pick];
            for (int i=0; i < divisor_array.length; i++){
                tmp.array[i] = new_tmp[i];
            }
        } else {
            double new_tmp_2[MAX_ARRAY_LENGTH]; // Assuming maximum length of the array
            double ones_array[MAX_ARRAY_LENGTH]; // Assuming maximum length of the array
            for(int i=0; i < pick;  i++){
                ones_array[i] = 1;
            }
            Array_Tuple ones_tuple = {ones_array, pick};
            Array_Tuple xor_out = CRC_xor(ones_tuple, tmp);
            int last_index;
            for (int i=0; i < xor_out.length; i++){
                new_tmp_2[i] = xor_out.array[i];
                last_index = i;
            }

            new_tmp_2[last_index + 1] = dividend_array.array[pick];
            for (int i=0; i < divisor_array.length; i++){
                tmp.array[i] = new_tmp_2[i];
            }
        }

        pick += 1;
    }
    
    if (tmp.array[0] == 1){
        Array_Tuple xor_out_2 = CRC_xor(divisor_array, tmp);
        Array_Tuple out_array = {xor_out_2.array, xor_out_2.length};
        return out_array;
    } else {
        double zeros_array[MAX_ARRAY_LENGTH]; // Assuming maximum length of the array
        for(int i=0; i < pick; i++){
            zeros_array[i] = 0;
        }
        Array_Tuple zeros_tuple = {zeros_array, pick};
        Array_Tuple xor_out_2 = CRC_xor(zeros_tuple, tmp);
        Array_Tuple out_array = {xor_out_2.array, xor_out_2.length};
        return out_array;
    }
}
// CRC_encodedData2
Array_Tuple CRC_encodeData(Array_Tuple data, Array_Tuple key){ //! Returns a calloc ptr
    double appended_data[data.length+key.length-1]; // Static allocation
    for (int i=0; i<data.length; i++){
        appended_data[i] = data.array[i];
    }
    for (int i=data.length; i<data.length+key.length-1; i++){
        appended_data[i] = 0;
    }
    Array_Tuple appended_data_tuple = {appended_data, data.length+key.length-1};
    Array_Tuple mod2div_out = CRC_mod2div(appended_data_tuple, key);
    double* codeword; // Static allocation
    codeword = (double*)calloc(data.length+mod2div_out.length, sizeof(double)); // need to free after
    for (int i=0; i<data.length; i++){
        codeword[i] = data.array[i];
    }
    for (int i=data.length; i<data.length+mod2div_out.length; i++){
        codeword[i] = mod2div_out.array[i - data.length];
    }
    Array_Tuple tuple = {codeword, data.length+mod2div_out.length};
    return tuple;
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
    return error;
}

