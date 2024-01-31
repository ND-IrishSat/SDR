#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include "standardArray.h"
#define MAX_ARRAY_LENGTH 1024


// CRC

// CRC_xor 
struct Array_Tuple CRC_xor(struct Array_Tuple a_array, struct Array_Tuple b_array){
    int length = b_array.length;
    static int out_array[MAX_ARRAY_LENGTH]; // Assuming maximum length of the array
    for (int i = 1; i < length; i++){
        int element1 = a_array.array[i];
        int element2 = b_array.array[i];
        if (element1 == element2){
            out_array[i-1] = 0;
        } else {
            out_array[i-1] = 1;
        }
    }
    struct Array_Tuple tuple = {out_array, length-1};
    return tuple;
}
// CRC_mod2div
struct Array_Tuple CRC_mod2div(struct Array_Tuple dividend_array, struct Array_Tuple divisor_array){
    int pick = divisor_array.length;
    static int temp[MAX_ARRAY_LENGTH]; // Assuming maximum length of the array
    
    for (int i=0; i < pick; i++){
        temp[i] = dividend_array.array[i];
    }
    struct Array_Tuple tmp = {temp, pick};

    while (pick < dividend_array.length){
        if (tmp.array[0] == 1){
            int new_tmp[MAX_ARRAY_LENGTH]; // Assuming maximum length of the array
            struct Array_Tuple x = CRC_xor(divisor_array, tmp);
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
            int new_tmp_2[MAX_ARRAY_LENGTH]; // Assuming maximum length of the array
            int ones_array[MAX_ARRAY_LENGTH]; // Assuming maximum length of the array
            for(int i=0; i < pick;  i++){
                ones_array[i] = 1;
            }
            struct Array_Tuple ones_tuple = {ones_array, pick};
            struct Array_Tuple xor_out = CRC_xor(ones_tuple, tmp);
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
        struct Array_Tuple xor_out_2 = CRC_xor(divisor_array, tmp);
        struct Array_Tuple out_array = {xor_out_2.array, xor_out_2.length};
        return out_array;
    } else {
        int zeros_array[MAX_ARRAY_LENGTH]; // Assuming maximum length of the array
        for(int i=0; i < pick; i++){
            zeros_array[i] = 0;
        }
        struct Array_Tuple zeros_tuple = {zeros_array, pick};
        struct Array_Tuple xor_out_2 = CRC_xor(zeros_tuple, tmp);
        struct Array_Tuple out_array = {xor_out_2.array, xor_out_2.length};
        return out_array;
    }
}
// CRC_encodedData2
struct Array_Tuple CRC_encodeData(struct Array_Tuple data, struct Array_Tuple key){
    int appended_data[data.length+key.length-1]; // Static allocation
    for (int i=0; i<data.length; i++){
        appended_data[i] = data.array[i];
    }
    for (int i=data.length; i<data.length+key.length-1; i++){
        appended_data[i] = 0;
    }
    struct Array_Tuple appended_data_tuple = {appended_data, data.length+key.length-1};
    struct Array_Tuple mod2div_out = CRC_mod2div(appended_data_tuple, key);
    int* codeword; // Static allocation
    codeword = (int*)calloc(data.length+mod2div_out.length, sizeof(int)); // need to free after
    for (int i=0; i<data.length; i++){
        codeword[i] = data.array[i];
    }
    for (int i=data.length; i<data.length+mod2div_out.length; i++){
        codeword[i] = mod2div_out.array[i - data.length];
    }
    struct Array_Tuple tuple = {codeword, data.length+mod2div_out.length};
    return tuple;
}

int CRC_check(struct Array_Tuple codeword, struct Array_Tuple key){
    struct Array_Tuple checkword = CRC_mod2div(codeword, key);
    int error = 0;
    for (int i = 0; i < checkword.length; i++) {
        if (checkword.array[i] == 0){
            error = 1;
            break;
        }
    }
    return error;
}

int * randomArray(int max_exclusive, int length){ // returns a pointer to an integer
    static int arr[MAX_ARRAY_LENGTH]; // Assuming maximum length of the array
    for (int i=0; i < length; i++){
        arr[i] = rand() % max_exclusive;
    }
    return arr; // to use this: int *array; array= randomArray(2,size);
}

int main(){
    // Declare Variables
    const int data_length = 10;
    double fs = 2.4e9; // carrier frequency 
    double Ts = 1.0 / fs; // sampling period in seconds
    double f0 = 0.0; // homodyne (0 HZ IF)
    double M = 8; // oversampling factor
    double L = 8; // pulse shaping filter length (Is this the ideal length?)
    char pulse_shape[] = "rrc"; // type of pulse shaping, also "rect"
    char scheme[] = "BPSK"; // Modulation scheme 'OOK', 'QPSK', or 'QAM'
    double alpha = 0.5; // roll-off factor of the RRC pulse-shaping filter

    struct Array_Tuple preamble = defineArray((int[]){0,1,0,0,0,0,1,1,0,0,0,1,0,1,0,0,1,1,1,1,0,1,0,0,0,1,1,1,0,0,1,0,0,1,0,1,1,0,1,1,1,0,1,1,0,0,1,1,0,1,0,1,0,1,1,1,1,1,1,0}, 60);// optimal periodic binary code for N = 63 https://ntrs.nasa.gov/citations/19800017860
    struct Array_Tuple data = {randomArray(2,data_length), data_length}; // this points to the random array generated;  creates 256 random bits

    struct Array_Tuple CRC_key = defineArray((int[]){1,0,0,1,1,0,0,0,0,1,1,1}, 12);// Best CRC polynomials: https://users.ece.cmu.edu/~koopman/crc/
    
    
    // Working Demo Code Here ----------------------------------------------------
    // calloc(length, size) is used in defineArray() and for the output variable of CRC_encodeData()

    struct Array_Tuple demoA = defineArray((int[]){1,0,0,1,0,0}, 6);
    struct Array_Tuple demoB = defineArray((int[]){1,1,0,1}, 4);
    struct Array_Tuple out = CRC_encodeData(demoA, demoB);
    
    printArray("Out", out.array, out.length);
    freeArrayMemory(demoA);
    freeArrayMemory(demoB);
    freeArrayMemory(out); // when done with the array, need to free its memory
    // ---------------------------------------------------------------------------
    //struct Array_Length_Tuple tuple = CRC_encodeData(data.array, CRC_key.array);
    return 0;
}