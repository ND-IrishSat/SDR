#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void printArray(char *string, int* arr, int length){
    for (int i = 0; i < strlen(string); i++)
    {
        printf("%c", string[i]);
        /* code */
    }
    printf("(%d): ", length);
    printf("[");
    for (int i=0; i < length; i++){
        if (i != 0){
            printf(", ");
        }
        printf("%d", arr[i]);
    }
    printf("]\n");
}
struct Array_Length_Tuple {
    int* array;
    int length;
};

// CRC
struct Array_Length_Tuple CRC_xor(struct Array_Length_Tuple a_array, struct Array_Length_Tuple b_array){
    int length = b_array.length;
    int* out_array;
    out_array = (int*)calloc(length-1, sizeof(int));
    for (int i = 1; i < length; i++){
        int element1 = a_array.array[i];
        int element2 = b_array.array[i];
        if (element1 == element2){
            out_array[i-1] = 0;
        } else {
            out_array[i-1] = 1;
        }
    }
    struct Array_Length_Tuple tuple = {out_array, length-1};
    return tuple;
}
struct Array_Length_Tuple CRC_mod2div(struct Array_Length_Tuple dividend_array, struct Array_Length_Tuple divisor_array){
    int pick = divisor_array.length;
    
    static int* temp;
    temp = (int*)calloc(pick, sizeof(int));
    
    for (int i=0; i < pick; i++){
        temp[i] = dividend_array.array[i];
    }
    struct Array_Length_Tuple tmp = {temp, pick};

    while (pick < dividend_array.length){
        if (tmp.array[0] == 1){
            // tmp = np.append(xor(divisor, tmp), dividend[pick])
            int* new_tmp;
            new_tmp = (int*)calloc(divisor_array.length, sizeof(int));
            struct Array_Length_Tuple x = CRC_xor(divisor_array, tmp); // comes out to 1 less length than tmp
            int last_index;
            for (int i=0; i < x.length; i++){
                new_tmp[i] = x.array[i];
                last_index = i;
            }
            new_tmp[last_index + 1] = dividend_array.array[pick];
            for (int i=0; i < divisor_array.length; i++){
                tmp.array[i] = new_tmp[i];
            }
        }else {
            int* new_tmp_2;
            new_tmp_2 = (int*)calloc(pick, sizeof(int));
            int* ones_array;
            ones_array = (int*)calloc(pick, sizeof(int));
            for(int i=0; i < pick;  i++){
                ones_array[i] = 1;
            }
            struct Array_Length_Tuple ones_tuple = {ones_array, pick};
            struct Array_Length_Tuple xor_out = CRC_xor(ones_tuple, tmp); // comes out to 1 less length than tmp
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
        struct Array_Length_Tuple xor_out_2 = CRC_xor(divisor_array, tmp);
        struct Array_Length_Tuple out_array = {xor_out_2.array, xor_out_2.length};
        return out_array;
    }
    else{
        int* zeros_array;
        for(int i=0; i < pick; i++){
            zeros_array[i] = 0;
        }
        struct Array_Length_Tuple zeros_tuple = {zeros_array, pick};
        struct Array_Length_Tuple xor_out_2 = CRC_xor(zeros_tuple, tmp);
        struct Array_Length_Tuple out_array = {xor_out_2.array, xor_out_2.length};
        return out_array;
    }
}
struct Array_Length_Tuple CRC_encodeData(struct Array_Length_Tuple data, struct Array_Length_Tuple key){
    // Appends n-1 zeroes at end of data
    int* appended_data;
    appended_data = (int*)calloc(data.length+key.length-1, sizeof(int)); // when do I need to and not need to define  the size?
    for (int i=0; i<data.length; i++){
        appended_data[i] = data.array[i];
    }
    for (int i=data.length; i<data.length+key.length-1; i++){
        appended_data[i] = 0;
    }
    printArray("appended", appended_data, data.length+key.length-1);
    struct Array_Length_Tuple appended_data_tuple = {appended_data, data.length+key.length-1};
    struct Array_Length_Tuple mod2div_out = CRC_mod2div(appended_data_tuple, key);
    // Append remainder in the original data
    int* codeword;
    codeword = (int*)calloc(data.length+mod2div_out.length, sizeof(int));
    for (int i=0; i<data.length; i++){
        codeword[i] = data.array[i];
    }
    for (int i=data.length; i<data.length+mod2div_out.length; i++){
        codeword[i] = mod2div_out.array[i - data.length];
    }
    struct Array_Length_Tuple tuple = {codeword, data.length+mod2div_out.length};
    return tuple;
}
int CRC_check(struct Array_Length_Tuple codeword, struct Array_Length_Tuple key){
    struct Array_Length_Tuple checkword = CRC_mod2div(codeword, key);
    int error = 0;
    for (int i = 0; i < checkword.length; i++)
    {
        if (checkword.array[i] == 0){
            error = 1;
            break;
        }
    }
    return error;
    
}
// CRC
struct Array_Length_Tuple defineArray(int array[], int length){
    int* arr;
    arr = (int*)calloc(length, sizeof(int));
    for (int i = 0; i < length; i++)
    {
        arr[i] = array[i];
    }
    struct Array_Length_Tuple tuple = {arr, length};
    return tuple;
}
int * randomArray(int max_exclusive, int length){ // returns a pointer to an integer
    static int* arr; // static integer pointer, static means it will persist throughout the programs runtime
    arr = (int*)calloc(length, sizeof(int)); // tells the program how much storage to save for this array
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

    int a[] = {0,1,0,0,0,0,1,1,0,0,0,1,0,1,0,0,1,1,1,1,0,1,0,0,0,1,1,1,0,0,1,0,0,1,0,1,1,0,1,1,1,0,1,1,0,0,1,1,0,1,0,1,0,1,1,1,1,1,1,0};
    struct Array_Length_Tuple preamble = defineArray(a, 60);// optimal periodic binary code for N = 63 https://ntrs.nasa.gov/citations/19800017860
    struct Array_Length_Tuple data = {randomArray(2,data_length), data_length}; // this points to the random array generated;  creates 256 random bits
    //printArray(data.array, data_length);
    int b[] = {1,0,0,1,1,0,0,0,0,1,1,1};
    struct Array_Length_Tuple CRC_key = defineArray(b, 12);// Best CRC polynomials: https://users.ece.cmu.edu/~koopman/crc/
    
    int dA[] = {1, 0, 0, 1, 0, 0, 1};
    int dB[] = {0, 1, 0, 1, 1};
    struct Array_Length_Tuple demoA = defineArray(dA, 7);
    struct Array_Length_Tuple demoB = defineArray(dB, 5);
    struct Array_Length_Tuple out = CRC_encodeData(demoA, demoB);
    printArray("Out", out.array, out.length);
    //struct Array_Length_Tuple tuple = CRC_encodeData(data.array, CRC_key.array);
    return 0;
}
