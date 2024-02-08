//standardArray.h
// Rylan Paul

// Comments:
// * highlight
//! Warning
//? Question
//TODO
//

#include <string.h>
#include <stdbool.h>
#include <stdlib.h>
#include <stdio.h>
// Prints an array printArray("Debug Text", ptr to int[], length) (prints out) --> Debug Text: [1, 0, ...]
void printArray(char *string, double* arr, int length){
    for (int i = 0; i < strlen(string); i++) {
        printf("%c", string[i]);
    }
    printf("(%d): ", length);
    printf("[");
    for (int i=0; i < length; i++){
        if (i != 0){
            printf(", ");
        }
        printf("%lf", arr[i]);
    }
    printf("]\n");
}
void printComplexArray(char *string, double* arr, double* imag_arr, int length){
    for (int i = 0; i < strlen(string); i++) {
        printf("%c", string[i]);
    }
    printf("(%d): ", length);
    printf("[");
    for (int i=0; i < length; i++){
        if (i != 0){
            printf(", ");
        }
        printf("%lf+%lfj", arr[i], imag_arr[i]);
    }
    printf("]\n");
}



// Array Tuple that keeps track of the length of the array and a pointer to the array
struct Array_Tuple {
    double* array;
    int length;
};
//
struct Complex_Array_Tuple {
    struct Array_Tuple real;
    struct Array_Tuple imaginary;
};

// Pass in an Array_Tuple and it will free the memory from that pointer, this only needs to be done on pointers that use the calloc() and malloc() functions
void freeArrayMemory(struct Array_Tuple array){
    free(array.array);
    return;
}
// Creates an Array_Tuple from a known array {1,1,1,0,...} and its known length
struct Array_Tuple defineArray(double array[], int length){ //! Returns a calloc ptr
    
    double* ptr;
    ptr = (double*)calloc(length, sizeof(double));
    if (ptr == NULL){
        printf("Not enough memory to allocate, what should happen?\n");
        exit(0);
    }
    for (int i = 0; i < length; i++)
    {
        ptr[i] = array[i];
    }
    
    struct Array_Tuple tuple = {ptr, length};
    return tuple;
}
double meanArray(double* array, int length){ // returns a pointer to an integer array
    double sum = 0;
    for (int i = 0; i < length; i++)
    {
        sum += array[i];
    }
    double avg = sum / (double)length;
    
    return avg; // to use this: int *array; array= randomArray(2,size);
}
// Generates a random pointer to an array
double * randomArray(int max_exclusive, int length){ //! Returns a calloc ptr
    double* arr; // Assuming maximum length of the array
    arr = (double*)calloc(length, sizeof(dobule))
    for (int i=0; i < length; i++){
        arr[i] = (double)(rand() % max_exclusive);
    }
    return arr; // to use this: int *array; array= randomArray(2,size);
}

struct Array_Tuple append_array(struct Array_Tuple a, struct Array_Tuple b){ //! Returns a calloc ptr
    int length = a.length + b.length;
    double* ptr;
    ptr = (double*)calloc(length, sizeof(double));
    for (int i = 0; i < a.length; i++)
    {
        ptr[i] = a.array[i];
    }
    for (int i = 0; i < b.length; i++)
    {
        ptr[a.length+i] = b.array[i];
    }
    struct Array_Tuple t = {ptr, length};
    return t; //once done, must free memory of ptr
}

// reverses the order of an array
struct Array_Tuple flip(struct Array_Tuple a){ //! Returns a calloc ptr
    double* ptr;
    ptr = (double*)calloc(a.length, sizeof(double));
    for (int i = 0; i < a.length; i++)
    {
        ptr[i] = a.array[a.length-1-i];
    }
    struct Array_Tuple t = {ptr, a.length};
    return t; //once done, must free memory of ptr
}

// Raises e to the element of every complex number in an array --> np.exp()
struct Complex_Array_Tuple exp_array(struct Complex_Array_Tuple array) //! Returns a calloc ptr
{
    // Euler's stuff: e to the power of a complex number -- cool!
    double e = M_E;
    double* real_ptr;
    real_ptr = (double*)calloc(array.real.length, sizeof(double));
    double* imag_ptr;
    imag_ptr = (double*)calloc(array.imaginary.length, sizeof(double));

    for (int i = 0; i < array.real.length; i++)
    {
        double power_num  = pow(e, array.real.array[i]);
        real_ptr[i] = power_num * cos(array.imaginary.array[i]); // technically array.imag.array[i] * ln(a), where a^(b+ci), but a=e, so ln(e)=1
        imag_ptr[i] = power_num * sin(array.imaginary.array[i]);
        /* code */
    }
    struct Array_Tuple real = {real_ptr, array.real.length};
    struct Array_Tuple imag = {imag_ptr, array.imaginary.length};
    struct Complex_Array_Tuple out = {real, imag};
    return out;
}

// returns a random number with a normal distribution
double rand_norm(double mu, double sigma){
    double U1, U2, W, mult;
    static double rand_norm_X1, rand_norm_X2;
    static int call_norm = 0;

    if (call_norm == 1)
    {
        call_norm = !call_norm;
        return (mu + sigma * (double)rand_norm_X2);
    }

    do
    {
        U1 = -1 + ((double)rand() / RAND_MAX) * 2;
        U2 = -1 + ((double)rand() / RAND_MAX) * 2;
        W = pow(U1, 2) + pow(U2, 2);
    } while (W >= 1 || W == 0);

    mult = sqrt((-2 * log(W)) / W);
    rand_norm_X1 = U1 * mult;
    rand_norm_X2 = U2 * mult;

    call_norm = !call_norm;

    return (mu + sigma * (double)rand_norm_X1);
}

// Takes a bit array and converts it to a pulse train from -1 to 1 with sps-1 zeros between each bit
struct Array_Tuple pulsetrain(struct Array_Tuple bits, int sps){ //! Returns a calloc ptr
    double* array_ptr;
    int length = bits.length * sps;
    array_ptr = (double*)calloc(length, sizeof(double));
    int offset = 0;
    for (int i = 0; i < bits.length; i++)
    {
        for (int x = 0; x < sps; x++)
        {
            if (x == 0){
                array_ptr[x+offset] = bits.array[i]*2 - 1; // sets 0 to -1 and 1 to 1
            } else {
                array_ptr[x+offset] = 0; // puts zeros between the bits to make it easier to read
            }
        }
        offset += sps;
    }
    struct Array_Tuple t = {array_ptr, length};
    return t;
    
}

struct Complex_Array_Tuple generateNoise(struct Complex_Array_Tuple testpacket){ //! Returns a calloc ptr
    // vars
    int max = 1;
    double mean = 0;
    double std_dev = 0.1;
    double noise_power = 0.2;
    int num_samples = testpacket.real.length;
    double phase_noise_strength = 0.1;

    // creating general normal distribution random arrays for real and complex
    double* awgn_comlpex_samples_real;
    double* awgn_comlpex_samples_imag;
    awgn_comlpex_samples_real = (double*)calloc(num_samples, sizeof(double));
    awgn_comlpex_samples_imag = (double*)calloc(num_samples, sizeof(double));
    for (int i = 0; i < num_samples; i++)
    {
        double phase_noise = rand_norm(mean, std_dev) * phase_noise_strength;
        awgn_comlpex_samples_real[i] = rand_norm(mean, std_dev) / sqrt(noise_power) * cos(phase_noise);
        awgn_comlpex_samples_imag[i] = rand_norm(mean, std_dev) / sqrt(noise_power) * sin(phase_noise);
    }

    struct Array_Tuple real = {awgn_comlpex_samples_real, num_samples};
    struct Array_Tuple imag = {awgn_comlpex_samples_imag, num_samples};
    struct Complex_Array_Tuple out = {real, imag};
    return out;
    //awgn_comlpex_samples_real = ;
    //awgn_comlpex_samples_imaginary = ;
}




















