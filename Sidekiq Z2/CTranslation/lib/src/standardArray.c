//standardArray.h
// Rylan Paul
#include <string.h>
#include <limits.h>
#include <stdbool.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <complex.h>
#include <time.h>

// checks if double is an integer
bool isInteger(double val)
{
    int truncated = (int)val;
    return (val == truncated);
}

// Array struct that keeps track of the length of the array and a pointer to the array
typedef struct Array_Tuple {
    double* array;
    int length;
}Array_Tuple;
// Two array structs to repesent array of complex numbers
typedef struct Complex_Array_Tuple {
    Array_Tuple real;
    Array_Tuple imaginary;
}Complex_Array_Tuple;
// Prints an array printArray("Debug Text", ptr to int[], length) (prints out) --> Debug Text: [1, 0, ...]
void printArray(char *string, Array_Tuple arr){
    printf("%s", string);
    printf("(%d): ", arr.length);
    printf("[");
    for (int i=0; i < arr.length; i++){
        if (i != 0){
            printf(", ");
        }
        printf("%lf", arr.array[i]);
    }
    printf("]\n");
}
void printComplexArray(char *string, Complex_Array_Tuple arr){
    for (int i = 0; i < strlen(string); i++) {
        printf("%c", string[i]);
    }
    printf("(%d): ", arr.real.length);
    printf("[ ");
    for (int i=0; i < arr.real.length; i++){
        if (i != 0){
            printf(", ");
        }
        printf("%lf%c%lfj", arr.real.array[i], (arr.imaginary.array[i] < 0 ? '\0' : '+'), arr.imaginary.array[i]);
    }
    printf(" ]\n");
}
// Pass in an Array_Tuple and it will free the memory from that pointer, this only needs to be done on pointers that use the calloc() and malloc() functions
void freeArrayMemory(Array_Tuple array){
    free(array.array);
    return;
}
void freeComplexArrayMemory(Complex_Array_Tuple array){
    free(array.real.array);
    free(array.imaginary.array);
    return;
}
// Creates an Array_Tuple from a known array {1,1,1,0,...} and its known length
Array_Tuple defineArray(double array[], int length){ 
    
    double* ptr;
    ptr = (double*)calloc(length, sizeof(double));
    if (ptr == NULL){
        printf("Not enough memory to allocate, what will happen?\n");
        exit(0);
    }
    for (int i = 0; i < length; i++)
    {
        ptr[i] = array[i];
    }
    
    Array_Tuple tuple = {ptr, length};
    return tuple;
}
// Creates an Complex_Array_Tuple length {0,0,0,0,0} and its known length
Complex_Array_Tuple zerosComplex(int length){ 
    
    double* real;
    real = (double*)calloc(length, sizeof(double));
    double* imaj;
    imaj = (double*)calloc(length, sizeof(double)); // initialized to zero
    
    Array_Tuple real_t = {real, length};
    Array_Tuple imaj_t = {imaj, length};
    Complex_Array_Tuple tuple = {real_t, imaj_t};
    return tuple;
}
// Creates an Array_Tuple given a length with each element set to the given value
Array_Tuple valueArray(int length, int value){ 
    
    double* real;
    real = (double*)calloc(length, sizeof(double));
    for (int i = 0; i < length; i++)
    {
        real[i] = value;
    }
    Array_Tuple real_t = {real, length};
    return real_t;
}
// Creates an Complex_Array_Tuple length {0,0,0,0,0} and its known length
Array_Tuple zerosArray(int length){ 
    
    double* real;
    real = (double*)calloc(length, sizeof(double));
    
    Array_Tuple real_t = {real, length};
    return real_t;
}

// Creates an Complex_Array_Tuple of complex conjugates of the input -> np.conj
Complex_Array_Tuple getConj(Complex_Array_Tuple input){ 
    int length = input.real.length;
    Complex_Array_Tuple out = zerosComplex(length);

    for (int i = 0; i < length; i++)
    {
        out.real.array[i] = input.real.array[i];
        out.imaginary.array[i] = -input.imaginary.array[i];
    }
    return out;
}
// gets the average of an array
double meanArray(double* array, int length){
    double sum = 0;
    for (int i = 0; i < length; i++)
    {
        sum += array[i];
    }
    double avg = sum / (double)length;
    
    return avg; // to use this: int *array; array= randomArray(2,size);
}
// gets the average of an array
double meanArrayTuple(Array_Tuple array){
    double sum = 0;
    for (int i = 0; i < array.length; i++)
    {
        sum += array.array[i];
    }
    double avg = sum / (double)array.length;
    
    return avg; // to use this: int *array; array= randomArray(2,size);
}
// Generates a random pointer to an array
double * randomArray(int max_exclusive, int length){ 
    static int b = 1;
    if (b == 1){
        b = 0;
        time_t t;
        srand((unsigned) time(&t));
    }
    double* arr; // Assuming maximum length of the array
    arr = (double*)calloc(length, sizeof(double));
    for (int i=0; i < length; i++){
        arr[i] = (double)(rand() % max_exclusive);
    }
    return arr; // to use this: int *array; array= randomArray(2,size);
}

Array_Tuple append_array(Array_Tuple a, Array_Tuple b){ 
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
    Array_Tuple t = {ptr, length};
    return t; //once done, must free memory of ptr
}

// reverses the order of an array
Array_Tuple flip(Array_Tuple a){ 
    Array_Tuple t = zerosArray(a.length);
    for (int i = 0; i < a.length; i++)
    {
        t.array[i] = a.array[a.length-1-i];
    }
    return t; //once done, must free memory of ptr
}

// Raises e to the element of every complex number in an array --> np.exp()
Complex_Array_Tuple expComplexArray(Complex_Array_Tuple array) 
{
    // Euler's stuff: e to the power of a complex number -- cool!
    double e = M_E;
    Complex_Array_Tuple out = zerosComplex(array.real.length);

    for (int i = 0; i < array.real.length; i++)
    {
        double power_num  = pow(e, array.real.array[i]);
        out.real.array[i] = power_num * cos(array.imaginary.array[i]); // technically array.imag.array[i] * ln(a), where a^(b+ci), but a=e, so ln(e)=1
        out.imaginary.array[i] = power_num * sin(array.imaginary.array[i]);
        /* code */
    }
    return out;
}
// Raises e to the element of every complex number in an array --> np.exp(), the input is 0+xj, where x is the input
Complex_Array_Tuple exp_imaginaryArray(Array_Tuple array) 
{
    // Euler's stuff: e to the power of a complex number -- cool!
    double e = M_E;
    Complex_Array_Tuple out = zerosComplex(array.length);
    for (int i = 0; i < array.length; i++)
    {
        double real = 0;
        double complex num = real + array.array[i] * I;
        double complex e_num = cexp(num);
        out.real.array[i] = creal(e_num); // technically array.imag.array[i] * ln(a), where a^(b+ci), but a=e, so ln(e)=1
        out.imaginary.array[i] = cimag(e_num);
        /* code */
    }
    return out;
}

//returns the maximum abs of complex array
double maxAbsoluteValue(Complex_Array_Tuple a){
    double max_sqr = -99;
    for (int i = 0; i < a.real.length; i++)
    {
        double x = a.real.array[i];
        double y = a.imaginary.array[i];
        double sqr = x*x + y*y;
        if (sqr > max_sqr){
            max_sqr = sqr;
        }
    }
    return sqrt(max_sqr);
    
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

//np.arange creates an array from start to end (exclusive) by a step
Array_Tuple arange(double start, double end, double step){ 
    int length = 0;
    double num = start;
    while (num <= end){
        num += step;
        length ++;
    }
    Array_Tuple out = zerosArray(length);
    num = start;
    for (int i = 0; i < length; i++)
    {
        out.array[i] = num;
        num += step;
    }
    return out;
}

//np.linspace creates an array from start to end with x number of points defined by length
Array_Tuple linspace(double start, double end, int length){ 
    Array_Tuple out = zerosArray(length);
    double step = (end - start) / (length-1);
    double num = start;
    for (int i = 0; i < length; i++)
    {
        out.array[i] = num;
        num += step;
    }
    return out;
}

Array_Tuple addArrays(Array_Tuple a, Array_Tuple b){ 
    int length = a.length;
    if (b.length > length){
        length = b.length;
    }
    Array_Tuple out = zerosArray(length);

    for (int i = 0; i < length; i++)
    {
        double sum = 0;
        if (i < a.length){
            sum += a.array[i];
        }
        if (i < b.length){
            sum += b.array[i];
        }
        out.array[i] = sum;
    }
    return out;
}
Array_Tuple subtractArrays(Array_Tuple a, Array_Tuple b){ 
    int length = a.length;
    if (b.length > length){
        length = b.length;
    }
    Array_Tuple out = zerosArray(length);

    for (int i = 0; i < length; i++)
    {
        double sum = 0;
        if (i < a.length){
            sum += a.array[i];
        }
        if (i < b.length){
            sum -= b.array[i];
        }
        out.array[i] = sum;
    }
    return out;
}
Array_Tuple subtractDoubleFromArray(Array_Tuple a, double b){ 
    int length = a.length;
    Array_Tuple out = zerosArray(length);

    for (int i = 0; i < length; i++)
    {
        out.array[i] = a.array[i] - b;
    }
    return out;
}
Array_Tuple divideDoubleFromArray(Array_Tuple a, double b){ 
    int length = a.length;
    Array_Tuple out = zerosArray(length);
    for (int i = 0; i < length; i++)
    {
        out.array[i] = a.array[i] / b;
    }
    return out;
}
Array_Tuple multiplyDoubleFromArray(Array_Tuple a, double b){ 
    int length = a.length;
    Array_Tuple out = zerosArray(length);

    for (int i = 0; i < length; i++)
    {
        out.array[i] = a.array[i] * b;
    }
    return out;
}
Array_Tuple multiplyArrays(Array_Tuple a, Array_Tuple b){ 
    int length = a.length;
    if (b.length < length){
        length = b.length;
    }
    Array_Tuple out = zerosArray(length);
    for (int i = 0; i < length; i++)
    {
        out.array[i] = a.array[i] * b.array[i];
    }
    return out;
}
Array_Tuple divideArrays(Array_Tuple a, Array_Tuple b){ 
    int length = a.length;
    if (b.length < length){
        length = b.length;
    }
    Array_Tuple out = zerosArray(length);

    for (int i = 0; i < length; i++)
    {
        if (b.array[i] != 0){
            out.array[i] = a.array[i] / b.array[i];
        } else {
            out.array[i] = INT_MAX;
        }
    }
    return out;
}
// multiplies two complex arrays
Complex_Array_Tuple multiplyComplexArrays(Complex_Array_Tuple x, Complex_Array_Tuple y){ 
    int length = x.real.length;
    if (y.real.length < length){
        length = y.real.length;
    }
    Complex_Array_Tuple out_tuple = zerosComplex(length);

    for (int i = 0; i < length; i++)
    {
        // (a+bi)*(c+di) = (ac-bd)+(ad+bc)i;
        double complex a = x.real.array[i] + I * x.imaginary.array[i];
        double complex b = y.real.array[i] + I * y.imaginary.array[i];
        double complex m = a * b;
        out_tuple.real.array[i] = creal(m);
        out_tuple.imaginary.array[i] = cimag(m);
    }
    return out_tuple;
}
// add two complex arrays together
Complex_Array_Tuple addComplexArrays(Complex_Array_Tuple a, Complex_Array_Tuple b){ 
    int length = a.real.length;
    if (b.real.length > length){
        length = b.real.length;
    }
    Complex_Array_Tuple out = zerosComplex(length);

    for (int i = 0; i < length; i++)
    {
        out.real.array[i] = a.real.array[i] + b.real.array[i];
        out.imaginary.array[i] = a.imaginary.array[i] + b.imaginary.array[i];
    }
    return out;
}

// subtract two complex arrays
Complex_Array_Tuple subtractComplexArrays(Complex_Array_Tuple a, Complex_Array_Tuple b){ 
    int length = a.real.length;
    if (b.real.length > length){
        length = b.real.length;
    }
    Complex_Array_Tuple out = zerosComplex(length);

    for (int i = 0; i < length; i++)
    {
        out.real.array[i] = a.real.array[i] - b.real.array[i];
        out.imaginary.array[i] = a.imaginary.array[i] - b.imaginary.array[i];
    }
    return out;
}

//same as np.sinc, does sin(pi * x) / (pi * x), but lim x-> 0 = 1, so if x=0 return 1
Array_Tuple sinc(Array_Tuple input){ 
    Array_Tuple out = zerosArray(input.length);

    for (int i = 0; i < input.length; i++)
    {
        double x = input.array[i];
        if (x != 0){ // dont wanna ÷ by zero
            out.array[i] = sin(M_PI * x) / (M_PI * x);
        } else {
            out.array[i] = 1; // lim x->0 = 1
        }
    }
    return out;
}
//np.sum Adds them all together
double sumArray(Array_Tuple input){
    double result = 0;
    for (int i = 0; i < input.length; i++)
    {
        result += input.array[i];
    }
    return result;
}

Complex_Array_Tuple everyOtherElement(Complex_Array_Tuple array, int offset){ 
    if (offset >= 1){
        offset = 1;
    } else {
        offset = 0;
    }
    int N = (int)array.real.length / (int)2;
    Complex_Array_Tuple out = zerosComplex(N);
    int l_index = 0;
    for (int i = 0; i+offset < array.real.length; i+=2)
    {
        out.real.array[l_index] = array.real.array[i+offset];
        out.imaginary.array[l_index] = array.imaginary.array[i+offset];
        l_index ++;
    }
    return out;
}

Array_Tuple absComplexArray(Complex_Array_Tuple array){ 
    
    // Store output in the provided struct
    Array_Tuple out = zerosArray(array.real.length);
    // Extract real and imaginary parts from output
    for (int i = 0; i < array.real.length; i++) {
        double real = array.real.array[i];
        double imaginary = array.imaginary.array[i];
        out.array[i] = sqrt(real * real + imaginary * imaginary);
    }
    return out;
}

int argMax(Array_Tuple input){
    int max_index = -1;
    int max_value = INT_MIN; // int min
    for (int i = 0; i < input.length; i++) {
        if (input.array[i] >= max_value){
            max_value = input.array[i];
            max_index = i;
        }
    }
    return max_index;
}
int argComplexMax(Complex_Array_Tuple input){
    int max_index = -1;
    int max_value = INT_MIN; // int min
    int max_value_imaj = INT_MIN; // int min
    for (int i = 0; i < input.real.length; i++) {
        if (input.real.array[i] > max_value){
            max_value = input.real.array[i];
            max_value_imaj = input.imaginary.array[i];
            max_index = i;
        } else if (input.real.array[i] == max_value && input.imaginary.array[i] > max_value_imaj){
            max_value = input.real.array[i];
            max_value_imaj = input.imaginary.array[i];
            max_index = i;
        }
    }
    return max_index;
}

//np.convolve: Circular Discrete Convolution --> https://en.wikipedia.org/wiki/Convolution
Complex_Array_Tuple convolve(Complex_Array_Tuple a, Array_Tuple v){ 
    int N = a.real.length + v.length - 1;
    Complex_Array_Tuple out = zerosComplex(N);

    // Perform convolution
    for (int n = 0; n < N; n++) {
        double sum_real = 0;
        double sum_imaj = 0;
        for (int M = 0; M < n+1; M++) {
            if (n-M >= 0 && n-M < v.length && M < a.real.length){
                sum_real += a.real.array[M] * v.array[n-M];
                sum_imaj += a.imaginary.array[M] * v.array[n-M];
            }
        }
        for (int M = n+1; M < N; M++)
        {
            if (n-M >= 0 && N+n-M < v.length && M < a.real.length){
                sum_real += a.real.array[M] * v.array[N+n-M];
                sum_imaj += a.imaginary.array[M] * v.array[N+n-M];
            }
        }
        
        out.real.array[n] = sum_real;
        out.imaginary.array[n] = sum_imaj;
    }
    return out;
}

Complex_Array_Tuple convolveSame(Complex_Array_Tuple a, Array_Tuple v){ 
    int N = a.real.length + v.length - 1;
    Complex_Array_Tuple long_out = zerosComplex(N);

    // Perform convolution
    for (int n = 0; n < N; n++) {
        double sum_real = 0;
        double sum_imaj = 0;
        for (int M = 0; M < n+1; M++) {
            if (n-M >= 0 && n-M < v.length && M < a.real.length){
                sum_real += a.real.array[M] * v.array[n-M];
                sum_imaj += a.imaginary.array[M] * v.array[n-M];
            }
        }
        for (int M = n+1; M < N; M++)
        {
            if (n-M >= 0 && N+n-M < v.length && M < a.real.length){
                sum_real += a.real.array[M] * v.array[N+n-M];
                sum_imaj += a.imaginary.array[M] * v.array[N+n-M];
            }
        }
        long_out.real.array[n] = sum_real;
        long_out.imaginary.array[n] = sum_imaj;
    }
    // remove difference/2 from both sides, this is what makes it different than np.convolve(mode='full')
    int output_length = a.real.length > v.length ? a.real.length : v.length;
    Complex_Array_Tuple shortened_out = zerosComplex(output_length);
    
    int difference = abs(N - output_length);
    double remove = (double)difference / 2.0;

    if (isInteger(remove)){
        for (int i = 0; i < output_length; i++)
        {
            shortened_out.real.array[i] = long_out.real.array[(int)remove + i];
            shortened_out.imaginary.array[i] = long_out.imaginary.array[(int)remove + i];
        }
    } else {
        for (int i = 0; i < output_length; i++)
        {
            shortened_out.real.array[i] = long_out.real.array[(int)floor(remove) + i];
            shortened_out.imaginary.array[i] = long_out.imaginary.array[(int)floor(remove) + i];
        }
    }
    freeComplexArrayMemory(long_out);
    return shortened_out;
}


void exportArray(Array_Tuple input, char filename[]){
    FILE *fpt;
    char export_name[1000] = "lib/graphs/";
    int start = strlen(export_name);
    for (int i = 0; i < strlen(filename); i++)
    {
        export_name[i+start] = filename[i];
    }
    export_name[strlen(export_name)] = '\0';
    fpt = fopen(export_name, "w+");
    for (int i = 0; i < input.length; i++)
    {
        fprintf(fpt, "%lf", input.array[i]);
        if (i != input.length-1){
            fprintf(fpt, "\n");
        }
    }
    fclose(fpt);
}

void exportComplexArray(Complex_Array_Tuple input, char filename[]){
    FILE *fpt;
    char export_name[1000] = "lib/graphs/";
    int start = strlen(export_name);
    for (int i = 0; i < strlen(filename); i++)
    {
        export_name[i+start] = filename[i];
    }
    export_name[strlen(export_name)] = '\0';
    fpt = fopen(export_name, "w+");
    for (int i = 0; i < input.real.length; i++)
    {
        if (input.imaginary.array[i] < 0){
            fprintf(fpt, "%lf%lfj", input.real.array[i], input.imaginary.array[i]);
        }else {
            fprintf(fpt, "%lf+%lfj", input.real.array[i], input.imaginary.array[i]);
        }
        
        if (i != input.real.length-1){
            fprintf(fpt, "\n");
        }
    }
    fclose(fpt);
}

void complexArrayToCharArray(Complex_Array_Tuple array) {
    char input[90000];
    int length = 0;
    for (int i = 0; i < array.real.length; i++)
    {
        char buffer_real[25];
        sprintf(buffer_real, "%.6f", array.real.array[i]);
        char buffer_imag[25];
        sprintf(buffer_imag, "%.6f", array.imaginary.array[i]);
        int x = 0;
        while(buffer_real[x] != '\0'){
            input[length] = buffer_real[x];
            length++;
            x++;
        }
        input[length] = ',';
        length++;
        x = 0;
        while(buffer_imag[x] != '\0'){
            input[length] = buffer_imag[x];
            length++;
            x++;
        }
        input[length] = ',';
        length++;
    }
    input[length] = '\0';
    length++;
    printf("%s", input);
    printf("\nLength: %d\n", length);
}

Complex_Array_Tuple charArrayToComplexArray(char *str){
    double *real = (double*)malloc(256 * sizeof(double));
    double *imag = (double*)malloc(256 * sizeof(double));

    for (int i = 0; i < strlen(str); i++)
    {
        /* code */
        printf("Hi");
    }
    

}