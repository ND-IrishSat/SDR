//standardArray.h
// Rylan Paul

// Comments:
// * highlight
//! Warning
//? Question
//TODO
//
#include <string.h>
#include <limits.h>
#include <stdbool.h>
#include <stdlib.h>
#include <stdio.h>
#include <fftw3.h>

// checks if double is an integer
bool isInteger(double val)
{
    int truncated = (int)val;
    return (val == truncated);
}

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
    printf("[ ");
    for (int i=0; i < length; i++){
        if (i != 0){
            printf(", ");
        }
        printf("%lf%c%lfj", arr[i], (imag_arr[i] < 0 ? '\0' : '+'), imag_arr[i]);
    }
    printf(" ]\n");
}

// Array Tuple that keeps track of the length of the array and a pointer to the array
typedef struct Array_Tuple {
    double* array;
    int length;
}Array_Tuple;
//
typedef struct Complex_Array_Tuple {
    Array_Tuple real;
    Array_Tuple imaginary;
}Complex_Array_Tuple;

// Pass in an Array_Tuple and it will free the memory from that pointer, this only needs to be done on pointers that use the calloc() and malloc() functions
void freeArrayMemory(Array_Tuple array){
    free(array.array);
    return;
}
// Creates an Array_Tuple from a known array {1,1,1,0,...} and its known length
Array_Tuple defineArray(double array[], int length){ //! Returns a calloc ptr
    
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
    
    Array_Tuple tuple = {ptr, length};
    return tuple;
}

// Creates an Complex_Array_Tuple length {0,0,0,0,0} and its known length
Complex_Array_Tuple zerosComplex(int length){ //! Returns a calloc ptr
    
    double* real;
    real = (double*)calloc(length, sizeof(double));
    double* imaj;
    imaj = (double*)calloc(length, sizeof(double)); // initialized to zero
    
    Array_Tuple real_t = {real, length};
    Array_Tuple imaj_t = {imaj, length};
    Complex_Array_Tuple tuple = {real_t, imaj_t};
    return tuple;
}
// Creates an Complex_Array_Tuple length {0,0,0,0,0} and its known length
Array_Tuple zerosArray(int length){ //! Returns a calloc ptr
    
    double* real;
    real = (double*)calloc(length, sizeof(double));
    
    Array_Tuple real_t = {real, length};
    return real_t;
}

// Creates an Complex_Array_Tuple of complex conjugates of the input -> np.conj
Complex_Array_Tuple getConj(Complex_Array_Tuple input){ //! Returns a calloc ptr
    int length = input.real.length;
    double* real;
    real = (double*)calloc(length, sizeof(double));
    double* imaj;
    imaj = (double*)calloc(length, sizeof(double)); // initialized to zero

    for (int i = 0; i < length; i++)
    {
        real[i] = input.real.array[i];
        imaj[i] = -input.imaginary.array[i];
    }
    
    
    Array_Tuple real_t = {real, length};
    Array_Tuple imaj_t = {imaj, length};
    Complex_Array_Tuple tuple = {real_t, imaj_t};
    return tuple;
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
double * randomArray(int max_exclusive, int length){ //! Returns a calloc ptr
    double* arr; // Assuming maximum length of the array
    arr = (double*)calloc(length, sizeof(double));
    for (int i=0; i < length; i++){
        arr[i] = (double)(rand() % max_exclusive);
    }
    return arr; // to use this: int *array; array= randomArray(2,size);
}

Array_Tuple append_array(Array_Tuple a, Array_Tuple b){ //! Returns a calloc ptr
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
Array_Tuple flip(Array_Tuple a){ //! Returns a calloc ptr
    double* ptr;
    ptr = (double*)calloc(a.length, sizeof(double));
    for (int i = 0; i < a.length; i++)
    {
        ptr[i] = a.array[a.length-1-i];
    }
    Array_Tuple t = {ptr, a.length};
    return t; //once done, must free memory of ptr
}


// Raises e to the element of every complex number in an array --> np.exp()
Complex_Array_Tuple expComplexArray(Complex_Array_Tuple array) //! Returns a calloc ptr
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
    Array_Tuple real = {real_ptr, array.real.length};
    Array_Tuple imag = {imag_ptr, array.imaginary.length};
    Complex_Array_Tuple out = {real, imag};
    return out;
}
// Raises e to the element of every complex number in an array --> np.exp(), the input is 0+xj, where x is the input
Complex_Array_Tuple exp_imaginaryArray(Array_Tuple array) //! Returns a calloc ptr
{
    // Euler's stuff: e to the power of a complex number -- cool!
    double e = M_E;
    double* real_ptr;
    real_ptr = (double*)calloc(array.length, sizeof(double));
    double* imag_ptr;
    imag_ptr = (double*)calloc(array.length, sizeof(double));
    for (int i = 0; i < array.length; i++)
    {
        double real = 0;
        double power_num  = pow(e, real);
        real_ptr[i] = power_num * cos(array.array[i]); // technically array.imag.array[i] * ln(a), where a^(b+ci), but a=e, so ln(e)=1
        imag_ptr[i] = power_num * sin(array.array[i]);
        /* code */
    }
    Array_Tuple real = {real_ptr, array.length};
    Array_Tuple imag = {imag_ptr, array.length};
    Complex_Array_Tuple out = {real, imag};
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
Array_Tuple arange(double start, double end, double step){ //! Returns a calloc ptr
    int length = 0;
    static double array[2048];
    double num = start;
    while (num < end){
        array[length] = num;
        num += step;
        length += 1;
    }
    double* out_array;
    out_array = (double*)calloc(length, sizeof(double));
    for (int i = 0; i < length; i++)
    {
        out_array[i] = array[i];
    }
    
    Array_Tuple out = {out_array, length};
    return out;
}

//np.linspace creates an array from start to end with x number of points defined by length
Array_Tuple linspace(double start, double end, int length){ //! Returns a calloc ptr
    double* out_array;
    out_array = (double*)calloc(length, sizeof(double));
    double step = (end - start) / (length-1);
    double num = start;
    for (int i = 0; i < length; i++)
    {
        out_array[i] = num;
        num += step;
    }
    
    Array_Tuple out = {out_array, length};
    return out;
}

Array_Tuple addArrays(Array_Tuple a, Array_Tuple b){ //! Returns a calloc ptr
    int length = a.length;
    if (b.length > length){
        length = b.length;
    }
    double* output;
    output = (double*) calloc(length, sizeof(double));

    for (int i = 0; i < length; i++)
    {
        double sum = 0;
        if (i < a.length){
            sum += a.array[i];
        }
        if (i < b.length){
            sum += b.array[i];
        }
        output[i] = sum;
    }
    Array_Tuple out = {output, length};
    return out;
}
Array_Tuple subtractArrays(Array_Tuple a, Array_Tuple b){ //! Returns a calloc ptr
    int length = a.length;
    if (b.length > length){
        length = b.length;
    }
    double* output;
    output = (double*) calloc(length, sizeof(double));

    for (int i = 0; i < length; i++)
    {
        double sum = 0;
        if (i < a.length){
            sum += a.array[i];
        }
        if (i < b.length){
            sum -= b.array[i];
        }
        output[i] = sum;
    }
    Array_Tuple out = {output, length};
    return out;
}
Array_Tuple subtractDoubleFromArray(Array_Tuple a, double b){ //! Returns a calloc ptr
    int length = a.length;
    double* output;
    output = (double*) calloc(length, sizeof(double));

    for (int i = 0; i < length; i++)
    {
        double sum = 0;
        sum += a.array[i] - b;
        output[i] = sum;
    }
    Array_Tuple out = {output, length};
    return out;
}
Array_Tuple divideDoubleFromArray(Array_Tuple a, double b){ //! Returns a calloc ptr
    int length = a.length;
    double* output;
    output = (double*) calloc(length, sizeof(double));

    for (int i = 0; i < length; i++)
    {
        output[i] = a.array[i] / b;
    }
    Array_Tuple out = {output, length};
    return out;
}
Array_Tuple multiplyDoubleFromArray(Array_Tuple a, double b){ //! Returns a calloc ptr
    int length = a.length;
    double* output;
    output = (double*) calloc(length, sizeof(double));

    for (int i = 0; i < length; i++)
    {
        output[i] = a.array[i] * b;
    }
    Array_Tuple out = {output, length};
    return out;
}
Array_Tuple multiplyArrays(Array_Tuple a, Array_Tuple b){ //! Returns a calloc ptr
    int length = a.length;
    if (b.length < length){
        length = b.length;
    }
    double* output;
    output = (double*) calloc(length, sizeof(double));

    for (int i = 0; i < length; i++)
    {
        output[i] = a.array[i] * b.array[i];
    }
    Array_Tuple out = {output, length};
    return out;
}
Array_Tuple divideArrays(Array_Tuple a, Array_Tuple b){ //! Returns a calloc ptr
    int length = a.length;
    if (b.length < length){
        length = b.length;
    }
    double* output;
    output = (double*) calloc(length, sizeof(double));

    for (int i = 0; i < length; i++)
    {
        if (b.array[i] != 0){
            output[i] = a.array[i] / b.array[i];
        } else {
            output[i] = INT_MAX;
        }
    }
    Array_Tuple out = {output, length};
    return out;
}
// multiplies two complex arrays
Complex_Array_Tuple multiplyComplexArrays(Complex_Array_Tuple x, Complex_Array_Tuple y){ //! Returns a calloc ptr
    int length = x.real.length;
    if (y.real.length < length){
        length = y.real.length;
    }
    double* output_real;
    output_real = (double*) calloc(length, sizeof(double));
    double* output_imaj;
    output_imaj = (double*) calloc(length, sizeof(double));

    for (int i = 0; i < length; i++)
    {
        // (a+bi)*(c+di) = (ac-bd)+(ad+bc)i;
        double a = x.real.array[i];
        double b = x.imaginary.array[i];
        double c = y.real.array[i];
        double d = y.imaginary.array[i];
        output_real[i] = a*c - b*d;
        output_imaj[i] = a*d + b*c;
    }
    Array_Tuple out_real = {output_real, length};
    Array_Tuple out_imaj = {output_imaj, length};
    Complex_Array_Tuple out_tuple = {out_real, out_imaj};
    return out_tuple;
}
// add two complex arrays together
Complex_Array_Tuple addComplexArrays(Complex_Array_Tuple a, Complex_Array_Tuple b){ //! Returns a calloc ptr
    int length = a.real.length;
    if (b.real.length > length){
        length = b.real.length;
    }
    double* output_real;
    output_real = (double*) calloc(length, sizeof(double));
    double* output_imaj;
    output_imaj = (double*) calloc(length, sizeof(double));

    for (int i = 0; i < length; i++)
    {
        output_real[i] = a.real.array[i] + b.real.array[i];
        output_imaj[i] = a.imaginary.array[i] + b.imaginary.array[i];
    }
    Array_Tuple out_real = {output_real, length};
    Array_Tuple out_imaj = {output_imaj, length};
    Complex_Array_Tuple out_tuple = {out_real, out_imaj};
    return out_tuple;
}

// subtract two complex arrays
Complex_Array_Tuple subtractComplexArrays(Complex_Array_Tuple a, Complex_Array_Tuple b){ //! Returns a calloc ptr
    int length = a.real.length;
    if (b.real.length > length){
        length = b.real.length;
    }
    double* output_real;
    output_real = (double*) calloc(length, sizeof(double));
    double* output_imaj;
    output_imaj = (double*) calloc(length, sizeof(double));

    for (int i = 0; i < length; i++)
    {
        output_real[i] = a.real.array[i] - b.real.array[i];
        output_imaj[i] = a.imaginary.array[i] - b.imaginary.array[i];
    }
    Array_Tuple out_real = {output_real, length};
    Array_Tuple out_imaj = {output_imaj, length};
    Complex_Array_Tuple out_tuple = {out_real, out_imaj};
    return out_tuple;
}

//same as np.sinc, does sin(pi * x) / (pi * x), but lim x-> 0 = 1, so if x=0 return 1
Array_Tuple sinc(Array_Tuple input){ //! Returns a calloc ptr
    double* ptr;
    ptr = (double*)calloc(input.length, sizeof(double));

    for (int i = 0; i < input.length; i++)
    {
        double x = input.array[i];
        if (x != 0){ // dont wanna ÷ by zero
            ptr[i] = sin(M_PI * x) / (M_PI * x);
        } else {
            ptr[i] = 1; // lim x->0 = 1
        }
    }
    
    Array_Tuple out = {ptr, input.length};
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
// Takes a bit array and converts it to a pulse train from -1 to 1 with sps-1 zeros between each bit
Array_Tuple pulsetrain(Array_Tuple bits, int sps){ //! Returns a calloc ptr
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
    Array_Tuple t = {array_ptr, length};
    return t;
    
}

Complex_Array_Tuple everyOtherElement(Complex_Array_Tuple array, int offset){ //! Returns a calloc ptr
    if (offset >= 1){
        offset = 1;
    } else {
        offset = 0;
    }
    double* real;
    real = (double*)calloc(array.real.length, sizeof(double));
    double* imaginary;
    imaginary = (double*)calloc(array.imaginary.length, sizeof(double));
    int l_index = 0;
    for (int i = 0; i+offset < array.real.length; i+=2)
    {
        real[l_index] = array.real.array[i+offset];
        imaginary[l_index] = array.imaginary.array[i+offset];
        l_index ++;
    }
    int N = (int)array.real.length / (int)2;
    Array_Tuple i_out = {real, N};
    Array_Tuple q_out = {imaginary, N};
    Complex_Array_Tuple complex_out = {i_out, q_out};
    return complex_out;
}
// Using the fftw3 library, calculates discrete fft on an array np.fft.fft()
Complex_Array_Tuple fft(Complex_Array_Tuple array){ //! Returns a calloc ptr
    int N = array.real.length;
    fftw_complex *in_fftw, *out_fftw;
    fftw_plan p;
    in_fftw = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
    out_fftw = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);

    p = fftw_plan_dft_1d(N, in_fftw, out_fftw, FFTW_FORWARD, FFTW_ESTIMATE);

    // Populate input array with real and imaginary parts
    for (int i = 0; i < N; i++) {
        double r = array.real.array[i];
        double im = array.imaginary.array[i];
        in_fftw[i][0] = r; // Real part
        in_fftw[i][1] = im; // Imaginary part
    }
    // Execute FFT
    fftw_execute(p);

    // Execute FFT
    fftw_execute(p);

    // Store output in the provided struct
    Complex_Array_Tuple out_struct = zerosComplex(N);

    // Extract real and imaginary parts from output
    for (int i = 0; i < N; i++) {
        double r = out_fftw[i][0];
        double im = out_fftw[i][1];
        out_struct.real.array[i] = r; // Real part
        out_struct.imaginary.array[i] = im; // Imaginary part
    }

    // Free memory and destroy plan
    fftw_destroy_plan(p);
    fftw_free(in_fftw);
    fftw_free(out_fftw);
    return out_struct;
}

Array_Tuple absComplexArray(Complex_Array_Tuple array){ //! Returns a calloc ptr
    
    // Store output in the provided struct
    double* result;
    result = (double*) calloc(array.real.length, sizeof(double));
    // Extract real and imaginary parts from output
    for (int i = 0; i < array.real.length; i++) {
        double real = array.real.array[i];
        double imaginary = array.imaginary.array[i];
        result[i] = sqrt(real * real + imaginary * imaginary);
    }

    Array_Tuple out = {result, array.real.length};
    return out;
}

int argMax(Array_Tuple input){
    int max_index = -1;
    int max_value = INT_MIN;
    // Extract real and imaginary parts from output
    for (int i = 0; i < input.length; i++) {
        if (input.array[i] > max_value){
            max_value = input.array[i];
            max_index = i;
        }
    }
    return max_index;
}

// Function to perform FFT shift along a single dimension
Array_Tuple fftshift(Array_Tuple data) { //! Returns a calloc ptr
    double temp;
    // Calculate the midpoint of the array
    int midpoint = data.length / 2;
    double* output;
    output = (double*) calloc(midpoint*2, sizeof(double));
    // Perform the shift
    for (int i = 0; i < midpoint; i++) {
        temp = data.array[i];
        output[i] = data.array[i + midpoint];
        output[i + midpoint] = temp;
    }

    Array_Tuple out = {output, midpoint*2};
    return out;
}

//fftshift for complex array
Complex_Array_Tuple complexfftshift(Complex_Array_Tuple input) { //! Returns two calloc ptrs from fft_shift
    // Perform FFT shift on both real and imaginary parts
    
    Array_Tuple real = fftshift(input.real);
    Array_Tuple imag = fftshift(input.imaginary);
    Complex_Array_Tuple out = {real, imag};
    return out;
}

//np.hamming
Array_Tuple hamming(int M){ //! Returns a calloc ptr
    double* output;
    output = (double*) calloc(M, sizeof(double));

    for (int i = 0; i < M; i++)
    {
        output[i] = 0.54 - 0.46 * cos((2 * M_PI * i)/(M - 1));
    }
    Array_Tuple out = {output, M};
    return out;
}

//np.convolve: Circular Discrete Convolution --> https://en.wikipedia.org/wiki/Convolution
Complex_Array_Tuple convolve(Complex_Array_Tuple a, Array_Tuple v){ //! Returns a calloc ptr
    int N = a.real.length + v.length - 1;
    double* output_real;
    output_real = (double*) calloc(N, sizeof(double));
    double* output_imaginary;
    output_imaginary = (double*) calloc(N, sizeof(double));

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
        
        output_real[n] = sum_real;
        output_imaginary[n] = sum_imaj;
    }

    Array_Tuple real_out = {output_real, N };
    Array_Tuple imaj_out = {output_imaginary, N};
    Complex_Array_Tuple out = {real_out, imaj_out};
    return out;
}

Complex_Array_Tuple convolveSame(Complex_Array_Tuple a, Array_Tuple v){ //! Returns a calloc ptr
    int N = a.real.length + v.length - 1;
    double* output_real;
    output_real = (double*) calloc(N, sizeof(double));
    double* output_imaginary;
    output_imaginary = (double*) calloc(N, sizeof(double));

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
        output_real[n] = sum_real;
        output_imaginary[n] = sum_imaj;
    }
    // remove difference/2 from both sides, this is what makes it different than np.convolve(mode='full')
    int output_length = a.real.length > v.length ? a.real.length : v.length;
    double* shortened_real;
    shortened_real = (double*) calloc(N, sizeof(double));
    double* shortened_imaj;
    shortened_imaj = (double*) calloc(N, sizeof(double));
    
    int difference = abs(N - output_length);
    double remove = (double)difference / 2.0;

    if (isInteger(remove)){
        for (int i = 0; i < output_length; i++)
        {
            shortened_real[i] = output_real[(int)remove + i];
            shortened_imaj[i] = output_imaginary[(int)remove + i];
        }
    } else {
        for (int i = 0; i < output_length; i++)
        {
            shortened_real[i] = output_real[(int)floor(remove) + i];
            shortened_imaj[i] = output_imaginary[(int)floor(remove) + i];
        }
    }

    

    Array_Tuple real_out = {shortened_real, output_length };
    Array_Tuple imaj_out = {shortened_imaj, output_length};
    Complex_Array_Tuple out = {real_out, imaj_out};
    free(output_real);
    free(output_imaginary);
    return out;
}
Array_Tuple firwin(int M, double cutoff){ //! Returns a calloc ptr

    double* output = (double*)calloc(M+1, sizeof(double));
    //scipy.signal.firwin
    for (int i = 0; i < M; i++)
    {
        double window = (0.54 - 0.46 * cos(2 * M_PI * i / M));
        
        double lowpassfilter = cutoff; // limit of sin(cutoff * M_PI * (i - M/2)) / (M_PI * (i - M/2)) as (i-M/2) --> 0 = cutoff
        if (i - M/2 != 0){
            lowpassfilter = sin(cutoff * M_PI * (i - M/2)) / (M_PI * (i - M/2));
        }
        output[i] = window * lowpassfilter;
    }
    
    Array_Tuple out = {output, M};
    return out;
}
//scipy.fftconvolve --> Ask isaac: almost same results as above convolve, just slightly off. By like e^-16. Worth rewriting?
// Complex_Array_Tuple fftconvolve(Complex_Array_Tuple a, Array_Tuple v){ //! Returns a calloc ptr

// }
// same as scipy.resample_poly which upsamples an array of data, can also specify a down sample
Complex_Array_Tuple resample_poly(Complex_Array_Tuple a, int up, int down){
    // upsample by factor of up
    double* upsampled_real = (double*)calloc(a.real.length * up, sizeof(double));
    double* upsampled_imaj = (double*)calloc(a.real.length * up, sizeof(double));
    int index = 0;
    // upsample the input
    for (int i = 0; i < a.real.length * up; i++)
    {
        if (i % up == 0){
            upsampled_real[i] = a.real.array[index];
            upsampled_imaj[i] = a.imaginary.array[index];
            index ++;
        }
        else {
            upsampled_real[i] = 0;
            upsampled_imaj[i] = 0;
        }
    }
    int greaterUpDown = (up >= down ? up : down);
    int filter_length = 10 * greaterUpDown;

    Array_Tuple filterCoeff = firwin(filter_length, 1.0 / (double)greaterUpDown);
    Array_Tuple up_real = {upsampled_real, a.real.length * up};
    Array_Tuple up_imaj = {upsampled_imaj, a.real.length * up};
    Complex_Array_Tuple up_complex = {up_real, up_imaj};
    

    Complex_Array_Tuple smoothed = convolveSame(up_complex, filterCoeff);
    //printArray("smoothed", smoothed.real.array, smoothed.real.length);
    int down_length = smoothed.real.length / down;
    double* down_real = (double*)calloc(down_length, sizeof(double)); //! can make these static
    double* down_imaj = (double*)calloc(down_length, sizeof(double));
    index = 0;
    double max_abs_value_sqr = -99;
    for (int i = 0; i < smoothed.real.length; i++)
    {
        if (i % down == 0){
            double x = smoothed.real.array[i];
            double y = smoothed.imaginary.array[i];
            double abs_value_sqr = x*x + y*y;
            down_real[index] = x;
            down_imaj[index] = y;
            index ++;
            if (abs_value_sqr > max_abs_value_sqr){
                max_abs_value_sqr = abs_value_sqr;
            }
        }
    }

    // normalize output then scale to input scale
    double* out_real = (double*)calloc(down_length, sizeof(double));
    double* out_imaj = (double*)calloc(down_length, sizeof(double));
    double max_abs_input = maxAbsoluteValue(a);
    for (int i = 0; i < down_length; i++)
    {
        out_real[i] = down_real[i] / sqrt(max_abs_value_sqr) * max_abs_input;
        out_imaj[i] = down_imaj[i] / sqrt(max_abs_value_sqr) * max_abs_input;
    }
    
    Array_Tuple real_out = {out_real, down_length};
    Array_Tuple imaj_out = {out_imaj, down_length};
    Complex_Array_Tuple complex_out = {real_out, imaj_out};

    // free all temp arrays
    free(upsampled_real);
    free(upsampled_imaj);
    free(down_real);
    free(down_imaj);
    free(filterCoeff.array);
    free(smoothed.real.array);
    free(smoothed.imaginary.array);

    return complex_out;
}

void exportArray(Array_Tuple input, char filename[]){
    FILE *fpt;
    fpt = fopen(filename, "w+");
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
    fpt = fopen(filename, "w+");
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

