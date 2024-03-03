#ifndef STANDARDARRAY_H
#define STANDARDARRAY_H
//standardArray.h
// Rylan Paul
#include "src/standardArray.c"
// checks if double is an integer
bool isInteger(double val);
// Array Tuple that keeps track of the length of the array and a pointer to the array
typedef struct Array_Tuple Array_Tuple;
typedef struct Complex_Array_Tuple Complex_Array_Tuple;
// Prints an array printArray("Debug Text", ptr to int[], length) (prints out) --> Debug Text: [1, 0, ...]
void printArray(char *string, Array_Tuple arr);
void printComplexArray(char *string, Complex_Array_Tuple arr);
// Pass in an Array_Tuple and it will free the memory from that pointer, this only needs to be done on pointers that use the calloc() and malloc() functions
void freeArrayMemory(Array_Tuple array);
void freeComplexArrayMemory(Complex_Array_Tuple array);
// Creates an Array_Tuple from a known array {1,1,1,0,...} and its known length
Array_Tuple defineArray(double array[], int length);
// Creates an Complex_Array_Tuple length {0,0,0,0,0} and its known length
Complex_Array_Tuple zerosComplex(int length);
// Creates an Array_Tuple given a length with each element set to the given value
Array_Tuple valueArray(int length, int value);
// Creates an Complex_Array_Tuple length {0,0,0,0,0} and its known length
Array_Tuple zerosArray(int length);
// Creates an Complex_Array_Tuple of complex conjugates of the input -> np.conj
Complex_Array_Tuple getConj(Complex_Array_Tuple input);
// gets the average of an array
double meanArray(double* array, int length);
// gets the average of an array
double meanArrayTuple(Array_Tuple array);
// Generates a random pointer to an array
double * randomArray(int max_exclusive, int length);
Array_Tuple append_array(Array_Tuple a, Array_Tuple b);
// reverses the order of an array
Array_Tuple flip(Array_Tuple a);
// Raises e to the element of every complex number in an array --> np.exp()
Complex_Array_Tuple expComplexArray(Complex_Array_Tuple array);
// Raises e to the element of every complex number in an array --> np.exp(), the input is 0+xj, where x is the input
Complex_Array_Tuple exp_imaginaryArray(Array_Tuple array);
//returns the maximum abs of complex array
double maxAbsoluteValue(Complex_Array_Tuple a);
// returns a random number with a normal distribution
double rand_norm(double mu, double sigma);
//np.arange creates an array from start to end (exclusive) by a step
Array_Tuple arange(double start, double end, double step);
//np.linspace creates an array from start to end with x number of points defined by length
Array_Tuple linspace(double start, double end, int length);
Array_Tuple addArrays(Array_Tuple a, Array_Tuple b);
Array_Tuple subtractArrays(Array_Tuple a, Array_Tuple b);
Array_Tuple subtractDoubleFromArray(Array_Tuple a, double b);
Array_Tuple divideDoubleFromArray(Array_Tuple a, double b);
Array_Tuple multiplyDoubleFromArray(Array_Tuple a, double b);
Array_Tuple multiplyArrays(Array_Tuple a, Array_Tuple b);
Array_Tuple divideArrays(Array_Tuple a, Array_Tuple b);
// multiplies two complex arrays
Complex_Array_Tuple multiplyComplexArrays(Complex_Array_Tuple x, Complex_Array_Tuple y);
// add two complex arrays together
Complex_Array_Tuple addComplexArrays(Complex_Array_Tuple a, Complex_Array_Tuple b);
// subtract two complex arrays
Complex_Array_Tuple subtractComplexArrays(Complex_Array_Tuple a, Complex_Array_Tuple b);
//same as np.sinc, does sin(pi * x) / (pi * x), but lim x-> 0 = 1, so if x=0 return 1
Array_Tuple sinc(Array_Tuple input);
//np.sum Adds them all together
double sumArray(Array_Tuple input);
Complex_Array_Tuple everyOtherElement(Complex_Array_Tuple array, int offset);
Array_Tuple absComplexArray(Complex_Array_Tuple array);
int argMax(Array_Tuple input);
int argComplexMax(Complex_Array_Tuple input);
//np.convolve: Circular Discrete Convolution --> https://en.wikipedia.org/wiki/Convolution
Complex_Array_Tuple convolve(Complex_Array_Tuple a, Array_Tuple v);
Complex_Array_Tuple convolveSame(Complex_Array_Tuple a, Array_Tuple v);
void exportArray(Array_Tuple input, char filename[]);
void exportComplexArray(Complex_Array_Tuple input, char filename[]);
#endif