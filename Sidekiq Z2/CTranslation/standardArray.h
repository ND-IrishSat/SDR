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
struct Array_Tuple defineArray(double array[], int length){
    
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
double * randomArray(int max_exclusive, int length){ // returns a pointer to an integer array
    static double arr[MAX_ARRAY_LENGTH]; // Assuming maximum length of the array
    for (int i=0; i < length; i++){
        arr[i] = (double)(rand() % max_exclusive);
    }
    return arr; // to use this: int *array; array= randomArray(2,size);
}

struct Array_Tuple append_array(struct Array_Tuple a, struct Array_Tuple b){
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

struct Array_Tuple flip(struct Array_Tuple a){
    double* ptr;
    ptr = (double*)calloc(a.length, sizeof(double));
    for (int i = 0; i < a.length; i++)
    {
        ptr[i] = a.array[a.length-1-i];
    }
    struct Array_Tuple t = {ptr, a.length};
    return t; //once done, must free memory of ptr
}
