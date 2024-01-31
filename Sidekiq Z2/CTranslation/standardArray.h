#include <string.h>
#include <stdbool.h>
#include <stdlib.h>
#include <stdio.h>
#define MAX_ARRAY_LENGTH 1024
// Prints an array printArray("Debug Text", ptr to int[], length) (prints out) --> Debug Text: [1, 0, ...]
void printArray(char *string, int* arr, int length){
    for (int i = 0; i < strlen(string); i++) {
        printf("%c", string[i]);
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


// Array Tuple that keeps track of the length of the array and a pointer to the array
struct Array_Tuple {
    int* array;
    int length;
};
// Pass in an Array_Tuple and it will free the memory from that pointer, this only needs to be done on pointers that use the calloc() and malloc() functions
void freeArrayMemory(struct Array_Tuple array){
    free(array.array);
    return;
}
// Takes in two arrays and appends them in the order of a then b. User can also opt to free the memory after the are appended
struct Array_Tuple append_array(struct Array_Tuple a, struct Array_Tuple b, bool free_memory_a, bool free_memory_b){
    int* ptr;
    ptr = (int*)calloc(a.length + b.length, sizeof(int));
    for (int i = 0; i < a.length; i++)
    {
        ptr[i] = a.array[i];
    }
    for (int i = 0; i < b.length; i++)
    {
        ptr[i+a.length] = b.array[i];
    }
    if (free_memory_a == true){
        freeArrayMemory(a);
    }
    if (free_memory_b == true){
        freeArrayMemory(b);
    }
    struct Array_Tuple out = {ptr, a.length + b.length};
    return out;
}
// Creates an Array_Tuple from a known array {1,1,1,0,...} and its known length
struct Array_Tuple defineArray(int array[], int length){
    
    int* ptr;
    ptr = (int*)calloc(length, sizeof(int));
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

// Generates a random pointer to an array
int * randomArray(int max_exclusive, int length){ // returns a pointer to an integer array
    static int arr[MAX_ARRAY_LENGTH]; // Assuming maximum length of the array
    for (int i=0; i < length; i++){
        arr[i] = rand() % max_exclusive;
    }
    return arr; // to use this: int *array; array= randomArray(2,size);
}