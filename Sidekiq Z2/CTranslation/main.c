#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include "standardArray.h"
#include "CRC.h"
#include "iq_imbalance.h"

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

    struct Array_Tuple preamble = defineArray((double[]){0,1,0,0,0,0,1,1,0,0,0,1,0,1,0,0,1,1,1,1,0,1,0,0,0,1,1,1,0,0,1,0,0,1,0,1,1,0,1,1,1,0,1,1,0,0,1,1,0,1,0,1,0,1,1,1,1,1,1,0}, 60);// optimal periodic binary code for N = 63 https://ntrs.nasa.gov/citations/19800017860
    struct Array_Tuple data = {randomArray(2,data_length), data_length}; // this points to the random array generated;  creates 256 random bits

    struct Array_Tuple CRC_key = defineArray((double[]){1,0,0,1,1,0,0,0,0,1,1,1}, 12);// Best CRC polynomials: https://users.ece.cmu.edu/~koopman/crc/
    // Working Demo Code Here ----------------------------------------------------
    // calloc(length, size) is used in defineArray() and for the output variable of CRC_encodeData()

    struct Array_Tuple demoA = defineArray((double[]){1.4, 5.3, 6.3, 1.2, -3.0, 9.3, 5.4, 1.2, -9.1, 8.2, 1.0, 4.5, -2.3}, 13);
    struct Array_Tuple out = averages(demoA, 10);
    printArray("Out", out.array, out.length);
    freeArrayMemory(demoA);
    freeArrayMemory(out); // when done with the array, need to free its memory
    // ---------------------------------------------------------------------------
    //struct Array_Length_Tuple tuple = CRC_encodeData(data.array, CRC_key.array);
    return 0;
}