/*
main.c
IrishSAT CLOVER SDR
February 2024
Rylan Paul
Notes:
    • link flag -lfftw3, run: brew install fftw // to install library for mac
            • for windows: http://www.fftw.org/install/windows.html
    • gcc main.c -o main_exec -lm -lfftw3
*/
#include <stdio.h>
#include <math.h>
#include <complex.h>
// Our libraries --------------------------------
#include "lib/irishsat_comms_lib.h" 
//-----------------------------------------------
// * Function Prototypes ------------------------
Complex_Array_Tuple pulse_Shaping_Main(Array_Tuple pulse_train, int sps, double fs, char pulse_shape[], double alpha, int L);

void DisplayOutput(Array_Tuple data, Array_Tuple demod_bits);
// ----------------------------------------------

int main(){
    // *DECLARE VARIABLES --------------
    const int data_length = 256;
    double fs = 2.38e9; // carrier frequency not 2.4e9!!!
    int L = 8; // pulse shaping filter length (Is this the ideal length?)
    char pulse_shape[] = "rrc"; // type of pulse shaping, also "rect"
    char scheme[] = "BPSK"; // Modulation scheme 'OOK', 'QPSK', or 'QAM'
    double alpha = 0.5; // roll-off factor of the RRC pulse-shaping filter
    int sps = 8; // samples per symbol, equal to oversampling factor
    // *Preamble, Data, and Matched Filter Coefficient Generation--------------
    Array_Tuple preamble = defineArray((double[]){0,1,0,0,0,0,1,1,0,0,0,1,0,1,0,0,1,1,1,1,0,1,0,0,0,1,1,1,0,0,1,0,0,1,0,1,1,0,1,1,1,0,1,1,0,0,1,1,0,1,0,1,0,1,1,1,1,1,1,0}, 60);// optimal periodic binary code for N = 63 https://ntrs.nasa.gov/citations/19800017860
    exportArray(preamble, "preamble.txt");
    Array_Tuple data = {randomArray(2,data_length), data_length}; // this points to the random array generated;  creates 256 random bits
    exportArray(data, "data.txt");
    Array_Tuple CRC_key = defineArray((double[]){1,0,0,1,1,0,0,0,0,1,1,1}, 12);// Best CRC polynomials: https://users.ece.cmu.edu/~koopman/crc/
    Array_Tuple data_encoded = CRC_encodeData(data, CRC_key);
    Array_Tuple bits = append_array(preamble, data_encoded);
    exportArray(bits, "bits.txt");
    
    // *Generation of Pulse Train
    Array_Tuple pulse_train = pulsetrain(bits, sps);
    exportArray(pulse_train, "pulsetrain.txt");
    freeArrayMemory(bits);


    // *Pulse Shaping
    Complex_Array_Tuple complexTestpacket = pulse_Shaping_Main(pulse_train, sps, fs, pulse_shape, alpha, L);
    exportComplexArray(complexTestpacket, "pulseshaping.txt");
    freeArrayMemory(pulse_train);

    //*############## Transmit ################
    // socket_complex_array(complexTestpacket)
    complexArrayToCharArray(complexTestpacket);
    //*############## complexTestpacket #######

    //printComplexArray("Transmitting", complexTestpacket);

    freeArrayMemory(preamble);
    freeArrayMemory(CRC_key);
    freeArrayMemory(data);
    freeArrayMemory(data_encoded);

    return 0;
}

// Functions --> Main Processes
Complex_Array_Tuple pulse_Shaping_Main(Array_Tuple pulse_train, int sps, double fs, char pulse_shape[], double alpha, int L){
    Array_Tuple testpacket_pulseshape = pulse_shaping(pulse_train, sps, fs, pulse_shape, alpha, L);
    //printArray("testpacket_pulseshape", testpacket_pulseshape.array, testpacket_pulseshape.length);
    Array_Tuple imajtestpacket = zerosArray(testpacket_pulseshape.length);
    Complex_Array_Tuple complexTestpacket = {testpacket_pulseshape, imajtestpacket};
    return complexTestpacket;
}

