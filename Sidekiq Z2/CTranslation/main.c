/*
main.c
IrishSAT CLOVER SDR
February 2024
Rylan Paul

Notes:
- For memory leaks, every function that returns a pointer using malloc() or calloc(), mark with //! Returns a calloc ptr
            - For those using vscode, install the better comments extension, then you will see comments starting with --> ! <-- in bright red


Funtcions to add:
- np.sinc
- np.hamming
- np.sum
- np.convolve
- signal.resample_poly
- signal.upfirdn
- np.mean
- signal.fftconvolve
- demod.symbol_demod
- CRC.CRCcheck (double check this works)
- pulse_shaping.pulse_shaping

*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdbool.h>

// Our libraries --------------------------------
#include "standardArray.h" //#include <fftw3.h> //link flag -lfftw3, run: brew install fftw // to install library
#include "CRC.h"
#include "iq_imbalance.h"
//-----------------------------------------------

void encodepacket();

int main(){
    encodepacket();
    return 0;
}

void encodepacket(){
    // Declare Variables
    const int data_length = 10;
    double fs = 2.4e9; // carrier frequency 
    double Ts = 1.0 / fs; // sampling period in seconds
    double f0 = 0.0; // homodyne (0 HZ IF)
    int M = 8; // oversampling factor
    double L = 8; // pulse shaping filter length (Is this the ideal length?)
    char pulse_shape[] = "rrc"; // type of pulse shaping, also "rect"
    char scheme[] = "BPSK"; // Modulation scheme 'OOK', 'QPSK', or 'QAM'
    double alpha = 0.5; // roll-off factor of the RRC pulse-shaping filter
    Array_Tuple preamble = defineArray((double[]){0,1,0,0,0,0,1,1,0,0,0,1,0,1,0,0,1,1,1,1,0,1,0,0,0,1,1,1,0,0,1,0,0,1,0,1,1,0,1,1,1,0,1,1,0,0,1,1,0,1,0,1,0,1,1,1,1,1,1,0}, 60);// optimal periodic binary code for N = 63 https://ntrs.nasa.gov/citations/19800017860
    Array_Tuple data = {randomArray(2,data_length), data_length}; // this points to the random array generated;  creates 256 random bits

    Array_Tuple CRC_key = defineArray((double[]){1,0,0,1,1,0,0,0,0,1,1,1}, 12);// Best CRC polynomials: https://users.ece.cmu.edu/~koopman/crc/
    Array_Tuple data_encoded = CRC_encodeData(data, CRC_key);
    freeArrayMemory(CRC_key); // if its not needed annymore
    Array_Tuple bits = append_array(preamble, data_encoded);

    int sps = M; // samples per symbol, equal to oversampling factor

    Array_Tuple matched_filter_coef = flip(preamble);

    Array_Tuple pulse_train = pulsetrain(bits, sps);
    
    Array_Tuple testpacket_pulseshape = pulse_shape(pulse_train, sps, fs, pulse_shape, alpha, L);
    double* imaj_test_packet = (double*)calloc(testpacket_pulseshape.length, sizeof(double));
    Array_Tuple imajtestpacket = {imaj_test_packet, testpacket_pulseshape.length};
    Complex_Array_Tuple complexTestpacket = {testpacket_pulseshape, imaj_test_packet};

    // Noise simulation ----------------------------
    Complex_Array_Tuple testpacket_noise = generateNoise(complexTestpacket);
    //printArray("T", testpacket_noise.array, testpacket_noise.length);
    // ---------------------------------------------

    ///////////////////////////////////
    // Add fractional delay

    // Create and apply fractional delay filter
    double delay = 0.4; // fractional delay, in samples
    int N = 21; // number of taps
    Array_Tuple n = arange(-N / 2, N / 2, 1); // ...-3,-2,-1,0,1,2,3...
    Array_Tuple temp = subtractDoubleFromArray(n, delay);
    Array_Tuple h = sinc(temp); // calc filter taps
    Array_Tuple ham = hamming(N);
    Array_Tuple temp2 = multiplyArrays(h, ham); // window the filter to make sure it decays to 0 on both sides
    Array_Tuple new_h = divideDoubleFromArray(temp2, sumArray(temp2)); // normalize to get unity gain, we don't want to change the amplitude/power
    Complex_Array_Tuple testpacket_filtered = convolve(testpacket_noise, new_h); // apply filter

    // apply a freq offset
    fs = 2.45e9; // arbitrary UHF frequency
    int fo = 61250; // Simulated frequency offset
    Ts = (double) 1 / fs; // calc sample period
    Array_Tuple t = arange(0, Ts*testpacket_filtered.real.length, Ts); // create time vector
    Array_Tuple t_new = multiplyDoubleFromArray(t, fo * M_PI * 2);
    Complex_Array_Tuple exp_t = exp_imaginaryArray(t_new);
    Complex_Array_Tuple testpacket_freq_shift = multiplyComplexArrays(testpacket_filtered, exp_t); // perform freq shift

    Complex_Array_Tuple samples_interpolated = resample_poly(testpacket_freq_shift, 16, 1);

    double mu = 0; // initial estimate of phase of sample
    Complex_Array_Tuple out = zerosComplex(testpacket_freq_shift.real.length + 10);
    Complex_Array_Tuple out_rail = zerosComplex(testpacket_freq_shift.real.length + 10);// stores values, each iteration we need the previous 2 values plus current value
    int i_in = 0; // input samples index
    int i_out = 2; // output index (let first two outputs be 0)
    //sps = 8; SPS is already 8

    while (i_out < testpacket_freq_shift.real.length && (i_in+16) < testpacket_freq_shift.real.length){
        //out[i_out] = samples_interpolated[i_in*16 + int(mu*16)]
        int value = (int)(i_in*16) + floor(mu*16);
        out.real.array[i_out] = samples_interpolated.real.array[value];
        out.imaginary.array[i_out] = samples_interpolated.imaginary.array[value];

        //out_rail[i_out] = int(np.real(out[i_out]) > 0) + 1j*int(np.imag(out[i_out]) > 0)
        out_rail.real.array[i_out] = (out.real.array[i_out] > 0);
        out_rail.imaginary.array[i_out] = (out.imaginary.array[i_out] > 0);
        
        //double x = (out_rail[i_out] - out_rail[i_out-2]) * out_conj[i_out-1];
        // (a+bi)*(c+di) = (ac-bd)+(ad+bc)i;
        Complex_Array_Tuple out_conj = conj(out);
        double x_sub_real = out_rail.real.array[i_out]- out_rail.real.array[i_out-2];
        double x_sub_imaj = out_rail.imaginary.array[i_out]- out_rail.imaginary.array[i_out-2];
        double x_conj_real = out_conj.real.array[i_out-1];
        double x_conj_imaj = out_conj.imaginary.array[i_out-1];
        
        double x_real = x_sub_real*x_conj_real - x_sub_imaj*x_conj_imaj;
        ////double x_imaj = x_sub_real*x_conj_imaj + x_sub_imaj*x_conj_real;

        //y = (out[i_out] - out[i_out-2]) * np.conj(out_rail[i_out-1])
        Complex_Array_Tuple out_rail_conj = conj(out_rail);
        double y_sub_real = out.real.array[i_out]- out.real.array[i_out-2];
        double y_sub_imaj = out.imaginary.array[i_out]- out.imaginary.array[i_out-2];
        double y_conj_real = out_rail_conj.real.array[i_out-1];
        double y_conj_imaj = out_rail_conj.imaginary.array[i_out-1];
        
        double y_real = y_sub_real*y_conj_real - y_sub_imaj*y_conj_imaj;
        ////double y_imaj = y_sub_real*y_conj_imaj + y_sub_imaj*y_conj_real;

        //mm_val = np.real(y - x)
        double mm_val = y_real - x_real;
        mu += sps + 0.3 * mm_val;
        i_in += (int)floor(mu);// round down to nearest int since we are using it as an index
        mu = mu - floor(mu); // remove the integer part of mu
        i_out++; // increment output index
    }
    //out = out[2:i_out] // remove the first two, and anything after i_out (that was never filled out)
    double* out_new_real_ptr = (double*)calloc(i_out-2, sizeof(double));;
    double* out_new_imaj_ptr = (double*)calloc(i_out-2, sizeof(double));
    for (int i = 2; i < i_out; i++)
    {
        out_new_real_ptr[i-2] = out.real.array[i];
        out_new_imaj_ptr[i-2] = out.imaginary.array[i];
    }
    Array_Tuple out_new_real = {out_new_real_ptr, i_out-2};
    Array_Tuple out_new_imaj = {out_new_imaj_ptr, i_out-2};
    Complex_Array_Tuple testpacket = {out_new_real, out_new_imaj}; // only include this line if you want to connect this code snippet with the Costas Loop later on

    //Coarse Frequency Detection and Correction
    //fft_samples = testpacket**2
    Complex_Array_Tuple fft_samples = multiplyComplexArrays(testpacket, testpacket);

    return;
}
