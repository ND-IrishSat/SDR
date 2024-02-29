/*
main.c
IrishSAT CLOVER SDR
February 2024
Rylan Paul

Notes:
- For memory leaks, every function that returns a pointer using malloc() or calloc(), mark with //! Returns a calloc ptr
            - For vscode, install the better comments extension, then you will see comments starting with ! in bright red
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdbool.h>
#include <complex.h>
// Our libraries --------------------------------
#include "lib/standardArray.h" //#include <fftw3.h> //link flag -lfftw3, run: brew install fftw // to install library
#include "lib/CRC.h"
#include "lib/iq_imbalance.h"
#include "lib/src/pulse_shaping.c"
#include "lib/src/symbol_demod_QAM.c"
//-----------------------------------------------


// * Function Prototypes ------------------------
Complex_Array_Tuple pulse_Shaping_Main(Array_Tuple pulse_train, int sps, double fs, char pulse_shape[], double alpha, int L);
Complex_Array_Tuple generateNoise(Complex_Array_Tuple testpacket);
Complex_Array_Tuple simulation_Channel(Complex_Array_Tuple testpacket_noise, double fs, double Ts);
Complex_Array_Tuple clock_Recovery(Complex_Array_Tuple testpacket, int sps);
Complex_Array_Tuple coarse_Frequency_Correction(Complex_Array_Tuple testpacket, double fs);
Complex_Array_Tuple fine_Frequency_Correction(Complex_Array_Tuple new_testpacket, double fs);
Complex_Array_Tuple frame_Sync(Complex_Array_Tuple testpacket, Array_Tuple matchedFilterCoef, Array_Tuple preamble, Array_Tuple data_encoded);
Array_Tuple demodulation(Complex_Array_Tuple recoveredData, char scheme[], Array_Tuple preamble, Array_Tuple CRC_key);
void DisplayOutput(Array_Tuple data, Array_Tuple demod_bits);
// ----------------------------------------------


int main(){
    // *DECLARE VARIABLES --------------
    const int data_length = 30;
    double fs = 2.4e9; // carrier frequency 
    double Ts = 1.0 / fs; // sampling period in seconds
    double f0 = 0.0; // homodyne (0 HZ IF)
    int M = 8; // oversampling factor
    int L = 8; // pulse shaping filter length (Is this the ideal length?)
    char pulse_shape[] = "rrc"; // type of pulse shaping, also "rect"
    char scheme[] = "BPSK"; // Modulation scheme 'OOK', 'QPSK', or 'QAM'
    double alpha = 0.5; // roll-off factor of the RRC pulse-shaping filter
    int sps = M; // samples per symbol, equal to oversampling factor

    // *Preamble, Data, and Matched Filter Coefficient Generation--------------
    Array_Tuple preamble = defineArray((double[]){0,1,0,0,0,0,1,1,0,0,0,1,0,1,0,0,1,1,1,1,0,1,0,0,0,1,1,1,0,0,1,0,0,1,0,1,1,0,1,1,1,0,1,1,0,0,1,1,0,1,0,1,0,1,1,1,1,1,1,0}, 60);// optimal periodic binary code for N = 63 https://ntrs.nasa.gov/citations/19800017860
    exportArray(preamble, "preamble.txt");
    Array_Tuple data = {randomArray(2,data_length), data_length}; // this points to the random array generated;  creates 256 random bits
    exportArray(data, "data.txt");
    Array_Tuple CRC_key = defineArray((double[]){1,0,0,1,1,0,0,0,0,1,1,1}, 12);// Best CRC polynomials: https://users.ece.cmu.edu/~koopman/crc/
    Array_Tuple data_encoded = CRC_encodeData(data, CRC_key);
    Array_Tuple bits = append_array(preamble, data_encoded);
    exportArray(bits, "bits.txt");
    Array_Tuple matched_filter_coef = flip(preamble);


    // *Generation of Pulse Train
    Array_Tuple pulse_train = pulsetrain(bits, sps);
    exportArray(pulse_train, "pulsetrain.txt");

    // *Pulse Shaping
    Complex_Array_Tuple complexTestpacket = pulse_Shaping_Main(pulse_train, sps, fs, pulse_shape, alpha, L);
    exportComplexArray(complexTestpacket, "pulseshaping.txt");

    //#############################################################################################
    // *----TRANSMISSION and Noise----
    Complex_Array_Tuple testpacket_noise = generateNoise(complexTestpacket);
    exportComplexArray(testpacket_noise, "noise.txt");
    //#############################################################################################

    // *Simulation of Channel
    fs = 2.45e9; // arbitrary UHF frequency
    Complex_Array_Tuple testpacket_freq_shift = simulation_Channel(testpacket_noise, fs, Ts);
    exportComplexArray(testpacket_freq_shift, "testpacketfreqshift.txt");

    // *Clock Recovery
    Complex_Array_Tuple testpacket = clock_Recovery(testpacket_freq_shift, sps);
    exportComplexArray(testpacket, "clockRecovery.txt");

    // *Coarse Frequency Correction
    Complex_Array_Tuple new_testpacket = coarse_Frequency_Correction(testpacket, fs);
    exportComplexArray(new_testpacket, "coarseFrequencyCorrection.txt");

    // *Fine Frequency Correction
    Complex_Array_Tuple costas_out = fine_Frequency_Correction(new_testpacket, fs);
    exportComplexArray(costas_out, "costasout.txt");

    // *IQ Imbalance Correction
    int mean_period = 100;
    Complex_Array_Tuple testpacket_IQ_Imbalance_Correct = IQImbalanceCorrect(costas_out, mean_period); //TODO clean function
    exportComplexArray(testpacket_IQ_Imbalance_Correct, "iqimbalanceout.txt");
    // *Frame Sync
    Complex_Array_Tuple recoveredData = frame_Sync(testpacket_IQ_Imbalance_Correct, matched_filter_coef, preamble, data_encoded); //TODO clean function
    // exports inside function
    // *Demodulation
    Array_Tuple demod_bits = demodulation(recoveredData, scheme, preamble, CRC_key); //TODO clean function
    exportArray(demod_bits, "demodbits.txt");
    // *Display Output
    DisplayOutput(data, demod_bits);

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
Complex_Array_Tuple generateNoise(Complex_Array_Tuple testpacket){
    // vars
    double mean = 0;
    double std_dev = 0.1;
    int num_samples = testpacket.real.length;

    double phase_noise_strength = 0.1;
    double noise_power = 10;

    // creating general normal distribution random arrays for real and complex
    Complex_Array_Tuple out = zerosComplex(num_samples);
    for (int i = 0; i < num_samples; i++)
    {
        double real = testpacket.real.array[i];
        double x = rand_norm(mean, std_dev);
        x /= sqrt(2.0);
        x /= sqrt(noise_power);
        real += x;

        double imaj = testpacket.real.array[i];
        double y = rand_norm(mean, std_dev);
        y /= sqrt(2.0);
        y /= sqrt(noise_power);
        imaj += y;

        double complex phase_noise = cexp(I * rand_norm(mean, std_dev) * phase_noise_strength);
        double complex packet_complex = real + I * imaj;
        double complex out_packet = packet_complex * phase_noise;

        out.real.array[i] = creal(out_packet);
        out.imaginary.array[i] = cimag(out_packet);
    }
    return out;
}
Complex_Array_Tuple simulation_Channel(Complex_Array_Tuple testpacket_noise, double fs, double Ts){
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
    Complex_Array_Tuple testpacket = convolve(testpacket_noise, new_h); // apply filter

    // apply a freq offset
    int fo = 61250; // Simulated frequency offset
    Ts = (double) 1 / fs; // calc sample period
    Array_Tuple t = arange(0, Ts*testpacket.real.length, Ts); // create time vector
    Array_Tuple t_new = multiplyDoubleFromArray(t, fo * M_PI * 2);
    Complex_Array_Tuple exp_t = exp_imaginaryArray(t_new);
    Complex_Array_Tuple testpacket_freq_shift = multiplyComplexArrays(testpacket, exp_t); // perform freq shift

    freeArrayMemory(n);
    freeArrayMemory(temp);
    freeArrayMemory(h);
    freeArrayMemory(ham);
    freeArrayMemory(temp2);
    freeArrayMemory(new_h);
    freeComplexArrayMemory(testpacket);
    freeArrayMemory(t);
    freeArrayMemory(t_new);
    freeComplexArrayMemory(exp_t);

    return testpacket_freq_shift;
    
}
Complex_Array_Tuple clock_Recovery(Complex_Array_Tuple testpacket, int sps){
    Complex_Array_Tuple samples_interpolated = resample_poly(testpacket, 16, 1);
    exportComplexArray(samples_interpolated, "samplesinterpolated.txt");
    double mu = 0; // initial estimate of phase of sample
    Complex_Array_Tuple out = zerosComplex(testpacket.real.length + 10);
    Complex_Array_Tuple out_rail = zerosComplex(testpacket.real.length + 10);// stores values, each iteration we need the previous 2 values plus current value
    int i_in = 0; // input samples index
    int i_out = 2; // output index (let first two outputs be 0)
    
    //int sps = 8; SPS is already 8

    while (i_out < testpacket.real.length && (i_in+16) < testpacket.real.length){
        //out[i_out] = samples_interpolated[i_in*16 + int(mu*16)]
        int value = (int)(i_in*16) + (int)floor(mu*16);
        out.real.array[i_out] = samples_interpolated.real.array[value];
        out.imaginary.array[i_out] = samples_interpolated.imaginary.array[value];

        //out_rail[i_out] = int(np.real(out[i_out]) > 0) + 1j*int(np.imag(out[i_out]) > 0)
        out_rail.real.array[i_out] = (int)(out.real.array[i_out] > 0);
        out_rail.imaginary.array[i_out] = (int)(out.imaginary.array[i_out] > 0);
        
        //double x = (out_rail[i_out] - out_rail[i_out-2]) * out_conj[i_out-1];
        // (a+bi)*(c+di) = (ac-bd)+(ad+bc)i;
        double complex out_num = out.real.array[i_out-1] + I * out.imaginary.array[i_out-1];
        double complex out_conj_complex = conj(out_num);
        double complex out_rail_num = out_rail.real.array[i_out-1] + I * out_rail.imaginary.array[i_out-1];
        double complex out_rail_conj_complex = conj(out_rail_num);

        double complex complex_out_rail = out_rail.real.array[i_out] + I * out_rail.imaginary.array[i_out];
        double complex complex_out_rail_2 = out_rail.real.array[i_out-2] + I * out_rail.imaginary.array[i_out-2];

        double complex complex_out = out.real.array[i_out] + I * out.imaginary.array[i_out];
        double complex complex_out_2 = out.real.array[i_out-2] + I * out.imaginary.array[i_out-2];

        double complex x = (complex_out_rail - complex_out_rail_2) * out_conj_complex;
        double complex y = (complex_out - complex_out_2) * out_rail_conj_complex;

        //mm_val = np.real(y - x)
        double mm_val = creal(y-x);
        mu += sps + 0.3 * mm_val;
        i_in += (int)floor(mu);// round down to nearest int since we are using it as an index
        mu = mu - floor(mu); // remove the integer part of mu
        i_out++; // increment output index
    }
    //out = out[2:i_out] // remove the first two, and anything after i_out (that was never filled out)
    Complex_Array_Tuple out_tuple = zerosComplex(i_out-2); // only include this line if you want to connect this code snippet with the Costas Loop later on
    for (int i = 0; i < (i_out-2); i++)
    {
        out_tuple.real.array[i] = out.real.array[i+2];
        out_tuple.imaginary.array[i] = out.imaginary.array[i+2];
    }
    freeComplexArrayMemory(samples_interpolated);
    freeComplexArrayMemory(out);
    freeComplexArrayMemory(out_rail);


    return out_tuple;
}
Complex_Array_Tuple coarse_Frequency_Correction(Complex_Array_Tuple testpacket, double fs){
    Complex_Array_Tuple fft_samples = multiplyComplexArrays(testpacket, testpacket);

    //psd = np.fft.fftshift(np.abs(np.fft.fft(fft_samples)))
    Complex_Array_Tuple tempfft = fft(fft_samples);
    Array_Tuple tempAbs = absComplexArray(tempfft);
    Array_Tuple psd = fftshift(tempAbs);
    exportArray(psd, "psd_coarseFrequencyCorrection.txt");

    //f = np.linspace(-fs/2.0, fs/2.0, len(psd))
    Array_Tuple f = linspace(-fs/2.0, fs/2.0, psd.length);
    exportArray(f, "f_coarseFrequencyCorrection.txt");

    //max_freq = f[np.argmax(psd)]
    int psd_argMax = argMax(psd);
    double max_freq = f.array[psd_argMax];

    double Ts = 1.0/fs; // calc sample period
    //t = np.arange(0, Ts*len(testpacket), Ts) # create time vector
    Array_Tuple t_coarse = arange(0, Ts * testpacket.real.length, Ts);
    //testpacket = testpacket * np.exp(-1j*2*np.pi*max_freq*t/2.0) # multiply by negative complex sinusoid at offset frequency
    double* empty_real = (double*)calloc(t_coarse.length, sizeof(double));
    Array_Tuple empty_real_tuple = {empty_real, t_coarse.length};
    for (int i = 0; i < t_coarse.length; i++)
    {
        double temp_value = 2 * M_PI * max_freq * t_coarse.array[i] / 2.0;
        t_coarse.array[i] = temp_value;
    }
    Complex_Array_Tuple complex_input_array = {empty_real_tuple, t_coarse}; // each freed separately
    Complex_Array_Tuple g = expComplexArray(complex_input_array);
    Complex_Array_Tuple new_testpacket = multiplyComplexArrays(testpacket, g);

    freeComplexArrayMemory(fft_samples);
    freeComplexArrayMemory(tempfft);
    freeArrayMemory(tempAbs);
    freeArrayMemory(psd);
    freeArrayMemory(f);
    freeArrayMemory(t_coarse);
    freeArrayMemory(empty_real_tuple);
    freeComplexArrayMemory(g);
    return new_testpacket;
}
Complex_Array_Tuple fine_Frequency_Correction(Complex_Array_Tuple new_testpacket, double fs){
    double alpha = 0.132; // how fast phase updates
    double beta = 0.00932; // how fast frequency updates
    //Costas Loop Fine Freq Correction
    int N = new_testpacket.real.length;
    double phase = 0;
    double freq = 0;
    // These next two params are what to adjust, to make the feedback loop faster or slower (which impacts stability) //! these can and should be adjusted to see what works
    //out = np.zeros(N, dtype=complex)
    Complex_Array_Tuple costas_out = zerosComplex(N);
    //freq_log = []
    double* freq_log = (double*)calloc(N, sizeof(double));

    for (int i = 0; i < N; i++)
    {
        // out[i] = testpacket[i] * np.exp(-1j*phase) # adjust the input sample by the inverse of the estimated phase offset
        double complex packet_complex = new_testpacket.real.array[i] + new_testpacket.imaginary.array[i] * I;
        double complex phase_complex = cexp(-I*phase);
        costas_out.real.array[i] = creal(packet_complex*phase_complex);
        costas_out.imaginary.array[i] = cimag(packet_complex*phase_complex);
        // error = np.real(out[i]) * np.imag(out[i]) # This is the error formula for 2nd order Costas Loop (e.g. for BPSK)
        double error = costas_out.real.array[i] * costas_out.imaginary.array[i];

        // # Advance the loop (recalc phase and freq offset)
        freq += (beta * error);
        // freq_log.append(freq * fs / (2*np.pi)) # convert from angular velocity to Hz for logging
        freq_log[i] = freq * fs / (2.0 * M_PI);
        // phase += freq + (alpha * error)
        phase += freq + alpha * error;
        // # Optional: Adjust phase so its always between 0 and 2pi, recall that phase wraps around every 2pi
        while (phase >= 2.0*M_PI){
            phase -= 2.0*M_PI;
        }
        while (phase < 0){
            phase += 2.0*M_PI;
        }
    }
    Array_Tuple freq_log_tuple = {freq_log, N};
    exportArray(freq_log_tuple, "costasFreqLog.txt");
    free(freq_log);
    //testpacket = costas_out

    // Testing frequency correction //! graphing only
    // fft_samples = testpacket**2
    Complex_Array_Tuple fft_samples_freq = multiplyComplexArrays(costas_out, costas_out);
    // psd = np.fft.fftshift(np.abs(np.fft.fft(fft_samples)))
    Complex_Array_Tuple fft_fft_samples = fft(fft_samples_freq);
    Array_Tuple fft_fft_samples_abs = absComplexArray(fft_fft_samples);
    Array_Tuple psd_freq_correction = fftshift(fft_fft_samples_abs);
    // f = np.linspace(-fs/2.0, fs/2.0, len(psd))
    Array_Tuple f_freq_correct = linspace(-fs/2.0, fs/2.0, psd_freq_correction.length);
    exportArray(psd_freq_correction, "psd_fine_freq_correct.txt");
    exportArray(f_freq_correct, "f_fine_freq_correct.txt");


    freeComplexArrayMemory(fft_samples_freq);
    freeComplexArrayMemory(fft_fft_samples);
    freeArrayMemory(fft_fft_samples_abs);
    freeArrayMemory(psd_freq_correction);
    freeArrayMemory(f_freq_correct);

    return costas_out;
}
Complex_Array_Tuple frame_Sync(Complex_Array_Tuple testpacket, Array_Tuple matchedFilterCoef, Array_Tuple preamble, Array_Tuple data_encoded){
    Array_Tuple abs_IQ_correct = absComplexArray(testpacket);
    double scale = meanArrayTuple(abs_IQ_correct);
    Complex_Array_Tuple out_frame_sync = zerosComplex(testpacket.real.length);
    for (int i = 0; i < testpacket.real.length; i++)
    {
        double real = testpacket.real.array[i];
        double imaj = testpacket.imaginary.array[i];
        out_frame_sync.real.array[i] = (real + scale) / 2 * scale;
        out_frame_sync.imaginary.array[i] = (imaj + scale) / 2 * scale;
    }
    //crosscorr = signal.fftconvolve(out,matched_filter_coef)
    Complex_Array_Tuple crosscorr = convolve(out_frame_sync, matchedFilterCoef);
    exportComplexArray(crosscorr, "crosscorr.txt");
    //idx = np.array(crosscorr).argmax()
    Array_Tuple crosscorr_abs = absComplexArray(crosscorr);
    int idx = argMax(crosscorr_abs);
    //recoveredPayload = testpacket[idx-len(preamble)+1:idx+len(data_encoded)+1] # Reconstruct original packet minus preamble
    int start = idx - preamble.length + 1;
    int end = idx + data_encoded.length + 1;
    int recoveredPayload_length = end - start;
    Complex_Array_Tuple recoveredPayload = zerosComplex(recoveredPayload_length);
    for (int i = start; i < end; i++)
    {
        int index = i - start;
        recoveredPayload.real.array[index] = out_frame_sync.real.array[i];
        recoveredPayload.imaginary.array[index] = out_frame_sync.imaginary.array[i];
    }
    exportComplexArray(recoveredPayload, "recoveredPayload.txt");
    //recoveredData = recoveredPayload[len(preamble):]
    int offset = preamble.length;
    int recoveredData_length = recoveredPayload.real.length - offset;
    Complex_Array_Tuple recoveredData = zerosComplex(recoveredData_length);
    for (int i = offset; i < recoveredPayload.real.length; i++)
    {
        int index = i-offset;
        double real = recoveredPayload.real.array[i];
        double imaj = recoveredPayload.imaginary.array[i];
        recoveredData.real.array[index] = real;
        recoveredData.imaginary.array[index] = imaj;
    }

    freeArrayMemory(abs_IQ_correct);
    freeComplexArrayMemory(crosscorr);
    freeArrayMemory(crosscorr_abs);
    freeComplexArrayMemory(out_frame_sync);
    freeComplexArrayMemory(recoveredPayload);
    return recoveredData; 
}
Array_Tuple demodulation(Complex_Array_Tuple recoveredData, char scheme[], Array_Tuple preamble, Array_Tuple CRC_key){
    Array_Tuple demod_bits = symbol_demod(recoveredData, scheme, 1, preamble.length); // gain has to be set to 1
    int error = CRC_check(demod_bits, CRC_key);
    printf("CRC error: %d\n", error);
    return demod_bits;
}
void DisplayOutput(Array_Tuple data, Array_Tuple demod_bits){
    // This just calculates how many bits were correctly received
    int num_correct = 0;
    int num = 0;
    printf("Tx [%3d]: |", data.length);
    for (int index = 0; index < data.length; index++)
    {
        num++;
        int value = (int)round(demod_bits.array[index]);
        int expected = (int)round(data.array[index]);
        if (value == expected){
            num_correct++;
        }
        printf("%1d|", expected);
    }
    printf("\n");
    printf("Rx [%3d]: |", demod_bits.length);
    for (int i = 0; i < data.length; i++)
    {
        int value = (int)round(demod_bits.array[i]);
        printf("%1d|", value);
    }
    printf("\n");
    
    printf("Received: %d / %d bits   |   %.1lf%%\n", num_correct, num, round((double)num_correct / (double)num * 100.0));
}