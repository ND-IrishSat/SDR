// signal_processing.c
#include <math.h>
#include "../fftw3/include/fftw3.h"
#include "../standardArray.h"

// Using the fftw3 library, calculates discrete fft on an array np.fft.fft()

// FFT function for Complex_Array_Tuple
Complex_Array_Tuple fft(Complex_Array_Tuple x) {
    int N = x.real.length;
    if (N <= 1) return x;  // base case for recursion

    // allocate memory for even and odd parts
    Complex_Array_Tuple even = zerosComplex(N / 2);
    Complex_Array_Tuple odd = zerosComplex(N / 2);

    // split the array into even and odd indices
    for (int i = 0; i < N / 2; i++) {
        even.real.array[i] = x.real.array[2 * i];
        even.imaginary.array[i] = x.imaginary.array[2 * i];
        odd.real.array[i] = x.real.array[2 * i + 1];
        odd.imaginary.array[i] = x.imaginary.array[2 * i + 1];
    }

    // Recursively apply FFT to both halves
    Complex_Array_Tuple fft_even = fft(even);
    Complex_Array_Tuple fft_odd = fft(odd);

    // Prepare output structure
    Complex_Array_Tuple result = zerosComplex(N);

    // Combine the results
    for (int k = 0; k < N / 2; k++) {
        double t = -2 * M_PI * k / N;
        double cos_t = cos(t);
        double sin_t = sin(t);

        double temp_real = cos_t * fft_odd.real.array[k] - sin_t * fft_odd.imaginary.array[k];
        double temp_imag = cos_t * fft_odd.imaginary.array[k] + sin_t * fft_odd.real.array[k];

        result.real.array[k] = fft_even.real.array[k] + temp_real;
        result.imaginary.array[k] = fft_even.imaginary.array[k] + temp_imag;

        result.real.array[k + N / 2] = fft_even.real.array[k] - temp_real;
        result.imaginary.array[k + N / 2] = fft_even.imaginary.array[k] - temp_imag;
    }

    // Free temporary structures
    free(even.real.array);
    free(even.imaginary.array);
    free(odd.real.array);
    free(odd.imaginary.array);

    return result;
}
//np.hamming
Array_Tuple hamming(int M){ 
    Array_Tuple out = zerosArray(M);

    for (int i = 0; i < M; i++)
    {
        out.array[i] = 0.54 - 0.46 * cos((2 * M_PI * i)/(M - 1));
    }
    return out;
}
// Function to perform FFT shift along a single dimension
Array_Tuple fftshift(Array_Tuple data) { 
    double temp;
    // Calculate the midpoint of the array
    int midpoint = data.length / 2;
    Array_Tuple out = zerosArray(midpoint*2);
    // Perform the shift
    for (int i = 0; i < midpoint; i++) {
        temp = data.array[i];
        out.array[i] = data.array[i + midpoint];
        out.array[i + midpoint] = temp;
    }
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
// Takes a bit array and converts it to a pulse train from -1 to 1 with sps-1 zeros between each bit
Array_Tuple pulsetrain(Array_Tuple bits, int sps){ 
    int length = bits.length * sps;
    Array_Tuple t = zerosArray(length);
    int offset = 0;
    for (int i = 0; i < bits.length; i++)
    {
        for (int x = 0; x < sps; x++)
        {
            if (x == 0){
                t.array[x+offset] = bits.array[i]*2 - 1; // sets 0 to -1 and 1 to 1
            } else {
                t.array[x+offset] = 0; // puts zeros between the bits to make it easier to read
            }
        }
        offset += sps;
    }
    return t;
    
}
Array_Tuple firwin(int M, double cutoff){ 
    Array_Tuple out = zerosArray(M+1);
    //scipy.signal.firwin
    for (int i = 0; i < M; i++)
    {
        double window = (0.54 - 0.46 * cos(2 * M_PI * i / M));
        
        double lowpassfilter = cutoff; // limit of sin(cutoff * M_PI * (i - M/2)) / (M_PI * (i - M/2)) as (i-M/2) --> 0 = cutoff
        if (i - M/2 != 0){
            lowpassfilter = sin(cutoff * M_PI * (i - M/2)) / (M_PI * (i - M/2));
        }
        out.array[i] = window * lowpassfilter;
    }
    
    return out;
}
// same as scipy.resample_poly which upsamples an array of data, can also specify a down sample
Complex_Array_Tuple resample_poly(Complex_Array_Tuple a, int up, int down){
    Complex_Array_Tuple up_complex = zerosComplex(a.real.length * up);
    // upsample by factor of up
    int index = 0;
    // upsample the input
    for (int i = 0; i < a.real.length * up; i++)
    {
        if (i % up == 0){
            up_complex.real.array[i] = a.real.array[index];
            up_complex.imaginary.array[i] = a.imaginary.array[index];
            index ++;
        }
        else {
            up_complex.real.array[i] = 0;
            up_complex.imaginary.array[i] = 0;
        }
    }
    int greaterUpDown = (up >= down ? up : down);
    int filter_length = 10 * greaterUpDown;

    Array_Tuple filterCoeff = firwin(filter_length, 1.0 / (double)greaterUpDown);
    

    Complex_Array_Tuple smoothed = convolveSame(up_complex, filterCoeff);

    int down_length = smoothed.real.length / down;
    Complex_Array_Tuple down_complex = zerosComplex(down_length);
    index = 0;
    double max_abs_value_sqr = -99;
    for (int i = 0; i < smoothed.real.length; i++)
    {
        if (i % down == 0){
            double x = smoothed.real.array[i];
            double y = smoothed.imaginary.array[i];
            double abs_value_sqr = x*x + y*y;
            down_complex.real.array[index] = x;
            down_complex.imaginary.array[index] = y;
            index ++;
            if (abs_value_sqr > max_abs_value_sqr){
                max_abs_value_sqr = abs_value_sqr;
            }
        }
    }

    // normalize output then scale to input scale
    Complex_Array_Tuple complex_out = zerosComplex(down_length);
    double max_abs_input = maxAbsoluteValue(a);
    for (int i = 0; i < down_length; i++)
    {
        complex_out.real.array[i] = down_complex.real.array[i] / sqrt(max_abs_value_sqr) * max_abs_input;
        complex_out.imaginary.array[i] = down_complex.imaginary.array[i] / sqrt(max_abs_value_sqr) * max_abs_input;
    }

    // free all temp arrays
    freeComplexArrayMemory(up_complex);
    freeArrayMemory(filterCoeff);
    freeComplexArrayMemory(smoothed);
    freeComplexArrayMemory(down_complex);

    return complex_out;
}

// returns complex noise on a complex array
Complex_Array_Tuple generateComplexNoise(Complex_Array_Tuple testpacket, double std_dev, double phase_noise_strength, double noise_power){
    // vars
    double mean = 0;
    int num_samples = testpacket.real.length;

    // creating general normal distribution random arrays for real and complex
    Complex_Array_Tuple out = zerosComplex(num_samples);
    for (int i = 0; i < num_samples; i++)
    {

        double real = testpacket.real.array[i];
        double imaj = testpacket.imaginary.array[i];

        double x = rand_norm(mean, std_dev);
        double x_2 = x / sqrt(2.0);
        double x_3 = x_2 / sqrt(noise_power);

        double y = rand_norm(mean, std_dev);
        double y_2 = y / sqrt(2.0);
        double y_3 = y_2 / sqrt(noise_power);

        double complex awgn_noise = x_3 + I*y_3;

        double complex phase_noise = cexp( I * (rand_norm(mean, std_dev) * phase_noise_strength) );
        double complex packet_complex = real + I * imaj;
        double complex out_packet = (packet_complex + awgn_noise) * phase_noise;
        out.real.array[i] = creal(out_packet);
        out.imaginary.array[i] = cimag(out_packet);
    }
    return out;
}


