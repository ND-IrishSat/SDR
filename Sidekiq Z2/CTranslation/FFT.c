#include <stdio.h>
#include <stdlib.h>
#include <math.h>

typedef struct {
    double real;
    double imag;
} Complex;

void fft(Complex *x, int N) {
    if (N <= 1) return;

    // Split array into even and odd indices
    Complex *even = (Complex *)malloc(N / 2 * sizeof(Complex));
    Complex *odd = (Complex *)malloc(N / 2 * sizeof(Complex));
    for (int i = 0; i < N / 2; i++) {
        even[i] = x[i * 2];
        odd[i] = x[i * 2 + 1];
    }

    // Recursively apply FFT on both halves
    fft(even, N / 2);
    fft(odd, N / 2);

    // Combine the results
    for (int k = 0; k < N / 2; k++) {
        double t = -2 * M_PI * k / N;
        Complex exp = {cos(t), sin(t)};
        Complex temp = {exp.real * odd[k].real - exp.imag * odd[k].imag,
                        exp.real * odd[k].imag + exp.imag * odd[k].real};
        
        x[k].real = even[k].real + temp.real;
        x[k].imag = even[k].imag + temp.imag;
        x[k + N / 2].real = even[k].real - temp.real;
        x[k + N / 2].imag = even[k].imag - temp.imag;
    }

    free(even);
    free(odd);
}

int main() {
    int N = 8; // Ensure N is a power of 2
    Complex x[N];

    // Example input (8 samples)
    for (int i = 0; i < N; i++) {
        x[i].real = i + 1; // Real part
        x[i].imag = 0;     // Imaginary part (0 for real input)
    }

    fft(x, N);

    // Output the FFT result
    printf("FFT result:\n");
    for (int i = 0; i < N; i++) {
        printf("x[%d] = %.5f + %.5fi\n", i, x[i].real, x[i].imag);
    }

    return 0;
}
