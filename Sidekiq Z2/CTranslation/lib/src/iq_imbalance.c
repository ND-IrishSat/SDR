// iq_imbalance.c
// Rylan Paul

#include "../standardArray.h"
// gets the average of every value within the period
Array_Tuple averages(Array_Tuple array, int period){
    Array_Tuple out = zerosArray(array.length);
    for (int index = 0; index < array.length; index++)
    {
        double num = 0;
        double sum = array.array[index];
        num ++;
        for (int i = 1; i < period+1; i++)
        {
            int end = 0;
            if (index-i >= 0){
                num += 1;
                sum += array.array[index-i];
            } else {
                end = 1;
            }
            if (index+i < array.length){
                num += 1;
                sum += array.array[index + i];
            } else if (end == 1){
                break;
            }
        }
        out.array[index] = sum / num;
    }
    return out;
}

// Code for IQImbalance Correction
Complex_Array_Tuple IQImbalanceCorrect(Complex_Array_Tuple packet, int mean_period)
{
    // This follows the formula shown by the matrix in this article:
    // https://www.faculty.ece.vt.edu/swe/argus/iqbal.pdf

    // divide into real and complex
    Array_Tuple I_tuple = {packet.real.array, packet.real.length}; // real
    Array_Tuple Q_tuple = {packet.imaginary.array, packet.imaginary.length}; // imag

    // I(t) = α cos (ωt) + βI
    // Q(t) = sin (ωt + ψ) + βQ
    // BI and BQ are simply the mean over X periods of I or Q, this value can be subtracted to remove
    Array_Tuple BI = averages(I_tuple, mean_period);
    Array_Tuple BQ = averages(Q_tuple, mean_period);
    double I_second[I_tuple.length];
    for (int i = 0; i < I_tuple.length; i++){
        I_second[i] = I_tuple.array[i] - BI.array[i];
    }
    double Q_second[Q_tuple.length];
    for (int i = 0; i < Q_tuple.length; i++){
        Q_second[i] = Q_tuple.array[i] - BQ.array[i];
    }
    Array_Tuple I_second_squared = zerosArray(I_tuple.length);
    for (int i = 0; i < I_tuple.length; i++){
        I_second_squared.array[i] = I_second[i] * I_second[i];
    }
    double a = sqrt(2.0 * meanArrayTuple(I_second_squared));
    Array_Tuple I_times_Q_second = zerosArray(I_tuple.length);
    for (int i = 0; i < I_tuple.length; i++){
        I_times_Q_second.array[i] = I_second[i] * Q_second[i];
    }

    double sin_psi = (2.0 / a) * meanArrayTuple(I_times_Q_second);

    double cos_psi = sqrt(1 - (sin_psi * sin_psi));
    double A = 1.0 / a;
    double C = -sin_psi / (a * cos_psi);
    double D = 1 / cos_psi;
    Complex_Array_Tuple out = zerosComplex(I_tuple.length);

    for (int i = 0; i < I_tuple.length; i++)
    {
        out.real.array[i] = A * I_tuple.array[i];
        out.imaginary.array[i] = C * I_tuple.array[i] + D * Q_tuple.array[i];
    }
    freeArrayMemory(BI);
    freeArrayMemory(BQ);
    freeArrayMemory(I_second_squared);
    freeArrayMemory(I_times_Q_second);
    return out;
}
