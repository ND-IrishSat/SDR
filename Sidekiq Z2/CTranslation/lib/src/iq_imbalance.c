// iq_imbalance.c
// Rylan Paul


// gets the average of every value within the period
Array_Tuple averages(Array_Tuple array, int period){ //! Returns a calloc ptr
    double* means;
    means = (double*)calloc(array.length, sizeof(double));
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
        means[index] = sum / num;
    }
    Array_Tuple out = {means, array.length};
    return out;
}

// Code for IQImbalance Correction
Complex_Array_Tuple IQImbalanceCorrect(Complex_Array_Tuple packet, int mean_period) //! Returns 2 calloc ptrs 
{
    // This follows the formula shown by the matrix in this article:
    // https://www.faculty.ece.vt.edu/swe/argus/iqbal.pdf

    // divide into real and complex
    /* data */
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
    double* I_second_squared;
    I_second_squared = (double*)calloc(I_tuple.length, sizeof(double));
    for (int i = 0; i < I_tuple.length; i++){
        I_second_squared[i] = I_second[i] * I_second[i];
    }
    double a = sqrt(2.0 * meanArray(I_second_squared, I_tuple.length));
    double* I_times_Q_second;
    I_times_Q_second = (double*)calloc(I_tuple.length, sizeof(double));
    for (int i = 0; i < I_tuple.length; i++){
        I_times_Q_second[i] = I_second[i] * Q_second[i];
    }

    double sin_psi = (2.0 / a) * meanArray(I_times_Q_second, I_tuple.length);

    double cos_psi = sqrt(1 - (sin_psi * sin_psi));
    double A = 1.0 / a;
    double C = -sin_psi / (a * cos_psi);
    double D = 1 / cos_psi;
    double* I_final;
    I_final = (double*)calloc(I_tuple.length, sizeof(double));
    double* Q_final;
    Q_final = (double*)calloc(Q_tuple.length, sizeof(double));

    for (int i = 0; i < I_tuple.length; i++)
    {
        I_final[i] = A * I_tuple.array[i];
        Q_final[i] = C * I_tuple.array[i] + D * Q_tuple.array[i];
    }
    freeArrayMemory(BI);
    freeArrayMemory(BQ);
    free(I_second_squared); // both temp arrays
    free(I_times_Q_second);
    Array_Tuple I_out = {I_final, I_tuple.length};
    Array_Tuple Q_out = {Q_final, Q_tuple.length};
    Complex_Array_Tuple out = {I_out, Q_out};
    return out;
}
