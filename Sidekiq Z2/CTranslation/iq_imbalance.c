// gets the average of every value within the period
struct Array_Tuple averages(struct Array_Tuple array, int period){ //! Returns a calloc ptr
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
    struct Array_Tuple out = {means, array.length};
    return out;
}

// Code for IQImbalance Correction
struct Complex_Array_Tuple IQImbalanceCorrect(struct Complex_Array_Tuple packet, int mean_period) //! Returns 2 calloc ptrs 
{
    // This follows the formula shown by the matrix in this article:
    // https://www.faculty.ece.vt.edu/swe/argus/iqbal.pdf

    // divide into real and complex
    /* data */
    struct Array_Tuple I = {packet.real.array, packet.real.length}; // real
    struct Array_Tuple Q = {packet.imaginary.array, packet.imaginary.length}; // imag

    // I(t) = α cos (ωt) + βI
    // Q(t) = sin (ωt + ψ) + βQ
    // BI and BQ are simply the mean over X periods of I or Q, this value can be subtracted to remove
    struct Array_Tuple BI = averages(I, mean_period);
    struct Array_Tuple BQ = averages(Q, mean_period);
    double I_second[I.length];
    for (int i = 0; i < I.length; i++){
        I_second[i] = I.array[i] - BI.array[i];
    }
    double Q_second[Q.length];
    for (int i = 0; i < Q.length; i++){
        Q_second[i] = Q.array[i] - BQ.array[i];
    }
    double* I_second_squared;
    I_second_squared = (double*)calloc(I.length, sizeof(double));
    for (int i = 0; i < I.length; i++){
        I_second_squared[i] = I_second[i] * I_second[i];
    }
    double a = sqrt(2.0 * meanArray(I_second_squared, I.length));
    double* I_times_Q_second;
    I_times_Q_second = (double*)calloc(I.length, sizeof(double));
    for (int i = 0; i < I.length; i++){
        I_times_Q_second[i] = I_second[i] * Q_second[i];
    }

    double sin_psi = (2.0 / a) * meanArray(I_times_Q_second, I.length);

    double cos_psi = sqrt(1 - (sin_psi * sin_psi));
    double A = 1.0 / a;
    double C = -sin_psi / (a * cos_psi);
    double D = 1 / cos_psi;
    double* I_final;
    I_final = (double*)calloc(I.length, sizeof(double));
    double* Q_final;
    Q_final = (double*)calloc(Q.length, sizeof(double));

    for (int i = 0; i < I.length; i++)
    {
        I_final[i] = A * I.array[i];
        Q_final[i] = C * I.array[i] + D * Q.array[i];
    }
    free(BI.array);
    free(BQ.array);
    free(I_second_squared); // both temp arrays
    free(I_times_Q_second);
    struct Array_Tuple I_out = {I_final, I.length};
    struct Array_Tuple Q_out = {Q_final, Q.length};
    struct Complex_Array_Tuple out = {I_out, Q_out};
    return out;
}
