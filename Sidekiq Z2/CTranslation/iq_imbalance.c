// gets the average of every value within the period, also be sure to free memory after use
struct Array_Tuple averages(struct Array_Tuple array, int period){
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

struct Array_Tuple IQImbalanceCorrect(struct Array_Tuple packet, int mean_period)
{

    // divide into real and complex
    double* I = packet.array; // real
    double* Q = packet.array; // imaginary
    /* data */
    struct Array_Tuple out = {packet.array, packet.length};
    return out;
}
