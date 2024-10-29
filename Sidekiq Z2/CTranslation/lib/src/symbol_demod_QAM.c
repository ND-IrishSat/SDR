
// new version of symbol demodulation
#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <math.h>
#include <string.h>
#include "../standardArray.h"

Array_Tuple symbol_demod(Complex_Array_Tuple baseband_symbols_complex, char scheme[], double channel_gain, double preamble_len)
{
  int l = baseband_symbols_complex.real.length;
  //baseband_symbols[0][0] = real
  //baseband_symbols[0][0] = imaginary

  int baseBandLen = baseband_symbols_complex.real.length;

  Array_Tuple a_demodulated = zerosArray(baseBandLen);
  // baseband_symbols_I = baseband_symbols_complex.real
  // baseband_symbols_Q = baseband_symbols_complex.imaginary

  // check this
  if (strcmp(scheme, "OOK") == 0)
  {
    double s_on = 1.0 * channel_gain;
    double s_off = 0 * channel_gain;

    Array_Tuple array_on = subtractDoubleFromArray(baseband_symbols_complex.real, s_on); // minus s_on
    Array_Tuple array_off = subtractDoubleFromArray(baseband_symbols_complex.real, s_off); // minus s_off
    Complex_Array_Tuple complex_array1 = {array_on, baseband_symbols_complex.imaginary};
    Complex_Array_Tuple complex_array2 = {array_off, baseband_symbols_complex.imaginary};
    Array_Tuple abs_baseband_symbols_1 = absComplexArray(complex_array1);
    Array_Tuple abs_baseband_symbols_2 = absComplexArray(complex_array2);
    for (int i = 0; i < baseband_symbols_complex.real.length; i++)
    {
      if (abs_baseband_symbols_1.array[i] < abs_baseband_symbols_2.array[i])
      {
          a_demodulated.array[i] = 1;
      }
      else
      {
        a_demodulated.array[i] = 0;
      }
    }
    freeArrayMemory(array_on);
    freeArrayMemory(array_off);
    freeArrayMemory(abs_baseband_symbols_1);
    freeArrayMemory(abs_baseband_symbols_2);
  } else if (strcmp(scheme, "BPSK") == 0)
  {
    double reference_plus = -1.0 * channel_gain;
    double reference_minus = -1.0 * reference_plus;
    
    Array_Tuple array_plus = subtractDoubleFromArray(baseband_symbols_complex.real, reference_plus); // minus s_on
    Array_Tuple array_minus = subtractDoubleFromArray(baseband_symbols_complex.real, reference_minus); // minus s_off
    Complex_Array_Tuple complex_array1 = {array_plus, baseband_symbols_complex.imaginary};
    Complex_Array_Tuple complex_array2 = {array_minus, baseband_symbols_complex.imaginary};
    Array_Tuple abs_baseband_symbols_1 = absComplexArray(complex_array1);
    Array_Tuple abs_baseband_symbols_2 = absComplexArray(complex_array2);
    
    for (int i = 0; i < baseband_symbols_complex.real.length; i++)
    {
      if (abs_baseband_symbols_1.array[i] < abs_baseband_symbols_2.array[i]) //not sure what goes here python code is kinda funky
        {
          a_demodulated.array[i] = 0;
        }
      else
      {
        a_demodulated.array[i] = 1;
      }
    }
    freeArrayMemory(array_plus);
    freeArrayMemory(array_minus);
    freeArrayMemory(abs_baseband_symbols_1);
    freeArrayMemory(abs_baseband_symbols_2);
  }
  else{
    printf("QPSK is not suppported yet, see commented code below");
  }
  return a_demodulated;
}
  /*
  if (strcmp(scheme, "QPSK") == 0)
  {
    float I_demodulated[] = {0};
    float Q_demodulated[] = {0};
    // construct the recieved symbols on the complex plane (signal space constellation))
    for (int i = 0; i < baseBandLen; i++)
    {
      baseband_symbols_complex[i] = baseband_symbols_I + 1 * j * baseband_symbols_Q;
    }
    // compute and define the reference signals in the (signal space constellation)
    float reference_00 = -1.0 * channel_gain - 1 * j * channel_gain;
    float reference_11 = 1.0 * channel_gain + 1 * j * channel_gain;
    float reference_01 = -1.0 * channel_gain + 1 * j * channel_gain;
    float reference_10 = 1.0 * channel_gain - 1 * j * channel_gain;

    // this part is kinda funky CHECK!!
    float symbol = 0;
    float comparison_string[4] = [ abs(symbol - reference_11), abs(symbol - reference_10), abs(symbol - reference_00), abs(symbol - reference_01) ];
    float min = comparison_string[0];
    int length = sizeof(I_demodulated) / sizeof(I_demodulated[0]);
    for (int i = 0; i < sizeof(baseband_symbols_complex) / sizeof(baseband_symbols_complex[0][0]); i++)
    {
      symbol = baseband_symbols_complex[i];
      for (int i = 1; i < 4; i++)
      {
        if (comparison_string[i] < min)
        {
          min = comparison_string[i];
        }
      }
      if (cabs(symbol - reference_11) == min)
      {
        I_demodulated[length] = 1;
        Q_demodulated[length] = 1;
      }
      else if (cabs(symbol - reference_10) == min)
      {
        I_demodulated[length] = 1;
        Q_demodulated[length] = 0;
      }
      else if (cabs(symbol - reference_00) == min)
      {
        I_demodulated[length] = 0;
        Q_demodulated[length] = 0;
      }
      else if (cabs(symbol - reference_01) == min)
      {
        I_demodulated[length] = 0;
        I_demodulated[length] = 1;
      }
      // in the python function we return the a_demodulated array from the function but unsure how this plays out in C?? maybe store in temp var
      for (int i = 0; i < sizeof(I_demodulated) / sizeof(I_demodulated[0]); i++)
      {
        a_demodulated[sizeof(a_demodulated) / sizeof(a_demodulated[0]) + 1] = I_demodulated[i];
      }
      for (int i = 0; i < sizeof(Q_demodulated) / sizeof(Q_demodulated[0]); i++)
      {
        a_demodulated[sizeof(a_demodulated) / sizeof(a_demodulated[0]) + 1] = Q_demodulated[i];
      }
    }
    if (strcmp(scheme, 'QAM') == 0)
    {
      static float preamble_symbols[preamble_len];
      for (int i = 0; i < preamble_len; i++)
      {
        preamble_symbols[i] = baseband_symbols[0][i];
      }
      preamble_demod[preamble_len];
      for (int i = 0; i < preamble_len; i++)
      {
        if (preamble_symbols[i] > 0.5)
        {
          preamble_demod[i] = 1;
        }
        else
        {
          preamble_demod[i] = 0;
        }
      }

      // variable for baseband symbols length
      newArrayLen = baseBandLen;
      newArray[newArrayLen];
      for (int i = 0; i < newArrayLen; i++)
      {
        newArray[0, i] = baseband_symbols[0, i + preamble_len];
      }
      newArray2[newArrayLen];
      for (int i = 0; i < newArrayLen; i++)
      {
        newArray[1, i] = baseband_symbols[1, i + preamble_len];
      }

      // help isaac
      int len_I = sizeof(baseband_symbols_I) / sizeof(baseband_symbols_I[0]);
      int len_Q = sizeof(baseband_symbols_Q) / sizeof(baseband_symbols_Q[0]);
      int delta = 0;
      float new_baseband_symbols_I[] = {0};
      float new_baseband_symbols_Q[] = {0};
      // need to remove the last delta elements from the end of the longer array so len_I = len_Q
      // thoughts = make a temp array to copy 0 to len_Q elements then copy those elements back into baseband_symbols_I
      // how to do without malloc? - function call, pointers, (unsure)

      if (len_I != len_Q)
      {
        if (len_I > len_Q)
        {
          delta = len_I - len_Q;
          // cutting off delta from the end of the baseband_symbols_I array so that I and Q are equal in length
          // rm last delta elements from I
          // baseband_symbols_I = baseband_symbols_I[...];

          for (int i = 0; i < len_Q; i++)
          {
            strcpy(new_baseband_symbols_I[i], baseband_symbols_I[i]);
            strcpy(new_baseband_symbols_Q[i], baseband_symbols_Q[i]);
          }
        }
        else if (len_Q > len_I)
        {
          delta = len_Q - len_I;
          for (int i = 0; i < len_I; i++)
          {
            strcpy(new_baseband_symbols_I[i], baseband_symbols_I[i]);
            strcpy(new_baseband_symbols_Q[i], baseband_symbols_Q[i]);
          }
        }
      }
      // change baseband_symbols_I[] in the following for loop to new_baseband_symbols_I[]

      // initialized I_demodulated and Q_demodulated again here but not sure way
      // Elizabeth - find way to write to include all values up to delta in the [...] arrays

      // construct part of python code --> baseband_symblos_complex = baseband_symbols_I + 1*j * baseband_symbols_Q
      for (int i = 0; i < strlen(baseband_symbols_complex); i++)
      {
        baseband_symbols_complex[i] = new_baseband_symbols_I[i] + 1 * j * new_baseband_symbols_Q[i];
      }
      // WHERE ALYSSA LEFT OFF --> compute and define the reference signals in the signal space (16 constellation points)
      float reference_0000 = -3.0 * channel_gain - 3 * I * channel_gain;
      float reference_0001 = -1.0 * channel_gain - 3 * I * channel_gain;
      float reference_0010 = -3.0 * channel_gain - 3 * I * channel_gain;
      float reference_0011 = 1.0 * channel_gain - 3 * I * channel_gain;
      float reference_0100 = -3.0 * channel_gain - 1 * I * channel_gain;
      float reference_0101 = -1.0 * channel_gain - 1 * I * channel_gain;
      float reference_0110 = 3.0 * channel_gain - 1 * I * channel_gain;
      float reference_0111 = 1.0 * channel_gain - 1 * I * channel_gain;
      float reference_1000 = -3.0 * channel_gain + 3 * I * channel_gain;
      //TODO float reference_1001 = -3.0 * channel_gain + 3 * I * channel_gain;
      float reference_1010 = 3.0 * channel_gain + 3 * I * channel_gain;
      float reference_1011 = 1.0 * channel_gain + 3 * I * channel_gain;
      float reference_1100 = -3.0 * channel_gain + 1 * I * channel_gain;
      float reference_1101 = -1.0 * channel_gain + 1 * I * channel_gain;
      float reference_1110 = 3.0 * channel_gain + 1 * I * channel_gain;
      float reference_1111 = 1.0 * channel_gain + 1 * I * channel_gain;
      // create a new array to use to find the minimum for the next large bit of code that compares the absolute values to populate a 1 or 0 in the Q_demodulated and I_demodulated arrays
      float new_comparison[16] = [ cabs(symbol - reference_1111), cabs(symbol - reference_1110), cabs(symbol - reference_1101), cabs(symbol - reference_1100), cabs(symbol - reference_1011), cabs(symbol - reference_1010), cabs(symbol - reference_1001), cabs(symbol - reference_1000), cabs(symbol - reference_0111), cabs(symbol - reference_0110), cabs(symbol - reference_0101), cabs(symbol - reference_0100), cabs(symbol - reference_0011), cabs(symbol - reference_0010), cabs(symbol - reference_0001), cabs(symbol - reference_0000) ];
      float newMin = new_comparison[0];
      for (int i = 1; i < 16; i++)
      {
        if (new_comparison[i] < newMin)
        {
          newMin = new_comparison[i];
        }
      }
      // Start a for-loop to iterate through all complex symbols and make a decision on 4-bits of data
      for (int i = 0; i < sizeof(baseband_symbols_complex) / sizeof(baseband_symbols_complex[0][0]); i++)
      {
        // ????
        int newLengthofQ = sizeof(Q_demodulated) / sizeof(Q_demodulated[0]);
        symbol = baseband_symbols_complex[i];
        if (cabs(symbol - reference_1111) == newMin)
        {
          Q_demodulated[newLengthofQ + 1] = 1;
          Q_demodulated[newLengthofQ + 2] = 1;
          I_demodulated[newLengthofQ + 1] = 1;
          I_demodulated[newLengthofQ + 2] = 1;
        }
        else if (cabs(symbol - reference_1110) == newMin)
        {
          Q_demodulated[newLengthofQ + 1] = 1;
          Q_demodulated[newLengthofQ + 2] = 1;
          I_demodulated[newLengthofQ + 1] = 1;
          I_demodulated[newLengthofQ + 2] = 0;
        }
        else if (cabs(symbol - reference_1101) == newMin)
        {
          Q_demodulated[newLengthofQ + 1] = 1;
          Q_demodulated[newLengthofQ + 2] = 1;
          I_demodulated[newLengthofQ + 1] = 0;
          I_demodulated[newLengthofQ + 2] = 1;
        }
        else if (cabs(symbol - reference_1100) == newMin)
        {
          Q_demodulated[newLengthofQ + 1] = 1;
          Q_demodulated[newLengthofQ + 2] = 1;
          I_demodulated[newLengthofQ + 1] = 1;
          I_demodulated[newLengthofQ + 2] = 0;
        }
        else if (cabs(symbol - reference_1011) == newMin)
        {
          Q_demodulated[newLengthofQ + 1] = 1;
          Q_demodulated[newLengthofQ + 2] = 0;
          I_demodulated[newLengthofQ + 1] = 1;
          I_demodulated[newLengthofQ + 2] = 1;
        }
        else if (cabs(symbol - reference_1010) == newMin)
        {
          Q_demodulated[newLengthofQ + 1] = 1;
          Q_demodulated[newLengthofQ + 2] = 0;
          I_demodulated[newLengthofQ + 1] = 1;
          I_demodulated[newLengthofQ + 2] = 0;
        }
        else if (cabs(symbol - reference_1001) == newMin)
        {
          Q_demodulated[newLengthofQ + 1] = 1;
          Q_demodulated[newLengthofQ + 2] = 0;
          I_demodulated[newLengthofQ + 1] = 0;
          I_demodulated[newLengthofQ + 2] = 1;
        }
        else if (cabs(symbol - reference_1000) == newMin)
        {
          Q_demodulated[newLengthofQ + 1] = 1;
          Q_demodulated[newLengthofQ + 2] = 0;
          I_demodulated[newLengthofQ + 1] = 0;
          I_demodulated[newLengthofQ + 2] = 0;
        }
        else if (cabs(symbol - reference_0111) == newMin)
        {
          Q_demodulated[newLengthofQ + 1] = 0;
          Q_demodulated[newLengthofQ + 2] = 1;
          I_demodulated[newLengthofQ + 1] = 1;
          I_demodulated[newLengthofQ + 2] = 1;
        }
        else if (cabs(symbol - reference_0110) == newMin)
        {
          Q_demodulated[newLengthofQ + 1] = 0;
          Q_demodulated[newLengthofQ + 2] = 1;
          I_demodulated[newLengthofQ + 1] = 1;
          I_demodulated[newLengthofQ + 2] = 0;
        }
        else if (cabs(symbol - reference_0101) == newMin)
        {
          Q_demodulated[newLengthofQ + 1] = 0;
          Q_demodulated[newLengthofQ + 2] = 1;
          I_demodulated[newLengthofQ + 1] = 0;
          I_demodulated[newLengthofQ + 2] = 1;
        }
        else if (cabs(symbol - reference_0100) == newMin)
        {
          Q_demodulated[newLengthofQ + 1] = 0;
          Q_demodulated[newLengthofQ + 2] = 1;
          I_demodulated[newLengthofQ + 1] = 0;
          I_demodulated[newLengthofQ + 2] = 0;
        }
        else if (cabs(symbol - reference_0011) == newMin)
        {
          Q_demodulated[newLengthofQ + 1] = 0;
          Q_demodulated[newLengthofQ + 2] = 0;
          I_demodulated[newLengthofQ + 1] = 1;
          I_demodulated[newLengthofQ + 2] = 1;
        }
        else if (cabs(symbol - reference_0010) == newMin)
        {
          Q_demodulated[newLengthofQ + 1] = 0;
          Q_demodulated[newLengthofQ + 2] = 0;
          I_demodulated[newLengthofQ + 1] = 1;
          I_demodulated[newLengthofQ + 2] = 0;
        }
        else if (cabs(symbol - reference_0001) == newMin)
        {
          Q_demodulated[newLengthofQ + 1] = 0;
          Q_demodulated[newLengthofQ + 2] = 0;
          I_demodulated[newLengthofQ + 1] = 0;
          I_demodulated[newLengthofQ + 2] = 1;
        }
        else if (cabs(symbol - reference_0000) == newMin)
        {
          Q_demodulated[newLengthofQ + 1] = 0;
          Q_demodulated[newLengthofQ + 2] = 0;
          I_demodulated[newLengthofQ + 1] = 0;
          I_demodulated[newLengthofQ + 2] = 0;
        }
      }
    }

    // this is the very last while loop after all the elif statements -- in line with the for loop that iterates through all complex symbols and make a decision on 4-bits of data (line 236 in py)
    int i = 0;
    int j = 0;
    a_demodulated[4];
    while (j < (sizeof(Q_demodulated) / sizeof(Q_demodulated[0])))
    {
      // 4 append statements
      // where to start a_demodulated
      a_demodulated[j] = Q_demodulated[i];
      a_demodulated[j + 1] = Q_demodulated[i + 1];
      a_demodulated[j + 2] = I_demodulated[i];
      a_demodulated[j + 3] = I_demodulated[i + 1];
      j = j + 4;
      i = i + 2;
    }
  }
  

*/
//}
