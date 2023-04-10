import numpy as np
import time
import matplotlib.pyplot as plt
from scipy import signal

# Define the demodulation function

# The demodulation function takes the following arguments as inputs:

# a_in1:                        This array contains the received samples to be demodulated, samples here include both the preamble and the payload.

# samples_perbit:               This is the oversampling factor on the receiver side (= 4/3*oversampling factor of the transmitter), each symbol you would like to demodulate
#                               is represented by "samples_perbit" bits

# payload_start:                This is the index of the received samples that indicates the beginning of the payload
#                               For example, if an array "a_in1" contains both the preamble and payload, the a_in1[payload_start] should be the first sample of the payload

# N                             This value is equal to the length of the payload in bits.

# Ts_in:                        This is the ADC input sampling period in sec, Ts_in = 1/fs_in


# The demodulation function returns the following arguments as outputs:


# baseband_symbols              This is the array of baseband symbols

# a_demodulated                 This is the bit array containing the demodulated payload


# Note: general (wireless) channel noncoherent demodulation requires both I,Q data as energy that was originally only in I could leak
# into Q channel as a result of the carrier freq & phase offset


def noncoherent_demodulation_wireline(a_in1, samples_perbit, payload_start, N, fs_in, Ts_in, preamble_length, showplot):


        payload = a_in1[payload_start:(payload_start + N*samples_perbit)]

        payload_energy = np.square(np.abs(payload))

        #Estimate the decision threshold for energy detection
        bits_used_for_estimation = 1000
        payload_energy_segment = payload_energy[0:int(bits_used_for_estimation*samples_perbit)]

        average_energy_per_symbol = np.sum(payload_energy_segment)/bits_used_for_estimation

        
        #frame sync error detected
        if(len(np.real(payload))!=int(N*samples_perbit)):
                print("error")
                print("payload length: ",len(np.real(payload)))
                print("N*samples_perbit", N*samples_perbit)
                return 0,0

        symbol_energy = np.mean(payload_energy.reshape(-1, samples_perbit), axis=1)

        bits_demodulated = np.zeros(N,dtype=int)
                
        for symbol_index in range(N):
                current_symbol_energy = samples_perbit*symbol_energy[symbol_index]
                if(current_symbol_energy>average_energy_per_symbol):
                        bits_demodulated[symbol_index] = 1
                        

               
        return bits_demodulated, None




