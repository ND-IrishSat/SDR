import numpy as np
import time
import matplotlib.pyplot as plt
from scipy import signal
import matched_filtering
import pyfftw
import phase_sync

# Define the demodulation function

# The demodulation function takes the following arguments as inputs:

# a_in1:                        This array contains the received samples to be demodulated, samples here include both the preamble and the payload.

# samples_perbit:               This is the oversampling factor on the receiver side (= 4/3*oversampling factor of the transmitter), each symbol you would like to demodulate
#                               is represented by "samples_perbit" bits

# carrier_frequency:            This is the estimated frequency of the IF carrier in Hz

# carrier_phase:                This is the estimated phase of the IF carrier in radian

# payload_start:                This is the index of the received samples that indicates the beginning of the payload
#                               For example, if an array "a_in1" contains both the preamble and payload, the a_in1[payload_start] should be the first sample of the payload

# N                             This value is equal to the length of the payload in bits.

# Ts_in:                        This is the ADC input sampling period in sec, Ts_in = 1/fs_in

# channel_gain                  This is the gain of the channel, channel impulse response is simply modeled as g(t)= channel_gain*delta(t)

# pulse_shaping                 This is a string indicating whether pulse shaping is applied on the transmitter side

# scheme                        Modulation Scheme

# The demodulation function returns the following arguments as outputs:


# (with pulse_shaping == "On" )

# baseband_symbols              This is the array of baseband symbols


# (with pulse_shaping == "Off" )

# a_demodulated                 This is the bit array containing the demodulated payload

def coherent_demodulation(a_in1, samples_perbit, carrier_frequency, payload_start, N, fs_in, Ts_in, pulse_shape, pulse_shaping, scheme, preamble_length, showplot):

        debugging = 1 # print out statistics for the developer

        t1 = time.time()

        if(scheme == "QPSK"):
                N = int(N/2)

        if(scheme == "QAM"):
                N = int(N/4)

        ### test to make sure this outputs the correct payload ###
        payload_before_correction = a_in1[payload_start:(payload_start + N*samples_perbit)]
        


        ones_length = preamble_length - 20
        
        print("payload_start-ones_length*samples_perbit:")
        print(payload_start-ones_length*samples_perbit)
        print("payload_start + N*samples_perbit:")
        print(payload_start + N*samples_perbit)
        
        print("a_in1:")
        print(a_in1)
        print("a_in1 length:")
        print(len(a_in1))
        
        payload_and_ones = a_in1[int(payload_start-ones_length*samples_perbit):(payload_start + N*samples_perbit)]
        
        
        print("payload_and_ones:")
        print(payload_and_ones)
        
        k = np.arange(len(payload_and_ones))
        
        Digital_LO = np.exp(-1j*2*np.pi*carrier_frequency*(k*Ts_in))
        print("Digital_LO: ")
        print(Digital_LO)

        #Correct the frequency and then extract & correct the phase

        payload_and_ones_corrected = phase_sync.phase_sync(payload_and_ones, Digital_LO, payload_start, ones_length, samples_perbit)
        print("payload_and_ones_corrected:")
        print(payload_and_ones_corrected)

# Save the frame sync result
        with open('phasesync_result.txt', 'w') as fs:
                fs.write('Phase Synchronization result:\n')
                fs.write('Packet Data: '+str(payload_and_ones)+'\n')
                fs.write('Digital LO: '+str(Digital_LO)+'\n')
                fs.write('Payload Start: '+str(payload_start)+'\n')
                fs.write('Preamble Length: '+str(ones_length)+'\n')
                fs.write('Oversampling Factor: '+str(samples_perbit)+'\n')
                fs.write('Phase Correct Length: '+str(payload_and_ones_corrected)+'\n')

        # save result to a numpy array as well
        np.savez('PhaseSync_data', payload_and_ones = payload_and_ones, Digital_LO = Digital_LO, payload_start = payload_start, \
                 ones_length = ones_length, samples_perbit = samples_perbit, payload_and_ones_corrected= payload_and_ones_corrected)

        payload_corrected = payload_and_ones_corrected[ones_length*samples_perbit:]
        print("len(payload_corrected)")
        print(len(payload_corrected))
        print("payload_corrected")
        print(payload_corrected[10000:10020])

        baseband_signal_I_new = np.real(payload_corrected)
        baseband_signal_Q_new = np.imag(payload_corrected)
        
        print("baseband_signal_I_new")
        print(baseband_signal_I_new[10000:10020])



        #correct for the DC offset through AC coupling
##        if(scheme == "QPSK"):
##                baseband_signal_I_new = baseband_signal_I_new - np.mean(baseband_signal_I_new)
##                baseband_signal_Q_new = baseband_signal_Q_new - np.mean(baseband_signal_Q_new)

#Matched filtering
        if(pulse_shaping == "On"):

                #Use Matched filter receiver realization when pulse shaping is applied on the transmitter side

                alpha = 0.9 #roll-off factor of the RRC matched-filter

                L = 8

                #frame sync error detected
                if(len(baseband_signal_I_new)==0):
                    print("frame sync eror")
                    return 0,0

                symbols_I = matched_filtering.matched_filtering(baseband_signal_I_new, samples_perbit, fs_in, alpha, L, pulse_shape)
                print("symbols_I")
                print(symbols_I[10000:10020])
                symbols_Q = matched_filtering.matched_filtering(baseband_signal_Q_new, samples_perbit, fs_in, alpha, L, pulse_shape)
            
                #Remove the L/2 samples at the beginning and L/2 samples at the end caused by group delay of the filter

                #symbols_I = symbols_I[int(L/2):(len(symbols_I)-int(L/2))]

                #symbols_Q = symbols_Q[int(L/2):(len(symbols_Q)-int(L/2))]
                
                if (pulse_shape == 'rect'):
                    symbols_I = symbols_I[1:]
                    symbols_Q = symbols_Q[1:]
                elif (pulse_shape == 'rrc'):
                    symbols_I = symbols_I[L-1:(len(symbols_I))-L]
                    symbols_Q = symbols_Q[L-1:(len(symbols_Q))-L]

                #Plot the received signal constellation
                if(showplot==True):
##                        plt.ylim((-0.8, 0.8))
##                        plt.xlim((-0.8, 0.8))
                        for i in range(len(symbols_I)):
                                plt.plot(symbols_I[i],symbols_Q[i], color='blue', marker='o', markersize=1)
                        plt.show()

                channel_gain = np.max(symbols_I)/1.0

                if(scheme == "QPSK"):
                        channel_gain = np.mean(np.abs(symbols_I))/(2/3)
                if(scheme == "QAM"):
                        channel_gain = np.mean(np.abs(symbols_I))/(3/2)
                    

                symbols = np.array([symbols_I,symbols_Q])

                return symbols, channel_gain
