import preamble_test
import symbol_mod_QAM
import symbol_demod_QAM
import numpy as np
import preamble_generator
import data_source
import mode_preconfiguration
import pulse_shaping
import matched_filtering
import time
import matplotlib.pyplot as plt

fs = 750000.0      
Ts = 1.0 / fs          # sampling period in seconds
f0 = 0.0               # homodyne (0 HZ IF)
M = int(input("Oversampling factor: "))                 # oversampling factor
T = M*Ts               # symbol period in seconds
Rs = 1/T               # symbol rate
segment_size = 1504    # One Transport Stream (TS) packet=188Bytes=1504 bits
R = int(input("Packet Ratio: "))                # Packet Ratio: number of segments contained in our larger OOK packet 
N = R*segment_size     # OOK Packet Length (equals R* segment_size)
L = 8                  # pulse shaping filter length
pulse_shape = input("Pulse Shape ('rect' for OOK or QPSK, 'rrc' for QAM): ") # type of pulse shaping
scheme = input("Modulation Scheme ('OOK', 'QPSK', or 'QAM'): ")        # Modulation scheme
alpha = 0.9 #roll-off factor of the RRC matched-filter
preamble_ones_length = int(input("Preamble ones length: "))
Mode = 4
b = "40k"             #bandwidth of Iperf client (UDP test)
test_packet_num = 20  #number of ADALM packets transmitted in Iperf test

#############################################
# Generate Packet

serverSock, generated_sequence, sequence_counter, l = mode_preconfiguration.tx_mode_preconfig(Mode, R, segment_size, N, b, test_packet_num)

preamble = preamble_generator.preamble_generator(preamble_ones_length)
#preamble = np.zeros(200)

Bits = data_source.data_source(Mode, serverSock, generated_sequence, sequence_counter, l)

packet_bits = np.append(preamble, Bits)

preamble_length = len(preamble)

total_samples, samples_perbit, fs_in, Ts_in =  mode_preconfiguration.rx_mode_preconfig(scheme,N,preamble_length,M,fs)

################################
# Generate Preamble check for matched filter
known_preamble_bits = preamble_generator.preamble_generator(preamble_ones_length)
known_preamble_symbols = symbol_mod_QAM.symbol_mod(known_preamble_bits, "OOK", len(known_preamble_bits))
known_preamble = np.abs(pulse_shaping.pulse_shaping(known_preamble_symbols, samples_perbit, fs_in, 'rect', 0.9, 8))

######################
# Generate symbols
start = time.time()
#print("1")
baseband_symbols = symbol_mod_QAM.symbol_mod(packet_bits, scheme, preamble_length)

######################
# Pulse Shaping
baseband = pulse_shaping.pulse_shaping(baseband_symbols, M, fs, pulse_shape, 0.9, 8)

bb_amplitude = 1
buffer = bb_amplitude*baseband

buffer_IQ = buffer

######################
# Separate Data Streams
I_data = np.real(buffer_IQ)
Q_data = np.imag(buffer_IQ)

a_in1 = buffer

#############################
# Simulate "Phase Sync"

oldN = N
if(scheme == "QPSK"):
        N = int(N/2)

if(scheme == "QAM"):
        N = int(N/4)

payload_start = len(known_preamble)

payload_before_correction = a_in1[payload_start:(payload_start + N*samples_perbit)]

ones_length = preamble_length - 20

payload_and_ones = a_in1[int(payload_start-ones_length*samples_perbit):(payload_start + N*samples_perbit)]

### In rf tx/rx scripts, phase sync goes here, but no need for phase sync in simulation

payload_corrected = payload_and_ones[ones_length*samples_perbit:]

baseband_signal_I_new = np.real(payload_corrected)
baseband_signal_Q_new = np.imag(payload_corrected)

N = oldN

######################################
# Matched filtering

symbols_I = matched_filtering.matched_filtering(baseband_signal_I_new, samples_perbit, fs_in, alpha, L, pulse_shape)

symbols_Q = matched_filtering.matched_filtering(baseband_signal_Q_new, samples_perbit, fs_in, alpha, L, pulse_shape)

#Remove the samples at the beginning and samples at the end caused by group delay of the filter

if (pulse_shape == 'rect'):
    symbols_I = symbols_I[1:]
    symbols_Q = symbols_Q[1:]
elif (pulse_shape == 'rrc'):
    symbols_I = symbols_I[L-1:(len(symbols_I))-L]
    symbols_Q = symbols_Q[L-1:(len(symbols_Q))-L]

symbols = np.array([symbols_I,symbols_Q])

######################################
# Calculate Channel Gain

if (scheme == 'QPSK'):  
    channel_gain = np.mean(np.abs(symbols_I))/(2/3)
elif (scheme == 'QAM'):  
    channel_gain = np.mean(np.abs(symbols_I))/(2)
else:
    channel_gain = np.mean(np.abs(symbols_I))

######################################
# Demodulation

a_demodulated = symbol_demod_QAM.symbol_demod(symbols, scheme, channel_gain, preamble_length)

#print("Demodulated bits match generated bits =")
if (pulse_shape == 'rrc'):
    if(np.array_equal(Bits[28:-32], a_demodulated, equal_nan=False)):
        print("BER:" + str(0.0))
        #BER = 0
    else:
        print("Non-zero BER")
        #BER = "non-zero"
else:
    if (np.array_equal(packet_bits[200:], a_demodulated, equal_nan=False)):
        print("BER:" + str(0.0))
        #BER = 0
    else:
        print("Non-zero BER")
        #BER = "non-zero"


######################################
# Throughput calculation

packet_kbits = len(Bits)
total_kbits = 0
total_kbits = total_kbits+packet_kbits

throughput = total_kbits/(time.time()-start)
#print("2")

#throughput = total_kbits

print("Throughput = ")
print(throughput)

showplot = True
if(showplot==True):
        #plt.ylim((-2, 1))
        #plt.xlim((1, 2))
        for i in range(len(symbols_I)):
                plt.plot(symbols_I[i]/channel_gain,symbols_Q[i]/channel_gain, color='blue', marker='o', markersize=1)
        plt.xlabel('I')
        plt.ylabel('Q')
        plt.title('IQ Plot for ' + scheme)
        plt.show()