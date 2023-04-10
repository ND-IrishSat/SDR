# +
# #%%writefile symbol_mod.py
# Run this cell to generate convert and save the script into a Python Program. 
# -

import numpy as np
import matplotlib.pyplot as plt
from scipy import signal

# Define the symbol_mod module

# # The symbol_mod module takes the following arguments as inputs:

# ## packet_bits                 The bit array to be mapped into symbols (including both the preamble bits and the
#                             payload bits)

# ## scheme                      A string indicating which scheme is adopted (e.g.: "OOK", "QPSK")

# ## preamble_length             Length of the preamble (in bits)

# # The symbol_mod function returns the following argument as output:

# ## baseband_symbols:           The baseband symbols obtained after mapping the bits

def symbol_mod(packet_bits, scheme, preamble_length): 

        if(scheme == 'OOK'):
                
                #print("hello")
                preamble = packet_bits[0:preamble_length]
                payload = packet_bits[preamble_length:len(packet_bits)]
                preamble_symbols = 1.0*preamble
                payload_symbols = 1.0*payload                
                baseband_symbols = np.append(preamble_symbols,payload_symbols)

        if(scheme == 'QPSK'):


################ if preamble different modulation scheme
                preamble = packet_bits[0:preamble_length]
                payload = packet_bits[preamble_length:len(packet_bits)]

# Setting the Preamble bits
                baseband_symbols_I = 1.0*preamble
                baseband_symbols_Q = np.zeros(preamble_length)
################
                #payload = packet_bits

# Imagine the two message bits arriving as being independently generated,\
# now split the payload into two to generate two message streams for I and Q channels.
                I_bits = payload[0:int(len(payload)/2)]
                Q_bits = payload[int(len(payload)/2):len(payload)]
# Now modulate the two streams separately according to the BPSK modulation scheme,\
# do not forget to scale the symbols appropriately to maintain the symbol energy.
                I_symbols = np.zeros(int(len(payload)/2))
                Q_symbols = np.zeros(int(len(payload)/2))
    
                for i in range(int(len(payload)/2)):
                    if I_bits[i] == 0:
                        I_symbols[i] = -1
                    elif I_bits[i] == 1:
                        I_symbols[i] = 1
                        
                #I_symbols = I_symbols/np.sqrt(2)


                for i in range(int(len(payload)/2)):
                    if Q_bits[i] == 0:
                        Q_symbols[i] = -1
                    elif Q_bits[i] == 1:
                        Q_symbols[i] = 1
                        
                #Q_symbols = Q_symbols/np.sqrt(2)
                preamble_symbols = baseband_symbols_I + 1j*baseband_symbols_Q
                print("1")
                data_symbols = I_symbols + 1j*Q_symbols
                print("2")
                
                # Similarly set the data symbols


# Scale QPSK payload to have the same per channel average signal energy as OOK
                data_symbols = data_symbols/np.sqrt(2)

# Add the preamble to the payload to create a packet               
                baseband_symbols = np.append(preamble_symbols, data_symbols)
                #baseband_symbols = data_symbols
                
        if(scheme == 'QAM'):

################ if preamble different modulation scheme
                preamble = packet_bits[0:preamble_length]
                payload = packet_bits[preamble_length:len(packet_bits)]

# Setting the Preamble bits
                baseband_symbols_I = 1.0*preamble
                baseband_symbols_Q = np.zeros(preamble_length)
################
                #payload = packet_bits[preamble_length:]

# #Imagine the two message bits arriving as being independently generated,\
# #now split the payload into two to generate two message streams for I and Q channels.

                I_bits = []
                Q_bits = []
                n = 0
                
                # ensure payload consists of full QAM symbols
                if ((len(payload) % 4) != 0):
                    n = len(payload) % 4
                    payload = payload[:-n]
                    
                #print("len payload after correction")
                #print(len(payload))
                
                i = 0
                while (i < len(payload)):
                    for j in range(2):
                        Q_bits.append(payload[i+j])
                    for k in [2,3]:
                        I_bits.append(payload[i+k])
                    i = i+4
                    
                #print("len I and Q bits after split: ")
                #print(len(I_bits))
                #print(len(Q_bits))
                               
                #Q_bits = payload[0:int(len(payload)/2)-1]
                #I_bits = payload[int(len(payload)/2):len(payload)]

# Now modulate the two streams separately according to the QAM modulation scheme,\
# do not forget to scale the symbols appropriately to maintain the symbol energy.
                I_symbols = np.zeros(int(len(I_bits)/2))
                Q_symbols = np.zeros(int(len(Q_bits)/2))
    
                for i in range(0,len(I_bits)-1,2):    
                    if I_bits[i] == 0 and I_bits[i+1] == 0:
                        I_symbols[int(i/2)] = -3
                    elif I_bits[i] == 0 and I_bits[i+1] == 1:
                        I_symbols[int(i/2)] = -1
                    elif I_bits[i] == 1 and I_bits[i+1] == 0:
                        I_symbols[int(i/2)] = 3
                    elif I_bits[i] == 1 and I_bits[i+1] == 1:
                        I_symbols[int(i/2)] = 1
                        
                #I_symbols = I_symbols/np.sqrt(2)

                for i in range(0,len(Q_bits)-1,2):
                    if Q_bits[i] == 0 and Q_bits[i+1] == 0:
                        Q_symbols[int(i/2)] = -3
                    elif Q_bits[i] == 0 and Q_bits[i+1] == 1:
                        Q_symbols[int(i/2)] = -1
                    elif Q_bits[i] == 1 and Q_bits[i+1] == 0:
                        Q_symbols[int(i/2)] = 3
                    elif Q_bits[i] == 1 and Q_bits[i+1] == 1:
                        Q_symbols[int(i/2)] = 1

                #Q_symbols = Q_symbols/np.sqrt(2)
                preamble_symbols = baseband_symbols_I + 1j*baseband_symbols_Q
                #print("len preamble symbols: ")
                #print(len(preamble_symbols))
                
                data_symbols = I_symbols + 1j*Q_symbols
                #print("len data symbols: ")
                #print(len(data_symbols))

# Scale QAM payload to have the same per channel average signal energy as OOK
                #data_symbols = data_symbols/np.sqrt(2)

# Add the preamble to the payload to create a packet               
                #baseband_symbols = data_symbols
                baseband_symbols = np.append(preamble_symbols, data_symbols)
                #baseband_symbols = np.reshape(baseband_symbols, (2, (int(len(baseband_symbols)/2))))
                #print("baseband_symbols: ")
                #print(baseband_symbols)
                
        return baseband_symbols
