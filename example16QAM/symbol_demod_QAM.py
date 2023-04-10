import numpy as np
import matplotlib.pyplot as plt
from scipy import signal
import time

# Define the symbol_mod module

# The symbol_demod module takes the following arguments as inputs:

# baseband symbols:            The symbols array to be mapped into bits

# scheme:                      A string indicating which scheme is adopted (e.g.: "OOK", "QPSK")

# channel_gain                 This is the gain of the channel, channel impulse response is simply modeled as g(t)= 
# channel_gain*delta(t)

# The symbol_demod module returns the following argument as output:

# a_demodulated:              An array containing the demodulated bits


def symbol_demod(baseband_symbols, scheme, channel_gain, preamble_len): #"a" is the bit array to be modulate

        a_demodulated = []
        
        if(scheme == 'OOK'):

                s_on = 1.0 * channel_gain

                s_off = 0* channel_gain

                baseband_symbols_I = baseband_symbols[0]

                baseband_symbols_Q = baseband_symbols[1]


                baseband_symbols_complex = baseband_symbols_I + 1j * baseband_symbols_Q


                for i in range( len(baseband_symbols_complex) ):

                        #Coherent: finding the minimum distance between the received symbol and all reference symbols 
                        # in the constellation plot.

                                if (np.abs(baseband_symbols_complex[i] - s_on) < np.abs(baseband_symbols_complex[i] - s_off)):
                                       a_demodulated.append(1)
                                else:
                                       a_demodulated.append(0)
 
        if(scheme == 'BPSK'):

                baseband_symbols_I = baseband_symbols[0]

                baseband_symbols_Q = baseband_symbols[1]

                reference_plus = 1.0*channel_gain

                reference_minus = -reference_plus

                baseband_symbols_complex = baseband_symbols_I + 1j * baseband_symbols_Q

                for i in range(len(baseband_symbols_complex)):

                    #Find the minimum distance between the received symbol and all reference symbols in the constellation plot.

                        if (np.abs(baseband_symbols_complex[i] - reference_plus) < np.abs(baseband_symbols_complex[i] - \
                                                                                          reference_minus)):
                               a_demodulated.append(0)
                        else:
                               a_demodulated.append(1)


        if(scheme == 'QPSK'):
                
                
                #if preamble not modulated
#                preamble_symbols = np.real(baseband_symbols[0,:preamble_len])
#                preamble_demod = preamble_symbols
#                
#                baseband_symbols = np.array([baseband_symbols[0,preamble_len:],baseband_symbols[1,preamble_len:]])
                
                
                
                #may need to move this line
                #baseband_symbols = np.reshape(baseband_symbols, (2, (int(len(baseband_symbols)/2))))

                baseband_symbols_I = baseband_symbols[0,:]

                baseband_symbols_Q = baseband_symbols[1,:]


                I_demodulated = []
                Q_demodulated = []

                #Construct the received symbols on the complex plane (signal space constellation)
                baseband_symbols_complex = baseband_symbols_I + 1j * baseband_symbols_Q

                #Compute and define the reference signals in the signal space (4 constellation points)
                reference_00 = -1.0*channel_gain  -1j* channel_gain

                reference_11 = 1.0*channel_gain + 1j* channel_gain

                reference_01 = -1.0*channel_gain + 1j* channel_gain

                reference_10 = 1.0*channel_gain  -1j* channel_gain

                #Start a for-loop to iterate through all complex symbols and make a decision on 2-bits of data 
                for i in range(len(baseband_symbols_complex)):

                        symbol = baseband_symbols_complex[i]

 #Find the minimum distance between the received symbol and all symbols in the constellation plot.

                        if (  np.abs(symbol - reference_11) == np.amin(  [np.abs(symbol - reference_11),  \
                                                                          np.abs(symbol - reference_10), \
                                                                          np.abs(symbol - reference_00), \
                                                                          np.abs(symbol - reference_01)]  )  ):
                               I_demodulated.append(1)
                               Q_demodulated.append(1)
                        elif(  np.abs(symbol - reference_10) == np.amin(  [np.abs(symbol - reference_11),  \
                                                                          np.abs(symbol - reference_10), \
                                                                          np.abs(symbol - reference_00), \
                                                                          np.abs(symbol - reference_01)]  )  ):
                               I_demodulated.append(1)
                               Q_demodulated.append(0)
                        elif(  np.abs(symbol - reference_00) == np.amin(  [np.abs(symbol - reference_11),  \
                                                                          np.abs(symbol - reference_10), \
                                                                          np.abs(symbol - reference_00), \
                                                                          np.abs(symbol - reference_01)]  )  ):
                               I_demodulated.append(0)
                               Q_demodulated.append(0)
                        elif(  np.abs(symbol - reference_01) == np.amin(  [np.abs(symbol - reference_11),  \
                                                                          np.abs(symbol - reference_10), \
                                                                          np.abs(symbol - reference_00), \
                                                                          np.abs(symbol - reference_01)]  )  ):
                               I_demodulated.append(0)
                               Q_demodulated.append(1)

                a_demodulated = np.append(I_demodulated, Q_demodulated)
                #return np.append(preamble_demod,a_demodulated)
                return a_demodulated
    
        if(scheme == 'QAM'):
                
                #if preamble not modulated
#                preamble_symbols = np.real(baseband_symbols[0,:preamble_len])
#                preamble_demod = preamble_symbols
#                
#                baseband_symbols = np.array([baseband_symbols[0,preamble_len:],baseband_symbols[1,preamble_len:]])
                
                #may need to move this line
                #baseband_symbols = np.reshape(baseband_symbols, (2, (int(len(baseband_symbols)/2))))

                baseband_symbols_I = baseband_symbols[0,:]
                #print("len(baseband_symbols[0])")
                #print(len(baseband_symbols[0]))

                baseband_symbols_Q = baseband_symbols[1,:]
                
                #baseband_symbols = np.reshape(baseband_symbols, (2, (int(len(baseband_symbols)/2))))
                ## Code to ensure that baseband symbols arrays are even and the same length
                #baseband_symbols_I = baseband_symbols[0,:]
                #baseband_symbols_I = baseband_symbols[0,int(len(baseband_symbols[0])/2)-1]
                #baseband_symbols_Q = baseband_symbols[1,:]

                #baseband_symbols_Q = baseband_symbols[int(len(baseband_symbols[1])/2),:]
                    
                if (len(baseband_symbols_I) != len(baseband_symbols_Q)):
                    
                    if (len(baseband_symbols_I) > len(baseband_symbols_Q)):
                        
                        delta = len(baseband_symbols_I) - len(baseband_symbols_Q)
                        baseband_symbols_I = baseband_symbols_I[:-delta]     
                        
                    elif (len(baseband_symbols_Q) > len(baseband_symbols_I)):
                        
                        delta = len(baseband_symbols_Q) - len(baseband_symbols_I)
                        baseband_symbols_Q = baseband_symbols_Q[:-delta]

                I_demodulated = []
                Q_demodulated = []

                #Construct the received symbols on the complex plane (signal space constellation)
                baseband_symbols_complex = baseband_symbols_I + 1j * baseband_symbols_Q
                
                #print("channel_gain: ")
                #print(channel_gain)

                #Compute and define the reference signals in the signal space (16 constellation points)
                reference_0000 = -3.0*channel_gain - 3j* channel_gain

                reference_0001 = -1.0*channel_gain - 3j* channel_gain

                reference_0010 = 3.0*channel_gain - 3j* channel_gain

                reference_0011 = 1.0*channel_gain - 3j* channel_gain
                
                reference_0100 = -3.0*channel_gain - 1j* channel_gain

                reference_0101 = -1.0*channel_gain - 1j* channel_gain

                reference_0110 = 3.0*channel_gain - 1j* channel_gain

                reference_0111 = 1.0*channel_gain - 1j* channel_gain
                
                reference_1000 = -3.0*channel_gain + 3j* channel_gain

                reference_1001 = -1.0*channel_gain + 3j* channel_gain

                reference_1010 = 3.0*channel_gain + 3j* channel_gain

                reference_1011 = 1.0*channel_gain + 3j* channel_gain
                
                reference_1100 = -3.0*channel_gain + 1j* channel_gain

                reference_1101 = -1.0*channel_gain + 1j* channel_gain

                reference_1110 = 3.0*channel_gain + 1j* channel_gain

                reference_1111 = 1.0*channel_gain + 1j* channel_gain

                #Start a for-loop to iterate through all complex symbols and make a decision on 4-bits of data 

                for i in range(len(baseband_symbols_complex)):

                        symbol = baseband_symbols_complex[i]


 #Find the minimum distance between the received symbol and all symbols in the constellation plot.

                        if (  np.abs(symbol - reference_1111) == np.amin(  [np.abs(symbol - reference_1111),  \
                                                                          np.abs(symbol - reference_1110), \
                                                                          np.abs(symbol - reference_1101), \
                                                                          np.abs(symbol - reference_1100), \
                                                                          \
                                                                          np.abs(symbol - reference_1011), \
                                                                          np.abs(symbol - reference_1010), \
                                                                          np.abs(symbol - reference_1001), \
                                                                          np.abs(symbol - reference_1000), \
                                                                          \
                                                                          np.abs(symbol - reference_0111), \
                                                                          np.abs(symbol - reference_0110), \
                                                                          np.abs(symbol - reference_0101), \
                                                                          np.abs(symbol - reference_0100), \
                                                                          \
                                                                          np.abs(symbol - reference_0011), \
                                                                          np.abs(symbol - reference_0010), \
                                                                          np.abs(symbol - reference_0001), \
                                                                          np.abs(symbol - reference_0000)]  )  ):
                               Q_demodulated.append(1)
                               Q_demodulated.append(1)
                               I_demodulated.append(1)
                               I_demodulated.append(1)
                        
                               
                        elif (  np.abs(symbol - reference_1110) == np.amin(  [np.abs(symbol - reference_1111),  \
                                                                          np.abs(symbol - reference_1110), \
                                                                          np.abs(symbol - reference_1101), \
                                                                          np.abs(symbol - reference_1100), \
                                                                          \
                                                                          np.abs(symbol - reference_1011), \
                                                                          np.abs(symbol - reference_1010), \
                                                                          np.abs(symbol - reference_1001), \
                                                                          np.abs(symbol - reference_1000), \
                                                                          \
                                                                          np.abs(symbol - reference_0111), \
                                                                          np.abs(symbol - reference_0110), \
                                                                          np.abs(symbol - reference_0101), \
                                                                          np.abs(symbol - reference_0100), \
                                                                          \
                                                                          np.abs(symbol - reference_0011), \
                                                                          np.abs(symbol - reference_0010), \
                                                                          np.abs(symbol - reference_0001), \
                                                                          np.abs(symbol - reference_0000)]  )  ):
                               Q_demodulated.append(1)
                               Q_demodulated.append(1)
                               I_demodulated.append(1)
                               I_demodulated.append(0)
                               
                        elif (  np.abs(symbol - reference_1101) == np.amin(  [np.abs(symbol - reference_1111),  \
                                                                          np.abs(symbol - reference_1110), \
                                                                          np.abs(symbol - reference_1101), \
                                                                          np.abs(symbol - reference_1100), \
                                                                          \
                                                                          np.abs(symbol - reference_1011), \
                                                                          np.abs(symbol - reference_1010), \
                                                                          np.abs(symbol - reference_1001), \
                                                                          np.abs(symbol - reference_1000), \
                                                                          \
                                                                          np.abs(symbol - reference_0111), \
                                                                          np.abs(symbol - reference_0110), \
                                                                          np.abs(symbol - reference_0101), \
                                                                          np.abs(symbol - reference_0100), \
                                                                          \
                                                                          np.abs(symbol - reference_0011), \
                                                                          np.abs(symbol - reference_0010), \
                                                                          np.abs(symbol - reference_0001), \
                                                                          np.abs(symbol - reference_0000)]  )  ):
                               Q_demodulated.append(1)
                               Q_demodulated.append(1)
                               I_demodulated.append(0)
                               I_demodulated.append(1)
                               
                        elif (  np.abs(symbol - reference_1100) == np.amin(  [np.abs(symbol - reference_1111),  \
                                                                          np.abs(symbol - reference_1110), \
                                                                          np.abs(symbol - reference_1101), \
                                                                          np.abs(symbol - reference_1100), \
                                                                          \
                                                                          np.abs(symbol - reference_1011), \
                                                                          np.abs(symbol - reference_1010), \
                                                                          np.abs(symbol - reference_1001), \
                                                                          np.abs(symbol - reference_1000), \
                                                                          \
                                                                          np.abs(symbol - reference_0111), \
                                                                          np.abs(symbol - reference_0110), \
                                                                          np.abs(symbol - reference_0101), \
                                                                          np.abs(symbol - reference_0100), \
                                                                          \
                                                                          np.abs(symbol - reference_0011), \
                                                                          np.abs(symbol - reference_0010), \
                                                                          np.abs(symbol - reference_0001), \
                                                                          np.abs(symbol - reference_0000)]  )  ):
                               Q_demodulated.append(1)
                               Q_demodulated.append(1)
                               I_demodulated.append(0)
                               I_demodulated.append(0)
                               
                        elif (  np.abs(symbol - reference_1011) == np.amin(  [np.abs(symbol - reference_1111),  \
                                                                          np.abs(symbol - reference_1110), \
                                                                          np.abs(symbol - reference_1101), \
                                                                          np.abs(symbol - reference_1100), \
                                                                          \
                                                                          np.abs(symbol - reference_1011), \
                                                                          np.abs(symbol - reference_1010), \
                                                                          np.abs(symbol - reference_1001), \
                                                                          np.abs(symbol - reference_1000), \
                                                                          \
                                                                          np.abs(symbol - reference_0111), \
                                                                          np.abs(symbol - reference_0110), \
                                                                          np.abs(symbol - reference_0101), \
                                                                          np.abs(symbol - reference_0100), \
                                                                          \
                                                                          np.abs(symbol - reference_0011), \
                                                                          np.abs(symbol - reference_0010), \
                                                                          np.abs(symbol - reference_0001), \
                                                                          np.abs(symbol - reference_0000)]  )  ):
                               Q_demodulated.append(1)
                               Q_demodulated.append(0)
                               I_demodulated.append(1)
                               I_demodulated.append(1)
                        
                               
                        elif (  np.abs(symbol - reference_1010) == np.amin(  [np.abs(symbol - reference_1111),  \
                                                                          np.abs(symbol - reference_1110), \
                                                                          np.abs(symbol - reference_1101), \
                                                                          np.abs(symbol - reference_1100), \
                                                                          \
                                                                          np.abs(symbol - reference_1011), \
                                                                          np.abs(symbol - reference_1010), \
                                                                          np.abs(symbol - reference_1001), \
                                                                          np.abs(symbol - reference_1000), \
                                                                          \
                                                                          np.abs(symbol - reference_0111), \
                                                                          np.abs(symbol - reference_0110), \
                                                                          np.abs(symbol - reference_0101), \
                                                                          np.abs(symbol - reference_0100), \
                                                                          \
                                                                          np.abs(symbol - reference_0011), \
                                                                          np.abs(symbol - reference_0010), \
                                                                          np.abs(symbol - reference_0001), \
                                                                          np.abs(symbol - reference_0000)]  )  ):
                               Q_demodulated.append(1)
                               Q_demodulated.append(0)
                               I_demodulated.append(1)
                               I_demodulated.append(0)
                               
                        elif (  np.abs(symbol - reference_1001) == np.amin(  [np.abs(symbol - reference_1111),  \
                                                                          np.abs(symbol - reference_1110), \
                                                                          np.abs(symbol - reference_1101), \
                                                                          np.abs(symbol - reference_1100), \
                                                                          \
                                                                          np.abs(symbol - reference_1011), \
                                                                          np.abs(symbol - reference_1010), \
                                                                          np.abs(symbol - reference_1001), \
                                                                          np.abs(symbol - reference_1000), \
                                                                          \
                                                                          np.abs(symbol - reference_0111), \
                                                                          np.abs(symbol - reference_0110), \
                                                                          np.abs(symbol - reference_0101), \
                                                                          np.abs(symbol - reference_0100), \
                                                                          \
                                                                          np.abs(symbol - reference_0011), \
                                                                          np.abs(symbol - reference_0010), \
                                                                          np.abs(symbol - reference_0001), \
                                                                          np.abs(symbol - reference_0000)]  )  ):
                               Q_demodulated.append(1)
                               Q_demodulated.append(0)
                               I_demodulated.append(0)
                               I_demodulated.append(1)
                               
                        elif (  np.abs(symbol - reference_1000) == np.amin(  [np.abs(symbol - reference_1111),  \
                                                                          np.abs(symbol - reference_1110), \
                                                                          np.abs(symbol - reference_1101), \
                                                                          np.abs(symbol - reference_1100), \
                                                                          \
                                                                          np.abs(symbol - reference_1011), \
                                                                          np.abs(symbol - reference_1010), \
                                                                          np.abs(symbol - reference_1001), \
                                                                          np.abs(symbol - reference_1000), \
                                                                          \
                                                                          np.abs(symbol - reference_0111), \
                                                                          np.abs(symbol - reference_0110), \
                                                                          np.abs(symbol - reference_0101), \
                                                                          np.abs(symbol - reference_0100), \
                                                                          \
                                                                          np.abs(symbol - reference_0011), \
                                                                          np.abs(symbol - reference_0010), \
                                                                          np.abs(symbol - reference_0001), \
                                                                          np.abs(symbol - reference_0000)]  )  ):
                               Q_demodulated.append(1)
                               Q_demodulated.append(0)
                               I_demodulated.append(0)
                               I_demodulated.append(0)
                               
                        elif (  np.abs(symbol - reference_0111) == np.amin(  [np.abs(symbol - reference_1111),  \
                                                                          np.abs(symbol - reference_1110), \
                                                                          np.abs(symbol - reference_1101), \
                                                                          np.abs(symbol - reference_1100), \
                                                                          \
                                                                          np.abs(symbol - reference_1011), \
                                                                          np.abs(symbol - reference_1010), \
                                                                          np.abs(symbol - reference_1001), \
                                                                          np.abs(symbol - reference_1000), \
                                                                          \
                                                                          np.abs(symbol - reference_0111), \
                                                                          np.abs(symbol - reference_0110), \
                                                                          np.abs(symbol - reference_0101), \
                                                                          np.abs(symbol - reference_0100), \
                                                                          \
                                                                          np.abs(symbol - reference_0011), \
                                                                          np.abs(symbol - reference_0010), \
                                                                          np.abs(symbol - reference_0001), \
                                                                          np.abs(symbol - reference_0000)]  )  ):
                               Q_demodulated.append(0)
                               Q_demodulated.append(1)
                               I_demodulated.append(1)
                               I_demodulated.append(1)
                        
                               
                        elif (  np.abs(symbol - reference_0110) == np.amin(  [np.abs(symbol - reference_1111),  \
                                                                          np.abs(symbol - reference_1110), \
                                                                          np.abs(symbol - reference_1101), \
                                                                          np.abs(symbol - reference_1100), \
                                                                          \
                                                                          np.abs(symbol - reference_1011), \
                                                                          np.abs(symbol - reference_1010), \
                                                                          np.abs(symbol - reference_1001), \
                                                                          np.abs(symbol - reference_1000), \
                                                                          \
                                                                          np.abs(symbol - reference_0111), \
                                                                          np.abs(symbol - reference_0110), \
                                                                          np.abs(symbol - reference_0101), \
                                                                          np.abs(symbol - reference_0100), \
                                                                          \
                                                                          np.abs(symbol - reference_0011), \
                                                                          np.abs(symbol - reference_0010), \
                                                                          np.abs(symbol - reference_0001), \
                                                                          np.abs(symbol - reference_0000)]  )  ):
                               Q_demodulated.append(0)
                               Q_demodulated.append(1)
                               I_demodulated.append(1)
                               I_demodulated.append(0)
                               
                        elif (  np.abs(symbol - reference_0101) == np.amin(  [np.abs(symbol - reference_1111),  \
                                                                          np.abs(symbol - reference_1110), \
                                                                          np.abs(symbol - reference_1101), \
                                                                          np.abs(symbol - reference_1100), \
                                                                          \
                                                                          np.abs(symbol - reference_1011), \
                                                                          np.abs(symbol - reference_1010), \
                                                                          np.abs(symbol - reference_1001), \
                                                                          np.abs(symbol - reference_1000), \
                                                                          \
                                                                          np.abs(symbol - reference_0111), \
                                                                          np.abs(symbol - reference_0110), \
                                                                          np.abs(symbol - reference_0101), \
                                                                          np.abs(symbol - reference_0100), \
                                                                          \
                                                                          np.abs(symbol - reference_0011), \
                                                                          np.abs(symbol - reference_0010), \
                                                                          np.abs(symbol - reference_0001), \
                                                                          np.abs(symbol - reference_0000)]  )  ):
                               Q_demodulated.append(0)
                               Q_demodulated.append(1)
                               I_demodulated.append(0)
                               I_demodulated.append(1)
                               
                        elif (  np.abs(symbol - reference_0100) == np.amin(  [np.abs(symbol - reference_1111),  \
                                                                          np.abs(symbol - reference_1110), \
                                                                          np.abs(symbol - reference_1101), \
                                                                          np.abs(symbol - reference_1100), \
                                                                          \
                                                                          np.abs(symbol - reference_1011), \
                                                                          np.abs(symbol - reference_1010), \
                                                                          np.abs(symbol - reference_1001), \
                                                                          np.abs(symbol - reference_1000), \
                                                                          \
                                                                          np.abs(symbol - reference_0111), \
                                                                          np.abs(symbol - reference_0110), \
                                                                          np.abs(symbol - reference_0101), \
                                                                          np.abs(symbol - reference_0100), \
                                                                          \
                                                                          np.abs(symbol - reference_0011), \
                                                                          np.abs(symbol - reference_0010), \
                                                                          np.abs(symbol - reference_0001), \
                                                                          np.abs(symbol - reference_0000)]  )  ):
                               Q_demodulated.append(0)
                               Q_demodulated.append(1)
                               I_demodulated.append(0)
                               I_demodulated.append(0)
                               
                        elif (  np.abs(symbol - reference_0011) == np.amin(  [np.abs(symbol - reference_1111),  \
                                                                          np.abs(symbol - reference_1110), \
                                                                          np.abs(symbol - reference_1101), \
                                                                          np.abs(symbol - reference_1100), \
                                                                          \
                                                                          np.abs(symbol - reference_1011), \
                                                                          np.abs(symbol - reference_1010), \
                                                                          np.abs(symbol - reference_1001), \
                                                                          np.abs(symbol - reference_1000), \
                                                                          \
                                                                          np.abs(symbol - reference_0111), \
                                                                          np.abs(symbol - reference_0110), \
                                                                          np.abs(symbol - reference_0101), \
                                                                          np.abs(symbol - reference_0100), \
                                                                          \
                                                                          np.abs(symbol - reference_0011), \
                                                                          np.abs(symbol - reference_0010), \
                                                                          np.abs(symbol - reference_0001), \
                                                                          np.abs(symbol - reference_0000)]  )  ):
                               Q_demodulated.append(0)
                               Q_demodulated.append(0)
                               I_demodulated.append(1)
                               I_demodulated.append(1)
                        
                               
                        elif (  np.abs(symbol - reference_0010) == np.amin(  [np.abs(symbol - reference_1111),  \
                                                                          np.abs(symbol - reference_1110), \
                                                                          np.abs(symbol - reference_1101), \
                                                                          np.abs(symbol - reference_1100), \
                                                                          \
                                                                          np.abs(symbol - reference_1011), \
                                                                          np.abs(symbol - reference_1010), \
                                                                          np.abs(symbol - reference_1001), \
                                                                          np.abs(symbol - reference_1000), \
                                                                          \
                                                                          np.abs(symbol - reference_0111), \
                                                                          np.abs(symbol - reference_0110), \
                                                                          np.abs(symbol - reference_0101), \
                                                                          np.abs(symbol - reference_0100), \
                                                                          \
                                                                          np.abs(symbol - reference_0011), \
                                                                          np.abs(symbol - reference_0010), \
                                                                          np.abs(symbol - reference_0001), \
                                                                          np.abs(symbol - reference_0000)]  )  ):
                               Q_demodulated.append(0)
                               Q_demodulated.append(0)
                               I_demodulated.append(1)
                               I_demodulated.append(0)
                               
                        elif (  np.abs(symbol - reference_0001) == np.amin(  [np.abs(symbol - reference_1111),  \
                                                                          np.abs(symbol - reference_1110), \
                                                                          np.abs(symbol - reference_1101), \
                                                                          np.abs(symbol - reference_1100), \
                                                                          \
                                                                          np.abs(symbol - reference_1011), \
                                                                          np.abs(symbol - reference_1010), \
                                                                          np.abs(symbol - reference_1001), \
                                                                          np.abs(symbol - reference_1000), \
                                                                          \
                                                                          np.abs(symbol - reference_0111), \
                                                                          np.abs(symbol - reference_0110), \
                                                                          np.abs(symbol - reference_0101), \
                                                                          np.abs(symbol - reference_0100), \
                                                                          \
                                                                          np.abs(symbol - reference_0011), \
                                                                          np.abs(symbol - reference_0010), \
                                                                          np.abs(symbol - reference_0001), \
                                                                          np.abs(symbol - reference_0000)]  )  ):
                               Q_demodulated.append(0)
                               Q_demodulated.append(0)
                               I_demodulated.append(0)
                               I_demodulated.append(1)
                               
                        elif (  np.abs(symbol - reference_0000) == np.amin(  [np.abs(symbol - reference_1111),  \
                                                                          np.abs(symbol - reference_1110), \
                                                                          np.abs(symbol - reference_1101), \
                                                                          np.abs(symbol - reference_1100), \
                                                                          \
                                                                          np.abs(symbol - reference_1011), \
                                                                          np.abs(symbol - reference_1010), \
                                                                          np.abs(symbol - reference_1001), \
                                                                          np.abs(symbol - reference_1000), \
                                                                          \
                                                                          np.abs(symbol - reference_0111), \
                                                                          np.abs(symbol - reference_0110), \
                                                                          np.abs(symbol - reference_0101), \
                                                                          np.abs(symbol - reference_0100), \
                                                                          \
                                                                          np.abs(symbol - reference_0011), \
                                                                          np.abs(symbol - reference_0010), \
                                                                          np.abs(symbol - reference_0001), \
                                                                          np.abs(symbol - reference_0000)]  )  ):
                               Q_demodulated.append(0)
                               Q_demodulated.append(0)
                               I_demodulated.append(0)
                               I_demodulated.append(0)

                a_demodulated = []  
                
                #print("len(I_demodulated)")
                #print(len(I_demodulated))
               
                i = 0
                while (i < (len(Q_demodulated))):
                    
                    a_demodulated.append(Q_demodulated[i])
                    a_demodulated.append(Q_demodulated[i+1])
                    a_demodulated.append(I_demodulated[i])
                    a_demodulated.append(I_demodulated[i+1])
                    i = i+2
                    
                #a_demodulated = np.append(I_demodulated, Q_demodulated)


        #return np.append(preamble_demod,a_demodulated)
        return a_demodulated


# Helper function

# Rotating a vector in constellation diagram:

# Take a complex number (2-D vector) as input (referenced through x, y coordinate) and return the coordinate of the 
# rotated vector (by angle radians)

def rotate(vector, angle):

        x = np.real(vector)

        y = np.imag(vector)

        x_new = np.cos(angle)*x - np.sin(angle)*y

        y_new = np.sin(angle)*x + np.cos(angle)*y

        vector_new = x_new + 1j*y_new
        

        return vector_new

# +
#qpskdemod_data = np.load('SymbolDemodResult.npz') #load

# +
#baseband_symbols = qpskdemod_data['baseband_symbols']
#scheme = qpskdemod_data['scheme']
#channel_gain = qpskdemod_data['channel_gain']
#demod_bits = qpskdemod_data['a_demodulated']

# +
#demodulated_bits = symbol_demod(baseband_symbols, scheme, channel_gain)

# +
# compare the obtained demodulated bits with the desired demodulated bits
#print(np.array_equal(demodulated_bits, demod_bits, equal_nan=False))
