import numpy as np
import CRC
import matplotlib.pyplot as plt
from scipy import signal

def preamble_generator(): #"a" is the bit array to be modulated


        #The preamble we see in previous labs consists of 20 ones followed by 20 zeros followed by 40 ones and then
        #followed by 20 zeros:

        preamble = np.array([1,0,1,0])
        preamble = np.append(preamble, np.ones(6))                       #np.ones(N) generates a numpy array containing N ones.
        preamble = np.append(preamble, np.zeros(10))   #np.append(a,b) appends array b after array a
        preamble = np.append(preamble, np.ones(40))
        preamble = np.append(preamble, np.zeros(20))
##        preamble = np.append(preamble, np.ones(60))
 
        
        return preamble

def packet_generator(data, key):
        
        preamble = preamble_generator()
        packet = np.append(preamble, data)
        packet = CRC.encodeData(packet, key)

        return packet

# Test
data = np.array([1,0,1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,1,1,0,1,1,0,1,1,1])
key = np.array([1,1,1,0,0,1,0,1])

packet = packet_generator(data, key)
print("Packet: ")
print(packet)

error = CRC.CRCcheck(packet, key)
print("Error: " + str(error))
