"""
Packet generation convenience functions

"""

import numpy as np
import CRC
import matplotlib.pyplot as plt
from scipy import signal

def preamble_generator(): #"a" is the bit array to be modulated

        preamble = np.array([0,1,0,0,0,0,1,1,0,0,0,1,0,1,0,0,1,1,1,1,0,1,0,0,0,1,1,1,0,0,1,0,0,1,0,1,1,0,1,1,1,0,1,1,0,0,1,1,0,1,0,1,0,1,1,1,1,1,1,0]).astype(float) # optimal periodic binary code for N = 63 https://ntrs.nasa.gov/citations/19800017860
        
        return preamble

def packet_generator(data, key):
        
        # preamble = preamble_generator(preamble_ones_length)
        # packet = np.append(preamble, data)
        # packet = CRC.encodeData(packet, key)

        preamble = preamble_generator()
        data = CRC.encodeData(data, key)
        packet = np.append(preamble, data)

        return packet

# Test
data = np.array([1,0,1,1,0,0,0,0,1,1,1,1,1,1,0,0,0,1,1,0,1,1,0,1,1,1])
key = np.array([1,1,1,0,0,1,0,1])
preamble_ones_length = 60

packet = packet_generator(data, key, preamble_ones_length)
# print("Packet: ")
# print(packet)

error = CRC.CRCcheck(packet, key)
# print("Error: " + str(error))