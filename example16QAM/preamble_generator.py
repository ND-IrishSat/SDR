import numpy as np
import matplotlib.pyplot as plt
from scipy import signal

# Define the preamble_generator module

# The modulation module takes the following arguments as inputs:


# The modulation function returns the following argument as output:

# preamble:                  The preamble sequence in bits

#preamble that worked for BER mode
def preamble_generator(preamble_ones_length): #"a" is the bit array to be modulated


        #The preamble we see in previous labs consists of 20 ones followed by 20 zeros followed by 40 ones and then
        #followed by 20 zeros:

        preamble = np.array([1,0,1,0])
        preamble = np.append(preamble, np.ones(6))                       #np.ones(N) generates a numpy array containing N ones.
        preamble = np.append(preamble, np.zeros(10))   #np.append(a,b) appends array b after array a
##        preamble = np.append(preamble, np.ones(40))
##        preamble = np.append(preamble, np.zeros(20))
        preamble = np.append(preamble, np.ones(preamble_ones_length))
 
        return preamble


# +
##new preamble under test
#def preamble_generator(): #"a" is the bit array to be modulated
#
#
###        preamble = np.array([1,0,1,0])
###        preamble = np.append(preamble, np.ones(6))                       #np.ones(N) generates a numpy array containing N ones.
###        preamble = np.append(preamble, np.zeros(10))   #np.append(a,b) appends array b after array a
###        preamble = np.append(preamble, np.ones(30))
###        preamble = np.append(preamble, np.zeros(30))
###        preamble = np.append(preamble, np.ones(30))
###        preamble = np.append(preamble, np.zeros(30))
###        preamble = np.append(preamble, np.ones(30))
###        preamble = np.append(preamble, np.zeros(30))
#
#        np.random.seed(2022)
#        preamble = np.random.randint(0,2,200) 
#        
#        return preamble
