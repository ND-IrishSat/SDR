import numpy as np
import matplotlib.pyplot as plt
from scipy import signal

# Define the pulse shaping module (Rectangular or Root-raised-cosine filters are used in the Tx and Rx side).
# (split the filters to TX & RX side so the receiver filter could be used to attenuate out-of-band channel noise)
# Source code for the commpy filters invoked here could be found in the following link:
# https://pydoc.net/scikit-commpy/0.3.0/commpy.filters/


# The pulse shaping module takes the following arguments as inputs:

# a:                        Base band symbols (containing both the preamble and payload). 

# M:                        This is the transmitter-side oversampling factor, each symbol (bit) being transmitted is represented
#                           by M samples

# fs:                       This is the sampling rate of the DAC

# pulse_shape:              "rect" (rectangular) or "rrc" (root-raised cosine)

# alpha:                    This is the roll-off factor for RRC (valid value [0,1])

# L:                        This is the length (span) of the pulse-shaping filter (in the unit of symbols)

# The pulse shaping module returns the following argument as output:

# baseband                   This is the baseband signal

#filter length: make it a parameter w/ default value
#TODO: introduce memory and make pulse shaping streaming
def pulse_shaping(a, M, fs, pulse_shape, alpha, L):
        
        #Upsample by a factor of M
        y = signal.upfirdn([1],a,M)

        if(pulse_shape == 'rrc'):

                #Root Raised-Cosine span
                N = L*M

                T_symbol = 1/(fs/M)

                ##Square root raised-cosine (SRRC)

                time, h = rrcosfilter(N,alpha, T_symbol, fs)


        if(pulse_shape == 'rect'):

                Ts = 1/fs

                #rectangular pulse
                h = np.ones(M)

        baseband = np.convolve(y,h)


        return baseband


def rrcosfilter(N, alpha, Ts, Fs):
    """
    Generates a root raised cosine (RRC) filter (FIR) impulse response.
 
    Parameters
    ----------
    N : int
        Length of the filter in samples.
 
    alpha : float
        Roll off factor (Valid values are [0, 1]).
 
    Ts : float
        Symbol period in seconds.
 
    Fs : float
        Sampling Rate in Hz.
 
    Returns
    ---------
 
    h_rrc : 1-D ndarray of floats
        Impulse response of the root raised cosine filter.
 
    time_idx : 1-D ndarray of floats
        Array containing the time indices, in seconds, for
        the impulse response.
    """
 
    T_delta = 1/float(Fs)
    time_idx = ((np.arange(N)-N/2))*T_delta
    sample_num = np.arange(N)
    h_rrc = np.zeros(N, dtype=float)
 
    for x in sample_num:
        t = (x-N/2)*T_delta
        if t == 0.0:
            h_rrc[x] = 1.0 - alpha + (4*alpha/np.pi)
        elif alpha != 0 and t == Ts/(4*alpha):
            h_rrc[x] = (alpha/np.sqrt(2))*(((1+2/np.pi)* \
                    (np.sin(np.pi/(4*alpha)))) + ((1-2/np.pi)*(np.cos(np.pi/(4*alpha)))))
        elif alpha != 0 and t == -Ts/(4*alpha):
            h_rrc[x] = (alpha/np.sqrt(2))*(((1+2/np.pi)* \
                    (np.sin(np.pi/(4*alpha)))) + ((1-2/np.pi)*(np.cos(np.pi/(4*alpha)))))
        else:
            h_rrc[x] = (np.sin(np.pi*t*(1-alpha)/Ts) +  \
                    4*alpha*(t/Ts)*np.cos(np.pi*t*(1+alpha)/Ts))/ \
                    (np.pi*t*(1-(4*alpha*t/Ts)*(4*alpha*t/Ts))/Ts)
 
    return time_idx, h_rrc





