import numpy as np
import matplotlib.pyplot as plt
from scipy import signal

# Define the matched filtering module (Root-raised-cosine filters are used in the Tx and Rx side to form an overall
# raised-cosine response). (split the filters to TX & RX side so the receiver filter could be used to attenuate out-of-band channel noise)
# Source code for the commpy filters invoked here could be found in the following link:
# https://pydoc.net/scikit-commpy/0.3.0/commpy.filters/


# The pulse shaping module takes the following arguments as inputs:

# baseband:                 This is the baseband signal

# M:                        This is the receiver-side oversampling factor, each symbol (bit) being transmitted is represented
#                           by M samples

# fs:                       This is the input sampling rate of the ADC

# alpha:                    This is the roll-off factor (valid value [0,1])

# L:                        This is the length (span) of the pulse-shaping filter (in the unit of symbols)


# The pulse shaping module returns the following argument as output:

# symbols                    The symbols are the results of the matched filtering, decimated.

def matched_filtering(baseband, M, fs, alpha, L, pulse_shape):

        if(pulse_shape == 'rrc'):
            #Root Raised-Cosine span: +/- 4 symbols
            N = L*M

            T_symbol = 1/(fs/M)

            ##Square root raised-cosine (SRRC)

            time, h = rrcosfilter(N,alpha, T_symbol, fs)

    # #        baseband_filtered = np.convolve(baseband,h)
    # #        plt.plot(baseband)
    # #        plt.title('Baseband Unfiltered')
    # #        plt.show()

            baseband_filtered = np.convolve(baseband,h)

    # #        plt.plot(baseband_filtered)
    # #        plt.title('Baseband Filtered')
    # #        plt.show()


            #Downsample by a factor of M
            baseband_filtered = baseband_filtered[(L*M):]
            symbols = signal.upfirdn([1],baseband_filtered,1, M)
        
        if(pulse_shape == 'rect'):

                Ts = 1/fs

                #rectangular pulse
                h = np.ones(M)
                
                baseband_filtered = np.convolve(baseband,h)
                baseband_filtered[(L*M):]
                symbols = signal.upfirdn([1],baseband_filtered,1, M)

        return symbols


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

