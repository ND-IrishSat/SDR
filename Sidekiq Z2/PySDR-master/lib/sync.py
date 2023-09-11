"""
Synchronization convenience functions

"""

import numpy as np
from scipy import signal

#####################################
# Muller and Muller Clock recovery
# Only works for BPSK, not for other modulation schemes

def mmcr(testpacket):
    samples_interpolated = signal.resample_poly(testpacket, 16, 1)

    mu = 0 # initial estimate of phase of sample
    out = np.zeros(len(testpacket) + 10, dtype=complex)
    out_rail = np.zeros(len(testpacket) + 10, dtype=complex) # stores values, each iteration we need the previous 2 values plus current value
    i_in = 0 # input samples index
    i_out = 2 # output index (let first two outputs be 0)
    sps=8

    while i_out < len(testpacket) and i_in+16 < len(testpacket):
        out[i_out] = samples_interpolated[i_in*16 + int(mu*16)]
        out_rail[i_out] = int(np.real(out[i_out]) > 0) + 1j*int(np.imag(out[i_out]) > 0)
        x = (out_rail[i_out] - out_rail[i_out-2]) * np.conj(out[i_out-1])
        y = (out[i_out] - out[i_out-2]) * np.conj(out_rail[i_out-1])
        mm_val = np.real(y - x)
        mu += sps + 0.3*mm_val
        i_in += int(np.floor(mu)) # round down to nearest int since we are using it as an index
        mu = mu - np.floor(mu) # remove the integer part of mu
        i_out += 1 # increment output index
    out = out[2:i_out] # remove the first two, and anything after i_out (that was never filled out)
    testpacket = out # only include this line if you want to connect this code snippet with the Costas Loop later on

    return testpacket

########################################
# Coarse Frequency Detection and Correction

def coarse_freq_correction(testpacket, fs):
    fft_samples = testpacket**2

    psd = np.fft.fftshift(np.abs(np.fft.fft(fft_samples)))
    f = np.linspace(-fs/2.0, fs/2.0, len(psd))

    max_freq = f[np.argmax(psd)]
    Ts = 1/fs # calc sample period
    t = np.arange(0, Ts*len(testpacket), Ts) # create time vector
    testpacket = testpacket * np.exp(-1j*2*np.pi*max_freq*t/2.0) # some error here with length of time vector.
    
    return testpacket

######################################
# Costas Loop Fine Freq Correction

def costas(testpacket):
    N = len(testpacket)
    phase = 0
    freq = 0
    # These next two params is what to adjust, to make the feedback loop faster or slower (which impacts stability)
    alpha = 0.132
    beta = 0.00932
    out = np.zeros(N, dtype=complex)
    freq_log = []
    for i in range(N):
        out[i] = testpacket[i] * np.exp(-1j*phase) # adjust the input sample by the inverse of the estimated phase offset
        error = np.real(out[i]) * np.imag(out[i]) # This is the error formula for 2nd order Costas Loop (e.g. for BPSK)

        # Advance the loop (recalc phase and freq offset)
        freq += (beta * error)
        freq_log.append(freq * fs / (2*np.pi)) # convert from angular velocity to Hz for logging
        phase += freq + (alpha * error)

        # Optional: Adjust phase so its always between 0 and 2pi, recall that phase wraps around every 2pi
        while phase >= 2*np.pi:
            phase -= 2*np.pi
        while phase < 0:
            phase += 2*np.pi

    testpacket = out

    return testpacket

#####################################
# Frame Sync

def frame_sync(testpacket, preamble, data_len):
    scale = np.mean(np.abs(testpacket))
    out = np.array([])

    matched_filter_coef = np.flip(preamble)

    for symbol in testpacket:
        symbol = (symbol + scale)/2*scale
        out = np.append(out, symbol)

    crosscorr = signal.fftconvolve(out,matched_filter_coef)

    peak_indices, _ = signal.find_peaks(crosscorr)

    # peak_indices, _ = signal.find_peaks(crosscorr, height =  peak_threshold, distance = int(0.8*total_samples))
    # this finds ALL packets contained in a given buffer

    idx = np.array(crosscorr).argmax()

    recoveredPayload = testpacket[idx-len(preamble)+1:idx+data_len+1]
    recoveredData = recoveredPayload[len(preamble):]
    
    return recoveredPayload, recoveredData
