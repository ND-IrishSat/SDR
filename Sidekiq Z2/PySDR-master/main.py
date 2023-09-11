"""

This is an end-to-end example of data generation, modulation, pulse shaping, 
channel simulation, synchronization, and demodulation. Much of this code was
informed by or taken from two sources: 

1. Source code from Dr. Laneman's RF Lab
2. Example code from PySDR (https://pysdr.org/content/rds.html)

An attempt has been made to provide explanation and description of each
step in the signal processing chain, along with plots to graphically illustrate
what is going on. Comments have also been inserted to indicate what future work 
needs to be done.

This framework can hopefully be used to build a functional SDR link. Work has
been started, but not finished, on getting this DSP framework to interface
with USRP SDR and Sidekiq Z2 SDR through their respective software APIs.

GENERAL OUTLINE OF CODE:
1. Declaration of Variables
2. Preamble, Data, and Matched Filter Coefficient Generation
3. Generation of Pulse Train // BPSK Modulation
4. Pulse Shaping
5. Simulation of Channel
6. Clock Recovery
7. Coarse Frequency Correction
8. Fine Frequency Correction
9. Frame Sync
10. Demodulation

"""

from lib import pulse_shaping
from lib import symbol_demod_QAM as demod
import numpy as np
from scipy import signal
import matplotlib.pyplot as plt
from lib import CRC

# Declare Variables
fs = 2.4e9             # carrier frequency 
Ts = 1.0 / fs          # sampling period in seconds
f0 = 0.0               # homodyne (0 HZ IF)
M = 8                  # oversampling factor
L = 8                  # pulse shaping filter length (Is this the ideal length?)
pulse_shape = 'rrc'    # type of pulse shaping, also "rect"
scheme = 'BPSK'        # Modulation scheme 'OOK', 'QPSK', or 'QAM'
alpha = 0.5            # roll-off factor of the RRC pulse-shaping filter

"""
Preamble, Data, and Matched Filter Coefficient Generation

The preamble is a special binary sequence used for its advantageous
autocorrelative properties. This means that when it is correlated with
itself, there is exactly one maximum peak, with the surrounding "peaks"
not exceeding a magnitude of 1

The data generated is a random binary sequence of length 256. This number was chosen
arbitrarily.

This data is then assigned a CRC polynomial based on the key below. This key is the "best"
12 bit key for a Hamming Distance of 2. Basically, the key is used as the divisor in a
modulo-2 division operation on the data. The remainder is another polynomial that is then
appended to the data. More explanation to come at the error check in the end.

The matched filter coefficients are used in frame detection. This is a flipped version
of the preamble sequence that is convolved with the received data to generate a crosscorrelation.

Resources:
https://ntrs.nasa.gov/citations/19800017860
https://www.geeksforgeeks.org/cyclic-redundancy-check-python/#

"""
preamble = np.array([0,1,0,0,0,0,1,1,0,0,0,1,0,1,0,0,1,1,1,1,0,1,0,0,0,1,1,1,0,0,1,0,0,1,0,1,1,0,1,1,1,0,1,1,0,0,1,1,0,1,0,1,0,1,1,1,1,1,1,0]).astype(float) # optimal periodic binary code for N = 63 https://ntrs.nasa.gov/citations/19800017860

data = np.random.randint(2, size=256).astype(float)

CRC_key = np.array([1,0,0,1,1,0,0,0,0,1,1,1]) # Best CRC polynomials: https://users.ece.cmu.edu/~koopman/crc/
data_encoded = CRC.encodeData(data, CRC_key)

bits = np.append(preamble,data_encoded)

sps = M # samples per symbol, equal to oversampling factor

matched_filter_coef = np.flip(preamble)

"""
Generation of Pulse Train, BPSK modulation

Bits are turned into a pulse train. Not exactly sure why, copied from PySDR.
Best guess is that it makes the following signal processing easier. It seems
that clock recovery doesn't work if I don't input a pulse train.

BPSK modulation is here done at the same time as the conversion of bits to
pulses. 0 is mapped to -1, 1 is mapped to 1 --> BPSK

Resources:
https://pysdr.org/content/digital_modulation.html

"""

pulse_train = np.array([])

for bit in bits:
    pulse = np.zeros(sps)
    pulse[0] = bit*2-1 # set the first value to either a 1 or -1
    pulse_train = np.concatenate((pulse_train, pulse)) # add the 8 samples to the signal

plt.stem(pulse_train, label="Pulse Train")
plt.stem(bits, 'ro', label="Original Bits")
plt.legend(loc="upper right")
plt.title("Pulse Train")
plt.show()

"""
Pulse Shaping (Raised Root Cosine)

Rasied Root Cosine filter applied. Basic idea is that it stretches out the 
symbols in time, thereby condensing them in frequency. Used to limit bandwidth
of symbol in frequency and also to avoid inter-symbol interference.

Resources:
https://pysdr.org/content/pulse_shaping.html
https://wirelesspi.com/pulse-shaping-filter/

"""

symbols_I = pulse_shaping.pulse_shaping(pulse_train, sps, fs, pulse_shape, alpha, L)

testpacket = symbols_I

plt.stem(symbols_I, 'ko')
plt.title("After pulse shaping")
plt.show()

#######################
# Transmission: Where transmission would occur in a real system

"""
Simulation of Channel

This bit of code simulates the non-idealities that would be present if this
data were to be transmitted over a real wireless channel.

1. Noise:
    Thermal noise in the channel would mean that the energy of each of your
    symbols would be slightly altered, resulting in slightly different sampled values
2. Fractional Delay:
    Since the local oscillators controlling the up/down conversion of your signal 
    (at transmitter and receiver, respectively) will have slightly different phase,
    this will be manifested as some fractional sample delay in your sampled bitstream
3. Frequency Offset:
    Since local oscillators are not perfect, even if the tx and rx LOs are tuned to the exact
    same frequency, there will be some small amount of error between them. In addition,
    since our receiver will be on a satellite moving at some high velocity relative to our
    ground station, there will be some doppler shift that must be accounted for.

Reference: 
https://pysdr.org/content/sync.html
https://wirelesspi.com/what-is-a-symbol-timing-offset-and-how-it-distorts-the-rx-signal/
https://wirelesspi.com/channel-estimation-in-wireless-communication/
https://wirelesspi.com/what-is-carrier-frequency-offset-cfo-and-how-it-distorts-the-rx-symbols/
https://wirelesspi.com/what-is-carrier-phase-offset-and-how-it-affects-the-symbol-detection/

"""

###########################
# Generate Noise

# Parameters for AWGN
mean = 0      # Mean of the Gaussian distribution (usually 0 for AWGN)
std_dev = 0.1   # Standard deviation of the Gaussian distribution
num_samples = len(testpacket)
awgn_samples = np.random.normal(mean, std_dev, num_samples)
testpacket += awgn_samples

#################################
# Add fractional delay

# Create and apply fractional delay filter
delay = 0.4 # fractional delay, in samples
N = 21 # number of taps
n = np.arange(-N//2, N//2) # ...-3,-2,-1,0,1,2,3...
h = np.sinc(n - delay) # calc filter taps
h *= np.hamming(N) # window the filter to make sure it decays to 0 on both sides
h /= np.sum(h) # normalize to get unity gain, we don't want to change the amplitude/power
testpacket = np.convolve(testpacket, h) # apply filter

###################################
# Add frequency offset - NOTE: Frequency offset of > 1% will result in inability of crosscorrelation operation to detect frame start

# apply a freq offset
fs = 2.45e9 # arbitrary UHF frequency
fo = 61250 # Simulated frequency offset
Ts = 1/fs # calc sample period
t = np.arange(0, Ts*(len(testpacket)), Ts) # create time vector
testpacket = testpacket * np.exp(1j*2*np.pi*fo*t) # perform freq shift

plt.stem(symbols_I, label="Original Pulse Shaped Waveform")
plt.stem(testpacket, 'ro', label="Non-ideal Waveform")
plt.stem(np.imag(testpacket), 'mo', label="Non-ideal Waveform")
plt.title("After fractional delay and frequency offset")
plt.legend(loc="upper left")
plt.show()

"""
Muller and Mueller Clock Recovery

Note: Only works for BPSK, not for other modulation schemes

This is the beginning of the receiver synchronization section
of this code. At this point, the signal would have been received and down-
converted by the receiver. The analog signal will have been converted to
digital IQ data. Present in this IQ data will be the non-idealities 
discussed previously.

This Muller and Mueller Clock Recovery algorithm acts to correct whatever fractional sampling
offset is present in the data. This is technically a digital implementation of a phase locked loop,
so it would be helpful to have some background information on these.

First, the data is interpolated to ensure that the algorithm is able to make a fractional 
delay determination. Next, the algorithm takes in data, continuously tuning a 'mu' parameter 
(representing the fractional delay) until it reaches some lower error threshold. As it operates, 
the algorithm pulls what it believes are the best samples out and into an output array, downsampling
and ensuring that the output contains one sample per symbol.

This algorithm acts in real time, which makes it useful for continuous operation in real systems.

References:
https://wirelesspi.com/mueller-and-muller-timing-synchronization-algorithm/
https://wirelesspi.com/phase-locked-loop-pll-in-a-software-defined-radio-sdr/
https://wirelesspi.com/phase-locked-loop-pll-for-symbol-timing-recovery/

"""

samples_interpolated = signal.resample_poly(testpacket, 16, 1)

plt.stem(np.real(samples_interpolated), 'bo', label="Non-ideal Interpolated Waveform")
plt.title("After interpolation")
plt.legend(loc="upper left")
plt.show()

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


plt.stem(signal.upfirdn([1],pulse_train,1, sps), 'ko', label="Original Pulse Train")
plt.stem(np.real(testpacket), 'ro', label="Real Part of Clock Recovered Waveform")
plt.stem(np.imag(testpacket), 'mo', label="Imaginary Part of Clock Recovered Waveform")
plt.title("After clock recovery")
plt.legend(loc="upper left")
plt.show()

"""
Coarse Frequency Detection and Correction

This coarse frequency detection and correction is for correcting large frequency offsets,
and results in a kind of "ballpark" corrected signal. First, the samples are squared in 
order to "remove the modulation" (covered in https://pysdr.org/content/sync.html Coarse 
Frequency Correction section). An FFT is then run on this data to find the frequency offset 
(which is manifested as a small sinusoid riding on top of the data signal: the magenta in the
plots you have seen). When this frequency offset is found (by taking the maximum argument of the FFT)
it is corrected by multiplying the sample data by the negative complex sinusoid at that frequency.

"""
fft_samples = testpacket**2

psd = np.fft.fftshift(np.abs(np.fft.fft(fft_samples)))
f = np.linspace(-fs/2.0, fs/2.0, len(psd))

max_freq = f[np.argmax(psd)]
Ts = 1/fs # calc sample period
t = np.arange(0, Ts*len(testpacket), Ts) # create time vector
testpacket = testpacket * np.exp(-1j*2*np.pi*max_freq*t/2.0) # multiply by negative complex sinusoid at offset frequency

plt.plot(f, psd)
plt.title("Frequency offset before correction")
plt.show()

# Plot
plt.stem(signal.upfirdn([1],pulse_train,1, sps), 'ko', label="Original Pulse Train")
plt.stem(np.real(testpacket), label="Real Part of Recovered Waveform")
plt.stem(np.imag(testpacket), 'ro', label="Imaginary Part of Recovered Waveform")
plt.legend(loc="upper left")
plt.title("After coarse frequency correction")
plt.show()

"""
Costas Loop Fine Freq Correction

The costas loop is a digital implementation of a Phase Locked Loop (like the
Muller and Mueller Clock Recovery Algorithm). The basic principle of a Costas Loop
is that it has two parameters, phase and frequency, that are continuously updated
through an error function. Like the M&M algorithm, it takes some samples to "lock on".
The Costas Loop will continuously track phase and frequency offset, making it suitable
for tracking and correcting for things like doppler shift. Since it can operate with 
more precision than the coarse frequency corrector, it is used for "fine frequency
correction"

References:
https://wirelesspi.com/costas-loop-for-carrier-phase-synchronization/
https://pysdr.org/content/sync.html

"""

N = len(testpacket)
phase = 0
freq = 0
# These next two params are what to adjust, to make the feedback loop faster or slower (which impacts stability)
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

# Plot freq over time to see how long it takes to hit the right offset
plt.plot(freq_log,'.-')
plt.title("Estimated frequency offset of Costas Loop vs Sample index")
plt.show()

testpacket = out

plt.stem(signal.upfirdn([1],pulse_train,1, sps), 'ko', label="Original Pulse Train")
plt.stem(np.real(testpacket), label="Real Part of Recovered Waveform")
plt.stem(np.imag(testpacket), 'ro', label="Imaginary Part of Recovered Waveform")
plt.legend(loc="upper left")
plt.title("After fine frequency correction")
plt.show()

# Testing frequency correction
fft_samples = testpacket**2
psd = np.fft.fftshift(np.abs(np.fft.fft(fft_samples)))
f = np.linspace(-fs/2.0, fs/2.0, len(psd))

plt.plot(f, psd)
plt.title("Frequency offset after correction")
plt.show()

"""
Frame Sync

Frame Sync relies on a known preamble sequence that is appended to 
the data at the transmit side and is used to correlate on the receive side.
Since the sequence that we chose has great autocorrelative properties, the
maximum argument of the resultant correlation of this known sequence with our
received data will almost always correspond to the starting index of our
payload. If you find that due to channel noise or other non-idealities, 
you are having trouble with frame sync, try increasing the proportion 
of preamble size to data size.

Reference: 
https://pysdr.org/content/sync.html
https://wirelesspi.com/correlation/
https://wirelesspi.com/convoluted-correlation-between-matched-filter-and-correlator/

"""

scale = np.mean(np.abs(testpacket))
out = np.array([])

for symbol in testpacket:
    symbol = (symbol + scale)/2*scale
    out = np.append(out, symbol)

crosscorr = signal.fftconvolve(out,matched_filter_coef)

plt.stem(np.real(testpacket), label="Recovered Waveform")
plt.stem(preamble, 'r', label="Preamble Sequence")
plt.stem(crosscorr, 'g', label="Crosscorrelation")
plt.title("Frame Synchronization: Crosscorrelation")
plt.legend(loc="upper right")
plt.show()

# peak_indices, _ = signal.find_peaks(crosscorr, height =  peak_threshold, distance = int(0.8*total_samples))
# this finds ALL packets contained in a given buffer
# You will need to do something like this for asynchronous processing

idx = np.array(crosscorr).argmax()

recoveredPayload = testpacket[idx-len(preamble)+1:idx+len(data_encoded)+1] # Reconstruct original packet minus preamble
recoveredData = recoveredPayload[len(preamble):]

# Plot
plt.stem(np.real(recoveredPayload), label="Recovered Payload")
plt.stem(np.imag(recoveredPayload), 'ro')
plt.stem(signal.upfirdn([1],pulse_train,1, sps), 'ko', label="Original Pulse Train")
plt.title("Data recovery")
plt.legend(loc="upper right")
plt.show()


"""
Demodulation // Error Checking

Demodulation is the process of mapping data symbols back to data bits. A
maximum likelihood detection scheme is implemented by making determinations
about decision boundaries and then classifying received symbols based
on these decision boundaries. For BPSK, the decision boundary is I=0.
Anything to the left of I=0 is interpreted as a -1, anything to the 
right is interpreted as a 1. 1 maps to 1, while -1 maps to 0. This
is how the original bits are reconstructed.

Error checking is done using Cyclic Redundancy Check. This basically does
modulo-2 division using a known polynomial key. When the data was initialized
and encoded, the remainder of the mod-2 division operation was appended to the
data. Now, when mod-2 division is done using the same key, the remainder should
be zero.

This is simply an error detecting operation. Forward error correction, which
involves sending redundant data, should be looked into in the future to 
minimize BER.

Reference: 
https://pysdr.org/content/digital_modulation.html
https://wirelesspi.com/packing-more-bits-in-one-symbol/
https://www.geeksforgeeks.org/cyclic-redundancy-check-python/#
http://www.sunshine2k.de/articles/coding/crc/understanding_crc.html

"""

demod_bits = demod.symbol_demod(recoveredData, scheme, 1, len(preamble)) # gain has to be set to 1
error = CRC.CRCcheck(demod_bits, CRC_key)
print("CRC error: " + str(error))

# Plot
plt.stem(np.real(demod_bits), label="Recovered Data")
plt.stem(0.5*np.real(data), 'ro', label="Original Data")
plt.title("Demodulated bit comparison")
plt.legend(loc="upper right")
plt.show()

####################################
# To do 

# 1. For packet detection, add threshold slightly less than number of ones in preamble
# 2. Look into equalizing channel gain
# 3. Look into IQ imbalance correction
# 4. Add error checking

# Look at the docs for further work and suggestions
