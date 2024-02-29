import numpy as np
from scipy import signal
import matplotlib.pyplot as plt

path = "Sidekiq Z2/CTranslation/lib/graphs/"

def ShowConstellationPlot(packet, mag=999, title="After IQ Imbalance Correction"):
    x = []
    y = []
    phase = []
    magnitudes = []
    for data in packet:
        phase.append(np.angle(data))
        if (mag == 999):
            magnitudes.append(np.abs(data))
        else:
            magnitudes.append(1)
    for index in range(0,len(phase)):
        magnitude = magnitudes[index]
        angle = phase[index]
        x_coord = magnitude * np.cos(angle)
        y_coord = magnitude * np.sin(angle)
        x.append(x_coord)
        y.append(y_coord)
    fig = plt.figure()
    ax = fig.add_subplot()
    ax.set_aspect('equal', adjustable='box')
    plt.title(title)
    plt.scatter(x, y, label="Constellation Diagram (Packet In)")
    plt.xlabel('Real (I)')
    plt.ylabel('Imaginary (Q)')
    plt.xlim(-2.5, 2.5)
    plt.ylim(-2.5, 2.5)
    plt.grid(True, which="both")
    plt.show()

def PlotWave(packet):
    I = np.real(packet)
    Q = np.imag(packet)
    plt.plot(np.arange(len(I)),I, '.-')
    plt.plot(np.arange(len(Q)),Q, '.-')
    plt.title("I & Q Packets")
    plt.grid()
    plt.show()
def documentToArray(filename):
    f = open(path + filename, "r")
    complex_num = False
    for line in f.readlines():
        if ("j" in line):
            complex_num = True
            break
    f.close()
    f = open(path + filename, "r")
    if (complex_num):
        arr = np.fromfile(f, dtype=complex, sep="\n")
        f.close()
        return arr
    arr = np.fromfile(f, dtype=float, sep="\n")
    f.close()
    return arr

if __name__ == "__main__":
    print("Graphing")
    preamble = documentToArray("preamble.txt")
    bits = documentToArray("bits.txt")
    pulse_train = documentToArray("pulsetrain.txt")
    print(f"Preamble: {len(preamble)}")
    print(f"Bits: {len(bits)}")
    # Pulse Train
    print(f"Pulse Train: {len(pulse_train)}")
    plt.stem(pulse_train, label="Pulse Train")
    plt.stem(bits, 'ro', label="Original Bits")
    plt.legend(loc="upper right")
    plt.title("Pulse Train")
    plt.show()
    #Pulse Shaping
    symbols_I = documentToArray("pulseshaping.txt")
    print(f"Pulse Shape: {len(symbols_I)}")
    plt.stem(np.real(symbols_I), 'ko')
    plt.title("After pulse shaping")
    plt.show()
    #Noise
    testpacket = documentToArray("noise.txt")
    print(f"Noise: {len(testpacket)}")
    #After Freq Offset
    testpacket = documentToArray("testpacketfreqshift.txt")
    print(f"After Freq Offset: {len(testpacket)}")
    plt.stem(np.real(symbols_I), label="Original Pulse Shaped Waveform")
    plt.stem(np.real(testpacket), 'ro', label="Non-ideal Waveform")
    plt.stem(np.imag(testpacket), 'mo', label="Non-ideal Waveform")
    plt.title("After fractional delay and frequency offset")
    plt.legend(loc="upper left")
    plt.show()
    #Clock Recovery — Interpolated
    samples_interpolated = documentToArray("samplesinterpolated.txt")
    print(f"samples_interpolated: {len(samples_interpolated)}")
    plt.stem(np.real(samples_interpolated), 'bo', label="Non-ideal Interpolated Waveform")
    plt.title("After interpolation")
    plt.legend(loc="upper left")
    plt.show()

    #After Clock Recovery
    sps = 8
    testpacket = documentToArray("clockRecovery.txt")
    print(f"After Clock Recovery: {len(testpacket)}")
    plt.stem(signal.upfirdn([1],pulse_train,1, sps), 'ko', label="Original Pulse Train")
    plt.stem(np.real(testpacket), 'ro', label="Real Part of Clock Recovered Waveform")
    plt.stem(np.imag(testpacket), 'mo', label="Imaginary Part of Clock Recovered Waveform")
    plt.title("After clock recovery")
    plt.legend(loc="upper left")
    plt.show()

    #Coarse Frequency Detection and Correction
    psd = documentToArray("psd_coarseFrequencyCorrection.txt")
    f = documentToArray("f_coarseFrequencyCorrection.txt")
    plt.plot(f, psd)
    plt.title("Frequency offset before correction")
    plt.show()
    testpacket = documentToArray("coarseFrequencyCorrection.txt")
    print(f"Coarse Frequency: {len(testpacket)}")
    # Plot
    plt.stem(signal.upfirdn([1],pulse_train,1, sps), 'ko', label="Original Pulse Train")
    plt.stem(np.real(testpacket), label="Real Part of Recovered Waveform")
    plt.stem(np.imag(testpacket), 'ro', label="Imaginary Part of Recovered Waveform")
    plt.legend(loc="upper left")
    plt.title("After coarse frequency correction")
    plt.show()

    #Costas Loop Fine Freq Correction
    freq_log = documentToArray("costasFreqLog.txt")
    print(f"freq_log: {len(freq_log)}")
    plt.plot(freq_log,'.-')
    plt.title("Estimated frequency offset of Costas Loop vs Sample index")
    plt.show()
    testpacket = documentToArray("costasout.txt")
    print(f"Costas: {len(testpacket)}")
    plt.stem(signal.upfirdn([1],pulse_train,1, sps), 'ko', label="Original Pulse Train")
    plt.stem(np.real(testpacket), label="Real Part of Recovered Waveform")
    plt.stem(np.imag(testpacket), 'ro', label="Imaginary Part of Recovered Waveform")
    plt.legend(loc="upper left")
    plt.title("After fine frequency correction")
    plt.show()
    #plot
    psd = documentToArray("psd_fine_freq_correct.txt")
    f = documentToArray("f_fine_freq_correct.txt")
    plt.plot(f, psd)
    plt.title("Frequency offset after correction")
    plt.show()

    # IQ Imbalance correction
    testpacket = documentToArray("iqimbalanceout.txt")
    print(f"IQ Imbalance: {len(testpacket)}")
    PlotWave(testpacket)
    ShowConstellationPlot(testpacket)

    #Frame Sync
    crosscorr = documentToArray("crosscorr.txt")
    print(f"Crosscorr: {len(crosscorr)}")
    plt.stem(np.real(testpacket), label="Recovered Waveform")
    plt.stem(preamble, 'r', label="Preamble Sequence")
    plt.stem(np.real(crosscorr), 'g', label="Crosscorrelation")
    plt.title("Frame Synchronization: Crosscorrelation")
    plt.legend(loc="upper right")
    plt.show()
    # Plot
    recoveredPayload = documentToArray("recoveredPayload.txt")
    print(f"Recovered Payload: {len(recoveredPayload)}")
    plt.stem(np.real(recoveredPayload), label="Recovered Payload")
    plt.stem(np.imag(recoveredPayload), 'ro')
    plt.stem(signal.upfirdn([1],pulse_train,1, sps), 'ko', label="Original Pulse Train")
    plt.title("Data recovery")
    plt.legend(loc="upper right")
    plt.show()
    
    demod_bits = documentToArray("demodbits.txt")
    data = documentToArray("data.txt")
    print(f"Demod Bits: {len(demod_bits)}")
    print(f"Data: {len(data)}")

    # This just calculates how many bits were correctly received
    num_correct = 0
    num = 0
    for index in range(0,len(data)):
        num += 1
        value = "NaN"
        if (index <= len(demod_bits)-1):
            value = demod_bits[index]
            if (value == data[index]):
                num_correct += 1
    print(f"Received: {num_correct} / {num} bits   |   {round(num_correct / num * 100)}%")

    #Check Receiving
    plt.stem(np.real(demod_bits), label="Recovered Data")
    plt.stem(0.5*np.real(data), 'ro', label="Original Data")
    plt.title("Demodulated bit comparison")
    plt.legend(loc="upper right")
    plt.show()