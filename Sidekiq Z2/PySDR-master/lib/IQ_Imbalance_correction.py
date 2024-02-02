#BPSK IQ Imbalance Correction
import numpy as np
import matplotlib.pyplot as plt
def Mean(array, period):
    means = []
    for index in range(0,len(array)):
        num = 0
        sum = array[index]
        num += 1
        for i in range(1, period+1):
            end = False
            if (index - i >= 0):
                num += 1
                sum += array[index - i]
            else:
                end = True
            if (index + i < len(array)):
                num += 1
                sum += array[index + i]
            elif end == True:
                break # speeds it up some
        means.append(sum / num)
    return means
def ShowConstellationPlot(packet, mag=999):
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
def IQ_Imbalance_Correct(packet, mean_period=10000):
    I = np.real(packet)
    Q = np.imag(packet)
    # I(t) = α cos (ωt) + βI
    # Q(t) = sin (ωt + ψ) + βQ
    # BI and BQ are simply the mean over X periods of I or Q, this value can be subtracted to remove
    BI = Mean(I, mean_period)
    BQ = Mean(Q, mean_period)
    I_second = I - BI
    Q_Second = Q - BQ
    a = np.sqrt(2 * np.mean(np.square(I_second)))
    sin_psi = (2 / 1) * np.mean(I_second * Q_Second)
    cos_psi = np.sqrt(np.subtract(1, np.square(sin_psi)))
    # Corrrection Matrix Parameters
    A = 1 / a
    C = - sin_psi / (a * cos_psi)
    D = 1 / cos_psi
    I_final = A * I + 0 * Q
    Q_final = C * I + D * Q
    crrected_packet = I_final + 1j * Q_final
    return crrected_packet

if __name__ == "__main__":
    a = np.array([1.4, 5.3, 6.3, 1.2, -3.0, 9.3, 5.4, 1.2, -9.1, 8.2, 1.0, 4.5, -2.3])
    print(Mean(a, 10))
    # testpacket = np.load('lib/IQImbalancetextpacket.npy')
    # ShowConstellationPlot(testpacket)
    # ShowConstellationPlot(testpacket,1)
    # corrected = IQ_Imbalance_Correct(testpacket)
    # ShowConstellationPlot(corrected, 1)
    # plt.scatter([-1,1], [0,0], c='y', label="Constellation Diagram (TestPacket)")
    # plt.scatter(np.real(testpacket), np.imag(testpacket), c='b', label="Constellation Diagram (TestPacket)")
    # plt.scatter(np.real(corrected), np.imag(corrected), c='r',label="Constellation Diagram (Corrected)")
    # plt.xlabel('Real (I)')
    # plt.ylabel('Imaginary (Q)')
    # plt.xlim(-1.5, 1.5)
    # plt.ylim(-1.5, 1.5)
    # plt.grid(True, which="both")
    # plt.show()



