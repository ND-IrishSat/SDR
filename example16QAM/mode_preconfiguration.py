import numpy as np
import matplotlib.pyplot as plt
from scipy import signal
import socket
import os

# This module takes the operation mode and system parameters as input and perform preconfiguration task.

def tx_mode_preconfig(Mode, R, segment_size, N, b, test_packet_num):
        
        #Length of each UDP datagram sent during the test, this formula ensures there are always 10 datagrams inside each OOK packet.
        l = int(R*(segment_size/8)/10)

        if(Mode == 1 or Mode == 2):
                #IP and Port used to receive datagrams containing the transport stream or Iperf test data
                UDP_IP_ADDRESS = "127.0.0.1"
                UDP_PORT_NO = 5010

                serverSock = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)
                serverSock.bind((UDP_IP_ADDRESS, UDP_PORT_NO))
        else:
                serverSock = None

        if(Mode == 1):
                #Start recoding the stream and encoding the video
                #note: you may want to modify the camera name after video= based on your PC camera setup
                os.system('echo Webcam recording starts...')
                os.system('start /b ffmpeg -loglevel panic -f dshow -i video="Integrated Camera" -b:v 90k -minrate 80k -maxrate 100k -vf scale=80:45 -r 8 -bufsize 3000k -flush_packets 0 udp://127.0.0.1:5010?pkt_size='+str(l)+'/stream.ts')

        if(Mode == 2):
                #Start an Iperf test, in each test we transmit "test_packet_num" OOK packets.
                os.system('start cmd /k iperf -u -c 127.0.0.1 -p 5010 -n '+str(l*10*test_packet_num)+' -l '+ str(l) + ' -i 1 -b '+b)

        if(Mode == 4):
                #Generate the random sequence used for error test
                #Receiver has a copy of the same sequence seeded by the same number (2020) so we can compare two sequences for BER measurement.
                np.random.seed(2021)
                generated_sequence = np.random.randint(0,2,N-18)  # 1-bit sequence number to detect packet loss, and 17-bit timestamp to measure the service time.
                generated_sequence[0] = 1
                generated_sequence[int(N*1/4)] = 1
                generated_sequence[int(N/2)] = 1
                generated_sequence[int(N*3/4)] = 1
                sequence_counter = 0
        else:
                generated_sequence = None
                sequence_counter = None
                
        return serverSock, generated_sequence, sequence_counter, l

def rx_mode_preconfig(scheme,N,preamble_length,M, fs):
        
        if(scheme == "OOK"):
                total_samples = int((N+preamble_length)*M*1.0)    # total number of samples corresponding to one OOK packet

        if(scheme == "BPSK"):
                total_samples = int((N+preamble_length)*M*1.0)    # total number of samples 
        
        if(scheme == "QPSK"):
                total_samples = int((int(N/2) + preamble_length)*M*1.0)    # total number of samples 

        if(scheme == "QAM"):
                total_samples = int((int(N/4)+preamble_length)*M*1.0)    # total number of samples

        samples_perbit = int(M)     
        fs_in = fs*1.0                                    # ADC input sampling rate
        Ts_in = 1/fs_in                                   # ADC input sampling period

        return total_samples, samples_perbit, fs_in, Ts_in
