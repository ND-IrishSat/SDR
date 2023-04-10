import numpy as np
import matplotlib.pyplot as plt
from scipy import signal
import time

# This module takes the operation mode and system parameters as input and generate the
# source data (bits) to feed the symbol modulator.

def data_source(Mode, serverSock, generated_sequence, sequence_counter, l):

        if(Mode == 1 or Mode == 2):
                
                #Hoarding 10 UDP datagrams to send in one OOK packet
                data, addr = serverSock.recvfrom(l)
                final_bytes = np.frombuffer(data,dtype=np.uint8)


                for i in range (9):
                        data, addr = serverSock.recvfrom(l)
                        temp_bytes = np.frombuffer(data,dtype=np.uint8)
                        final_bytes = np.append(final_bytes,temp_bytes)
##                        print("i: ",i)

                        
                Bits = np.unpackbits(final_bytes) #unpacking bytes into bits for modulation

# #                if(Bits.size != N):
# #                        continue

                print ("10 datagrams received and ready to be transimitted!")
                print ("Size of each datagram: ", l)

        if(Mode == 4):
                Bits = generated_sequence

                #Add the one-bit field length sequence number

##                sequence_number = sequence_counter%2
                sequence_number = sequence_counter%4
                if sequence_number == 0:
                        sequence_code = np.array([0,0])
                elif sequence_number == 1:
                        sequence_code = np.array([0,1])
                elif sequence_number == 2:
                        sequence_code = np.array([1,0])
                else:
                        sequence_code = np.array([1,1])
                        
                Bits = np.append(Bits, sequence_code)
##                Bits = np.append(Bits, sequence_number)

# #                #Add the seventeen-bit timestamp
# #
# #                time_stamp_int = int(time.time()*100)%100000
# #                time_stamp_bin = '{0:017b}'.format(time_stamp_int)
# #
# #                for i in range (17):
# #
# #                        Bits = np.append(Bits, int(time_stamp_bin[i]))

                #Add the timestamp

                time_stamp_int = int(time.time()*100)%100000
                time_stamp_bin = '{0:016b}'.format(time_stamp_int)

                for i in range (16):

                        Bits = np.append(Bits, int(time_stamp_bin[i]))
                        
        #Bits = np.append(np.append(Bits, Bits),np.append(Bits,Bits))
        #Bits = np.append(Bits, Bits) 
        return Bits



