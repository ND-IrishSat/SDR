import numpy as np
import matplotlib.pyplot as plt
from scipy import signal
import time
import os
import socket

                
                        
#Define the data_sink / postprocessing function
def data_sink(Mode,q4, R, N, debugging, q3, visualize, q8, segment_size):

        if(Mode == 1 or Mode == 2):
                #Length of each UDP datagram sent during the test or video streaming
                l = int(R*(segment_size/8)/10)
                
                #IP and Port where the datagrams containing iperf test data or video stream are sent
                UDP_IP2 = "127.0.0.1"
                UDP_PORT2 = 5005

                clientSock = socket.socket(socket.AF_INET,socket.SOCK_DGRAM)

                streaming = 1 #streaming flag: 1->streaming ready to begin, 0->streaming ongoing

        if(Mode==2):
                #Start an Iperf test
                os.system('start cmd /k iperf  -u -s -p 5005 -l '+str(l)+' -i 1')
        
        if(Mode==4):
                #Generate the random sequence used for error test
                np.random.seed(2021)
                generated_sequence = np.random.randint(0,2,N-1)
                generated_sequence[0] = 1
                generated_sequence[int(N*1/4)] = 1
                generated_sequence[int(N/2)] = 1
                generated_sequence[int(N*3/4)] = 1
##                print("generated_sequence: ",generated_sequence[0:1000])

                
                expected_sequence_number = 0
                #Counting the packets received and packets with errors for the frame error rate calculation
                packet_errors = 0.0
                packet_loss = 0.0
                short_term_packet_loss = 0.0


        #Parameters for throughput calculation
        packet_kbits = N/1024
        total_kbits = 0
        window_counter = 0
        long_term_BER = 0
        long_term_BER_counter = 0
        long_term_BER_sum = 0

        #temp for debugging
        counter = 0

        if(visualize == 1):
                frame_rate_factor = 30
                visualize_counter = 0


        while(Mode == 1 or Mode == 2 or Mode == 4):

                if(window_counter == 0):
                        start = time.time()
                        total_kbits = 0

                a_demodulated, baseband_symbols, data = q4.get()

                if(Mode == 1 and streaming == 1):
                        os.system('start /b ffplay -loglevel panic -x 320 -y 180 -left 800 -top 60 udp://127.0.0.1:5005') #follow the stream
                        streaming = 1   #not resetting the streaming flag in wireless case for quick recovery/restart from burst interference

                if(Mode == 1 or Mode == 2):
          
                        #received packet containing 10 UDP datagrams
                        received_packets = np.packbits(a_demodulated)

                        received_datagrams = np.array_split(received_packets,10)

                        for i in range(10):
                                udp_packet = received_datagrams[i].tobytes()
                                clientSock.sendto(udp_packet, (UDP_IP2, UDP_PORT2))

                        print ("10 Datagrams received and sent!")
                        print ("Size of one datagram: ", len(udp_packet))

                if(Mode==4):

                        #Compare a_demodulated with the generated_sequence
                        packet_flag = 0
                        packet_error_bits = 0 #report number of error bits in each packet

                        #Record the histogram (distribution) of error bits in a packet
                        error_bin = [0]*5

##                        print("first 20 bits (demod): ",  a_demodulated[0:20])
##                        print("first 20 bits (actual): ",  generated_sequence[0:20])
                        
                        for i in range(len(a_demodulated)-18):
                                if(a_demodulated[i] != generated_sequence[i]):
                                        packet_error_bits = packet_error_bits+1
                                        packet_flag = 1
                                        if(i < int(0.2*len(a_demodulated))):
                                                error_bin[0] = error_bin[0]+1
                                        elif(i < int(0.4*len(a_demodulated))):
                                                error_bin[1] = error_bin[1]+1
                                        elif(i < int(0.6*len(a_demodulated))):
                                                error_bin[2] = error_bin[2]+1
                                        elif(i < int(0.8*len(a_demodulated))):
                                                error_bin[3] = error_bin[3]+1
                                        else:
                                                error_bin[4] = error_bin[4]+1
                        if(debugging == 1):                       
                                print()
                                print("bin 1:",error_bin[0])
                                print("bin 2:",error_bin[1])
                                print("bin 3:",error_bin[2])
                                print("bin 4:",error_bin[3])
                                print("bin 5:",error_bin[4])
                                print()
                       
                        if (packet_flag):
                                packet_errors = packet_errors + 1.0


                        #Check the 1-bit sequence number

##                        detected_sequence_number = a_demodulated[len(a_demodulated)-18]
##                        print("detected sequence number: ", detected_sequence_number)

                        #Check the 2-bit sequence number

                        detected_sequence_number = a_demodulated[(len(a_demodulated)-18):(len(a_demodulated)-16)]
                        if detected_sequence_number[0] == 0 and  detected_sequence_number[1] == 0:
                                sequence_number = 0
                        elif detected_sequence_number[0] == 0 and  detected_sequence_number[1] == 1:
                                sequence_number = 1
                        elif detected_sequence_number[0] == 1 and  detected_sequence_number[1] == 0:
                                sequence_number = 2
                        else:
                                sequence_number = 3
                        
                        print("detected sequence number: ", sequence_number)
                        
                        expected_sequence_number = (expected_sequence_number + 1)%4

                        if(expected_sequence_number != sequence_number):
                                print("packet loss detected!")
                                long_term_BER_counter = long_term_BER_counter + 1
                                long_term_BER_sum = long_term_BER_sum + 0.5
##                                print("expected sequence number: ", expected_sequence_number)
##                                print("actual sequence detected: ", sequence_number)
                                #cannot trust sequence number of a 0.5 BER "wrong packet"
                        else:
                                pass
##                                print("sequence number as expected")

##                        if (detected_sequence_number != expected_sequence_number):
##                                packet_loss = packet_loss + 1.0
##                                short_term_packet_loss = short_term_packet_loss + 1.0
####                                print("Packet loss detected!")

##                        expected_sequence_number = 1 - detected_sequence_number


                        #Extract the timestamp

                        timestamp_bin = a_demodulated[(len(a_demodulated)-17):(len(a_demodulated))]

                        timestamp_int = 0
                                                    
                        for i in range (17):
                                if(timestamp_bin[i] == 1):
                                        timestamp_int = timestamp_int+pow(2,16-i)
                                                    
                                

                        current_time = (time.time()*100)%100000

                        time_passed = (current_time-timestamp_int)/100.0

##                        if(debugging == 1):
##                                print("overall service time for this packet: ",time_passed)
                                
                        #Print the number of error bits in every packet
                        if(debugging == 1):
                                print()
                                print("Number of bit errors in this packet: ", packet_error_bits)
                                print()
                                print("Total number of bits in this packet: ", len(a_demodulated))
                                print()
                        packet_BER = packet_error_bits/len(a_demodulated)
                        print("packet BER: ", packet_BER)
##                        if(packet_BER !=0):
##                                plt.plot(np.real(data))
##                                plt.title('raw data for packet with error')
##                                plt.show()
##                                plt.plot(few_before)
##                                plt.title('few before the packet with error')
##                                plt.show()
##                        else:
##                                counter = counter + 1
##                                if(counter%5==4):
##                                        plt.plot(np.real(data))
##                                        plt.title('raw data for correct packet')
##                                        plt.show()
                                        
                        long_term_BER_counter = long_term_BER_counter + 1
                        long_term_BER_sum = long_term_BER_sum + packet_BER
                        #long term BER for the latest 100 packets
                        if (long_term_BER_counter%100 == 0):
                                long_term_BER = long_term_BER_sum/long_term_BER_counter
                                long_term_BER_counter = 0
                                long_term_BER_sum = 0
                                print("long term BER (past 100 packets, nonoverlapping): ", long_term_BER)
                                if(long_term_BER!=0):
                                        print("Long term BER nonzero!!!")


                        #Visualization
                        if(visualize == 1):
                                visualize_counter = visualize_counter + 1
                                #send the symbols for visualize in the constellation diagram
                                if(visualize_counter == frame_rate_factor):
                                        #cast packet_BER and long term BER to scientific notation
                                        packet_BER_str = np.format_float_scientific(packet_BER,precision=4)
                                        long_term_BER_str = np.format_float_scientific(long_term_BER,precision=4)
                                        q3.put([baseband_symbols,packet_BER_str,long_term_BER_str])
                                        visualize_counter = 0


                #Print out throughput information
   
                window_counter = window_counter+1
                total_kbits = total_kbits+packet_kbits


                if(window_counter == 20):#another x packets received
                        window_throughput = total_kbits/(time.time()-start)
##                        if(q8.qsize()<2):
##                            q8.put(window_throughput)
##                        print("queue size at source: ", q8.qsize())
                        if(debugging == 0):
                                print("Average Rx Throughput in kbps:",window_throughput)
##                        print("Average Rx Throughput in kbps:",window_throughput)
                        window_counter = 0
