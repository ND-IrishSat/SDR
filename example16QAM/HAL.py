#Import signal processing libraries
import numpy as np
import matplotlib.pyplot as plt
from scipy import signal
import math

#Import ADALM Python libraries
import adi
import iio
import libm2k

#Debugging lib tools
import time

#PLL SPI helper functions (reused from Radiohound)        
def PLL_TX(byte1,byte2,byte3,byte4):
    data = bytearray([byte1,byte2,byte3,byte4])
    return data

def vcoSubsys2(k):
        return (0x6010 | (k << 7))


def identify_device():
        #Automatically determine the connected device using iio scan context API
        scanned_contexts = iio.scan_contexts()
        m2k_device = "ADALM-2000" #default device: ADALM-2000
        for value in scanned_contexts.values():
                if m2k_device in value:
                        device = 'm2k'
                        print("ADALM 2000 detected!")
                else:
                        device = 'pluto'
                        print("ADALM pluto detected!")
                break
        return device

def setup_tx_device_control(device, same_machine, ftuned, tx_gain, fs):

        aout = None
        spi_desc = None
        sdr = None
                
        if(device == 'm2k'):
        
                #Initialize the transmitter                   
                if(same_machine == True):                #running in the same machine, specified via the command line or the GUI
                        print("Running programs in the same machine!")
                        uri = libm2k.getAllContexts()
                        sorted_uri = sorted(uri)
                        if(len(uri) == 2):
                                ctx = libm2k.m2kOpen(sorted_uri[0]) #two devices available, the first plugged-in device is treated as the TX
                        else:
                                ctx = libm2k.m2kOpen()   #only one device (the Tx) available, the Rx is already running             
                else:                           
                        ctx = libm2k.m2kOpen()           #running in different machines, choose the only available device

                if ctx is None:
                        print("Connection Error: No ADALM2000 device available/connected to your PC.")
                        
                ctx.calibrateADC()
                ctx.calibrateDAC()

                ps = ctx.getPowerSupply()
                ps.enableChannel(0, True)

                aout = ctx.getAnalogOut()

                if(ftuned != None):

                    # setup SPI
                    m2k_spi_init = libm2k.m2k_spi_init()
                    m2k_spi_init.clock = 6 # SCLK for HMC833
                    m2k_spi_init.mosi = 14  # SDI for HMC833
                    m2k_spi_init.miso = 5  # SDO for HMC833
                    m2k_spi_init.bit_numbering = libm2k.MSB
                    m2k_spi_init.cs_polarity = libm2k.ACTIVE_HIGH 
                    m2k_spi_init.context = ctx

                    spi_init_param = libm2k.spi_init_param()
                    #Specify the SPI interface clock speed
                    # maximum of 4 MHz clock used for HMC833 as in RH main.c
                    # 6 MHz used for the USB interfaceboard
                    spi_init_param.max_speed_hz = 1000000

                    spi_init_param.mode = libm2k.SPI_MODE_0
                    spi_init_param.chip_select = 15 #SEN for HMC833
                    spi_init_param.extra = m2k_spi_init

                    spi_desc = libm2k.spi_init(spi_init_param)

                    ps.pushChannel(0, 3.3)

                    print("spi interface initialization completed!")
                    print("Starting PLL initialization...")

                    libm2k.spi_write_and_read(spi_desc, PLL_TX(0b00001110,0b00000000,0b00000010,0b10011010))       #&H7 &H14D
                    libm2k.spi_write_and_read(spi_desc, PLL_TX(0b00011000,0b00000000,0b00000000,0b00000000))       #&HC &H0
                    libm2k.spi_write_and_read(spi_desc, PLL_TX(0b00011000,0b00000000,0b00000000,0b00000000))       #&HC &H0
                    libm2k.spi_write_and_read(spi_desc, PLL_TX(0b00000000,0b00010100,0b11110010,0b11101010))       #&H0 &HA7975
                    libm2k.spi_write_and_read(spi_desc, PLL_TX(0b00000010,0b00000000,0b00000000,0b00000100))       #&H1 &H2
                    libm2k.spi_write_and_read(spi_desc, PLL_TX(0b00000100,0b00000000,0b00000000,0b00000010))       #&H2 &H5 reference Divider = 1
                    libm2k.spi_write_and_read(spi_desc, PLL_TX(0b00000110,0b00000000,0b00000000,0b01111000))       #&H3 &H3C
                    libm2k.spi_write_and_read(spi_desc, PLL_TX(0b00001010,0b00000000,0b00101100,0b01010000))       #&H5 &H1628
                    libm2k.spi_write_and_read(spi_desc, PLL_TX(0b00001010,0b00000000,0b11000001,0b01000000))       #&H5 &H60A0
                    libm2k.spi_write_and_read(spi_desc, PLL_TX(0b00001010,0b00000000,0b00000000,0b00000000))       #&H5 &H0
                    libm2k.spi_write_and_read(spi_desc, PLL_TX(0b00001100,0b00000110,0b00011110,0b10010100))       #&H6 &H30F4A
                    libm2k.spi_write_and_read(spi_desc, PLL_TX(0b00001110,0b00000000,0b01000010,0b10011010))       #&H7 &H429A Will try relocking if Lock Detect fails (only tries once)
                    libm2k.spi_write_and_read(spi_desc, PLL_TX(0b00010001,0b10000011,0b01111101,0b11111110))       #&H8 &HC1BEFF
                    libm2k.spi_write_and_read(spi_desc, PLL_TX(0b00010010,0b10101000,0b11111111,0b11111110))       #&H9 &H547FFF
                    libm2k.spi_write_and_read(spi_desc, PLL_TX(0b00010100,0b00000000,0b01000000,0b10001110))       #&HA &H22047 FSM = 50MHz (Ref clock), 256 measurement periods for Vtune accuracy of +/-98KHz
                    libm2k.spi_write_and_read(spi_desc, PLL_TX(0b00010110,0b00001111,0b10000000,0b01000010))       #&HB &H7C021
                    libm2k.spi_write_and_read(spi_desc, PLL_TX(0b00011000,0b00000000,0b00000000,0b00000000))       #&HC &H0
                    libm2k.spi_write_and_read(spi_desc, PLL_TX(0b00011010,0b00000000,0b00000000,0b00000000))       #&HD &H0
                    libm2k.spi_write_and_read(spi_desc, PLL_TX(0b00011100,0b00000000,0b00000000,0b00000000))       #&HE &H0
                    libm2k.spi_write_and_read(spi_desc, PLL_TX(0b00011110,0b00000000,0b00000001,0b00000010))       #&HF &H102
                    libm2k.spi_write_and_read(spi_desc, PLL_TX(0b00100000,0b00000000,0b00000001,0b01100100))       #&H10 &HB2
                    libm2k.spi_write_and_read(spi_desc, PLL_TX(0b00100010,0b00010000,0b00000000,0b00000100))       #&H11 &H80002
                    libm2k.spi_write_and_read(spi_desc, PLL_TX(0b00100100,0b00000000,0b00000000,0b00000110))       #&H12 &H3
                    libm2k.spi_write_and_read(spi_desc, PLL_TX(0b00100110,0b00000000,0b00100100,0b10110010))       #&H13 &H1259
                    libm2k.spi_write_and_read(spi_desc, PLL_TX(0b00001000,0b00110011,0b00110011,0b00110100))       #&H4 &H19999A
                    libm2k.spi_write_and_read(spi_desc, PLL_TX(0b00011000,0b00000000,0b00000000,0b00000000))       #&HC &H0
                    libm2k.spi_write_and_read(spi_desc, PLL_TX(0b00011000,0b00000000,0b00000000,0b00000000))       #&HC &H0
                    libm2k.spi_write_and_read(spi_desc, PLL_TX(0b00011000,0b00000000,0b00000000,0b00000000))       #&HC &H0

                    print("PLL initialization completed!") 
                    print("Start frequency tuning...")
                    
                    freq = ftuned
                    f_osc = 50     #MHz
                    R_div = 1
                    f_pd = f_osc/R_div
                    R_div_0 = (R_div << 1) & 0b11111111    
                    R_div_1 = (R_div >> 7) & 0b11111111    
                    R_div_2 = (R_div >> 15) & 0b11111111
                    R_div_3 = (R_div >> 23) & 0b00000001
                    R_div_3 = R_div_3 | 0b00000100

                    if(freq>3000):
                        k1=0.5
                        Nint = math.floor(k1*freq/f_pd)
                        Nfrac = round((k1*freq/f_pd-Nint)*16777216)

                    else:
                        k = math.floor(3000/freq)
                        if(k>1):
                            k = k - (k%2)

                        Nint = math.floor(k*freq/f_pd)
                        Nfrac = round((k*freq/f_pd-Nint)*16777216)

                        Nfrac_0 = (Nfrac << 1) & 0b11111111
                        Nfrac_1 = (Nfrac >> 7) & 0b11111111
                        Nfrac_2 = (Nfrac >> 15) & 0b11111111
                        Nfrac_3 = (Nfrac >> 23) & 0b00000001
                        Nfrac_3 = Nfrac_3 | 0b00001000
                        Nint_0 = (Nint << 1) & 0b11111111
                        Nint_1 = (Nint >> 7) & 0b11111111
                        Nint_2 = (Nint >> 15) & 0b00001111

                    if(freq > 3000):
                        pll_freq = (f_pd)*(Nint+Nfrac/16777216.00)/k1

                        PLL_TX(0b00011000,0b00000000,0b00000000,0b00000000)                #&HC &H0
                        PLL_TX(0b00001010,0b00000000,0b11000001,0b00100000)                #&H5 &H6090 set divider accordingly
                        PLL_TX(0b00001010,0b00000001,0b01011010,0b00000000)                #&H5 &HAD00 Using High freq sub band
                        PLL_TX(0b00010010,0b10101000,0b11111111,0b11111110)                #&H9 &H547FFF //Charge Pump

                        PLL_TX(0b00011000,0b00000000,0b00000000,0b00000000)                #&HC &H0
                        PLL_TX(0b00001100,0b00000110,0b00011110,0b10010100)                #&H6 &H30F4A
                        PLL_TX(0b00000110,0b00000110,0b00011110,0b10010100)                #&H6 &H30F4A
                        PLL_TX(R_div_3,R_div_2,R_div_1,R_div_0)                            #&H2 &H1 // set reference divider to R_div

                        PLL_TX(0b00001010,0b00000000,0b01000000,0b00110000)                #&H5 &H2018 Set to doubler mode
                        PLL_TX(0b00001010,0b00000000,0b00000000,0b00000000)                #&H5 &H0
                        PLL_TX(0b00000110,Nint_2,Nint_1,Nint_0)                            # &H3 Set integer reg to freq
                        PLL_TX(0b00001010,0b00000000,0b00000000,0b00000000)                #&H5 &H0

                        PLL_TX(Nfrac_3,Nfrac_2,Nfrac_1,Nfrac_0)                            #&H4 Set Frequency Fractional Part
                        PLL_TX(0b00011000,0b00000000,0b00000000,0b00000000)                #&HC &H0

                    else:

                        pll_freq = f_pd*(Nint+Nfrac/16777216.00)/k
                        vcoSub = vcoSubsys2(k)
                        vcoSub_0 = (vcoSub << 1) & 0b11111111
                        vcoSub_1 = (vcoSub >> 7) & 0b11111111

                        libm2k.spi_write_and_read(spi_desc, PLL_TX(0b00011000,0b00000000,0b00000000,0b00000000))                #&HC &H0
                        libm2k.spi_write_and_read(spi_desc, PLL_TX(0b00001010,0b00000000,vcoSub_1,vcoSub_0))                    #&H5 &H6090 set divider accordingly
                        libm2k.spi_write_and_read(spi_desc, PLL_TX(0b00001010,0b00000000,0b01000010,0b00000000))                #&H5 &H2100 tuning default
                        libm2k.spi_write_and_read(spi_desc, PLL_TX(0b00010010,0b10101000,0b11111111,0b11111110))                #&H9 &H547FFF //Charge Pump
                        libm2k.spi_write_and_read(spi_desc, PLL_TX(0b00011000,0b00000000,0b00000000,0b00000000))                #&HC &H0
                        libm2k.spi_write_and_read(spi_desc, PLL_TX(0b00001100,0b00000110,0b00011110,0b10010100))                #&H6 &H30F4A //Sigma Delta Bypass for Integer mode
                        libm2k.spi_write_and_read(spi_desc, PLL_TX(0b00001100,0b00000110,0b00011110,0b10010100))                #&H6 &H30F4A
                        libm2k.spi_write_and_read(spi_desc, PLL_TX(R_div_3,R_div_2,R_div_1,R_div_0))                            #&H2 &H1 set reference divider to R_div
                        libm2k.spi_write_and_read(spi_desc,  PLL_TX(0b00001010,0b00000000,0b01010001,0b00110000))                #&H5 &H2898 Set to fundamental Mode
                        libm2k.spi_write_and_read(spi_desc, PLL_TX(0b00001010,0b00000000,0b00000000,0b00000000))                #&H5 &H0
                        libm2k.spi_write_and_read(spi_desc, PLL_TX(0b00000110,Nint_2,Nint_1,Nint_0))                            #&H3 Set integer reg to freq
                        libm2k.spi_write_and_read(spi_desc, PLL_TX(0b00001010,0b00000000,0b00000000,0b00000000))                #&H5 &H0
                        libm2k.spi_write_and_read(spi_desc, PLL_TX(Nfrac_3,Nfrac_2,Nfrac_1,Nfrac_0))                            #&H4 Set Frequency Fractional Part
                        libm2k.spi_write_and_read(spi_desc, PLL_TX(0b00011000,0b00011000,0b00000000,0b00000000))                #&HC &H0
                    
                    print("initial frequency tuned to (MHz) : ",ftuned)


                aout.setCyclic(False)
                aout.setSampleRate(0, fs)
                aout.enableChannel(0, True)             #sending data through channel 0
                aout.setSampleRate(1, fs)
                aout.enableChannel(1, True)             #sending data through channel 1


                #Set the kernel buffer count for the DAC, default is 4
                aout.setKernelBuffersCount(0,64)
                aout.setKernelBuffersCount(1,64)

                return aout, spi_desc, sdr

        if(device == 'pluto'):
                
                #Initialize the transmitter   
                if(same_machine == True):           #running in the same machine, specified via the command line or the GUI
                        print("Running programs in the same machine!")
                        scanned_contexts = iio.scan_contexts()
                        uri = list(scanned_contexts.keys())
                        sorted_uri = sorted(uri)

                        print("sorted uri: ",sorted_uri)
                        if(len(uri) == 2):
                                sdr = adi.Pluto(uri = sorted_uri[0]) #two devices available, the first plugged-in device is treated as the TX
                                print("my uri: ", sorted_uri[0])
                        else:
                                sdr = adi.Pluto()   #only one device (the Tx) available, the Rx is already running             
                else:                           
                        sdr = adi.Pluto()           #running in different machines, choose the only available device

                if sdr is None:
                        print("Connection Error: No Pluto device available/connected to your PC.")
    
                sdr.tx_lo = ftuned*1000000 #MHz to Hz
                sdr.sample_rate = fs
                sdr.tx_cyclic_buffer = False
                sdr.tx_rf_bandwidth = int(sdr.sample_rate*0.5) 
                sdr.gain_control_mode_chan0 = "manual"
                #Specify the Pluto tx gain, available values are: [-89.750000 0.250000 0.000000]
                sdr.tx_hardwaregain_chan0 = tx_gain

                return aout, spi_desc, sdr



def setup_rx_device_control(device, same_machine, ftuned, fs_in, total_samples,rx_gain):
    
        if(device == 'm2k'):
                #Initialize the receiver

                if(same_machine == True):                     #running in the same machine, specified via the command line or the GUI
                        print("Running programs in the same machine!")
                        uri = libm2k.getAllContexts()
                        sorted_uri = sorted(uri)
                        if(len(uri) == 2):                    #two available devices, indicating the Tx has not started yet
                                ctx = libm2k.m2kOpen(sorted_uri[1])  #the second plugged-in device is treated as the RX
                        else:
                                ctx = libm2k.m2kOpen()        #only one device available (since the Tx has already started), use the only available device as Rx                        
                else:
                        ctx = libm2k.m2kOpen()                #running in different machines, choose the only available device

                if ctx is None:
                        print("Connection Error: No ADALM2000 device available/connected to your PC.")

                ctx.calibrateADC()
                ctx.calibrateDAC()

                ps = ctx.getPowerSupply()
                ps.enableChannel(0, True)
                #ps.pushChannel(0, Vtune)

                #Initialize the receiver

                if(same_machine == True):                     #running in the same machine, specified via the command line or the GUI
                        print("Running programs in the same machine!")
                        uri = libm2k.getAllContexts()
                        sorted_uri = sorted(uri)
                        if(len(uri) == 2):                    #two available devices, indicating the Tx has not started yet
                                ctx = libm2k.m2kOpen(sorted_uri[1])  #the seco
                ain = ctx.getAnalogIn()
                ain.enableChannel(0,True)
                ain.enableChannel(1,True)
                ain.setSampleRate(fs_in)                   #input sampling rate (Hz) available value: 1000 10000 100000 1000000 10000000 100000000)
                ain.setRange(0,-2,2)
                ain.setRange(1,-2,2)
                ain.setKernelBuffersCount(48)

                aout = ctx.getAnalogOut()
                trig = ain.getTrigger()
                digital = ctx.getDigital()

                trig.setAnalogSource(0) 
                trig.setAnalogDelay(0)                      #trigger is centered
                trig.setAnalogMode(0, libm2k.ALWAYS)        #streaming mode

                if(ftuned != None):
                    # setup SPI
                    m2k_spi_init = libm2k.m2k_spi_init()
                    m2k_spi_init.clock = 6 # SCLK for HMC833
                    ##    digital.setOutputMode(1,1) #open drain: 0
                    m2k_spi_init.mosi = 14  # SDI for HMC833
                    m2k_spi_init.miso = 5  # SDO
                    m2k_spi_init.bit_numbering = libm2k.MSB
                    m2k_spi_init.cs_polarity = libm2k.ACTIVE_HIGH #ENABLE PIN for HMC833
                    m2k_spi_init.context = ctx
                    ##
                    spi_init_param = libm2k.spi_init_param()
                    ##    #Maximum of 4 MHz for HMC833 as in RH main.c
                    #6 MHz used for the USB interfaceboard
                    spi_init_param.max_speed_hz = 1000000

                    spi_init_param.mode = libm2k.SPI_MODE_0
                    spi_init_param.chip_select = 15
                    spi_init_param.extra = m2k_spi_init

                    spi_desc = libm2k.spi_init(spi_init_param)


                    ps.pushChannel(0, 3.3)

                    print("spi interface initialization completed!")
                    print("starting the PLL board initialization..")

                    libm2k.spi_write_and_read(spi_desc, PLL_TX(0b00001110,0b00000000,0b00000010,0b10011010))       #&H7 &H14D
                    libm2k.spi_write_and_read(spi_desc, PLL_TX(0b00011000,0b00000000,0b00000000,0b00000000))       #&HC &H0
                    libm2k.spi_write_and_read(spi_desc, PLL_TX(0b00011000,0b00000000,0b00000000,0b00000000))       #&HC &H0
                    libm2k.spi_write_and_read(spi_desc, PLL_TX(0b00000000,0b00010100,0b11110010,0b11101010))       #&H0 &HA7975
                    libm2k.spi_write_and_read(spi_desc, PLL_TX(0b00000010,0b00000000,0b00000000,0b00000100))       #&H1 &H2
                    libm2k.spi_write_and_read(spi_desc, PLL_TX(0b00000100,0b00000000,0b00000000,0b00000010))       #&H2 &H5 reference Divider = 1
                    libm2k.spi_write_and_read(spi_desc, PLL_TX(0b00000110,0b00000000,0b00000000,0b01111000))       #&H3 &H3C
                    libm2k.spi_write_and_read(spi_desc, PLL_TX(0b00001010,0b00000000,0b00101100,0b01010000))       #&H5 &H1628
                    libm2k.spi_write_and_read(spi_desc, PLL_TX(0b00001010,0b00000000,0b11000001,0b01000000))       #&H5 &H60A0
                    libm2k.spi_write_and_read(spi_desc, PLL_TX(0b00001010,0b00000000,0b00000000,0b00000000))       #&H5 &H0
                    libm2k.spi_write_and_read(spi_desc, PLL_TX(0b00001100,0b00000110,0b00011110,0b10010100))       #&H6 &H30F4A
                    libm2k.spi_write_and_read(spi_desc, PLL_TX(0b00001110,0b00000000,0b01000010,0b10011010))       #&H7 &H429A Will try relocking if Lock Detect fails (only tries once)
                    libm2k.spi_write_and_read(spi_desc, PLL_TX(0b00010001,0b10000011,0b01111101,0b11111110))       #&H8 &HC1BEFF
                    libm2k.spi_write_and_read(spi_desc, PLL_TX(0b00010010,0b10101000,0b11111111,0b11111110))       #&H9 &H547FFF
                    libm2k.spi_write_and_read(spi_desc, PLL_TX(0b00010100,0b00000000,0b01000000,0b10001110))       #&HA &H22047 FSM = 50MHz (Ref clock), 256 measurement periods for Vtune accuracy of +/-98KHz
                    libm2k.spi_write_and_read(spi_desc, PLL_TX(0b00010110,0b00001111,0b10000000,0b01000010))       #&HB &H7C021
                    libm2k.spi_write_and_read(spi_desc, PLL_TX(0b00011000,0b00000000,0b00000000,0b00000000))       #&HC &H0
                    libm2k.spi_write_and_read(spi_desc, PLL_TX(0b00011010,0b00000000,0b00000000,0b00000000))       #&HD &H0
                    libm2k.spi_write_and_read(spi_desc, PLL_TX(0b00011100,0b00000000,0b00000000,0b00000000))       #&HE &H0
                    libm2k.spi_write_and_read(spi_desc, PLL_TX(0b00011110,0b00000000,0b00000001,0b00000010))       #&HF &H102
                    libm2k.spi_write_and_read(spi_desc, PLL_TX(0b00100000,0b00000000,0b00000001,0b01100100))       #&H10 &HB2
                    libm2k.spi_write_and_read(spi_desc, PLL_TX(0b00100010,0b00010000,0b00000000,0b00000100))       #&H11 &H80002
                    libm2k.spi_write_and_read(spi_desc, PLL_TX(0b00100100,0b00000000,0b00000000,0b00000110))       #&H12 &H3
                    libm2k.spi_write_and_read(spi_desc, PLL_TX(0b00100110,0b00000000,0b00100100,0b10110010))       #&H13 &H1259
                    libm2k.spi_write_and_read(spi_desc, PLL_TX(0b00001000,0b00110011,0b00110011,0b00110100))       #&H4 &H19999A
                    libm2k.spi_write_and_read(spi_desc, PLL_TX(0b00011000,0b00000000,0b00000000,0b00000000))       #&HC &H0
                    libm2k.spi_write_and_read(spi_desc, PLL_TX(0b00011000,0b00000000,0b00000000,0b00000000))       #&HC &H0
                    libm2k.spi_write_and_read(spi_desc, PLL_TX(0b00011000,0b00000000,0b00000000,0b00000000))       #&HC &H0

                    print("PLL initialization completed!")
                    

                    #LO tuning and I/Q imbalance estimate (recalibration starts here)        
                    print("Starting the frequency tuning...")

                    freq = ftuned
                    ##    freq = freq - 110                              # Subtract IF frequency of 110 MHz. This is LO frequency of PLL
                    f_osc = 50     #MHz
                    R_div = 1
                    f_pd = f_osc/R_div
                    R_div_0 = (R_div << 1) & 0b11111111    # should be equal to 0b00000010
                    R_div_1 = (R_div >> 7) & 0b11111111    # should be equal to 0b10000000
                    R_div_2 = (R_div >> 15) & 0b11111111
                    R_div_3 = (R_div >> 23) & 0b00000001
                    R_div_3 = R_div_3 | 0b00000100

                    if(freq>3000):
                        k1=0.5
                        Nint = math.floor(k1*freq/f_pd)
                        Nfrac = round((k1*freq/f_pd-Nint)*16777216)

                    else:
                        k = math.floor(3000/freq)
                        if(k>1):
                            k = k - (k%2)

                        Nint = math.floor(k*freq/f_pd)
                        Nfrac = round((k*freq/f_pd-Nint)*16777216)

                        Nfrac_0 = (Nfrac << 1) & 0b11111111
                        Nfrac_1 = (Nfrac >> 7) & 0b11111111
                        Nfrac_2 = (Nfrac >> 15) & 0b11111111
                        Nfrac_3 = (Nfrac >> 23) & 0b00000001
                        Nfrac_3 = Nfrac_3 | 0b00001000
                        Nint_0 = (Nint << 1) & 0b11111111
                        Nint_1 = (Nint >> 7) & 0b11111111
                        Nint_2 = (Nint >> 15) & 0b00001111

                    if(freq > 3000):
                        pll_freq = (f_pd)*(Nint+Nfrac/16777216.00)/k1

                        PLL_TX(0b00011000,0b00000000,0b00000000,0b00000000)                #&HC &H0
                        PLL_TX(0b00001010,0b00000000,0b11000001,0b00100000)                #&H5 &H6090 set divider accordingly
                        PLL_TX(0b00001010,0b00000001,0b01011010,0b00000000)                #&H5 &HAD00 Using High freq sub band
                        PLL_TX(0b00010010,0b10101000,0b11111111,0b11111110)                #&H9 &H547FFF //Charge Pump

                        PLL_TX(0b00011000,0b00000000,0b00000000,0b00000000)                #&HC &H0
                        PLL_TX(0b00001100,0b00000110,0b00011110,0b10010100)                #&H6 &H30F4A
                        PLL_TX(0b00000110,0b00000110,0b00011110,0b10010100)                #&H6 &H30F4A
                        PLL_TX(R_div_3,R_div_2,R_div_1,R_div_0)                            #&H2 &H1 // set reference divider to R_div

                        PLL_TX(0b00001010,0b00000000,0b01000000,0b00110000)                #&H5 &H2018 Set to doubler mode
                        PLL_TX(0b00001010,0b00000000,0b00000000,0b00000000)                #&H5 &H0
                        PLL_TX(0b00000110,Nint_2,Nint_1,Nint_0)                            # &H3 Set integer reg to freq
                        PLL_TX(0b00001010,0b00000000,0b00000000,0b00000000)                #&H5 &H0

                        PLL_TX(Nfrac_3,Nfrac_2,Nfrac_1,Nfrac_0)                            #&H4 Set Frequency Fractional Part
                        PLL_TX(0b00011000,0b00000000,0b00000000,0b00000000)                #&HC &H0

                    else:

                        pll_freq = f_pd*(Nint+Nfrac/16777216.00)/k
                        vcoSub = vcoSubsys2(k)
                        vcoSub_0 = (vcoSub << 1) & 0b11111111
                        vcoSub_1 = (vcoSub >> 7) & 0b11111111

                        libm2k.spi_write_and_read(spi_desc, PLL_TX(0b00011000,0b00000000,0b00000000,0b00000000))                #&HC &H0
                        libm2k.spi_write_and_read(spi_desc, PLL_TX(0b00001010,0b00000000,vcoSub_1,vcoSub_0))                    #&H5 &H6090 set divider accordingly
                        libm2k.spi_write_and_read(spi_desc, PLL_TX(0b00001010,0b00000000,0b01000010,0b00000000))                #&H5 &H2100 tuning default
                        libm2k.spi_write_and_read(spi_desc, PLL_TX(0b00010010,0b10101000,0b11111111,0b11111110))                #&H9 &H547FFF //Charge Pump
                        libm2k.spi_write_and_read(spi_desc, PLL_TX(0b00011000,0b00000000,0b00000000,0b00000000))                #&HC &H0
                        libm2k.spi_write_and_read(spi_desc, PLL_TX(0b00001100,0b00000110,0b00011110,0b10010100))                #&H6 &H30F4A //Sigma Delta Bypass for Integer mode
                        libm2k.spi_write_and_read(spi_desc, PLL_TX(0b00001100,0b00000110,0b00011110,0b10010100))                #&H6 &H30F4A
                        libm2k.spi_write_and_read(spi_desc, PLL_TX(R_div_3,R_div_2,R_div_1,R_div_0))                            #&H2 &H1 set reftunedference divider to R_div
                        libm2k.spi_write_and_read(spi_desc,  PLL_TX(0b00001010,0b00000000,0b01010001,0b00110000))                #&H5 &H2898 Set to fundamental Mode
                        libm2k.spi_write_and_read(spi_desc, PLL_TX(0b00001010,0b00000000,0b00000000,0b00000000))                #&H5 &H0
                        libm2k.spi_write_and_read(spi_desc, PLL_TX(0b00000110,Nint_2,Nint_1,Nint_0))                            #&H3 Set integer reg to freq
                        libm2k.spi_write_and_read(spi_desc, PLL_TX(0b00001010,0b00000000,0b00000000,0b00000000))                #&H5 &H0
                        libm2k.spi_write_and_read(spi_desc, PLL_TX(Nfrac_3,Nfrac_2,Nfrac_1,Nfrac_0))                            #&H4 Set Frequency Fractional Part
                        libm2k.spi_write_and_read(spi_desc, PLL_TX(0b00011000,0b00011000,0b00000000,0b00000000))                #&HC &H0

                    print("frequency tuned to (MHz) : ",ftuned)

                    #Measure to estimate the DC offset at I and Q input
                    print("measuring I,Q input DC offsetl...")
                    dc_offset_data = ain.getSamples(1000000)
                    dc_offset_data_I = dc_offset_data[0]
                    dc_offset_data_Q = dc_offset_data[1]

                else: #wireline
                    spi_desc = None
                    dc_offset_data_I = None
                    dc_offset_data_Q = None
                
                sdr = None


        if(device == 'pluto'):
                #Initialize the receiver (Pluto)

                if(same_machine == True):                     #running in the same machine, specified via the command line or the GUI
                        print("Running programs in the same machine!")
                        scanned_contexts = iio.scan_contexts()
                        uri = list(scanned_contexts.keys())
                        print("uri list: ",uri)
                        sorted_uri = sorted(uri)
                        if(len(uri) == 2):                    #two available devices, indicating the Tx has not started yet
                                sdr = adi.Pluto(uri = sorted_uri[1])  #the second plugged-in device is treated as the RX
                                print("my uri: ",sorted_uri[1])
                        else:
                                sdr = adi.Pluto()             #only one device available (since the Tx has already started), use the only available device as Rx                        
                else:
                        sdr = adi.Pluto()                     #running in different machines, choose the only available device

                if sdr is None:
                        print("Connection Error: No Pluto device available/connected to your PC.")

                #kernel buffer count is 4 by default for ADALM (also can be implicitly verifid through this program)

                sdr._rxadc.set_kernel_buffers_count(8)

                sdr.rx_lo = 2375000000

                sdr.sample_rate = fs_in

                sdr.rx_rf_bandwidth = int(0.5*sdr.sample_rate)

                sdr.gain_control_mode_chan0 = "manual"

                #hardware gain rx available values [-3 1 71]
                sdr.rx_hardwaregain_chan0 = rx_gain
                
                
                #Measure to estimate the DC offset at I and Q input
                print("measuring I,Q input DC offsetl...")

                #adalm-pluto seems not allowing changing the buffer size within a program
                number_packets_in_buffer = 8

                sdr.rx_buffer_size = int(total_samples*number_packets_in_buffer) #to be adjusted
##                print("current buffer size: ", int(total_samples*number_packets_in_buffer))
                print("current buffer size: ", sdr.rx_buffer_size)
                
                dc_offset_data = sdr.rx()
                dc_offset_data_I = np.real(dc_offset_data)
                dc_offset_data_Q = np.imag(dc_offset_data)

                ain = None
                spi_desc = None

        return ain, spi_desc, sdr, dc_offset_data_I, dc_offset_data_Q



def freq_tune(device, spi_desc, sdr, freq):
        if(device == 'm2k'):
                f_osc = 50     #MHz
                R_div = 1
                f_pd = f_osc/R_div
                R_div_0 = (R_div << 1) & 0b11111111    
                R_div_1 = (R_div >> 7) & 0b11111111    
                R_div_2 = (R_div >> 15) & 0b11111111
                R_div_3 = (R_div >> 23) & 0b00000001
                R_div_3 = R_div_3 | 0b00000100

                if(freq>3000):
                    k1=0.5
                    Nint = math.floor(k1*freq/f_pd)
                    Nfrac = round((k1*freq/f_pd-Nint)*16777216)

                else:
                    k = math.floor(3000/freq)
                    if(k>1):
                        k = k - (k%2)

                    Nint = math.floor(k*freq/f_pd)
                    Nfrac = round((k*freq/f_pd-Nint)*16777216)

                    Nfrac_0 = (Nfrac << 1) & 0b11111111
                    Nfrac_1 = (Nfrac >> 7) & 0b11111111
                    Nfrac_2 = (Nfrac >> 15) & 0b11111111
                    Nfrac_3 = (Nfrac >> 23) & 0b00000001
                    Nfrac_3 = Nfrac_3 | 0b00001000
                    Nint_0 = (Nint << 1) & 0b11111111
                    Nint_1 = (Nint >> 7) & 0b11111111
                    Nint_2 = (Nint >> 15) & 0b00001111

                if(freq > 3000):
                    pll_freq = (f_pd)*(Nint+Nfrac/16777216.00)/k1

                    PLL_TX(0b00011000,0b00000000,0b00000000,0b00000000)                #&HC &H0
                    PLL_TX(0b00001010,0b00000000,0b11000001,0b00100000)                #&H5 &H6090 set divider accordingly
                    PLL_TX(0b00001010,0b00000001,0b01011010,0b00000000)                #&H5 &HAD00 Using High freq sub band
                    PLL_TX(0b00010010,0b10101000,0b11111111,0b11111110)                #&H9 &H547FFF //Charge Pump

                    PLL_TX(0b00011000,0b00000000,0b00000000,0b00000000)                #&HC &H0
                    PLL_TX(0b00001100,0b00000110,0b00011110,0b10010100)                #&H6 &H30F4A
                    PLL_TX(0b00000110,0b00000110,0b00011110,0b10010100)                #&H6 &H30F4A
                    PLL_TX(R_div_3,R_div_2,R_div_1,R_div_0)                            #&H2 &H1 // set reference divider to R_div

                    PLL_TX(0b00001010,0b00000000,0b01000000,0b00110000)                #&H5 &H2018 Set to doubler mode
                    PLL_TX(0b00001010,0b00000000,0b00000000,0b00000000)                #&H5 &H0
                    PLL_TX(0b00000110,Nint_2,Nint_1,Nint_0)                            # &H3 Set integer reg to freq
                    PLL_TX(0b00001010,0b00000000,0b00000000,0b00000000)                #&H5 &H0

                    PLL_TX(Nfrac_3,Nfrac_2,Nfrac_1,Nfrac_0)                            #&H4 Set Frequency Fractional Part
                    PLL_TX(0b00011000,0b00000000,0b00000000,0b00000000)                #&HC &H0

                else:

                    pll_freq = f_pd*(Nint+Nfrac/16777216.00)/k
                    vcoSub = vcoSubsys2(k)
                    vcoSub_0 = (vcoSub << 1) & 0b11111111
                    vcoSub_1 = (vcoSub >> 7) & 0b11111111

                    libm2k.spi_write_and_read(spi_desc, PLL_TX(0b00011000,0b00000000,0b00000000,0b00000000))                #&HC &H0
                    libm2k.spi_write_and_read(spi_desc, PLL_TX(0b00001010,0b00000000,vcoSub_1,vcoSub_0))                    #&H5 &H6090 set divider accordingly
                    libm2k.spi_write_and_read(spi_desc, PLL_TX(0b00001010,0b00000000,0b01000010,0b00000000))                #&H5 &H2100 tuning default
                    libm2k.spi_write_and_read(spi_desc, PLL_TX(0b00010010,0b10101000,0b11111111,0b11111110))                #&H9 &H547FFF //Charge Pump
                    libm2k.spi_write_and_read(spi_desc, PLL_TX(0b00011000,0b00000000,0b00000000,0b00000000))                #&HC &H0
                    libm2k.spi_write_and_read(spi_desc, PLL_TX(0b00001100,0b00000110,0b00011110,0b10010100))                #&H6 &H30F4A //Sigma Delta Bypass for Integer mode
                    libm2k.spi_write_and_read(spi_desc, PLL_TX(0b00001100,0b00000110,0b00011110,0b10010100))                #&H6 &H30F4A
                    libm2k.spi_write_and_read(spi_desc, PLL_TX(R_div_3,R_div_2,R_div_1,R_div_0))                            #&H2 &H1 set reference divider to R_div
                    libm2k.spi_write_and_read(spi_desc,  PLL_TX(0b00001010,0b00000000,0b01010001,0b00110000))                #&H5 &H2898 Set to fundamental Mode
                    libm2k.spi_write_and_read(spi_desc, PLL_TX(0b00001010,0b00000000,0b00000000,0b00000000))                #&H5 &H0
                    libm2k.spi_write_and_read(spi_desc, PLL_TX(0b00000110,Nint_2,Nint_1,Nint_0))                            #&H3 Set integer reg to freq
                    libm2k.spi_write_and_read(spi_desc, PLL_TX(0b00001010,0b00000000,0b00000000,0b00000000))                #&H5 &H0
                    libm2k.spi_write_and_read(spi_desc, PLL_TX(Nfrac_3,Nfrac_2,Nfrac_1,Nfrac_0))                            #&H4 Set Frequency Fractional Part
                    libm2k.spi_write_and_read(spi_desc, PLL_TX(0b00011000,0b00011000,0b00000000,0b00000000))                #&HC &H0

        else:
                sdr.tx_lo = int(freq)
                sdr.rx_lo = int(freq)
        return 0


def get_buffer(device,ain,sdr,buffer_size):
        if(device == 'm2k'):

##                print("no issue before get samples")
##                print("buffer_size", buffer_size)
                buffer = ain.getSamples(int(buffer_size))
##                print("no issue after get sample")
                return buffer

        else:
                pluto_data = sdr.rx()
                buffer = []
                I_data = np.real(pluto_data)
                buffer.append(I_data)
                Q_data = np.imag(pluto_data)
                buffer.append(Q_data)
                
                return buffer

def push_buffer(device,aout,sdr,buffer_data,scheme):
        if(device == 'm2k'):
                if(scheme=="OOK" or scheme=="BPSK"):
                        return aout.push(0,buffer_data)
                else:
                        return aout.push(buffer_data)
        else:
                buffer_IQ = buffer_data*2 ** 14
                return sdr.tx(buffer_IQ)

