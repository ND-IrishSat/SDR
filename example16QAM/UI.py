import argparse
from argparse import RawTextHelpFormatter

#Create TX CLI parser
def create_tx_CLI_parser():
        parser = argparse.ArgumentParser(formatter_class=RawTextHelpFormatter)

        parser.add_argument("-m","--mode", type=int, help="Specify which mode you would like to run the transmitter program in:\n"
                            "Mode 1:Video streaming mode\n"
                            "Mode 2:Iperf test mode\n"
                            "Mode 3:Overhead test mode(outdated) \n"
                            "Mode 4:BER test mode\n"
                            "Press Crtl+Break to exit the program")

        parser.add_argument("-r","--ratio", type=int, help="Specify the size of each packet, the number of bits contained in one "
                            "packet size = ratio*segment_size. The ratios chosen need to be multiples of 10.")

        parser.add_argument("-o","--oversampling_factor", type=int, help="Each symbol transmitted is represented by OVERSAMPLING_FACTOR samples.")

        parser.add_argument("-fs","--sampling_rate", type=float, help="Specify the sampling rate of the DAC (output sampling rate)")

        parser.add_argument("-b","--bandwidth", help="Specify the speed (bandwidth) that iperf is sending UDP datagrams during the throughput test")

        parser.add_argument("-s","--same_machine", action='store_true', help="Run both the transmitter and receiver programs on the same machine. Make sure you plugged in the Tx first "
                            "and then the Rx so the programs know which device they should talk to")

        parser.add_argument("-p","--payload_length", type=int, help="Parameter for the overhead testing mode. Specify the length in bits of the packet payload.")

        parser.add_argument("-d","--debugging", action='store_true', help="Used by the developer to dubug and troubleshoot the system, additional stats info are printed.")

        parser.add_argument("-gain","--tx_gain", type=float, help="Specify the pluto TX gain in dB")
        
        return parser

#Create RX CLI parser
def create_rx_CLI_parser():
        parser = argparse.ArgumentParser(formatter_class=RawTextHelpFormatter)

        parser.add_argument("-m","--mode", type=int, help="Specify which mode you would like to run the OOK transmitter in:\n"
                            "Mode 1:Streaming mode\n"
                            "Mode 2:Iperf mode\n"
                            "Mode 3:Overhead test mode (outdated)\n"
                            "Mode 4:Error test mode\n"
                            "Press Crtl+Break to exit the program")

        parser.add_argument("-r","--ratio", type=int, help="Specify the size of each OOK packet, the number of bits contained in one "
                            "OOK packet length = ratio*segment_size. The ratios chosen need to be multiples of 10.")

        parser.add_argument("-o","--oversampling_factor", type=int, help="Each bit transmitted is represented by OVERSAMPLING_FACTOR samples of the OOK signals. The "
                            "oversamplig factors specified need to be multiples of 3")

        parser.add_argument("-fs","--sampling_rate", type=float, help="Specify the sampling rate of the ADC (Output sampling rate)")

        parser.add_argument("-s","--same_machine", action='store_true', help="Run both the transmitter and receiver programs on the same machine. Make sure you plugged in the Tx first "
                            "and then the Rx so the programs know which device they should talk to")

        parser.add_argument("-p","--payload_length", type=int, help="Parameter for the overhead testing mode. Specify the length in bits of the OOK packet payload.")

        parser.add_argument("-v","--visualized_bits", type=int, help="Visualize the received baseband signal on signal space constellation")

        parser.add_argument("-d","--debugging", action='store_true', help="Used by the developer to dubug and troubleshoot the system, additional stats info are printed.")

        parser.add_argument("-show","--showplot", action='store_true', help="Show constellation plot")

        parser.add_argument("-gain","--rx_gain", type=float, help="Specify the pluto RX gain in dB")

        return parser


#Create TX runtime GUI (showing system parameters and allowing manual LO tuning)
from tkinter import *
tx_current_freq = 2375  #default initial frequency of the LO
step_freq = 5 #default stepsize for LO frequency tuning
def tx_runtime_GUI(q7, fs, M, N, ftuned, preamble_length):
        def up():
                global tx_current_freq
                tx_current_freq = tx_current_freq + step_freq
                label3.config(text="Current Carrier Frequency:"+str(tx_current_freq),font = 'Helvetica 12 bold')
                q7.put("up")
        def down():
                global tx_current_freq
                tx_current_freq = tx_current_freq - step_freq
                label3.config(text="Current Carrier Frequency:"+str(tx_current_freq),font = 'Helvetica 12 bold')
                q7.put("down")
        base = Tk()
        base.title('RadioWare TX GUI')
        base.geometry('300x480')
        #Display the WI logo
        img = PhotoImage(file='wi.png')
        label4 = Label(image = img)
        label4.pack()

        label0 = Label(text="RadioWare Tx Configuration\n", font='Helvetica 12 bold', foreground="#070636")
        label0.pack()
        
        label1 = Label(
        text="DAC Sampling Rate (MHz): "+str(fs/1000000)+'\n'+"Initial TX Carrier Frequency (MHz): "+str(ftuned)+'\n'\
        +"Symbol rate (kBd): "+str(round(fs/M/1000,3))+'\n'+"Preamble length (bits): "+str(preamble_length)+'\n'+"Payload length (bits): "+str(N)+'\n',
        font = 4,
        foreground="black",  
        background="white"  
        )
        label1.pack()
        
        button1 = Button(
        text="Turn up TX carrier frequency",
        command = up,
        width=28,
        height=5,
        bg="blue",
        fg="white",
        )
        button1.pack()

        label2 = Label(
        text="\n"
        )
        label2.pack()

        button2 = Button(
        text="Turn down TX carrier frequency",
        command = down,
        width=28,
        height=5,
        bg="blue",
        fg="white",
        )
        button2.pack()

        label3 = Label(
        text="Current Carrier Frequency (MHz): "+str(ftuned),
        font = 'Helvetica 12 bold'
        )
        label3.pack()
        
        mainloop()



rx_current_freq = 2375
step_freq = 5

def rx_runtime_GUI(q7, q8, fs, M, N, ftuned, preamble_length):
        def up():
                global rx_current_freq
                rx_current_freq = rx_current_freq + step_freq
                label3.config(text="Current Carrier Frequency:"+str(rx_current_freq),font = 'Helvetica 12 bold')
                q7.put("up")
        def down():
                global rx_current_freq
                rx_current_freq = rx_current_freq - step_freq
                label3.config(text="Current Carrier Frequency:"+str(rx_current_freq),font = 'Helvetica 12 bold')
                q7.put("down")
                
        def update_throughput():
                if(q8.qsize()!=0):
                        current_throughput = q8.get()
                        label5.config(text="Current Throughpur (kbps): "+str(round(current_throughput,3)),font = 'Helvetica 12 bold')
    
        base = Tk()
        base.after(1000, update_throughput)
        base.title('RadioWare RX GUI')
        base.geometry('300x480')
        #Add the WI logo
        img = PhotoImage(file='wi.png')
        label4 = Label(image = img)
        label4.pack()

        label0 = Label(text="RadioWare Rx Configuration\n", font='Helvetica 12 bold', foreground="#070636")
        label0.pack()

        label1 = Label(
        text="ADC Sampling Rate (MHz): "+str(fs/1000000)+'\n'+"Initial RX Carrier Frequency (MHz): "+str(ftuned)+'\n'\
        +"Symbol rate (kBd): "+str(round(fs/M/1000,3))+'\n'+"Preamble length (bits): "+str(preamble_length)+'\n'+"Payload length (bits): "+str(N)+'\n'\
        +"Maximum throughput (kbps): "+str(round(fs/M/1000*4*(N/(N+preamble_length)),3))+'\n',
        font = 4,
        foreground="black",  # Set the text color to white
        background="white"  # Set the background color to black
        )
        label1.pack()

        label5 = Label(
        text="Current Throughput (kbps): ",
        font = 'Helvetica 12 bold'
        )
        label5.pack()

        button1 = Button(
        text="Turn up RX carrier frequency",
        command = up,
        width=25,
        height=5,
        bg="blue",
        fg="white",
        )
        button1.pack()

        label2 = Label(
        text="\n"
        )
        label2.pack()

        button2 = Button(
        text="Turn down RX carrier frequency",
        command = down,
        width=25,
        height=5,
        bg="blue",
        fg="white",
        )
        button2.pack()

        label3 = Label(
        text="Current Carrier Frequency (MHz):"+str(ftuned),
        font = 'Helvetica 12 bold'
        )
        label3.pack()

        #Listen to event (blocking) if there is no new item from the throughput statistics queue    
        update_throughput()
        mainloop()


#Create TX system parameter initialization GUI
import PySimpleGUI as sg
def create_tx_GUI_parser():
        sg.change_look_and_feel('DefaultNoMoreNagging')
        
        #Implement GUI
        layout = [[sg.Text('Operating Mode:')],      
                  [sg.InputText('1')],
                  [sg.Text('Packet ratio:')],
                  [sg.InputText('40')],
                  [sg.Text('Oversampling factor:')],
                  [sg.InputText('12')],
                  [sg.Text('Sampling rate:(Hz)')],
                  [sg.Radio('75000000',"radio 1", default=True),sg.Radio('7500000',"radio 1"),sg.Radio('750000',"radio 1"),
                   sg.Radio('75000',"radio 1"),sg.Radio('7500',"radio 1"),sg.Radio('750',"radio 1")],
                  [sg.Text('Iperf Bandwidth:')],
                  [sg.InputText('100k')],
                  [sg.Checkbox('Run programs in the same machine',default=True)],
                  [sg.Submit(), sg.Cancel()]]      

        window = sg.Window('OOK Transmitter', layout)
        return window

#Create RX system parameter initialization GUI
def create_rx_GUI_parser():
        sg.change_look_and_feel('DefaultNoMoreNagging')

        #Implement GUI
        layout = [[sg.Text('Operating Mode:')],      
                  [sg.InputText('1')],
                  [sg.Text('Packet ratio:')],
                  [sg.InputText('40')],
                  [sg.Text('Oversampling factor:')],
                  [sg.InputText('12')],
                  [sg.Text('Sampling rate:(Hz)')],
                  [sg.Radio('75000000',"radio 1", default=True),sg.Radio('7500000',"radio 1"),sg.Radio('750000',"radio 1"),
                   sg.Radio('75000',"radio 1"),sg.Radio('7500',"radio 1"),sg.Radio('750',"radio 1")],
                  [sg.Checkbox('Run programs in the same machine',default=True)],
                  [sg.Text('Packet payload length for the overhead test mode:')],
                  [sg.InputText('40')],
                  [sg.Checkbox('Visualize the received signals')],
                  [sg.Text('Number of bits to be visualized:')],
                  [sg.InputText('100')],

                  [sg.Submit(), sg.Cancel()]]      

        window = sg.Window('OOK Receiver', layout)
        return window
