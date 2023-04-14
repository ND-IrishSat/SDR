This is a project I did last year for the Communications Lab. It's kind of spaghetti code - most of it was adapted by and not written by me - but it can be used as a
reference for our project. The best way to get an idea of what's going on is to look at modemTest.py first and then all the other files. Don't bother with the .ipynb
files -- these are jupyter notebook files, so unless you have that set up you won't be able to open them.

INSTRUCTIONS FOR VERIFYING BEHAVIOR OF SYSTEM

To test modulation/demodulation by themselves:
1. Open and run modem_test.ipynb
2. User can input their own bit stream and select their desired modulation scheme
3. It will display the generated bits along with the demodulated bits, and will tell you if they match

To test modulation/demodulation with pulse shaping and matched filtering:
1. Open and run FinalModemTest.ipynb
2. User can input their own bit stream, modulation scheme, pulse shape, upsampling factor, packet ratio
3. It will display the generated bits along with the demodulated bits, and will tell you if they match

Easiest way to test full system:
1. Open ThroughputAvg.ipynb
2. Running this code will execute a battery of tests that will plot averaged throughput against various parameters
3. User can also call ModemTest function that takes scheme, length of preamble ones, upsamplign factor, packet ratio, pulse shaping filter length, packet length multiplier, and pulse shape as inputs and returns the throughput of the system.

If you do not want to open a jupyter notebook:
1. Run FullModemTest.py in terminal
2. It will take user inputs for scheme, length of preamble ones, upsamplign factor, packet ratio, pulse shaping filter length, packet length multiplier, and pulse shape
3. It will then return the throughput of the system and an IQ plot

Note* Group delay issue for rrc filter not fixed for QPSK or OOK, resulting in non-zero bit error. IQ plot will have intersymbol interference with rect filter, but there will be zero bit error. If you wish to see a clean IQ plot, but with non-zero bit error, then use rrc filter.*
