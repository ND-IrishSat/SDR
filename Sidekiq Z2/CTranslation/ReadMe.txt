ReadMe.txt


Report: 1/28/2024
The main.c file containes the py CRC library translation. Currently, the script needs revised such that calloc() is not used, this way we can avoid memory mismanagement.
- Rylan Paul

Report: 2/8/2024
- Reorganization and abstraction in main.c is definitely needed
- Functions needed
        - np.fft.fft
        - np.fft.fftshift
        - scipy.signal.fftconvolve
        - np.arange
        - np.linspace --> 
        - int np.argmax --> gets the index of the max value in an array or real numbers
        - double* np.abs --> absolute value of an array, with complex number support
- pulse_shaping.c needs finished
- Working on exporting data as CSV to be graphed in python with matplotlib
- Rylan Paul

Report: 2/10/2024
- I've now added np.linspace, np.arange, np.fft.fftshift, np.fft.fft, np.abs
- Reorganization and abstraction in everything is still needed
  DONE  - np.sinc
  DONE  - np.hamming
  DONE  - np.sum
  DONE  - add arrays
  DONE  - subtract arrays
  DONE  - multiply arrays
  DONE  - divide arrays
  DONE  - np.mean
  DONE  - np.convolve(mode="full")
  DONE  - np.convolve(mode="same")
  DONE  - signal.resample_poly
  DONE  - signal.firwin
  np.convolve instead - scipy.signal.fftconvolve
        - demod.symbol_demod
        - CRC.CRCcheck (double check this works)
        - pulse_shaping.pulse_shaping

- Working on exporting data as CSV to be graphed in python with matplotlib
- Rylan Paul

Funtcions to add:

