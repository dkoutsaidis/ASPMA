import numpy as np
import sys
sys.path.append('../../software/models/')
from scipy.signal import blackmanharris, triang
from scipy.fftpack import ifft
import math
import dftModel as DFT
import utilFunctions as UF

def sineModel_MultiRes(x, fs, w1 , w2, w3, N1, N2, N3, t, B1, B2, B3):
    """
    MultiResolution Analysis/synthesis of a sound using the sinusoidal model, without sine tracking
    x: input array sound, [w1,w2,w3]: 3 analysis windows, [N1,N2,N3]: 3 sizes of complex spectrum,
    t: threshold in negative dB, [B1,B2,B3]: 3 frequency bands
    returns y: output array sound
    """
	
    h1M1 = int(math.floor((w1.size+1)/2))                   # half analysis window 1 size by rounding
    h1M2 = int(math.floor(w1.size/2))                       # half analysis window 1 size by floor
    h2M1 = int(math.floor((w2.size+1)/2))                   # half analysis window 2 size by rounding  
    h2M2 = int(math.floor(w2.size/2))                       # half analysis window 2 size by floor
    h3M1 = int(math.floor((w3.size+1)/2))                   # half analysis window 3 size by rounding  
    h3M2 = int(math.floor(w3.size/2))                       # half analysis window 3 size by floor

    Ns = 512                                                # FFT size for synthesis (even)
    H = Ns/4                                                # Hop size used for analysis and synthesis
    hNs = Ns/2                                              # half of synthesis FFT size

    pin = max(hNs, h1M1, h2M1, h3M1)                        # init sound pointer in middle of biggest anal window       
    pend = x.size - pin                                     # last sample to start a frame

    fftbuffer = np.zeros(max(N1,N2,N3))                     # initialize buffer for FFT

    yw = np.zeros(Ns)                                       # initialize output sound frame
    y = np.zeros(x.size)                                    # initialize output array

    w1 = w1 / sum(w1)                                       # normalize analysis window 1
    w2 = w2 / sum(w2)                                       # normalize analysis window 2
    w3 = w3 / sum(w3)                                       # normalize analysis window 3

    sw = np.zeros(Ns)                                       # initialize synthesis window
    ow = triang(2*H)                                        # triangular window
    sw[hNs-H:hNs+H] = ow                                    # add triangular window
    bh = blackmanharris(Ns)                                 # blackmanharris window
    bh = bh / sum(bh)                                       # normalized blackmanharris window
    sw[hNs-H:hNs+H] = sw[hNs-H:hNs+H] / bh[hNs-H:hNs+H]     # normalized synthesis window

    while pin<pend:                                         # while input sound pointer is within sound 
    #-----analysis-----

        #same frames with different window sizes centered to pin.
	x1 = x[pin-h1M1:pin+h1M2]                                   # select frame 1
	x2 = x[pin-h2M1:pin+h2M2]                                   # select frame 2
	x3 = x[pin-h3M1:pin+h3M2]                                   # select frame 3

	mX1, pX1 = DFT.dftAnal(x1, w1, N1)                          # compute dft of frame 1
        mX2, pX2 = DFT.dftAnal(x2, w2, N2)                          # compute dft of frame 2
        mX3, pX3 = DFT.dftAnal(x3, w3, N3)                          # compute dft of frame 3

	ploc1 = UF.peakDetection(mX1, t)                            # detect locations of peaks of frame 1
        ploc2 = UF.peakDetection(mX2, t)                            # detect locations of peaks of frame 2
	ploc3 = UF.peakDetection(mX3, t)                            # detect locations of peaks of frame 3

	iploc1, ipmag1, ipphase1 = UF.peakInterp(mX1, pX1, ploc1)   # refine peak values of frame 1 by interpolation
	iploc2, ipmag2, ipphase2 = UF.peakInterp(mX2, pX2, ploc2)   # refine peak values of frame 2 by interpolation
	iploc3, ipmag3, ipphase3 = UF.peakInterp(mX3, pX3, ploc3)   # refine peak values of frame 3 by interpolation

	ipfreq1 = fs*iploc1/float(N1)                               # convert peak locations of frame 1 to Hertz
        ipfreq2 = fs*iploc2/float(N2)                               # convert peak locations of frame 2 to Hertz
        ipfreq3 = fs*iploc3/float(N3)                               # convert peak locations of frame 3 to Hertz
                
        #constracting final arrays according to frequency bands.
        finalfreq = []
        finalmag = []
        finalphase = []
        for i in range(ipfreq1.size):
            if (ipfreq1[i]>=0 and ipfreq1[i]<=B1):
                finalfreq.append(ipfreq1[i])
                finalmag.append(ipmag1[i])
                finalphase.append(ipphase1[i])
        for i in range(ipfreq2.size):
            if (ipfreq2[i]>B1 and ipfreq2[i]<=B2):
                finalfreq.append(ipfreq2[i])
                finalmag.append(ipmag2[i])
                finalphase.append(ipphase2[i])
        for i in range(ipfreq3.size):
            if (ipfreq3[i]>B2 and ipfreq3[i]<=B3):
                finalfreq.append(ipfreq3[i])
                finalmag.append(ipmag3[i])
                finalphase.append(ipphase3[i])

        finalfreq = np.array(finalfreq)
        finalmag = np.array(finalmag)
        finalphase = np.array(finalphase)
                
    #-----synthesis-----
	Y = UF.genSpecSines(finalfreq, finalmag, finalphase, Ns, fs)   # generate sines in the spectrum         
	fftbuffer = np.real(ifft(Y))                                   # compute inverse FFT
	yw[:hNs-1] = fftbuffer[hNs+1:]                                 # undo zero-phase window
	yw[hNs-1:] = fftbuffer[:hNs+1] 
	y[pin-hNs:pin+hNs] += sw*yw                                    # overlap-add and apply a synthesis window
	pin += H                                                       # advance sound pointer

    return y
