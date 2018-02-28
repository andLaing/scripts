#!/usr/bin/env python
import matplotlib
import numpy as np
from scipy import signal
import matplotlib.pyplot as plt
import sys
import tables as tb


pms = [ 1, 4, 5, 8, 9, 18, 19, 22, 23, 30, 31 ]

def main():
    """ a look at the spectrogram (sFFT) of PMT data """

    meanAbs = []
    ## Filter definition
    fSample = 40E6
    freqLPF = 100E3
    freqLPFd = 2*freqLPF / fSample
    b, a = signal.butter(1, freqLPFd, 'low', analog=False)
    ##
    runs = []
    nfile = len(sys.argv)-1
    for i in xrange(nfile):
        with load_input(i) as dataF:
            runs.append(sys.argv[1+i][13:-8])
            lenWf = get_wf_len(dataF)
            meanAbs.append([ np.zeros(lenWf/2+1, np.float64) for j in pms ])
            if i == 0:
                frq_plot = np.empty(lenWf/2+1)

            npm = len(pms)#len(dataF.root.RD.pmtrwf[0])
            nevt = len(dataF.root.RD.pmtrwf)
            for ipm in range(npm):
                jpm = getPMid(dataF, ipm)    
                for ievt in range(nevt):
                    
                    sg = getWF(dataF, ipm, ievt)
                    
                    sg = sg - np.mean(sg)
                    ## Filtering
                    #sg = signal.lfilter(b, a, sg)
                    ##

                    ft = np.fft.rfft(sg)
                    freq = np.fft.rfftfreq(len(sg), d=25E-9)
                    if i == 0 and ipm == 0 and ievt == 0:
                        frq_plot = freq

                    meanAbs[i][jpm] += np.absolute(ft)/nevt
                    
                    #fq, t, sp = signal.spectrogram(sg, fs=1.0/25E-9, nperseg=32000, window='blackmanharris', noverlap=16000)

    fig, axes = plt.subplots(nrows=3, ncols=4)
    ## 
    for i in xrange(nfile):
        for j, (ax, row) in enumerate(zip(axes.flatten(), meanAbs[i])):
            ax.plot(frq_plot, row, label='run '+runs[i])
            if i == nfile-1:
                ax.set_xlim(0, 100000)
                ax.legend(loc=0)
                ax.set_title('Mean FFT abs. val, pmt '+str(pms[j]))

    fig.show()
    raw_input("Happy?")
                
    

def load_input(i):
    """ load and return the input file """

    fileN = 'hdf5/'+sys.argv[1+i]+'/dst_waves.gdcsnext.'+sys.argv[-1].zfill(3)+'_'+sys.argv[1+i]+'.root.h5'
    print(fileN)
    return tb.open_file( fileN, 'r' )


def get_wf_len(inF):
    """ check the length of the wf in samples """

    return len(inF.root.RD.pmtrwf[0][0])


def getWF(inF, ipm, ievt):
    """ Get a specific waveform from file """

    return inF.root.RD.pmtrwf[ievt][ipm]


def getPMid(inF, i):
    """ return the index for the pm for identification """

    return pms.index(inF.root.Sensors.DataPMT[i][0])

if __name__ == '__main__':
    main()
