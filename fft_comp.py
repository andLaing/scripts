#!/usr/bin/env python
import sys
import numpy as np
import tables as tb
import matplotlib.pyplot as plt
from scipy import signal
from invisible_cities.database import load_db



def main():
    """ version of fft_tests.py stripped down and adapted for MC wf """
    

    ## Filter definition
    fSample = 40E6
    freqLPF = 100E3
    freqLPFd = 2*freqLPF / fSample
    b, a = signal.butter(1, freqLPFd, 'low', analog=False)
    ##

    DataPMT = load_db.DataPMT(4723)
    
    defsDone = False
    with load_input(sys.argv[1]) as file0:

        npmt = len(DataPMT.ChannelID)
        if not defsDone:
            wf_len = get_wf_len(file0)
            meanAbs = [ np.zeros(int(wf_len/2+1), np.float64) for i in range(npmt) ]
            frq_plot = np.empty(int(wf_len/2+1))
            defsDone = True

        for i in range(npmt):
            for ievt in range(len(file0.root.RD.pmtrwf)):

                sg = file0.root.RD.pmtrwf[ievt][i]
                sg = sg - np.mean(sg)
                ## Filter wf
                sg = signal.lfilter(b, a, sg)

                ft = np.fft.rfft(sg)
                if i == 0 and ievt == 0:
                    frq_plot = np.fft.rfftfreq(len(sg), d=25E-9)
                meanAbs[i] += np.abs(ft)

            meanAbs[i] /= len(file0.root.RD.pmtrwf)

    fig, axes = plt.subplots(nrows=4,ncols=3, figsize=(20,6))
    for i, (ax, row) in enumerate(zip(axes.flatten(), meanAbs)):
        ax.plot(frq_plot, row)
        ax.set_title( 'Mean FFT magnitude, Eid '+str(DataPMT.ChannelID[i])+', Sensor '+DataPMT.PmtID[i])
        ax.set_xlim(300, 25000)
    plt.tight_layout()
    fig.show()
    fig.savefig('avFFTMag_simOsc2.png')

    input('happy with plots?')


def load_input(fName):
    """ load and return the input file """

    return tb.open_file( fName, 'r' )


def get_wf_len(inF):
    """ check the length of the wf in samples """

    return len(inF.root.RD.pmtrwf[0][0])


if __name__ == '__main__':
    main()
