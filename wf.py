#!/usr/bin/env python
import sys
import numpy as np
from scipy import signal
import matplotlib.pyplot as plt
from specto import load_input, getWF, get_wf_len


def main():
    """ Waveform plotter with low pass filter """


    ## Groups showing similar noise profile
    grp1 = [ 1, 4, 5, 8, 9 ]
    grp2 = [ 18, 19, 22, 23, 30, 31 ]
    #grp3 = [ 18, 19, 22, 23, 26, 27 ]
    with load_input(0) as dataF:

        npm = len(dataF.root.Sensors.DataPMT)#len(dataF.root.RD.pmtrwf[0])
        nevt = len(dataF.root.RD.pmtrwf)

        ## Filter definition
        fSample = 40E6
        freqLPF = 100E3
        freqLPFd = 2*freqLPF / fSample
        b, a = signal.butter(1, freqLPFd, 'low', analog=False)
        ##
        fig, axes = plt.subplots(nrows=3, ncols=2, figsize=(20,6))
        #fig.tight_layout()
        fig.show()
        wf_len = get_wf_len(dataF)
        if len(sys.argv) > 3:
            wf_len = wf_len/2+1      
        elif len(sys.argv) == 3:
            g1_first = np.zeros(wf_len, np.float64)
            g2_first = np.zeros(wf_len, np.float64)
            g3_first = np.zeros(wf_len, np.float64)
            mean_first = np.zeros(wf_len, np.float64)
        ##
        for ievt in range(nevt):
            ## clear the axies
            for ax in axes.flatten():
                ax.cla()
            plt_frq = np.zeros(wf_len, np.float64)
            fwf_mean = np.zeros(wf_len, np.float64)
            wf_mean = np.zeros(wf_len, np.float64) # No filter
            g1_mean = np.zeros(wf_len, np.float64)
            g2_mean = np.zeros(wf_len, np.float64)
            g3_mean = np.zeros(wf_len, np.float64)
            for ipm in range(npm):

                sg = getWF(dataF, ipm, ievt)
                sg = sg - np.mean(sg)

                sgf = signal.lfilter(b, a, sg)
                ## remove mean again just in case
                sgf = sgf - np.mean(sgf)
                #sgf = sg

                pmID = getPMid(dataF, ipm)

                if len(sys.argv) == 3:
                    axes[0][0].plot(sgf, label='pmt '+str(pmID))
                    fwf_mean += sgf/npm
                    wf_mean += sg/npm
                    if pmID in grp1:
                        g1_mean += sgf/len(grp1)
                    elif pmID in grp2:
                        g2_mean += sgf/len(grp2)
                    elif pmID in grp3:
                        g3_mean += sgf/len(grp3)
                else:
                    ft = np.fft.rfft(sgf)
                    freq = np.fft.rfftfreq(len(sgf), d=25E-9)
                    if ipm == 0:
                        plt_frq = freq
                    if sys.argv[2] == 'mag':
                        ft_mag = np.absolute(ft)
                        axes[0][0].plot(freq, ft_mag, label='pmt '+str(pmID))
                        fwf_mean += ft_mag/npm
                        if pmID in grp1:
                            g1_mean += ft_mag/len(grp1)
                        elif pmID in grp2:
                            g2_mean += ft_mag/len(grp2)
                        elif pmID in grp3:
                            g3_mean += ft_mag/len(grp3)
                    elif sys.argv[2] == 'phase':
                        ft_pha = np.angle(ft)
                        axes[0][0].plot(freq, ft_pha, label='pmt '+str(pmID))
                        fwf_mean += ft_pha/npm
                        if pmID in grp1:
                            g1_mean += ft_pha/len(grp1)
                        elif pmID in grp2:
                            g2_mean += ft_pha/len(grp2)
                        elif pmID in grp3:
                            g3_mean += ft_pha/len(grp3)
            
            
            ## The axes not set
            if len(sys.argv) == 3:
                axes[0][1].plot(g1_mean)
                axes[0][1].set_title('Group 1 mean waveform')
                axes[1][0].plot(g2_mean)
                axes[1][0].set_title('Group 2 mean waveform')
                axes[1][1].plot(g3_mean)
                axes[1][1].set_title('Group 3 mean waveform')
                axes[2][0].plot(fwf_mean)
                axes[2][0].set_title('Mean waveform')
                if ievt == 0:
                    g1_first = g1_mean
                    g2_first = g2_mean
                    g3_first = g3_mean
                    mean_first = fwf_mean
                else:
                    axes[0][1].plot(g1_first)
                    axes[1][0].plot(g2_first)
                    axes[1][1].plot(g3_first)
                    axes[2][0].plot(mean_first)
                axes[2][1].plot(wf_mean)
                axes[2][1].set_title('Mean waveform and corrected')
                axes[2][1].plot(wf_mean-fwf_mean)
                axes[2][1].set_xlim(0, 1000)
            else:
                axes[0][0].set_xlim(0,50000)
                axes[0][1].plot(plt_frq, g1_mean)
                axes[0][1].set_title('Group 1 mean '+sys.argv[2])
                axes[0][1].set_xlim(0,50000)
                axes[1][0].plot(plt_frq, g2_mean)
                axes[1][0].set_title('Group 2 mean '+sys.argv[2])
                axes[1][0].set_xlim(0,50000)
                axes[1][1].plot(plt_frq, g3_mean)
                axes[1][1].set_title('Group 3 mean '+sys.argv[2])
                axes[1][1].set_xlim(0,50000)
                axes[2][0].plot(plt_frq, fwf_mean)
                axes[2][0].set_title('Mean '+sys.argv[2])
                axes[2][0].set_xlim(0,50000)
            plt.draw()
            #fig.legend(loc=0)
            catcher = input("next plot?")
            if catcher == 'q':
                exit()
            plt.cla()


def getPMid(inF, i):
    """ return the pmt elecid for identification """

    return inF.root.Sensors.DataPMT[i][0]


if __name__ == '__main__':
    main()
