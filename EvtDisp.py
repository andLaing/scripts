#!/usr/bin/env python
import sys
import numpy as np
from scipy import signal
import matplotlib.pyplot as plt
import tables as tb


#grp1 = [ 1, 4, 5, 8, 9 ]
#grp2 = [ 18, 19, 22, 23, 30, 31 ]
#grp1 = [ 0, 1, 4, 5, 6, 7 ]
#grp2 = [ 12, 13, 16, 17, 18, 19 ]
#grp1 = [ 0, 1, 12, 13, 16, 17 ]
#grp2 = [ 6, 7, 18, 19 ]
#grp3 = [ 4, 5 ]
grp1 = [ 0,  1, 4, 16, 17, 24 ]
grp2 = [ 12, 13, 28, 29 ]
grp3 = [ 8, 9 ]

def main():

    """ Display pmt waveforms (filtered if required) and FFTs """

    #runNo = sys.argv[1]
    #fileNo = sys.argv[2]
    #fName = 'hdf5/'+runNo+'/dst_waves.gdcsnext.'+fileNo.zfill(3)+'_'+runNo+'.root.h5'
    fName = sys.argv[1]

    ## Filter definition
    fSample = 40E6
    freqLPF = 100E3
    freqLPFd = 2*freqLPF / fSample
    b, a = signal.butter(1, freqLPFd, 'low', analog=False)
    ##
    ftAbsSave = np.empty(1)
    ##
    with tb.open_file( fName, 'r' ) as dataF:
    
        npm = len(dataF.root.Sensors.DataPMT)
        nevt = len(dataF.root.RD.pmtrwf)

        fig, axes = plt.subplots(nrows=3, ncols=2, figsize=(20,6))
        fig.show()

        for ievt in range(nevt):
            ## clear the axies
            for ax in axes.flatten():
                ax.cla()

            grp1Setter = False
            grp2Setter = False
            grp3Setter = False
            for ipm in range(npm):

                rwf = dataF.root.RD.pmtrwf[ievt][ipm]
                rwf = rwf - np.mean(rwf)
                # Filter
                fwf = signal.lfilter(b, a, rwf)
                ## remove mean again just in case
                fwf = fwf - np.mean(fwf)

                ## FFT
                ft = np.fft.rfft(fwf)
                freq = np.fft.rfftfreq(len(fwf), d=25E-9)
                ft_mag = np.absolute(ft)
                ##
                
                pmID = dataF.root.Sensors.DataPMT[ipm][0]

                if not grp1Setter and  pmID in grp1:
                    ftAbsSave = ft_mag
                elif not grp2Setter and pmID in grp2:
                    ftAbsSave2 = ft_mag
                elif not grp3Setter and pmID in grp3:
                    ftAbsSave3 = ft_mag

                if pmID in grp1:
                    axes[0][0].plot(fwf)
                    axes[1][1].plot(freq, ft_mag, label='pmt '+str(pmID))
                    #print([ x for x in np.select([ft_mag>=4000],[freq]) if x != 0.])
                    print('Mags: ', pmID, [ (x, y) for (x, y) in zip(freq, ft_mag) if x >= 15000 and x < 25000 ])
                    #axes[2][1].plot(freq, ft_mag/ftAbsSave)
                    if not grp1Setter:
                        axes[0][0].set_title('Group 1 filtered waveforms')
                        axes[1][1].set_title('Group 1 filtered FFT magnitudes')
                        axes[1][1].set_xlim(500, 25000)
                        #axes[2][1].set_title('Group 1 filtered FFT magnitudes normalised')
                        #axes[2][1].set_xlim(500, 25000)
                        #axes[2][1].set_ylim(0, 2)
                        grp1Setter = True
                    #else:
                        #print(pmID, [ x for (x, y) in zip(freq, ft_mag/ftAbsSave) if x > 4000 and x < 10000 and y >= 0.9 and y <= 1.1])
                elif pmID in grp2:
                    axes[0][1].plot(fwf)
                    axes[2][0].plot(freq, ft_mag, label='pmt '+str(pmID))
                    ## print('Mags: ', pmID, [ (x, y) for (x, y) in zip(freq, ft_mag) if x < 10000 ])
                    if not grp2Setter:
                        axes[0][1].set_title('Group 2 filtered waveforms')
                        axes[2][0].set_title('Group 2 filtered FFT magnitudes')
                        axes[2][0].set_xlim(500, 25000)
                        grp2Setter = True
                    #else:
                        #print(pmID, [ x for (x, y) in zip(freq, ft_mag/ftAbsSave2) if x > 500 and x < 10000 and y >= 0.9 and y <= 1.1])
                elif pmID in grp3:
                    axes[1][0].plot(fwf)
                    axes[2][1].plot(freq, ft_mag, label='pmt '+str(pmID))
                    if not grp3Setter:
                        axes[1][0].set_title('Group 3 filtered waveforms')
                        axes[2][1].set_title('Group 3 filtered FFT magnitudes')
                        axes[2][1].set_xlim(500, 25000)
                        grp2Setter = True
                #if pmID == grp1[0] or pmID == grp2[0]:
                    #axes[1][0].plot(fwf)
                    #if pmID == grp1[0]:
                     #   axes[1][0].set_title('WF from each group')

            plt.tight_layout()
            plt.draw()
            catcher = input("next plot?")
            if catcher == 'q':
                exit()
            elif catcher == 'save':
                plt.savefig('eventDisplay_'+str(ievt)+'.pdf')
            plt.cla()

if __name__ == '__main__':
    main()
