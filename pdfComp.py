#!/usr/bin/env python
import sys
import numpy as np
import tables as tb
import matplotlib.pyplot as plt
from glob import iglob
from invisible_cities.database import load_db as DB
from invisible_cities.core.random_sampling    import NoiseSampler as SiPMsNoiseSampler


def main():

    """ Plot some pdfs """

    in_base = sys.argv[1]#, sys.argv[2] ]#, sys.argv[3], sys.argv[4] ]
    ## assumes order mean data run, mean DC run, mode data run, mode DC run
    lbls = [ 'r5906 median ped', 'r5906 mode ped' ]#, 'r4821 mean ped', 'r4821 mode ped' ]

    sipmDats = DB.DataSiPM(5906)
    noise_sampler = SiPMsNoiseSampler(5906)

    binsS = 0
    sipmHist = [ 0, 0 ]#, 0, 0 ]
    defsDone = [ False, False ]#, False, False ]
    #for i in range(4):
    #for i in range(2):
    for fN in iglob(in_base+'*.h5'):
        with tb.open_file(fN, 'r') as data:
            if not defsDone[i]:
                    #if i == 0:
                binsS = np.array(data.root.HIST.sipm_bins)
                sipmHist[i] = np.zeros((len(data.root.HIST.sipm[0]),
                                        len(data.root.HIST.sipm[0][0])))
                defsDone[i] = True
                
                sipmHist[i] += data.root.HIST.sipm[0]

    ## comparisons between plots and saved pdfs
    badCh = [ id for id, act in zip(sipmDats.SensorID,sipmDats.Active) if act==0]
    channel = input('Which channel do you want to see? ')
    if channel in badCh:
        channel = input('Channel bad, pick another ')
    while not 'q' in channel:
        inx = sipmDats.SensorID[sipmDats.SensorID==int(channel)].index[0]

        for j, hst in enumerate(sipmHist):
            plt.bar(binsS, hst[inx], width=0.1, log=True, label=lbls[j])
        plt.title('PDF comparison channel '+channel)
        plt.xlabel('pe')
        plt.ylabel('AU')
        plt.legend()
        plt.show(block=False)
        #plt.savefig('SipmSpecCompCh'+channel+'.png')
        channel = input('Move to another channel or quit? ')
        plt.clf()
        plt.close()
    exit()


if __name__ == '__main__':
    main()
