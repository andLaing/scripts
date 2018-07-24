#!/usr/bin/env python
import sys
import matplotlib.pyplot as plt
import numpy as np
import tables as tb
from glob import glob

#import invisible_cities.core.fit_functions as fitf


def main(file_list=None):
    """ Get the pedestal and rms event by event and histogram """

    pedestals = [ [] for pm in range(12) ]
    rmss      = [ [] for pm in range(12) ]
    ids = []
    defsDone = False

    seven_hundred_mus = 700 / 0.025
    ten_mus           = 10 / 0.025
    if not file_list:
        file_list = glob(sys.argv[1]+'*_waveforms.h5')
    for fileName in file_list:
        with tb.open_file(fileName, 'r') as dataF:
            if not defsDone:
                ids = np.array(dataF.root.Sensors.DataPMT[:])
                defsDone = True
            npmt = len(ids)
            for rwf in dataF.root.RD.pmtrwf:
                for pm in range(npmt):

                    ## Try to emulate deconvolve_pmts
                    pedestals[pm].append( np.mean(rwf[pm]) )
                    rms = 0
                    for i in range(ten_mus):
                        rms += np.power(pedestals[pm][-1] - rwf[pm][i], 2)
                    rmss[pm].append( np.sqrt(rms / ten_mus) )
                    #pedestals[pm].append( np.mean(rwf[pm]) )
                    #rmss[pm].append( np.std(rwf[pm], ddof=1) )

    fig0, axes0 = plt.subplots(nrows=4,ncols=3, figsize=(20,6))
    for i, (ax, means) in enumerate(zip(axes0.flatten(), pedestals)):
        ax.hist(means)
        ax.set_title('Pedestal for SensorID '+str(i)+', ElecID '+str(ids[i][0]))
        ax.set_xlabel('Baseline mean (ADC)')
        ax.set_ylabel('AU')
        print('SensorID '+str(i)+' average pedestal = '+str(round(np.mean(means), 2))+'+/-'+str(round(np.std(means, ddof=1), 2)))
    plt.tight_layout()
    fig0.show()

    fig1, axes1 = plt.subplots(nrows=4,ncols=3, figsize=(20,6))
    for j, (ax, stds) in enumerate(zip(axes1.flatten(), rmss)):
        ax.hist(stds)
        ax.set_title('Pedestal RMS for SensorID '+str(j)+', ElecID '+str(ids[j][0]))
        ax.set_xlabel('Baseline RMS (ADC)')
        ax.set_ylabel('AU')
        print('SensorID '+str(j)+' average rms = '+str(round(np.mean(stds), 2))+'+/-'+str(round(np.std(stds, ddof=1), 2)))
    plt.tight_layout()
    fig1.show()

    input("ready to move on?")

    mean_mean = np.array(pedestals).mean(axis=1)
    std_mean  = np.array(pedestals).std(ddof=1, axis=1)
    mean_rms  = np.array(rmss).mean(axis=1)
    std_rms   = np.array(rmss).std(ddof=1, axis=1)

    return mean_mean, std_mean, mean_rms, std_rms

if __name__ == '__main__':
    main()
