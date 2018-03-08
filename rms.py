#!/usr/bin/env python
import sys
import matplotlib.pyplot as plt
import numpy as np
import tables as tb
from glob import iglob


def main():
    """ Get the pedestal and rms event by event and histogram """

    pedestals = [ [] for pm in range(12) ]
    rmss      = [ [] for pm in range(12) ]
    ids = []
    defsDone = False
    for fileName in iglob(sys.argv[1]+'*_waveforms.h5'):
        with tb.open_file(fileName, 'r') as dataF:
            if not defsDone:
                ids = np.array(dataF.root.Sensors.DataPMT[:])
                defsDone = True
            npmt = len(ids)
            for rwf in dataF.root.RD.pmtrwf:
                for pm in range(npmt):
                    
                    pedestals[pm].append( np.mean(rwf[pm]) )
                    rmss[pm].append( np.std(rwf[pm], ddof=1) )

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

if __name__ == '__main__':
    main()
