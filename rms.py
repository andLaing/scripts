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
    for fileName in iglob(sys.argv[1]+'*_waveforms.h5'):
        with tb.open_file(fileName, 'r') as dataF:

            npmt = len(dataF.root.Sensors.DataPMT)
            for rwf in data.root.RD.pmtrwf:
                for pm in range(npmt):
                    sensID = int(dataF.root.Sensors.DataPMT[pm][1])
                    
                    pedestals[sensID].append( np.mean(rwf) )
                    rmss[sensID].append( np.std(rwf, ddof=1) )

    fig0, axes0 = plt.subplots(nrows=4,ncols=3, figsize=(20,6))
    for i, (ax, means) in enumerate(zip(axes0.flatten(), pedestals)):
        ax.hist(means)
        ax.set_title('Pedestal for SensorID '+str(i))
        ax.set_xlabel('Baseline mean (ADC)')
        ax.set_ylabel('AU')
    fig0.show()

    fig1, axes1 = plt.subplots(nrows=4,ncols=3, figsize=(20,6))
    for j, (ax, stds) in enumerate(zip(axes1.flatten(), rmss)):
        ax.hist(stds)
        ax.set_title('Pedestal RMS for SensorID '+str(i))
        ax.set_xlabel('Baseline RMS (ADC)')
        ax.set_ylabel('AU')
    fig1.show()

if __name__ == '__main__':
    main()
