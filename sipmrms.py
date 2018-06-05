#!/usr/bin/env python
import sys
import matplotlib.pyplot as plt
import numpy as np
import tables as tb
from glob import iglob


def sipm_rms():
    """
    Calculate rms for each sipm DICE
    """

    file_base = sys.argv[1]

    rmss = []

    evt_no = 0
    for file_name in iglob(file_base+'*.h5'):
        print('File: ', file_name)
        with tb.open_file(file_name) as dataIn:

            try:
                sipmrwf = dataIn.root.RD.sipmrwf
            except tb.NoSuchNodeError:
                print('No rw in file: ', file_name)
                continue

            if len(rmss) == 0:
                rmss = [ [] for i in range(sipmrwf.shape[1]//64) ]

            chNos = dataIn.root.Sensors.DataSiPM[:]['sensorID']
            for evtrwf in sipmrwf:
                evt_no += 1
                if evt_no % 100 == 0:
                    print('event count: ', evt_no)
                for chNo, rwf in zip(chNos, evtrwf):

                    rms = np.std(rwf, ddof=1)

                    rmss[chNo // 1000 - 1].append(rms)

    ## Plot distributions
    figs = []
    for idef in range(7):
        figs.append(plt.subplots(nrows=2,ncols=2, figsize=(10,10)))

        read_start = idef * 4
        read_end   = read_start + len(figs[-1][1])
        for j, (rms, ax) in enumerate(zip(rmss[read_start:read_end], figs[-1][1])):
            ax.hist(rms)
            ax.set_title('Pedestal RMS for DICE '+str(read_start + j + 1))
            ax.set_xlabel('Baseline RMS (ADC)')
            ax.set_ylabel('AU')

        plt.tight_layout()
        figs[-1][0].show()

    input('Hit enter when done with plots')
    

if __name__ == '__main__':
    sipm_rms()
