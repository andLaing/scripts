#!/usr/bin/env python
import sys
import matplotlib.pyplot as plt
import numpy as np
import tables as tb
from glob import iglob
from datetime import datetime

from invisible_cities.icaro.hst_functions import shift_to_bin_centers


def main():
    """ Basic plot of event rate """

    # bin width in seconds
    bin_width = int(sys.argv[2])
    max_evt = np.inf
    if len(sys.argv) > 3:
        max_evt = int(sys.argv[3])
    hist_end = int(datetime.now().timestamp())
    bins_def = False
    bins = None
    stamp_hist = None
    nevt = 0
    true_last = 0
    for fileName in iglob(sys.argv[1]+'*_waveforms.h5'):
        if nevt < max_evt:
            with tb.open_file(fileName, 'r') as dataF:
                timestamps = np.fromiter((evt[1]/1000 for evt in dataF.root.Run.events), np.int)
                nevt += len(timestamps)
                if not bins_def:
                #start time in seconds
                    start_time = int(timestamps[0])
                    bins = np.arange(start_time, hist_end, bin_width)
                    stamp_hist = np.zeros(len(bins)-1, dtype=np.int)
                    bins_def = True
           
                # Ignore last bin as probably not equal width
                true_last = max(true_last, np.argwhere(bins>timestamps[-1])[0][0])-1
                stamp_hist += np.histogram(timestamps, bins)[0]

    print('check: ', datetime.fromtimestamp(shift_to_bin_centers(bins)[0]))
    ## Convert to human readable
    bins = shift_to_bin_centers(bins)[:true_last]
    dates = [ datetime.fromtimestamp(t) for t in bins ]
    plt.plot(dates, stamp_hist[:true_last]/bin_width)
    plt.title('Trigger rate in bins of '+str(bin_width/60)+' minutes')
    plt.xlabel('Timestamp (s)')
    plt.ylabel('events per second')
    plt.show()
            

if __name__ == '__main__':
    main()
