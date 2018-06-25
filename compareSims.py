import os
import sys
import random
import tables as tb
import numpy  as np
from glob import iglob as gl

import matplotlib.pyplot as plt

import invisible_cities.reco.paolina_functions as plf
import invisible_cities.reco.peak_functions    as pf
import invisible_cities.reco.xy_algorithms     as xy
from   invisible_cities.reco.dst_functions     import load_xy_corrections
from   invisible_cities.reco.corrections       import LifetimeCorrection

from   invisible_cities.database import load_db as DB

from   invisible_cities.core.configure  import configure
from   invisible_cities.core.exceptions import SipmEmptyList
from   invisible_cities.core.exceptions import ClusterEmptyList
from   invisible_cities.core.exceptions import SipmZeroCharge

from   invisible_cities.io.pmaps_io     import load_pmaps_as_df

import invisible_cities.core.fit_functions     as fitf
from   invisible_cities.icaro.hst_functions    import hist


def compare_sims():
    """
    Just a copy (more or less) of Paola's notebook for
    simulation low level plots adjusted to compare simulation
    with and without the low frequency noise.
    """

    run_number      = -4735
    DataSiPM        = DB.DataSiPM(run_number)
    data_xs         = DataSiPM.X.values
    data_ys         = DataSiPM.Y.values

    ## conf parameters
    ## n_rebin = 1
    n_files = 10
    drift_v = 1.

    ## Containers for plots ##
    s1Q_NO       = []
    s1H_NO       = []
    s1W_NO       = []
    s1per_evt_NO = []
    s1Q_O        = []
    s1H_O        = []
    s1W_O        = []
    s1per_evt_O  = []
    #maxPMTtot_ratio = []

    s2per_evt_NO = []
    s2Q_NO       = []
    s2H_NO       = []
    s2W_NO       = []
    s2per_evt_O  = []
    s2Q_O        = []
    s2H_O        = []
    s2W_O        = []

    no_iter = gl("/Users/laingandrew/Documents/NEXT/IC_stuff/simdats/noOsc_paola/Tl_DISK_CLOSE_TO_ANODE_7bar_noiseSiPM_pmaps.*.h5")
    o_iter  = gl("/Users/laingandrew/Documents/NEXT/IC_stuff/simdats/Tl_DISK_CLOSE_TO_ANODE_7bar_noiseSiPM_pmaps.*.h5")
    for no_file, o_file in zip(no_iter, o_iter):
    ##for ifile in range(n_files):
     
        ##PATH_IN_NOOSC = "/Users/laingandrew/Documents/NEXT/IC_stuff/simdats/noOsc_paola/Tl_DISK_CLOSE_TO_ANODE_7bar_noiseSiPM_pmaps.{}.h5".format(ifile)
        ##PATH_IN_OSC   = "/Users/laingandrew/Documents/NEXT/IC_stuff/simdats/Tl_DISK_CLOSE_TO_ANODE_7bar_noiseSiPM_pmaps.{}.h5".format(ifile)
        print(no_file)
        #pmapNO_dict = load_pmaps(no_file)
        s1sNO, s2sNO, _, _, _ = load_pmaps_as_df(no_file)
        print('NO pmaps got')
        #pmapO_dict  = load_pmaps(o_file)
        s1sO, s2sO, _, _, _ = load_pmaps_as_df(o_file)
        print('O pmaps got')
        
        ## event_numbers_NO, timestamps_NO = get_event_numbers_and_timestamps_from_file_name(PATH_IN_NOOSC)
        ## event_numbers_O, timestamps_O   = get_event_numbers_and_timestamps_from_file_name(PATH_IN_OSC)
        event_numbers_NO = s2sNO['event'].unique()
        event_numbers_O  = s2sO['event'].unique()
        print(len(event_numbers_NO), len(event_numbers_NO))
        for evtNO, evtO in zip(event_numbers_NO, event_numbers_O):
            ## pmapNO = pmapNO_dict.get(evtNO, None)
            ## pmapO  = pmapO_dict.get(evtO, None)
            s1s_NO = s1sNO[s1sNO['event']==evtNO]
            s1s_O  = s1sO[s1sO['event']==evtO]
            s2s_NO = s2sNO[s2sNO['event']==evtNO]
            s2s_O  = s2sO[s2sO['event']==evtO]
            #print(len(s1s_NO), len(s1s_O), len(s2s_NO), len(s2s_O))
            #print('evtNO = ', evtNO, ' evtO = ', evtO,' pmaps got')
            if len(s1s_NO) == 0 or len(s1s_O) == 0 or len(s2s_NO) == 0 or len(s2s_O) == 0:
                continue
            ## if pmapNO:
            ##     s1s_NO    = pmapNO.s1s
            ##     s2sraw_NO = pmapNO.s2s
            ## else:
            ##     continue
            ## if pmapO:
            ##     s1s_O    = pmapO.s1s
            ##     s2sraw_O = pmapO.s2s
            ## else:
            ##     continue

            #if not s1s_NO or not s2sraw_NO or not s1s_O or not s2sraw_O: continue

            s1per_evt_NO.append(s1s_NO['peak'].nunique())
            s1per_evt_O .append(s1s_O['peak'].nunique())
            #print('n_s1_nosc = ', s1s_NO['peak'].nunique(), ', n_s1_osc = ', s1s_O['peak'].nunique())
            if s1s_NO['peak'].nunique() != 1: continue
            if s1s_O['peak'].nunique()  != 1: continue

            ## s2sNO = [rebin_peak(s2raw_NO, n_rebin) for s2raw_NO in s2sraw_NO]
            ## s2sO  = [rebin_peak(s2raw_O, n_rebin) for s2raw_O in s2sraw_O]

            ## s1tNO = s1s_NO[0].times
            ## s1eNO = s1s_NO[0].pmts
            ## S1tNO = s1t_NO[np.argmax(s1eNO)]
            ## s1tO  = s1s_O[0].times
            ## s1eO  = s1s_O[0].pmts
            ## S1tO  = s1t_O[np.argmax(s1eO)]
            #print('Aqui?')
            #fill histograms
            s1H_NO.append(s1s_NO['ene'].max())
            s1W_NO.append(s1s_NO['time'].iloc[-1] - s1s_NO['time'].iloc[0])
            s1Q_NO.append(s1s_NO['ene'].sum())
            s1H_O.append(s1s_O['ene'].max())
            s1W_O.append(s1s_O['time'].iloc[-1] - s1s_O['time'].iloc[0])
            s1Q_O.append(s1s_O['ene'].sum())

            s2per_evt_NO.append(s2s_NO['peak'].nunique())
            s2per_evt_O.append(s2s_O['peak'].nunique())

            for is2 in range(s2s_NO['peak'].nunique()):
                s2 = s2s_NO[s2s_NO['peak'] == is2]
                s2Q_NO.append(s2['ene'].sum())
                ## if s2['time'].iloc[-1] - s2['time'].iloc[0] < 0:
                ##     print(evtNO, ' has -ve width S2 in NO')
                s2W_NO.append(s2['time'].iloc[-1] - s2['time'].iloc[0])
                s2H_NO.append(s2['ene'].max())
            for js2 in range(s2s_O['peak'].nunique()):
                s2 = s2s_O[s2s_O['peak'] == js2]
                s2Q_O.append(s2['ene'].sum())
                ## if s2['time'].iloc[-1] - s2['time'].iloc[0] < 0:
                ##     print(evtO, ' has -ve width S2 in O')
                s2W_O.append(s2['time'].iloc[-1] - s2['time'].iloc[0])
                s2H_O.append(s2['ene'].max())

    print(min(s1per_evt_NO), max(s1per_evt_NO), min(s1per_evt_O), max(s1per_evt_O))
    print(min(s2per_evt_NO), max(s2per_evt_NO), min(s2per_evt_O), max(s2per_evt_O))
    print(min(s1Q_NO), max(s1Q_NO), min(s2Q_O), max(s2Q_O))
    print(min(s1H_NO), max(s1H_NO), min(s2H_O), max(s2H_O))
    print(min(s1W_NO), max(s1W_NO), min(s2W_NO), max(s2W_NO))
    print(min(s1W_O), max(s1W_O), min(s2W_O), max(s2W_O))

    fig, axes = plt.subplots(nrows=2, ncols=2, figsize=(20,6))
    axes[0][0].hist(s1per_evt_NO, bins=10, range=(0, 10), label='No. S1s no low frequency noise', log=True, histtype='step')
    axes[0][0].hist(s1per_evt_O , bins=10, range=(0, 10), label='No. S1s with low frequency noise', log=True, histtype='step')
    axes[0][0].legend()
    axes[0][0].set_xlabel('No. s1 per event')
    axes[0][0].set_ylabel('AU')
    
    axes[0][1].hist(s1H_NO, bins=200, range=(0, 200), label='S1 height no low frequency noise', log=True, histtype='step')
    axes[0][1].hist(s1H_O , bins=200, range=(0, 200), label='S1 height with low frequency noise', log=True, histtype='step')
    axes[0][1].legend()
    axes[0][1].set_xlabel('S1 height (pes)')
    axes[0][1].set_ylabel('AU')

    axes[1][0].hist(s1Q_NO, bins=100, range=(0, 1500), label='S1 charge no low frequency noise', log=True, histtype='step')
    axes[1][0].hist(s1Q_O , bins=100, range=(0, 1500), label='S1 charge with low frequency noise', log=True, histtype='step')
    axes[1][0].legend()
    axes[1][0].set_xlabel('S1 charge (pes)')
    axes[1][0].set_ylabel('AU')

    axes[1][1].hist(s1W_NO, bins=40, range=(0, 1000), label='S1 width no low frequency noise', log=True, histtype='step')
    axes[1][1].hist(s1W_O , bins=40, range=(0, 1000), label='S1 width with low frequency noise', log=True, histtype='step')
    axes[1][1].legend()
    axes[1][1].set_xlabel('S1 width (ns)')
    axes[1][1].set_ylabel('AU')
    fig.show()
    fig.savefig('s1_plots_pmaps_oscNoosc.png')

    fig2, axes2 = plt.subplots(nrows=2, ncols=2, figsize=(20,6))
    axes2[0][0].hist(s2per_evt_NO, bins=10, range=(0, 10), label='No. S2s no low frequency noise', log=True, histtype='step')
    axes2[0][0].hist(s2per_evt_O , bins=10, range=(0, 10), label='No. S2s with low frequency noise', log=True, histtype='step')
    axes2[0][0].legend()
    axes2[0][0].set_xlabel('No. s2 per event')
    axes2[0][0].set_ylabel('AU')
    
    axes2[0][1].hist(s2H_NO, bins=1000, range=(10, 65000), label='S2 height no low frequency noise', log=True, histtype='step')
    axes2[0][1].hist(s2H_O , bins=1000, range=(10, 65000), label='S2 height with low frequency noise', log=True, histtype='step')
    axes2[0][1].legend()
    axes2[0][1].set_xlabel('S2 height (pes)')
    axes2[0][1].set_ylabel('AU')

    axes2[1][0].hist(s2Q_NO, bins=1000, range=(10, 750000), label='S2 charge no low frequency noise', log=True, histtype='step')
    axes2[1][0].hist(s2Q_O , bins=1000, range=(10, 750000), label='S2 charge with low frequency noise', log=True, histtype='step')
    axes2[1][0].legend()
    axes2[1][0].set_xlabel('S2 charge (pes)')
    axes2[1][0].set_ylabel('AU')

    axes2[1][1].hist(s2W_NO, bins=235, range=(1500, 350000), label='S2 width no low frequency noise', log=True, histtype='step')
    axes2[1][1].hist(s2W_O , bins=235, range=(1500, 350000), label='S2 width with low frequency noise', log=True, histtype='step')
    axes2[1][1].legend()
    axes2[1][1].set_xlabel('S2 width (ns)')
    axes2[1][1].set_ylabel('AU')
    fig2.show()
    fig2.savefig('s2_plots_pmaps_oscNoosc.png')
    input('Done with plots?')


if __name__ == '__main__':
    compare_sims()
