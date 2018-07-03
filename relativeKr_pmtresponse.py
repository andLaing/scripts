import sys

import numpy as np

import matplotlib.pyplot as plt

from glob import iglob

from invisible_cities.io.pmaps_io import load_pmaps_as_df


def relative_pmt_response():
    """
    Script which uses pmaps (will be generalised in future to check XY dependence)
    to look at the relative response of the PMTs in Kr events and
    compares to the results on Poisson mu from calibrations
    """

    pmap_file_base = sys.argv[1]

    s1hists = {x : [] for x in range(12)}
    s2hists = {x : [] for x in range(12)}

    for fn in iglob(pmap_file_base + '*.h5'):

        ## This version just using pmt databases
        _, _, _, s1pmtdf, s2pmtdf = load_pmaps_as_df(fn)

        for evt in s1pmtdf['event'].unique():
            s1evt = s1pmtdf[s1pmtdf['event'] == evt]
            s2evt = s2pmtdf[s2pmtdf['event'] == evt]
            for peak in s1evt['peak'].unique():
                s1peak = s1evt[s1evt['peak'] == peak]
                pmt1Q = s1peak[s1peak['npmt'] == 1]['ene'].sum()
                for pmt in s1peak['npmt']:
                    if pmt != 1:
                        s1hists[pmt].append(s1peak[s1peak['npmt'] == pmt]['ene'].sum()/pmt1Q)
                    else:
                        s1hists[pmt].append(pmt1Q)
            for peak in s2evt['peak'].unique():
                s2peak = s2evt[s2evt['peak'] == peak]
                pmt1Q = s2peak[s2peak['npmt'] == 1]['ene'].sum()
                for pmt in s2peak['npmt']:
                    if pmt != 1:
                        s2hists[pmt].append(s2peak[s2peak['npmt'] == pmt]['ene'].sum()/pmt1Q)
                    else:
                        s2hists[pmt].append(pmt1Q)

    ## Make the plots
    s1bins = np.arange(-2, 10, 0.1)
    s2bins = np.arange(0, 2000, 10)
    figs1, axess1 = plt.subplots(nrows=3, ncols=4, figsize=(20,6))
    for (key, val), ax in zip(s1hists.items(), axess1.flatten()):
        if key == 1:
            ax.title('PMT 1 S1 charge distribution')
            ax.xlabel('integrated charge (pe)')
            ax.ylabel('AU')
            ax.hist(val)
        else:
            ax.title('PMT '+str(key)+' relative charge distribution')
            ax.xlabel('pmt q / pmt1 q')
            ax.ylabel('AU')
            ax.hist(val, bins=s1bins)
    plt.tight_layout()
    figs1.show()
    figs2, axess2 = plt.subplots(nrows=3, ncols=4, figsize=(20,6))
    for (key, val), ax in zip(s2hists.items(), axess2.flatten()):
        if key == 1:
            ax.title('PMT 1 S1 charge distribution')
            ax.xlabel('integrated charge (pe)')
            ax.ylabel('AU')
            ax.hist(val)
        else:
            ax.title('PMT '+str(key)+' S2 relative charge distribution')
            ax.xlabel('pmt q / pmt1 q')
            ax.ylabel('AU')
            ax.hist(val, bins=s1bins)
    plt.tight_layout()
    figs2.show()
    input('plots good?')

                        
if __name__ == '__main__':
    relative_pmt_response()
