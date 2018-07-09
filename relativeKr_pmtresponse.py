import sys

import numpy  as np
import tables as tb

import matplotlib.pyplot as plt

from glob      import iglob
from functools import partial

from invisible_cities.io.pmaps_io          import load_pmaps_as_df
from invisible_cities.io.channel_param_io  import subset_param_reader as spr
import invisible_cities.core.fit_functions as fitf
from invisible_cities.icaro.hst_functions  import shift_to_bin_centers


def relative_pmt_response():
    """
    Script which uses pmaps (will be generalised in future to check XY dependence)
    to look at the relative response of the PMTs in Kr events and
    compares to the results on Poisson mu from calibrations
    """

    pmap_file_base = sys.argv[1]

    run_number = pmap_file_base.split('/')[2]

    s1hists = {x : [] for x in range(12)}
    s2hists = {x : [] for x in range(12)}
    s1sumh  = []
    s2sumh  = []

    for fn in iglob(pmap_file_base + '*.h5'):

        ## This version just using pmt databases
        s1df, s2df, _, s1pmtdf, s2pmtdf = load_pmaps_as_df(fn)

        for evt in s1pmtdf['event'].unique():
            s1evt  = s1pmtdf[s1pmtdf['event'] == evt]
            s2evt  = s2pmtdf[s2pmtdf['event'] == evt]
            s1sevt = s1df[s1df['event'] == evt]
            s2sevt = s2df[s2df['event'] == evt]
            for peak in s1evt['peak'].unique():
                s1peak = s1evt[s1evt['peak'] == peak]
                s1sumh.append(s1sevt[s1sevt['peak'] == peak]['ene'].sum())
                pmt1Q = s1peak[s1peak['npmt'] == 1]['ene'].sum()
                for pmt in s1peak['npmt'].unique():
                    if pmt != 1:
                        s1hists[pmt].append(s1peak[s1peak['npmt'] == pmt]['ene'].sum()/pmt1Q)
                    else:
                        s1hists[pmt].append(pmt1Q)
                        
            for peak in s2evt['peak'].unique():
                s2peak = s2evt[s2evt['peak'] == peak]
                s2sumh.append(s2sevt[s2sevt['peak'] == peak]['ene'].sum())
                pmt1Q = s2peak[s2peak['npmt'] == 1]['ene'].sum()
                for pmt in s2peak['npmt'].unique():
                    if pmt != 1:
                        s2hists[pmt].append(s2peak[s2peak['npmt'] == pmt]['ene'].sum()/pmt1Q)
                    else:
                        s2hists[pmt].append(pmt1Q)

    ## Make the plots
    s1bins = np.arange(-2, 10, 0.1)
    s2bins = np.arange(0, 1.2, 0.02)
    figs1, axess1 = plt.subplots(nrows=3, ncols=4, figsize=(20,6))
    s1pmt1 = np.array(s1hists[1])
    s1sumh = np.array(s1sumh)
    for (key, val), ax in zip(s1hists.items(), axess1.flatten()):
        if key == 1:
            #ax.hist(val, bins=100)
            ax.scatter(s1sumh[(s1sumh>2) & (s1sumh<15)], np.array(val)[(s1sumh>2) & (s1sumh<15)])
            ax.set_title('PMT 1 S1 charge distribution vs s1 sum charge')
            ax.set_xlabel('integrated charge in PMT sum (pe)')
            ax.set_ylabel('integrated charge in PMT1 (pe)')
            #ax.set_ylabel('AU)')
            #ax.set_xlabel('integrated charge in PMT1 (pe)')
        else:
            #ax.hist(val, bins=s1bins)
            ## ax.scatter(s1pmt1[np.abs(val) < 10], np.array(val)[np.abs(val) < 10])
            ax.scatter(s1sumh[(s1sumh>2) & (s1sumh<15)], np.array(val)[(s1sumh>2) & (s1sumh<15)])
            ax.set_title('PMT '+str(key)+' S1 relative charge distribution vs PMT sum charge')
            ax.set_xlabel('integrated charge in PMT sum (pe)')
            ax.set_ylabel('pmt q / pmt1 q')
            #ax.set_ylabel('AU')
            #ax.set_xlabel('pmt q / pmt1 q')
    plt.tight_layout()
    figs1.show()
    figs1.savefig('s1relativechargeCORRs1SUMFILT_R'+run_number+'.png')

    fitVals = {}
    figs2, axess2 = plt.subplots(nrows=3, ncols=4, figsize=(20,6))
    s2pmt1 = np.array(s2hists[1])
    s2sumh = np.array(s2sumh)
    for (key, val), ax in zip(s2hists.items(), axess2.flatten()):
        if key == 1:
            ax.set_title('PMT 1 S2 charge distribution vs s2 sum charge')
            ax.set_xlabel('integrated charge in PMT sum (pe)')
            ax.set_ylabel('integrated charge in PMT1 (pe)')
            #ax.set_ylabel('AU')
            #ax.set_xlabel('integrated charge in PMT1 (pe)')
            #ax.hist(val, bins=100)
            ax.scatter(s2sumh[(s2sumh>4000) & (s2sumh<12000)], np.array(val)[(s2sumh>4000) & (s2sumh<12000)])
        else:
            #ax.set_title('PMT '+str(key)+' S2 relative charge distribution vs PMT sum charge')
            #ax.set_xlabel('pmt1 q (pe)')
            ax.set_xlabel('integrated charge in PMT sum (pe)')
            ax.set_ylabel('pmt q / pmt1 q')
            #ax.set_ylabel('AU')
            #ax.scatter(s2pmt1[np.abs(val) < 10], np.array(val)[np.abs(val) < 10])
            ax.scatter(s2sumh[(s2sumh>4000) & (s2sumh<12000)], np.array(val)[(s2sumh>4000) & (s2sumh<12000)])
            #vals, bins, _ = ax.hist(val, bins=s2bins)
            ## limit fit to region with stat error <= 10% Poisson
            #useful_bins = np.argwhere(vals>=100)
            #b1 = useful_bins[0][0]
            #b2 = useful_bins[-1][0]
            #errs = np.sqrt(vals[b1:b2])
            #fvals = fitf.fit(fitf.gauss, shift_to_bin_centers(bins)[b1:b2], vals[b1:b2],
            #                 seed=(vals.sum(), bins[vals.argmax()], 0.1),
            #                 sigma=errs)
            #ax.plot(shift_to_bin_centers(bins),
            #        fitf.gauss(shift_to_bin_centers(bins), *fvals.values))
            #fitVals[key] = (fvals.values[1], fvals.values[2])
            #print('Fit PMT '+str(key), fvals.values, fvals.errors, fvals.chi2)
    plt.tight_layout()
    figs2.show()
    figs2.savefig('s2relativechargeCORRs2SUM_R'+run_number+'.png')

    ## figcal, axcal = plt.subplots()
    ## axcal.errorbar(list(fitVals.keys()),
    ##                np.fromiter((x[0] for x in fitVals.values()), np.float),
    ##                yerr=np.fromiter((x[1] for x in fitVals.values()), np.float),
    ##                label='Average response of PMTs to Kr relative to PMT 1')
    ## ## Get the calibration info for comparison.
    ## cal_files = [ fname for fname in sys.argv[2:] ]
    ## read_params = partial(spr, table_name='FIT_pmt_scaled_dark_pedestal',
    ##                       param_names=['poisson_mu'])
    ## ## Assumes ordering, ok?
    ## for i, fn in enumerate(cal_files):
    ##     cal_run = fn.split('_')[1]
    ##     with tb.open_file(fn) as cal_in:
    ##         pmt1Val = 0
    ##         pmt1Err = 0
    ##         cVals = []
    ##         cErrs = []
    ##         for sens, (pars, errs) in read_params(cal_in):
    ##             if sens != 1:
    ##                 cVals.append(pars['poisson_mu'])
    ##                 cErrs.append(errs['poisson_mu'])
    ##             else:
    ##                 pmt1Val = pars['poisson_mu']
    ##                 pmt1Err = errs['poisson_mu']
    ##         normVals = np.array(cVals) / pmt1Val
    ##         normErrs = normVals * np.sqrt(np.power(np.array(cErrs)/np.array(cVals), 2) +
    ##                                       np.power(pmt1Err/pmt1Val, 2))
    ##         axcal.errorbar(list(fitVals.keys()), normVals,
    ##                        yerr=normErrs, label='Calibration '+cal_run)
    ## axcal.legend()
    ## axcal.set_xlabel('PMT sensor ID')
    ## axcal.set_ylabel('Response relative to that of PMT 1')
    ## figcal.show()
    ## figcal.savefig('calPoisKrRelCompStats.png')
    input('plots good?')

                        
if __name__ == '__main__':
    relative_pmt_response()
