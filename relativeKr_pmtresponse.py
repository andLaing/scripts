import sys

import numpy  as np
import tables as tb

import matplotlib.pyplot as plt

from glob      import glob
from functools import partial

from vetoedPDFs import sorter_func

from invisible_cities.io.pmaps_io          import load_pmaps_as_df
from invisible_cities.io.dst_io            import load_dsts
from invisible_cities.io.channel_param_io  import subset_param_reader as spr
import invisible_cities.core.fit_functions as fitf
from invisible_cities.icaro.hst_functions  import shift_to_bin_centers
from invisible_cities.database             import load_db as DB


def relative_pmt_response():
    """
    Script which uses pmaps (will be generalised in future to check XY dependence)
    to look at the relative response of the PMTs in Kr events and
    compares to the results on Poisson mu from calibrations
    """

    pmap_file_base = sys.argv[1]
    dst_file_base  = sys.argv[2]

    run_number = pmap_file_base.split('/')[2][1:]

    pmt_dats = DB.DataPMT(int(run_number))

    s1hists = {x : [] for x in range(12)}
    s2hists = {x : [] for x in range(12)}
    s1sumh  = []
    s2sumh  = []
    hitPMTdist = {x : [] for x in range(12)}
    hitPMTZpos = {x : [] for x in range(12)}

    pmap_sorter = sorter_func(pmap_file_base)
    pmap_file_list = sorted(glob(pmap_file_base + '*.h5'), key=pmap_sorter)

    dst_sorter = sorter_func(dst_file_base)
    dst_file_list = sorted(glob(dst_file_base + '*.h5'), key=dst_sorter)

    ## dst_frame = load_dsts(dst_file_list, 'DST', 'Events')
    dst_frame = load_dsts(dst_file_list, 'RECO', 'Events')

    dst_evt_list = dst_frame['event'].unique()

    #for fn in iglob(pmap_file_base + '*.h5'):
    for fn in pmap_file_list:

        ## This version just using pmt databases
        s1df, s2df, _, s1pmtdf, s2pmtdf = load_pmaps_as_df(fn)

        common_evts = np.intersect1d(s1pmtdf['event'].unique(), dst_evt_list)

        for evt in common_evts:
        #for evt in s1pmtdf['event'].unique():
            #evt    = dst_evt_iter[0]
            s1evt  = s1pmtdf[s1pmtdf['event'] == evt]
            s2evt  = s2pmtdf[s2pmtdf['event'] == evt]
            s1sevt = s1df[s1df['event'] == evt]
            s2sevt = s2df[s2df['event'] == evt]
            hit_evt = dst_frame[dst_frame['event'] == evt]
            ## if hit_evt['nS2'].iloc[0] == 1 and len(s2evt['peak'].unique()) == 1 and len(s1evt['peak'].unique()) == 1:
            if hit_evt['npeak'].nunique() == 1 and len(s2evt['peak'].unique()) == 1 and len(s1evt['peak'].unique()) == 1:
                ## Not well defined for multi-S2 events
                hit_x = hit_evt['X'].iloc[0]
                hit_y = hit_evt['Y'].iloc[0]
                hit_z = hit_evt['Z'].iloc[0]
                for peak in s1evt['peak'].unique():
                    s1peak = s1evt[s1evt['peak'] == peak]
                    s1sumh.append(s1sevt[s1sevt['peak'] == peak]['ene'].sum())
                    pmt1Q = s1peak[s1peak['npmt'] == 1]['ene'].sum()
                    for pmt in s1peak['npmt'].unique():
                        hitPMTdist[pmt].append(np.sqrt(np.power(hit_x-pmt_dats[pmt_dats['SensorID'] == pmt].X.values, 2)+np.power(hit_y-pmt_dats[pmt_dats['SensorID'] == pmt].Y.values, 2)))
                        hitPMTZpos[pmt].append(hit_z)
                        if pmt != 1:
                            s1hists[pmt].append(s1peak[s1peak['npmt'] == pmt]['ene'].sum()/pmt1Q)
                        else:
                            s1hists[pmt].append(pmt1Q)
                        
                for peak in s2evt['peak'].unique():
                    s2peak = s2evt[s2evt['peak'] == peak]
                    s2sumh.append(s2sevt[s2sevt['peak'] == peak]['ene'].sum())
                    if s2sumh[-1] > 4000:# and s2sumh[-1] < 12000:
                        ## pmt1Q = s2peak[s2peak['npmt'] == 1]['ene'].values[5:-5]
                        pmt1Q = s2peak[s2peak['npmt'] == 1]['ene'].sum()
                        for pmt in s2peak['npmt'].unique():
                            if pmt != 1:
                                ## s2hists[pmt].append(s2peak[s2peak['npmt'] == pmt]['ene'].values[5:-5]/pmt1Q)
                                s2hists[pmt].append(s2peak[s2peak['npmt'] == pmt]['ene'].sum()/pmt1Q)
                            else:
                                s2hists[pmt].append(pmt1Q)

            #dst_evt_iter.iternext()

    ## Make the plots
    s1sumh = np.array(s1sumh)
    s2sumh = np.array(s2sumh)
    figs0, axes0 = plt.subplots(nrows=1, ncols=2)
    axes0[0].hist(s1sumh)
    axes0[0].set_title('PMT sum S1 distribution')
    axes0[1].hist(s2sumh)
    axes0[1].set_title('PMT sum S2 distribution')
    plt.tight_layout()
    figs0.show()
    figs0.savefig('SumChargescharge_R'+run_number+'.png')
    figs1, axess1 = plt.subplots(nrows=3, ncols=4, figsize=(20,6))
    s1pmt1 = np.array(s1hists[1])
    s1bins = np.arange(-2, 4, 0.1)
    s2bins = np.arange(0.4, 1.1, 0.005)
    s1select = (s1sumh>2) & (s1sumh<150)
    for (key, val), ax in zip(s1hists.items(), axess1.flatten()):
        if key == 1:
            ax.hist(np.array(val)[s1select], bins=100)
            #ax.scatter(s1sumh[s1select], np.array(val)[s1select])
            ## ax.scatter(np.array(hitPMTdist[key])[s1select], np.array(val)[s1select])
            #ax.scatter(np.array(hitPMTZpos[key])[s1select], np.array(val)[s1select])
            ax.set_title('PMT 1 S1 charge')
            #ax.set_xlabel('integrated charge in PMT sum (pe)')
            #ax.set_xlabel('z pos.')
            #ax.set_ylabel('integrated charge in PMT1 (pe)')
            ax.set_ylabel('AU)')
            ax.set_xlabel('integrated charge in PMT1 (pe)')
            sh_hits = np.array(hitPMTZpos[key])[s1select].shape
            sh_val = np.array(val)[s1select].shape
            covar = np.cov(np.array(hitPMTZpos[key])[s1select].reshape(1, sh_hits[0]), np.array(val)[s1select].reshape(1, sh_val[0]))[0, 1]
            corr_coef = covar / (np.std(np.array(val)[s1select], ddof=1)*np.std(np.array(hitPMTZpos[key])[s1select], ddof=1))
            print('Sensor ', key, ' correlation coefficient = ', corr_coef)
        else:
            vals, bins, _ = ax.hist(np.array(val)[s1select], bins=s1bins)
            ## ax.scatter(s1pmt1[np.abs(val) < 10], np.array(val)[np.abs(val) < 10])
            #ax.scatter(s1sumh[s1select], np.array(val)[s1select])
            ## ax.scatter(np.array(hitPMTdist[key])[s1select & (np.abs(val) < 10)], np.array(val)[s1select & (np.abs(val) < 10)])
            #ax.scatter(np.array(hitPMTZpos[key])[s1select & (np.abs(val) < 10)], np.array(val)[s1select & (np.abs(val) < 10)])
            ax.set_title('PMT '+str(key)+' S1 relative charge')
            #ax.set_xlabel('integrated charge in PMT sum (pe)')
            ## ax.set_xlabel('PMT-hit dist. (mm)')
            #ax.set_xlabel('hit Z pos')
            #ax.set_ylabel('pmt q / pmt1 q')
            ax.set_ylabel('AU')
            ax.set_xlabel('pmt q / pmt1 q')
            ## sh_hits = np.array(hitPMTZpos[key])[s1select].shape
            ## sh_val = np.array(val)[s1select].shape
            ## covar = np.cov(np.array(hitPMTZpos[key])[s1select].reshape(1, sh_hits[0]), np.array(val)[s1select].reshape(1, sh_val[0]))[0, 1]
            ## corr_coef = covar / (np.std(np.array(val)[s1select], ddof=1)*np.std(np.array(hitPMTZpos[key])[s1select], ddof=1))
            s1select2 = s1select & (np.array(val) > 0) & (np.array(val) <= 2)
            print('Sensor ', key, ' mean = ', np.mean(np.array(val)[s1select2]))
            useful_bins = np.argwhere(vals>=100)
            b1 = useful_bins[0][0]
            b2 = useful_bins[-1][0]
            errs = np.sqrt(vals[b1:b2])
            fvals = fitf.fit(fitf.gauss, shift_to_bin_centers(bins)[b1:b2],
                             vals[b1:b2],
                             seed=(vals.sum(), bins[vals.argmax()], 0.1),
                             sigma=errs)
            ax.plot(shift_to_bin_centers(bins)[b1:b2],
                    fvals.fn(shift_to_bin_centers(bins)[b1:b2]))
            print('Fit S1 '+str(key), fvals.values, fvals.errors, fvals.chi2)
    plt.tight_layout()
    figs1.show()
    figs1.savefig('s1relativechargeThzoom_R'+run_number+'.png')

    fitVals = {}
    figs2, axess2 = plt.subplots(nrows=3, ncols=4, figsize=(20,6))
    s2pmt1 = np.array(s2hists[1])
    for (key, val), ax in zip(s2hists.items(), axess2.flatten()):
        if key == 1:
            ## ax.set_title('PMT 1 S2 charge')
            ax.set_title('PMT 1 S2 charge vs s2 sum charge')
            ax.set_xlabel('S2pmt sum charge (pe)')
            #ax.set_xlabel('integrated charge in PMT sum (pe)')
            #ax.set_xlabel('z pos')
            ax.set_ylabel('integrated charge in PMT1 (pe)')
            #ax.set_ylabel('AU')
            ax.set_xlabel('integrated charge in PMT1 (pe)')
            ## ax.hist(np.array(val)[(s2sumh>4000) & (s2sumh<12000)], bins=100)
            ##ax.hist(np.concatenate(val), bins=100)
            ## ax.hist(np.array(val)[s2sumh>4000], bins=100)
            ax.scatter(s2sumh[s2sumh>4000], np.array(val)[s2sumh>4000])
            #ax.scatter(s2sumh[(s2sumh>4000) & (s2sumh<12000)], np.array(val)[(s2sumh>4000) & (s2sumh<12000)])
            ## ax.scatter(np.array(hitPMTdist[key])[(s2sumh>4000) & (s2sumh<12000)], np.array(val)[(s2sumh>4000) & (s2sumh<12000)])
            #ax.scatter(np.array(hitPMTZpos[key])[(s2sumh>4000) & (s2sumh<12000)], np.array(val)[(s2sumh>4000) & (s2sumh<12000)])
            ## sh_hits = np.array(hitPMTdist[key])[(s2sumh>4000) & (s2sumh<12000)].shape
            ## sh_val = np.array(val)[(s2sumh>4000) & (s2sumh<12000)].shape
            ## covar = np.cov(np.array(hitPMTdist[key])[(s2sumh>4000) & (s2sumh<12000)].reshape(1, sh_hits[0]), np.array(val)[(s2sumh>4000) & (s2sumh<12000)].reshape(1, sh_val[0]))[0, 1]
            ## corr_coef = covar / (np.std(np.array(val)[(s2sumh>4000) & (s2sumh<12000)], ddof=1)*np.std(np.array(hitPMTdist[key])[(s2sumh>4000) & (s2sumh<12000)], ddof=1))
            ## print('Sensor ', key, ' correlation coefficient = ', corr_coef)
        else:
            ## ax.set_title('PMT '+str(key)+' S2 relative charge')
            ax.set_title('PMT '+str(key)+' S2 relative charge vs pmt sum')
            ax.set_ylabel('pmt q / pmt1 q')
            #ax.set_xlabel('zpos')
            ax.set_xlabel('integrated charge in PMT sum (pe)')
            #ax.set_xlabel('pmt q / pmt1 q')
            #ax.set_ylabel('AU')
            #ax.scatter(s2pmt1[np.abs(val) < 10], np.array(val)[np.abs(val) < 10])
            #ax.scatter(s2sumh[(s2sumh>4000) & (s2sumh<12000)], np.array(val)[(s2sumh>4000) & (s2sumh<12000)])
            ## vals, bins, _ = ax.hist(np.array(val)[(s2sumh>4000) & (s2sumh<12000)], bins=s2bins)
            ## vals, bins, _ = ax.hist(np.concatenate(val), bins=s2bins)
            ## vals, bins, _ = ax.hist(np.array(val)[s2sumh>4000], bins=s2bins)
            ax.scatter(s2sumh[s2sumh>4000], np.array(val)[s2sumh>4000])
            ## ax.scatter(np.array(hitPMTdist[key])[(s2sumh>4000) & (s2sumh<12000)], np.array(val)[(s2sumh>4000) & (s2sumh<12000)])
            #ax.scatter(np.array(hitPMTZpos[key])[(s2sumh>4000) & (s2sumh<12000)], np.array(val)[(s2sumh>4000) & (s2sumh<12000)])
            sh_sum = s2sumh[s2sumh>4000].shape
            sh_vals = np.array(val)[s2sumh>4000].shape
            covar = np.cov(s2sumh[s2sumh>4000].reshape(1, sh_sum[0]), np.array(val)[s2sumh>4000].reshape(1, sh_vals[0]))[0, 1]
            corr_coef = covar / (np.std(np.array(val)[s2sumh>4000], ddof=1)*np.std(s2sumh[s2sumh>4000], ddof=1))
            ## sh_hits = np.array(hitPMTdist[key])[(s2sumh>4000) & (s2sumh<12000)].shape
            ## sh_val = np.array(val)[(s2sumh>4000) & (s2sumh<12000)].shape
            ## covar = np.cov(np.array(hitPMTdist[key])[(s2sumh>4000) & (s2sumh<12000)].reshape(1, sh_hits[0]), np.array(val)[(s2sumh>4000) & (s2sumh<12000)].reshape(1, sh_val[0]))[0, 1]
            ## corr_coef = covar / (np.std(np.array(val)[(s2sumh>4000) & (s2sumh<12000)], ddof=1)*np.std(np.array(hitPMTdist[key])[(s2sumh>4000) & (s2sumh<12000)], ddof=1))
            print('Sensor ', key, ' correlation coefficient = ', corr_coef)
            ## limit fit to region with stat error <= 10% Poisson
            ## useful_bins = np.argwhere(vals>=200)
            ## b1 = useful_bins[0][0]
            ## b2 = useful_bins[-1][0]
            ## errs = np.sqrt(vals[b1:b2])
            ## print('Seed check: ', (vals.sum(), bins[vals.argmax()], 0.02))
            ## fvals = fitf.fit(fitf.gauss, shift_to_bin_centers(bins)[b1:b2], vals[b1:b2],
            ##                  seed=(vals.sum(), bins[vals.argmax()], 0.02),
            ##                  sigma=errs, bounds=[(0, 0, 0.00001), (1e10, 2, 3)])
            ## ax.plot(shift_to_bin_centers(bins),
            ##         fitf.gauss(shift_to_bin_centers(bins), *fvals.values))
            ## fitVals[key] = (fvals.values[1], fvals.values[2])
            ## print('Fit PMT '+str(key), fvals.values, fvals.errors, fvals.chi2)
    plt.tight_layout()
    figs2.show()
    figs2.savefig('s2relativechargeThvsSum_R'+run_number+'.png')

    ## figcal, axcal = plt.subplots()
    ## axcal.errorbar(list(fitVals.keys()),
    ##                np.fromiter((x[0] for x in fitVals.values()), np.float),
    ##                yerr=np.fromiter((x[1] for x in fitVals.values()), np.float),
    ##                label='Average response of PMTs to Kr relative to PMT 1')
    ## ## Get the calibration info for comparison.
    ## cal_files = [ fname for fname in sys.argv[3:] ]
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
    ## figcal.savefig('calPoisKrRelCompStatsFILT.png')
    input('plots good?')

                        
if __name__ == '__main__':
    relative_pmt_response()
