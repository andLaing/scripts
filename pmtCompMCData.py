import sys

import numpy  as np
import tables as tb
import pandas as pd

from scipy.optimize import leastsq

import matplotlib.pyplot as plt

from glob      import glob

from vetoedPDFs import sorter_func

from invisible_cities.io.pmaps_io import load_pmaps
import invisible_cities.core.fit_functions as fitf
from invisible_cities.icaro.hst_functions  import shift_to_bin_centers
from invisible_cities.io.dst_io import load_dsts
from invisible_cities.reco.dst_functions import load_lifetime_xy_corrections
from invisible_cities.reco               import pmaps_functions          as pmf

lt_corr = load_lifetime_xy_corrections('corrs/corrections_run6198.h5', group='XYcorrections', node='Elifetime')


def compare_mc():
    """
    Looks at MC and data for the PMTs and checks scale levels
    at the level of individual PMTs and the sum.
    Attempt to check values for relative scaling for PMTs.
    run as python pmtCompMCData.py <MC data file base> <Data data file base>
    """

    mc_file_base = sys.argv[1]
    da_file_base = sys.argv[2]
    dst_mc_base  = sys.argv[3]
    dst_da_base  = sys.argv[4]

    run_number = da_file_base[da_file_base.find('/r')+2:da_file_base.find('/r')+6]

    mc_sorter = sorter_func(mc_file_base)
    mc_file_list = sorted(glob(mc_file_base + '*.h5'), key=mc_sorter)
    da_sorter = sorter_func(da_file_base)
    da_file_list = sorted(glob(da_file_base + '*.h5'), key=da_sorter)

    mc_hit = load_dsts(glob(dst_mc_base + '*.h5'), 'RECO', 'Events')
    da_hit = load_dsts(glob(dst_da_base + '*.h5'), 'RECO', 'Events')
    
    dfcols = ['wf_sum', 'p0', 'cp0', 'p1', 'cp1', 'p2', 'cp2', 'p3', 'cp3', 'p4', 'cp4', 'p5', 'cp5', 'p6', 'cp6', 'p7', 'cp7', 'p8', 'cp8', 'p9', 'cp9', 'p10', 'cp10', 'p11', 'cp11']
    pmt_scales = [1, 0.79, 1, 0.80, 0.72, 1.11, 1.03, 0.82, 0.82, 1.03, 0.89, 0.95, 0.82]
    mc_sums = pd.DataFrame(columns=dfcols)
    for fn in mc_file_list:
        print('Reading mc file ', fn)
        pmaps = load_pmaps(fn)
        print('...data got')
        for evt, pmap in pmaps.items():

            if len(pmap.s2s) == 1:
                try:
                    mc_hit[mc_hit.event==evt].X.values[0]
                    hx = mc_hit[mc_hit.event==evt].X.values
                    hy = mc_hit[mc_hit.event==evt].Y.values
                    hz = mc_hit[mc_hit.event==evt].Z.values
                    #for s2 in pmap.s2s:
                    s2 = pmap.s2s[0]
                    rs2 = pmf.rebin_peak(s2, 2)
                    if hz.shape[0] == len(rs2.times):
                        new_row = [s2.pmts.waveform(x).sum() for x in range(12)]
                        cn_row  = [life_correction(hx, hy, hz, rs2.pmts.waveform(x)) for x in range(12)]
                        new_row = np.column_stack((new_row, cn_row)).flatten()
                        ## new_row.insert(0, s2.total_energy)
                        new_row = np.insert(new_row, 0, s2.total_energy)
                    
                        mc_sums.loc[len(mc_sums)] = list(new_row)
                except IndexError:
                    continue

    da_sums = pd.DataFrame(columns=dfcols)
    for fn in da_file_list:
        print('Reading data file ', fn)
        pmaps = load_pmaps(fn)
        print('...data got')
        for evt, pmap in pmaps.items():

            if len(pmap.s2s) == 1:
                try:
                    da_hit[da_hit.event==evt].X.values[0]
                    hx = da_hit[da_hit.event==evt].X.values
                    hy = da_hit[da_hit.event==evt].Y.values
                    hz = da_hit[da_hit.event==evt].Z.values
                    #for s2 in pmap.s2s:
                    s2 = pmap.s2s[0]
                    rs2 = pmf.rebin_peak(s2, 1)
                    #print('Check: ', hz.shape[0], len(rs2.times))
                    if hz.shape[0] == len(rs2.times):
                        new_row = [s2.pmts.waveform(x).sum() for x in range(12)]
                        cn_row  = [life_correction(hx, hy, hz, rs2.pmts.waveform(x)) for x in range(12)]
                        new_row = np.column_stack((new_row, cn_row)).flatten()
                        #new_row.insert(0, s2.total_energy)
                        new_row = np.insert(new_row, 0, s2.total_energy)
                    
                        da_sums.loc[len(da_sums)] = list(new_row)
                except IndexError:
                    continue

    trg0 = mc_sums['p0'] * pmt_scales[1] > 8835
    trg2 = mc_sums['p2'] * pmt_scales[3] > 7836
    ## Make some plots
    mc_sums[trg0 & trg2].wf_sum.plot.hist(bins=np.linspace(0, 1.2e6, 100),
                                          label='MC',
                                          density=True, histtype='step')
    da_sums.wf_sum.plot.hist(bins=np.linspace(0, 1.2e6, 100), label='data',
                             density=True, histtype='step')
    plt.title('PMT sum')
    plt.xlabel('Summed PMT charge (pe)')
    plt.yscale('log')
    plt.show()

    ## Attempt big fit. (only lifetime corrected [1::2] done in function)
    efunc = general_chi2(mc_sums.drop('wf_sum', axis=1).values.T)
    ## full_dats = np.apply_along_axis(np.histogram, 1,
    ##                                 da_sums.drop('wf_sum', axis=1).values.T,
    ##                                 bins=np.linspace(0, 120000, 100),
    ##                                 density=True)[:, 0]
    full_dats = np.apply_along_axis(np.histogram, 1,
                                    da_sums.drop('wf_sum', axis=1).values.T[1::2],
                                    bins=np.linspace(0, 120000, 100))[:, 0]
    dat_norms = np.fromiter((s.sum() for s in full_dats), np.int)
    full_dats = np.concatenate(full_dats)
    errs = np.sqrt(full_dats)
    errs[errs<=0] = 3
    par_seed = pmt_scales[1:]
    pfit, cov, infodict, msg, ier = leastsq(efunc, par_seed,
                                            args=(full_dats, errs, dat_norms),
                                            full_output=True)
    print('Fit res: ', pfit, ier, infodict, msg)
    trg0 = mc_sums['p0'] * pfit[1] * pfit[0] > 8835
    trg2 = mc_sums['p2'] * pfit[1] * pfit[2] > 7836
    fig, axes = plt.subplots(nrows=3, ncols=4, figsize=(20,6))
    #mc_sums['new_sum'] = mc_sums.drop('wf_sum', axis=1).sum(axis=1)
    for cname, ax, p in zip(dfcols[1::2], axes.flatten(), pfit):
        ax.set_title('PMT '+cname[2:]+' pe distribution')
        ax.set_xlabel('Photoelectrons')
        ax.set_ylabel('AU')
        mc_sums[trg0 & trg2][cname].multiply(p).plot.hist(ax=ax,
                                              bins=np.linspace(0, 120000, 100),
                                              label='MC', density=True,
                                              histtype='step', log=True)
        da_sums[cname].plot.hist(ax=ax, bins=np.linspace(0, 120000, 100),
                                 label='data', density=True, histtype='step', log=True)
        ## if 'p1' == cname:
        ##     mc_vals = mc_sums[trg0 & trg2][cname].values
        ##     da_vals = da_sums[cname].values
        ##     ffunc = simple_pmt1_fit(mc_vals)
        ##     dcv, hbins = np.histogram(da_vals, density=True,
        ##                               bins=np.linspace(0, 120000, 100))
        ##     hbins = shift_to_bin_centers(hbins)
        ##     errs = np.sqrt(dcv)
        ##     errs[errs==0] = 3
        ##     fvals = fitf.fit(ffunc, hbins, dcv, seed=(1), sigma=errs)
        ##     ax.plot(hbins, fvals.fn(hbins), label='fit attempt')
        ##     print('fit result: ', fvals.values, fvals.errors)
            
        ax.legend()
    plt.tight_layout()
    fig.show()
    plt.show()

    ## mc_sums.new_sum.plot.hist(bins=100, label='MC',
    ##                           density=True, histtype='step')
    ## da_sums.wf_sum.plot.hist(bins=100, label='data',
    ##                          density=True, histtype='step')
    ## plt.title('PMT sum')
    ## plt.xlabel('Summed PMT charge (pe)')
    ## plt.yscale('log')
    ## plt.show()
        


def life_correction(x, y, z, pmt_zspec):
    """
    take in the spectrum (pmap) for a PMT
    and apply the xy binned lifetime corrections.
    """

    corrected_val = np.sum(pmt_zspec * lt_corr(z, x, y).value)
    return corrected_val


def simple_pmt1_fit(mc_vals):
    
    def pmt1_fit(x, gscale):
        mc_spec, _ = np.histogram(gscale * mc_vals, density=True,
                                  bins=np.linspace(0, 120000, 100))
        return mc_spec
    return pmt1_fit


def general_chi2(mc):
    
    indcs = [0, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11]
    mc_raw = mc[::2]
    mc_corr = mc[1::2]
    #p_vals = [0.79, 0.80, 0.72, 1.11, 1.03, 0.82, 0.82, 1.03, 0.89, 0.95, 0.82]
    #p_errs = 0.01 ## Dummy for now
    def errorfunc(p, y, yerr, ynorms):

        trg0 = mc_raw[0] * p[1] * p[0] > 8835
        trg2 = mc_raw[2] * p[1] * p[2] > 7836

        ## spectra = [np.histogram(p[1] * p[i] * mc[i][trg0 & trg2],
        ##                         bins=np.linspace(0, 120000, 100),
        ##                         density = True)[0]
        ##            for i in indcs                                ]
        ## spectra.insert(1, np.histogram(p[1] * mc[1][trg0 & trg2]       ,
        ##                                bins=np.linspace(0, 120000, 100),
        ##                                density=True)[0]                 )
        spectra = [np.histogram(p[1] * p[i] * mc_corr[i][trg0 & trg2],
                                bins=np.linspace(0, 120000, 100))[0] for i in indcs]
        spectra.insert(1, np.histogram(p[1] * mc_corr[1][trg0 & trg2]       ,
                                       bins=np.linspace(0, 120000, 100))[0])
        mc_norms = np.fromiter((s.sum() for s in spectra), np.int)
        spectra = ynorms.reshape(ynorms.shape[0], 1) * spectra / mc_norms.reshape(mc_norms.shape[0], 1)
        spectra = np.concatenate(spectra)
        spec_err = np.sqrt(spectra)
        spec_err[spec_err<=0] = 3

        errfunc = (spectra - y) / np.sqrt(yerr**2 + spec_err**2)
        return errfunc
    return errorfunc

        

    

    
        

    

if __name__ == '__main__':
    compare_mc()
