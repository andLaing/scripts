
import numpy as np
import tables as tb
import matplotlib.pyplot as plt
from functools import partial

from calutils import weighted_av_std

import invisible_cities.core.fit_functions as fitf
import invisible_cities.reco.spe_response  as speR


GainSeeds = [21.3, 23.4, 26.0, 25.7, 30.0, 22.7, 25.1, 32.7, 23.1, 25.5, 20.8, 22.0]
SigSeeds  = [11.3, 11.5, 10.6, 11.9, 13.1, 9.9, 11.0, 14.7, 10.6, 10.4, 9.3, 10.0]

ffuncs = {'ngau':speR.poisson_scaled_gaussians(n_gaussians=7),
          'intgau':speR.poisson_scaled_gaussians(min_integral=100),
          'dfunc':partial(speR.scaled_dark_pedestal, min_integral=100),
          'conv':partial(speR.dark_convolution, min_integral=100)}


darr = np.zeros(3)
def scaler(x, mu):
    global darr
    return mu * darr


def seeds_and_bounds(indx, func, bins, spec, ped_vals, ped_errs, lim_ped):

    norm_seed = spec.sum()
    
    ped_seed = ped_vals[1]
    ped_min  = ped_seed - lim_ped * ped_errs[1]
    ped_max  = ped_seed + lim_ped * ped_errs[1]

    ped_sig_seed = ped_vals[2]
    ped_sig_min  = max(0.001, ped_sig_seed - lim_ped * ped_errs[2])
    ped_sig_max  = ped_sig_seed + lim_ped * ped_errs[2]

    ## Remove the ped prediction and check try to get seeds for 1pe
    # first scale the dark pedestal
    dscale = spec[bins<0].sum() / fitf.gauss(bins[bins<0], *ped_vals).sum()
    GSeed  = GainSeeds[indx]
    GSSeed = SigSeeds[indx]
        
    ## Test scale
    ftest = fitf.fit(scaler, bins[bins<0], spec[bins<0], (dscale))

    if 'gau' in func:
        # There are 6 variables: normalization, pedestal pos., spe mean, poisson mean, pedestal sigma, 1pe sigma
        sd0 = (norm_seed, -np.log(ftest.values[0]), ped_seed, ped_sig_seed, GSeed, GSSeed)
        bd0 = [(0, 0, ped_min, ped_sig_min, 0, 0.001), (1e10, 10000, ped_max, ped_sig_max, 10000, 10000)]
        return sd0, bd0
    ## The other functions only have four parameters: normalization, spe mean, poisson mean, 1pe sigma
    sd0 = (norm_seed, -np.log(ftest.values[0]), GSeed, GSSeed)
    bd0 = [(0, 0, 0, 0.001), (1e10, 10000, 10000, 10000)]
    return sd0, bd0


def fit_dataset(dataF_table, funcName, min_stat, limit_ped):

    bins = np.array(dataF_table.root.HIST.pmt_dark_bins)
    specsD = np.array(dataF_table.root.HIST.pmt_dark).sum(axis=0)
    specsL = np.array(dataF_table.root.HIST.pmt_spe).sum(axis=0)

    ## pedSig err poissonMu err gain err gSig err chi2
    fitVals = np.zeros((12, 9), dtype=np.float)
    for i, (dspec, lspec) in enumerate(zip(specsD, specsL)):

        b1 = 0
        b2 = len(dspec)
        if min_stat != 0:
            valid_bins = np.argwhere(lspec>=min_stat)
            b1 = valid_bins[0][0]
            b2 = valid_bins[-1][0]

        ## Fit the dark spectrum with a Gaussian (not really necessary for the conv option)
        gb0 = [(0, -100, 0), (1e99, 100, 10000)]
        av, rms = weighted_av_std(bins, dspec)
        sd0 = (dspec.sum(), av, rms)
        errs = np.sqrt(dspec)
        errs[errs==0] = 0.0001
        gfitRes = fitf.fit(fitf.gauss, bins, dspec, sd0, sigma=errs, bounds=gb0)

        fitVals[i,0] = gfitRes.values[2]
        fitVals[i,1] = gfitRes.errors[2]

        scale = lspec.sum() / dspec.sum()
        
        if 'dfunc' in funcName:
            respF = ffuncs[funcName](dark_spectrum=dspec[b1:b2] * scale,
                                     pedestal_mean=gfitRes.values[1],
                                     pedestal_sigma=gfitRes.values[2])
        elif 'conv' in funcName:
            respF = ffuncs[funcName](dark_spectrum=dspec[b1:b2] * scale,
                                     bins=bins[b1:b2])
        else:
            respF = ffuncs[funcName]

        ped_vals = np.array([gfitRes.values[0] * scale, gfitRes.values[1], gfitRes.values[2]])

        binR = bins[b1:b2]
        global darr
        darr = dspec[b1:b2] * scale
        darr = darr[binR<0]
        seeds, bounds = seeds_and_bounds(i, funcName, bins[b1:b2], lspec[b1:b2],
                                         ped_vals, gfitRes.errors, limit_ped)

        ## The fit
        errs = np.sqrt(lspec[b1:b2])
        if not 'gau' in funcName:
            errs = np.sqrt(errs**2 + np.exp(-2 * seeds[1]) * dspec[b1:b2])
        errs[errs==0] = 1

        rfit = fitf.fit(respF, bins[b1:b2], lspec[b1:b2], seeds, sigma=errs, bounds=bounds)

        fitVals[i,2] = rfit.values[1]
        fitVals[i,3] = rfit.errors[1]
        fitVals[i,8] = rfit.chi2
        if 'gau' in funcName:
            fitVals[i,4] = rfit.values[4]
            fitVals[i,5] = rfit.errors[4]
            fitVals[i,6] = rfit.values[5]
            fitVals[i,7] = rfit.errors[4]
        else:
            fitVals[i,4] = rfit.values[2]
            fitVals[i,5] = rfit.errors[2]
            fitVals[i,6] = rfit.values[3]
            fitVals[i,7] = rfit.errors[3]

    return fitVals


def optPMTCal(fileNames, intWidths, funcName, min_stat, limit_ped):

    fResults = []
    for i in range(len(fileNames)):

        with tb.open_file(fileNames[i], 'r') as dataF:
            fResults.append(fit_dataset(dataF, funcName, min_stat, limit_ped))

    fig, axes = plt.subplots(nrows=3, ncols=2, figsize=(20,6))
    fig.show()
    axis_titles = ['Pedestal sigma', 'Poisson mu', 'Gain', 'Gain sigma', 'chi2']
    for j in range(12):
        ## clear the axes first
        for k, ax in enumerate(axes.flatten()):
            ax.cla()

            if k < 5:
                vals = np.fromiter((pars[j][k*2] for pars in fResults), np.float)
                if k < 4:
                    errs = np.fromiter((pars[j][k*2+1] for pars in fResults), np.float)
                ax.errorbar(intWidths, vals, yerr=errs, fmt='r.', ecolor='r')
                ax.set_title(axis_titles[k]+' vs integral width for PMT '+str(j))
        plt.tight_layout()
        plt.draw()
        catcher = input("next plot? q to stop ")
        if catcher == 'q':
            exit()
        plt.cla()
