import sys
import numpy as np
import tables as tb
from scipy.signal import find_peaks_cwt
from pandas import DataFrame
import matplotlib.pyplot as plt

import invisible_cities.core.fit_functions as fitf

"""
run as:
python pmtCalFit.py <input file name> <function name: ngau, intgau, dfunc, conv> opt: <minimum stats per bin to set fit limits> <number of sigmas for bound on pedestal parameters>

function name meanings:
ngau : 7 gaussians fitted, can be changed adding line respF.nGau = x
intgau : Gaussians until integral < min_integ=100, can be changed adding line respF.min_integ = x
dfunc : scaled dark spectrum + Gaussians up to integ < 100
conv : expicit convolution of dark spectrum with Gaussians up to integ < 100

"""

## this should probably be somehwere useful if it doesn't already exist
def weighted_av_std(values, weights):

    avg = np.average(values, weights=weights)

    var = np.average((values-avg)**2, weights=weights)
    # renormalize
    var = values.sum() * var / (values.sum()-1)

    return avg, np.sqrt(var)


## Probably shite
darr = np.zeros(3)
def scaler(x, mu):
    global darr
    return mu * darr


def seeds_and_bounds(func, bins, spec, ped_vals, ped_errs, lim_ped):

    norm_seed = spec.sum()
    
    ped_seed = ped_vals[1]
    ped_min  = ped_seed - lim_ped * ped_errs[1]
    ped_max  = ped_seed + lim_ped * ped_errs[1]
    #print('ped check: ', ped_seed, ped_min, ped_max)

    ped_sig_seed = ped_vals[2]
    ped_sig_min  = max(0.001, ped_sig_seed - lim_ped * ped_errs[2])
    ped_sig_max  = ped_sig_seed + lim_ped * ped_errs[2]
    #print('rms check: ', ped_sig_seed, ped_sig_min, ped_sig_max)

    ## Remove the ped prediction and check try to get seeds for 1pe
    # first scale the dark pedestal
    dscale = spec[bins<0].sum() / fitf.gauss(bins[bins<0], *ped_vals).sum()
    l_subtract_d = spec - fitf.gauss(bins, *ped_vals) * dscale
    pDL = find_peaks_cwt(l_subtract_d, np.arange(10, 20), min_snr=1, noise_perc=5)
    p1pe = pDL[(bins[pDL]>15) & (bins[pDL]<50)]
    p1pe = p1pe[spec[p1pe].argmax()]
    ## Now fit a Gaussian
    fgaus = fitf.fit(fitf.gauss, bins[p1pe-10:p1pe+10], l_subtract_d[p1pe-10:p1pe+10],
                     (l_subtract_d[p1pe-10:p1pe+10].max(), bins[p1pe], 7), sigma=np.sqrt(l_subtract_d[p1pe-10:p1pe+10]))
    #print('1pe fit check: ', fgaus.values, fgaus.errors)
    ## Test scale
    ftest = fitf.fit(scaler, bins[bins<0], spec[bins<0], (dscale))
    #print('ftest par = ', ftest.values[0], -np.log(ftest.values[0]))

    if 'gau' in func:
        # There are 6 variables: normalization, pedestal pos., spe mean, poisson mean, pedestal sigma, 1pe sigma
        sd0 = (norm_seed, ped_seed, fgaus.values[1]-ped_vals[1], -np.log(ftest.values[0]), ped_sig_seed, np.sqrt(fgaus.values[2]**2 - ped_vals[2]**2))
        bd0 = [(0, ped_min, 0, 0, ped_sig_min, 0.001), (1e99, ped_max, 10000, 10000, ped_sig_max, 10000)]
        #print('Seed check: ', sd0)
        return sd0, bd0
    ## The other functions only have four parameters: normalization, spe mean, poisson mean, 1pe sigma
    sd0 = (norm_seed, fgaus.values[1]-ped_vals[1], -np.log(ftest.values[0]), np.sqrt(fgaus.values[2]**2 - ped_vals[2]**2))
    bd0 = [(0, 0, 0, 0.001), (1e99, 10000, 10000, 10000)]
    #print('Seed check: ', sd0)
    return sd0, bd0
    

def main():
    """ Fitting for pmt response to ~spe led pulses """

    fileName = sys.argv[1]
    funcName = sys.argv[2]
    min_stat  = 0
    limit_ped = 10000.
    if len(sys.argv) > 3:
        min_stat = int(sys.argv[3])
        limit_ped = int(sys.argv[4])

    dats = tb.open_file(fileName, 'r')
    bins = np.array(dats.root.HIST.pmtdar_bins)
    specsD = np.array(dats.root.HIST.pmtdar).sum(axis=0)
    specsL = np.array(dats.root.HIST.pmtspe).sum(axis=0)

    respF = fitf.SensorSpeResponse(bins)
    ffuncs = {'ngau':respF.set_gaussians, 'intgau':respF.min_integ_gaussians, 'dfunc': respF.scaled_dark_pedestal, 'conv':respF.dark_convolution}

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

        ## Scale just in case we lost a different amount of integrals in dark and led
        scale = lspec.sum() / dspec.sum()
        #print('Scale check: ', scale)
        respF.set_dark_func(dspec[b1:b2], cent=gfitRes.values[1], sig=gfitRes.values[2], scale=scale)
        respF.redefine_bins(bins[b1:b2])

        ## Take into account the scale in seed finding (could affect Poisson mu)
        ped_vals = np.array([gfitRes.values[0] * scale, gfitRes.values[1], gfitRes.values[2]])

        binR = bins[b1:b2]
        global darr
        darr = dspec[b1:b2] * scale
        darr = darr[binR<0]
        seeds, bounds = seeds_and_bounds(funcName, bins[b1:b2], lspec[b1:b2],
                                         ped_vals, gfitRes.errors, limit_ped)

        ## The fit
        errs = np.sqrt(lspec[b1:b2])
        if not 'gau' in funcName:
            errs = np.sqrt(errs**2 + np.exp(-2 * seeds[2]) * dspec[b1:b2])
        errs[errs==0] = 1#0.001
        rfit = fitf.fit(ffuncs[funcName], bins[b1:b2], lspec[b1:b2], seeds, sigma=errs, bounds=bounds)
        ## plot the result
        plt.errorbar(bins, lspec, xerr=0.5*np.diff(bins)[0], yerr=np.sqrt(lspec), fmt='b.')
        plt.plot(bins[b1:b2], rfit.fn(bins[b1:b2]), 'r')
        plt.plot(bins[b1:b2], ffuncs[funcName](bins[b1:b2], *seeds), 'g')
        plt.title('Spe response fit to channel')
        plt.xlabel('ADC')
        plt.ylabel('AU')
        print('Fit values: ', rfit.values)
        print('Fit errors: ', rfit.errors)
        print('Number of Gaussians: ', respF.nGau)
        print('Fit chi2: ', rfit.chi2)
        plt.show(block=False)
        next_plot = input('press enter to move to next fit')
        plt.clf()
        plt.close()

    
if __name__ == '__main__':
    main()
