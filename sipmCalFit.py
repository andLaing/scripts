import sys
import numpy as np
import tables as tb
from scipy.signal import find_peaks_cwt
from pandas import DataFrame
import matplotlib.pyplot as plt
from functools import partial

#from JA_NEWCalibration.calib import responseWithDataPed
import invisible_cities.reco.spe_response as speR
import invisible_cities.core.fit_functions as fitf
from invisible_cities.database import load_db as DB

from pmtCalFit import weighted_av_std

useSavedSeeds = True
GainSeeds = []
SigSeeds  = []


## Probably shite
darr = np.zeros(3)
def scaler(x, mu):
    global darr
    return mu * darr


def seeds_and_bounds(indx, func, bins, spec, ped_vals, ped_errs, lim_ped):

    global GainSeeds, SigSeeds
    
    norm_seed = spec.sum()

    GSeed  = 0
    GSSeed = 0
    if useSavedSeeds:
        GSeed  = GainSeeds[indx]
        GSSeed = SigSeeds[indx]
    else:
        pDL = find_peaks_cwt(led, np.arange(4, 20), min_snr=1, noise_perc=5)
        p1pe = pDL[(bins[pDL]>10) & (bins[pDL]<20)]
        p1pe = p1pe[spec[p1pe].argmax()]
        p1 = fitf.fit(fitf.gauss, bins[p1pe-5:p1pe+5], spec[p1pe-5:p1pe+5], seed=(spec[p1pe], bins[p1pe], 3.))
        GSeed = p1.values[1]
        GSSeed = p1.values[2]

    dscale = spec[bins<5].sum() / fitf.gauss(bins[bins<5], *ped_vals).sum()
    fscale = fitf.fit(scaler, bins[bins<5], spec[bins<5], (dscale))
    muSeed = -np.log(fscale.values[0])
    if muSeed < 0: muSeed = 0.001

    if 'gau' in func:
        ped_seed = ped_vals[1]
        ped_min  = ped_seed - lim_ped * ped_errs[1]
        ped_max  = ped_seed + lim_ped * ped_errs[1]
        ped_sig_seed = ped_vals[2]
        ped_sig_min  = max(0.001, ped_sig_seed - lim_ped * ped_errs[2])
        ped_sig_max  = ped_sig_seed + lim_ped * ped_errs[2]
        sd0 = (norm_seed, muSeed, ped_seed, ped_sig_seed, GSeed, GSSeed)
        bd0 = [(0, 0, ped_min, ped_sig_min, 0, 0.001), (1e10, 10000, ped_max, ped_sig_max, 10000, 10000)]
        #print('Seed check: ', sd0)
        return sd0, bd0

    sd0 = (norm_seed, muSeed, GSeed, GSSeed)
    bd0 = [(0, 0, 0, 0.001), (1e10, 10000, 10000, 10000)]
    return sd0, bd0


def fit_dataset(dataF=None, funcName=None, minStat=None, limitPed=None):

    """ Check new fit function on SiPM spectra """
    global useSavedSeeds, GainSeeds, SigSeeds

    file_name = dataF
    func_name = funcName
    min_stat = minStat
    limit_ped = limitPed
    run_no = file_name[file_name.find('R')+1:file_name.find('R')+5]
    optimise = True
    if not file_name:
        optimise = False
        file_name = sys.argv[1]
        func_name = sys.argv[2]
        min_stat  = 0
        limit_ped = 10000.
        if len(sys.argv) > 3:
            useSavedSeeds = True if 'true' in sys.argv[3] else False
            min_stat = int(sys.argv[4])
            limit_ped = int(sys.argv[5])

    run_no = int(run_no)
    chNos = DB.DataSiPM(run_no).SensorID.values
    if useSavedSeeds:
        dodgy = DB.DataSiPM(run_no).index[DB.DataSiPM(run_no).Active==0].values
        GainSeeds = DB.DataSiPM(run_no).adc_to_pes.values
        SigSeeds  = DB.DataSiPM(run_no).Sigma.values
        ## Give generic values to previously dead or dodgy channels
        GainSeeds[dodgy] = 15
        SigSeeds[dodgy] = 2

    sipmIn = tb.open_file(file_name, 'r')

    ## Bins are the same for dark and light, just use light for now
    bins = np.array(sipmIn.root.HIST.sipm_spe_bins)
    ## LED correlated and anticorrelated spectra:
    specsL = np.array(sipmIn.root.HIST.sipm_spe).sum(axis=0)
    specsD = np.array(sipmIn.root.HIST.sipm_dark).sum(axis=0)
    
    ffuncs = {'ngau':speR.poisson_scaled_gaussians(n_gaussians=7),
              'intgau':speR.poisson_scaled_gaussians(min_integral=100),
              'dfunc':partial(speR.scaled_dark_pedestal, min_integral=100),
              'conv':partial(speR.dark_convolution, min_integral=100)}

    ## Loop over the specra:
    outData = []
    ## Extra protection since 3065 is weird
    knownDead = [ 3056, 8056, 14010, 25049 ]
    specialCheck = [1006, 1007, 3000, 3001, 7000, 22029, 28056, 28057]
    for ich, (led, dar) in enumerate(zip(specsL, specsD)):
        if chNos[ich] in knownDead:
            outData.append([chNos[ich], [0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0], 0, 0])
            print('no peaks in dark spectrum, spec ', ich)
            continue
        ## Limits for safe fit
        b1 = 0
        b2 = len(dar)
        if min_stat != 0:
            valid_bins = np.argwhere(led>=min_stat)
            b1 = valid_bins[0][0]
            b2 = valid_bins[-1][0]
        # Seed finding
        pD = find_peaks_cwt(dar, np.arange(2, 20), min_snr=2)
        if len(pD) == 0:
            ## Try to salvage in case not a masked channel
            ## Masked channels have al entries in one bin.
            if led[led>0].size == 1:
                outData.append([0., 0., 0., 0., 0., 0., 0.])
                print('no peaks in dark spectrum, spec ', ich)
                continue
            else:
                pD = np.array([dar.argmax()])
        ## Fit the dark spectrum with a Gaussian (not really necessary for the conv option)
        gb0 = [(0, -100, 0), (1e99, 100, 10000)]
        sd0 = (dar.sum(), 0, 2)
        errs = np.sqrt(dar[pD[0]-5:pD[0]+5])
        errs[errs==0] = 0.1
        gfitRes = fitf.fit(fitf.gauss, bins[pD[0]-5:pD[0]+5], dar[pD[0]-5:pD[0]+5], sd0, sigma=errs, bounds=gb0)

        ## Scale just in case we lost a different amount of integrals in dark and led
        scale = led.sum() / dar.sum()
        if 'dfunc' in func_name:
            respF = ffuncs[func_name](dark_spectrum=dar[b1:b2] * scale,
                                     pedestal_mean=gfitRes.values[1],
                                     pedestal_sigma=gfitRes.values[2])
        elif 'conv' in func_name:
            respF = ffuncs[func_name](dark_spectrum=dar[b1:b2] * scale,
                                     bins=bins[b1:b2])
        else:
            respF = ffuncs[func_name]

        ## Take into account the scale in seed finding (could affect Poisson mu)????
        ped_vals = np.array([gfitRes.values[0] * scale, gfitRes.values[1], gfitRes.values[2]])

        binR = bins[b1:b2]
        global darr
        darr = dar[b1:b2] * scale
        darr = darr[binR<5]
        seeds, bounds = seeds_and_bounds(ich, func_name, bins[b1:b2], led[b1:b2],
                                         ped_vals, gfitRes.errors, limit_ped)
        
        ## The fit
        errs = np.sqrt(led)
        if not 'gau' in func_name:
            errs = np.sqrt(errs**2 + np.exp(-2 * seeds[1]) * dar)
        errs[errs==0] = 0.001
        print('About to fit channel ', chNos[ich])
        rfit = fitf.fit(respF, bins[b1:b2], led[b1:b2], seeds, sigma=errs[b1:b2], bounds=bounds)
        chi = rfit.chi2
        if not optimise:
            if chNos[ich] in specialCheck or chi >= 10:
                if chNos[ich] in specialCheck: print('Special check channel '+str(chNos[ich]))
                print('Channel fit: ', rfit.values, chi)
                plt.errorbar(bins, led, xerr=0.5*np.diff(bins)[0], yerr=errs, fmt='b.')
                plt.plot(bins[b1:b2], respF(bins[b1:b2], *rfit.values), 'r')
                plt.title('Spe response fit to channel '+str(chNos[ich]))
                plt.xlabel('ADC')
                plt.ylabel('AU')
                plt.show()
        outData.append([chNos[ich], rfit.values, rfit.errors, respF.n_gaussians, chi])

    ## Couple of plots
    gainIndx = 2
    if 'gau' in func_name:
        gainIndx = 4
   
    pVals = [np.fromiter((ch[1][gainIndx] for ch in outData), np.float),
             np.fromiter((ch[1][gainIndx+1] for ch in outData), np.float),
             np.fromiter((ch[1][1] for ch in outData), np.float),
             np.fromiter((ch[4] for ch in outData), np.float)]
    if optimise:
        sipmIn.close()
        return pVals

    fig, axes = plt.subplots(nrows=2, ncols=2, figsize=(20,6))
    for ax, val in zip(axes.flatten(), pVals):
        ax.hist(val, bins=100)
    fig.show()
    input('finished with plots?')

    with open(file_name[:-3]+'_Fit_'+func_name+'.dat', 'w') as dbF:
        dbF.write('Minimum statistics: '+str(min_stat)+'\n\n')
        for vals in outData:
            dbF.write(str(vals)+'\n')
        

if __name__ == '__main__':
    fit_dataset()
