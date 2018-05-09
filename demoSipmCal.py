import sys
import numpy as np
import tables as tb
from scipy.signal import find_peaks_cwt
import matplotlib.pyplot as plt
from functools import partial
from pandas import DataFrame

import invisible_cities.reco.spe_response as speR
import invisible_cities.core.fit_functions as fitf
from invisible_cities.core.core_functions import weighted_mean_and_var
from invisible_cities.database import load_db as DB

## Probably shite
darr = np.zeros(3)
def scaler(x, mu):
    global darr
    return mu * darr


def seeds_and_bounds(bins, spec, ped_vals):

    norm_seed = spec.sum()

    GSeed  = 0
    GSSeed = 0
    pDL = find_peaks_cwt(spec, np.arange(4, 20), min_snr=1, noise_perc=5)
    p1pe = pDL[(bins[pDL]>10) & (bins[pDL]<20)]
    if len(p1pe) == 0:
        p1pe = np.argwhere(bins==15)[0][0]
    else:
        p1pe = p1pe[spec[p1pe].argmax()]
    p1 = fitf.fit(fitf.gauss, bins[p1pe-5:p1pe+5], spec[p1pe-5:p1pe+5], seed=(spec[p1pe], bins[p1pe], 3.))
    GSeed = p1.values[1] - ped_vals[1]
    if p1.values[2] <= ped_vals[2]:
        GSSeed = 0.5
    else:
        GSSeed = np.sqrt(p1.values[2]**2 - ped_vals[2]**2)

    dscale = spec[bins<5].sum() / fitf.gauss(bins[bins<5], *ped_vals).sum()
    errs = np.sqrt(spec[bins<5])
    errs[errs==0] = 1
    fscale = fitf.fit(scaler, bins[bins<5], spec[bins<5], (dscale), sigma=errs)
    muSeed = -np.log(fscale.values[0])
    if muSeed < 0: muSeed = 0.001

    sd0 = (norm_seed, muSeed, GSeed, GSSeed)
    bd0 = [(0, 0, 0, 0.001), (1e10, 10000, 10000, 10000)]
    return sd0, bd0


def fit_dataset():

    file_name = sys.argv[1]
    min_stat  = 0
    if len(sys.argv) > 2:
        min_stat = int(sys.argv[2])

    run_no = file_name[file_name.find('R')+1:file_name.find('R')+5]

    sipmIn = tb.open_file(file_name, 'r')

    ## Bins are the same for dark and light, just use light for now
    bins = np.array(sipmIn.root.HIST.sipm_spe_bins)
    ## LED correlated and anticorrelated spectra:
    specsL = np.array(sipmIn.root.HIST.sipm_spe).sum(axis=0)
    specsD = np.array(sipmIn.root.HIST.sipm_dark).sum(axis=0)

    ffunc = partial(speR.scaled_dark_pedestal, min_integral=100)

    ## This needs to be from a DB or some tempory mapping!!
    if 'Sensors' in sipmIn.root:
        atcaNos = np.fromiter((x['channel'] for x in sipmIn.root.Sensors.DataSiPM), np.int)
        chNos   = np.fromiter((x['sensorID'] for x in sipmIn.root.Sensors.DataSiPM), np.int)
        if np.all(chNos == -1):
            dats = DB.DataSiPM(5166)
            chns = dats.ChannelID
            chNos = np.fromiter((dats.loc[chns == x, 'SensorID'] for x in atcaNos), np.int)
    else:
        print('No sensor info in file, assuming DB ordering for run 5166')
        atcaNos = DB.DataSiPM(5166).ChannelID.values
        chNos = DB.DataSiPM(5166).SensorID.values
    n0 = ['Normalization', 'norm error', 'poisson mu', 'poisson error', 'Pedestal', 'pedestal error', 'Pedestal sigma', 'ped sig error', 'gain', 'gain error', 'gain sigma', 'gain sig error', 'chi2']
    outData = []

    for ich, (led, dar) in enumerate(zip(specsL, specsD)):
        print('channel index = ', ich, ', sensor: ', chNos[ich], ', atca: ', atcaNos[ich])
        channelRes = [0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.]

        plt.cla()
        plt.ion()
        plt.show()

        #protections
        avL, _ = weighted_mean_and_var(bins, led, True)
        avD, _ = weighted_mean_and_var(bins, dar, True)
        print(avL - avD)
        if avL - avD < 1.0:
            plt.bar(bins, led, width=1)
            plt.bar(bins, dar, width=1)
            plt.title('Spe dark and led spectra for sensor '+str(chNos[ich])+', atca '+str(atcaNos[ich]))
            plt.xlabel('ADC')
            plt.ylabel('AU')
            plt.draw()
            plt.pause(0.001)
            status = input('Potential dead or zero light channel, do fit? [yes/no] ')
            if 'no' in status:
                outData.append(channelRes)
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
                outData.append(channelRes)
                print('no peaks in dark spectrum, spec ', ich)
                continue
            else:
                pD = np.array([dar.argmax()])
        ## Fit the dark spectrum with a Gaussian
        gb0 = [(0, -100, 0), (1e99, 100, 10000)]
        sd0 = (dar.sum(), 0, 2)
        errs = np.sqrt(dar[pD[0]-5:pD[0]+5])
        errs[errs==0] = 0.1
        gfitRes = fitf.fit(fitf.gauss, bins[pD[0]-5:pD[0]+5], dar[pD[0]-5:pD[0]+5], sd0, sigma=errs, bounds=gb0)

        channelRes[5] = gfitRes.values[1]
        channelRes[6] = gfitRes.errors[1]
        channelRes[7] = gfitRes.values[2]
        channelRes[8] = gfitRes.errors[2]

        respF = ffunc(dark_spectrum=dar[b1:b2],
                      pedestal_mean=gfitRes.values[1],
                      pedestal_sigma=gfitRes.values[2])

        ped_vals = np.array([gfitRes.values[0], gfitRes.values[1], gfitRes.values[2]])

        binR = bins[b1:b2]
        global darr
        darr = dar[b1:b2]
        darr = darr[binR<5]
        seeds, bounds = seeds_and_bounds(bins[b1:b2], led[b1:b2], ped_vals)

        ## The fit
        errs = np.sqrt(led + np.exp(-2 * seeds[1]) * dar)
        errs[errs==0] = 0.001
        rfit = fitf.fit(respF, bins[b1:b2], led[b1:b2], seeds, sigma=errs[b1:b2], bounds=bounds)
        chi = rfit.chi2
        if chi >= 7 or rfit.values[3] >= 2.5 or rfit.values[3] <= 1:
            ## The offending parameter seems to be the sigma in most cases
            nseed = rfit.values
            nseed[3] = 1.7
            nbound = [(bounds[0][0], bounds[0][1], bounds[0][2], 1),
                      (bounds[1][0], bounds[1][1], bounds[1][2], 2.5)]
            rfit = fitf.fit(respF, bins[b1:b2], led[b1:b2],
                            nseed, sigma=errs[b1:b2], bounds=nbound)
            chi = rfit.chi2

        if chi >= 10 or rfit.values[2] < 12 or rfit.values[3] > 3:
            plt.errorbar(bins, led, xerr=0.5*np.diff(bins)[0], yerr=errs, fmt='b.')
            plt.plot(bins[b1:b2], respF(bins[b1:b2], *rfit.values), 'r')
            plt.plot(bins[b1:b2], respF(bins[b1:b2], *seeds), 'g')
            plt.title('Spe response fit to sensor '+str(chNos[ich])+', atca '+str(atcaNos[ich]))
            plt.xlabel('ADC')
            plt.ylabel('AU')
            plt.draw()
            plt.pause(0.001)
            input('bad fit, just so you know')

        channelRes[0] = rfit.values[0]
        channelRes[1] = rfit.errors[0]
        channelRes[2] = rfit.values[1]
        channelRes[3] = rfit.errors[1]
        channelRes[8] = rfit.values[2]
        channelRes[9] = rfit.errors[2]
        channelRes[10] = rfit.values[3]
        channelRes[11] = rfit.errors[3]
        channelRes[12] = chi
        outData.append(channelRes)

    dfDat = DataFrame(outData, columns=n0, index=chNos)

    dfDat.chi2.plot.hist(bins=250)
    input('Chi^2 dist ok?')

    plt.cla()
    dfDat.gain.plot.hist(bins=250)
    input('Gain dist ok?')
    
    dfDat.to_csv('SiPMFits_file-'+file_name.split('/')[-1]+'.txt')

if __name__ == '__main__':
    fit_dataset()
        
