import numpy as np
import tables as tb
from scipy.optimize import curve_fit
from scipy.signal import find_peaks_cwt
from pandas import DataFrame
import matplotlib.pyplot as plt

#from JA_NEWCalibration.calib import responseWithDataPed
import invisible_cities.reco.spe_response as speR
import invisible_cities.core.fit_functions as fitf
from fitCheck import gauFit
from invisible_cities.database import load_db as DB


## Probably shite
darr = np.zeros(3)
def scaler(x, mu):
    global darr
    return mu * darr


def main():

    """ Check new fit function on SiPM spectra """

    saved_g = np.array(DB.DataSiPM(4832).adc_to_pes, np.float)
    chNos = np.array(DB.DataSiPM(4832).SensorID)

    ## sipmIn = tb.open_file('../outdats/sipmGainTest_R4832.h5', 'r')
    sipmIn = tb.open_file('../outdats/sipmSpeSpecs_gainmode_R5354.h5', 'r')

    ## Bins are the same for dark and light, just use light for now
    bins = np.array(sipmIn.root.HIST.sipm_spe_bins)
    ## LED correlated and anticorrelated spectra:
    specsL = np.array(sipmIn.root.HIST.sipm_spe).sum(axis=0)
    specsD = np.array(sipmIn.root.HIST.sipm_dark).sum(axis=0)

    ## General definition of fit class.
    #respF = SensorSpeResponse(bins, specsL[0].sum())

    ## The bounds on the parameters (norm, 1pe mean, poisson mean, 1pe sig)
    ## Are loose and general
    b0 = [(0., 0., 0., 0.001),(1e20, 100., 100., 10000.)]
    ## Names for output
    n0 = ['poisson', 'gain', 'sigma', 'chi2', 'saved gain', 'pull_gain', 'sigPe']

    ## Loop over the specra:
    outData = []
    dbDat = np.zeros((1792, 4), dtype=np.float)
    ## Extra protection since 3065 is weird
    knownDead = [ 3056, 8056, 14010, 25049 ]
    specialCheck = [1006, 1007, 3000, 3001, 7000, 22029, 28056, 28057]
    for ich, (led, dar, gn) in enumerate(zip(specsL, specsD, saved_g)):
        if chNos[ich] in knownDead:
            outData.append([0., 0., 0., 0., 0., 0., 0.])
            print('no peaks in dark spectrum, spec ', ich)
            continue
        ## Limits for safe fit
        b1 = np.argwhere(led>=50)[0][0]
        b2 = np.argwhere(led>=50)[-1][0]
        # Seed finding
        ## Dark:
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
        ## ped, pE = curve_fit(gauFit, bins[pD[0]-5:pD[0]+5], dar[pD[0]-5:pD[0]+5], p0=[dar[bins==0][0], 0., 2.])
        gfunc = fitf.fit(fitf.gauss, bins[pD[0]-5:pD[0]+5], dar[pD[0]-5:pD[0]+5], seed=(dar[bins==0][0], 0., 2.))
        ## Set the dark values to the fit function (still not ideal)
        #respF.set_dark_func(dar[b1:b2], ped[1], np.abs(ped[2]))
        respF = speR.scaled_dark_pedestal(dar[b1:b2], gfunc.values[1], np.abs(gfunc.values[2]), 100)
        global darr
        binR = bins[b1:b2]
        darr = dar[b1:b2]
        darr = darr[binR<0]
        ## LED
        pDL = find_peaks_cwt(led, np.arange(4, 20), min_snr=1, noise_perc=5)
        if len(pDL) <= 1 or abs(pDL[1]-pDL[0]) > 25:
            ## default values
            ledR = led[b1:b2]
            fscaler = fitf.fit(scaler, binR[binR<0], ledR[binR<0], (led[bins<0].max() / dar[bins<10].max()))
            muSeed = -np.log(fscaler.values[0])
            if muSeed < 0: muSeed = 0.001
            ## sd0 = (led.sum(), 15., muSeed, 2.)
            sd0 = (led.sum(), muSeed, 15., 2.)
        else:
            z0 = pDL[bins[pDL]<10]
            print('z0: ', z0)
            muSeed = 0.
            if z0.size == 0:
                print('check1: ', led[bins==0])
                ledR = led[b1:b2]
                fscaler = fitf.fit(scaler, binR[binR<0], ledR[binR<0], (led[bins<0].max() / dar[bins<10].max()))
                muSeed = -np.log(fscaler.values[0])
                if muSeed < 0: muSeed = 0.001
            else:
                ## ped, pE = curve_fit(gauFit, bins[z0[0]-5:z0[0]+5], led[z0[0]-5:z0[0]+5], p0=[led[z0[0]], ped[1], np.abs(ped[2])])
                vals = fitf.fit(fitf.gauss, bins[z0[0]-5:z0[0]+5], led[z0[0]-5:z0[0]+5], seed=(led[z0[0]], gfunc.values[1], np.abs(gfunc.values[2])))
                ped = vals.values
                print('Ped sees: ', ped)
            z1 = pDL[(bins[pDL]>10) & (bins[pDL]<30)]
            print('z1: ', z1)
            if z1.size == 0:
                print('check2: ', led[(bins>10) & (bins<30)].max())
                p1 = [led[(bins>10) & (bins<30)].max(), 15., 2.]
                if muSeed == 0.:
                    ledR = led[b1:b2]
                    fscaler = fitf.fit(scaler, binR[binR<0], ledR[binR<0], (led[bins<0].max() / dar[bins<10].max()))
                    muSeed = -np.log(fscaler.values[0])
                    if muSeed < 0: muSeed = 0.001
            else:
                z1 = z1[led[z1].argmax()]
                ## p1, p1E = curve_fit(gauFit, bins[z1-6:z1+6], led[z1-6:z1+6], p0=[led[z1], bins[z1], 3.])
                p1 = fitf.fit(fitf.gauss, bins[z1-6:z1+6], led[z1-6:z1+6], seed=(led[z1], bins[z1], 3.))
                p1 = p1.values
                if muSeed == 0.:
                    ledR = led[b1:b2]
                    fscaler = fitf.fit(scaler, binR[binR<0], ledR[binR<0], (led[bins<0].max() / dar[bins<10].max()))
                    muSeed = -np.log(fscaler.values[0])
                    if muSeed < 0: muSeed = 0.001
                print('1pe seed fit: ', p1)
            print(p1[0], ped[0], muSeed)
            # improve 1pe sigma seed?
            sigSeed = np.abs(p1[2])
            sd0 = (led.sum(), muSeed, p1[1], sigSeed)
        print('Pre fit seed check: ', sd0)
        ##
        ## Fit function
        ## vals, errs = curve_fit(respF.scaled_dark_pedestal,
        ##                            bins[b1:b2], led[b1:b2],
        ##                            sigma=np.sqrt(led[b1:b2]+dar[b1:b2]),
        ##                            p0=sd0, bounds=b0)
        rfun = fitf.fit(respF, bins[b1:b2], led[b1:b2], seed=sd0, sigma=np.sqrt(led[b1:b2]+dar[b1:b2]), bounds=b0)
        chDic = {}
        e1 = rfun.errors#np.sqrt(np.diagonal(rfun.errors))
        #chi = np.sum(((led[b1:b2]-respF.scaled_dark_pedestal(bins[b1:b2], *vals))/np.sqrt(led[b1:b2]+dar[b1:b2]))**2)
        #chi = chi/(len(bins[b1:b2])-len(vals))
        chi = rfun.chi2
        for i in range(1, len(rfun.values)):
            chDic[n0[i-1]] = (rfun.values[i], e1[i])
        #print(ich, chDic)
        print('indx = ', ich, 'channel = ', chNos[ich], ' seeds=', sd0)
        if chNos[ich] in specialCheck:
        ## if chi > 10:
            plt.bar(bins, led, width=1)
             ## plt.plot(bins[b1:b2], respF.scaled_dark_pedestal(bins[b1:b2], *vals), 'r')
            plt.plot(bins[b1:b2], respF(bins[b1:b2], *rfun.values), 'r')
            plt.show()
        print('nGau = ', respF.n_gaussians)
        dbDat[ich][0] = chDic['gain'][0]
        dbDat[ich][1] = chDic['gain'][1]
        dbDat[ich][2] = chDic['sigma'][0]
        dbDat[ich][3] = chDic['sigma'][1]
        outData.append([chDic['poisson'], chDic['gain'], chDic['sigma'], chi, gn, (gn-chDic['gain'][0])/chDic['gain'][1], chDic['sigma'][0]/chDic['gain'][0]])
    dfDat = DataFrame(outData, columns=n0, index=chNos)
    #print('Check: ', dfDat[(dfDat.chi2<10.) & (dfDat.chi2!=0)].pull_gain)
    dfDat.chi2.plot.hist(bins=250)
    plt.show()
    #dfDat[(dfDat.chi2!=0) & (abs(dfDat.pull_gain)<50)].pull_gain.hist(bins=500)
    ents, bins = np.histogram(np.array(dfDat[(dfDat.chi2!=0) & (abs(dfDat.pull_gain)<50)].pull_gain), bins=np.arange(-10, 10., 0.1))
    #print(len(ents), len(bins))
    print('Check entries: ', ents.sum())
    plt.bar(bins[:-1], ents, width=0.1)
    b1 = np.argwhere(ents>=10)[0][0]
    b2 = np.argwhere(ents>=10)[-1][0]
    pv, pvE = curve_fit(gauFit, bins[b1:b2], ents[b1:b2], sigma=np.sqrt(ents[b1:b2]), p0=[ents.max(), -1, 1])
    print('pull fit: ', pv, pvE)
    plt.plot(bins[:-1], gauFit(bins[:-1], *pv), 'r')
    print(dfDat[(dfDat.chi2!=0) & (abs(dfDat.pull_gain)>=50)].pull_gain)
    plt.show()
    dfDat[(dfDat.chi2<10.) & (dfDat.chi2!=0) & (dfDat.sigPe<10)].sigPe.plot.hist(bins=100)
    plt.show()
    #print(dfDat.sigma)
    dfDat.to_csv('SiPMDatTest.dat')

    with open('SiPMGain_database.dat', 'w') as dbF:
        ## for i in range(len(n0)):
        ##     dbF.write('5166, 100000, '+str(n0[i])+', '+str(dbDat[i][0])+', '+str(dbDat[i][1])+', '+str(dbDat[i][2])+', '+str(dbDat[i][3])+'\n')
        for ch, (gain, gE, sig, sigE) in zip(chNos, dbDat):
            dbF.write('5166, 100000, '+str(ch)+', '+str(gain)+', '+str(gE)+', '+str(sig)+', '+str(sigE)+'\n')
        

if __name__ == '__main__':
    main()
