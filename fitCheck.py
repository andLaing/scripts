#!/usr/bin/env python


import math
import sys
import numpy as np
import tables as tb
from scipy.signal import find_peaks_cwt
from scipy.optimize import curve_fit


def gauFit(val, norm, mu, sig):
    efacto = 1./math.sqrt(2.*math.pi*sig**2)
    
    return norm * efacto*np.exp(-(val-mu)**2/(2*sig**2))

def main():
    """ Check some fits """

    file_in = sys.argv[1]

    mcg = tb.open_file(file_in, 'r')

    bins = np.array(mcg.root.HIST.pmtspe_bins)
    specs = np.array(mcg.root.HIST.pmtspe).sum(axis=0)

    allMean = 0.
    allErr = 0.
    for indx, sipm in enumerate(specs):
        fitVals = []
        ErrVals = []
        peaks = find_peaks_cwt(sipm, np.arange(4,10), min_snr=3)
        for indx2, p in enumerate(peaks):
            if indx2 < 5:
                vals, errs = curve_fit(gauFit, bins[p-5:p+5],
                                        sipm[p-5:p+5], p0=[sipm[p], bins[p], 5],
                                        sigma=np.sqrt(sipm[p-5:p+5]),
                                        bounds=([0., bins[p-5], 0.],[1000000, bins[p+5], 100000.]))
                fitVals.append(vals)
                ErrVals.append(np.diagonal(errs))
        stds = np.fromiter((v[2] for v in fitVals), np.float)
        #print(indx, fitVals)
        print(indx, 'approx Poisson mu = ',fitVals[1][0]/fitVals[0][0],'mean rms = ', stds.mean(), ' +/- ', stds.std(ddof=1), 'approx gain = ',fitVals[1][1]-fitVals[0][1])
        meanG = 0.
        meanE = 0.
        for indx3, (fvals, evals, std) in enumerate(zip(fitVals, ErrVals, stds)):
            if indx3 != 0:
                g1 = (fvals[1]-fitVals[0][1]) / (indx3)
                e1 = np.sqrt(evals[1]+ErrVals[0][1]) / indx3
                meanG += g1
                meanE += e1*e1
                print('Gain from peak ', indx3, ' = ', g1, '+/-', e1)
                print('Peak sigma: ', std)
        allMean += meanG/(len(fitVals)-1)
        allErr += (np.sqrt(meanE)/(len(fitVals)-1))**2
        print('Mean gain = ', meanG/(len(fitVals)-1), ' +/-',np.sqrt(meanE)/(len(fitVals)-1))
    print('Mean of all: ', allMean/len(specs), '+/-', np.sqrt(allErr)/len(specs))
    

if __name__ == '__main__':
    main()
