#!/usr/bin/env python
import sys
import numpy as np
import tables as tb
import matplotlib.pyplot as plt
from glob import iglob
from invisible_cities.database import load_db as DB
from invisible_cities.core.random_sampling    import NoiseSampler as SiPMsNoiseSampler


def main():

    """ Plot some pdfs """

    in_base = sys.argv[1]
    runNo = int(sys.argv[2])

    sipmDats = DB.DataSiPM(runNo)
    noise_sampler = SiPMsNoiseSampler(runNo)

    binsS = 0
    binsM = 0
    sipmHist = 0
    mauHist = 0
    defsDone = False
    for fN in iglob(in_base+'*.h5'):
        with tb.open_file(fN, 'r') as data:
            if not defsDone:
                binsS = np.array(data.root.HIST.sipm_bins)
                binsM = np.array(data.root.HIST.sipmMAU_bins)
                sipmHist = np.zeros((len(data.root.HIST.sipm[0]),
                                     len(data.root.HIST.sipm[0][0])))
                mauHist = np.zeros((len(data.root.HIST.sipmMAU[0]),
                                    len(data.root.HIST.sipmMAU[0][0])))
                defsDone = True
                
            sipmHist += data.root.HIST.sipm[0]
            mauHist += data.root.HIST.sipmMAU[0]

    ## comparisons between plots and saved pdfs
    badCh = [ id for id, act in zip(sipmDats.SensorID,sipmDats.Active) if act==0]
    for i, (ped, mau, pdf) in enumerate(zip(sipmHist, mauHist, noise_sampler.probs)):
        sensor_id = sipmDats.SensorID[i]
        if not sensor_id in badCh:
            mean1, rms1 = wav(binsS, ped)
            mean2, rms2 = wav(binsM, mau)
            mean3, rms3 = wav2(noise_sampler.xbins, pdf)
            mnBins = binsS[np.nonzero(ped)[0][0]], binsM[np.nonzero(mau)[0][0]], noise_sampler.xbins[np.nonzero(pdf)[0][0]]
            mxBins = binsS[np.nonzero(ped)[0][-1]], binsM[np.nonzero(mau)[0][-1]], noise_sampler.xbins[np.nonzero(pdf)[0][-1]]
            ints = ped.sum(), mau.sum()
            maxPos = binsS[np.argmax(ped)], binsM[np.argmax(mau)], noise_sampler.xbins[np.argmax(pdf)]
            print('Sipm ', sensor_id)
            print('Means: ', mean1, mean2, mean3)
            print('RMS: ', rms1, rms2, rms3)
            print('Min bin: ', mnBins)
            print('Max bin: ', mxBins)
            print('Integral: ', ints)#Check if there was overflow
            print('Max val at: ', maxPos)
            ratS = np.where(pdf[:-1] != 0, ped/(ped.sum()*pdf[:-1]), 0.).sum()
            ratM = np.where(pdf[:-1] != 0, mau/(mau.sum()*pdf[:-1]), 0.).sum()
            print('Ratios: ', ratS, ratM)
            #if ratS > 100:
            if int(sensor_id)>=2000:## and int(sensor_id) < 11064 == 0:
            #if sensor_id == 4044 or np.abs(mean1-mean3) > 0.5 or np.abs(mean2-mean3) > 0.5 or np.abs(rms1-rms3) > 0.5 or np.abs(rms2-rms3) > 0.5 or np.abs(mnBins[0]-mnBins[2]) > 5 or np.abs(mnBins[1]-mnBins[2]) > 5 or np.abs(mxBins[0]-mxBins[2]) > 10 or np.abs(mxBins[1]-mxBins[2]) > 10:
                fig, axes = plt.subplots(nrows=1, ncols=3, figsize=(20,6))
                axes[0].bar(binsS, ped, width=0.1, log=True)
                axes[0].set_title('Spectrum rwf-mean, sipm '+str(sensor_id))
                axes[0].set_xlabel('pe')
                axes[0].set_ylabel('probability')
                axes[1].bar(binsM, mau, width=0.1, log=True)
                axes[1].set_title('Spectrum rwf-mean-mau, sipm '+str(sensor_id))
                axes[1].set_xlabel('pe')
                axes[1].set_ylabel('probability')
                axes[2].bar(noise_sampler.xbins, pdf, width=0.1, log=True)
                axes[2].set_title('PDF '+str(sensor_id))
                axes[2].set_xlabel('pe')
                axes[2].set_ylabel('probability')
                plt.tight_layout()
                fig.show()
                catch = input('test good?')
                if catch == 'q':
                    exit()
                elif catch == 's':
                    fig.savefig('SpecsCh'+str(sensor_id)+'.png')
                else:
                    plt.close(fig)
                    plt.clf()


def wav(bins, weis):

    mean = np.average(bins, weights=weis)
    var = np.average((bins-mean)**2, weights=weis)
    rms = np.sqrt(weis.sum()*var/(weis.sum()-1))

    return mean, rms

def wav2(bins, weis):
    """ Second version that gives the unbiased rms
    in the case where the weights are not frequencies (PDFs already normalised) """

    mean = np.average(bins, weights=weis)
    var = np.average((bins-mean)**2, weights=weis)
    ##
    sum1 = weis.sum()
    sum2 = (weis**2).sum()
    rms = np.sqrt(sum1*var/(sum1-sum2/sum1))

    return mean, rms


if __name__ == '__main__':
    main()
