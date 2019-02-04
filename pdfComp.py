#!/usr/bin/env python
import sys
import numpy as np
import tables as tb
import matplotlib.pyplot as plt
from invisible_cities.database             import load_db as DB
from invisible_cities.core.random_sampling import NoiseSampler
from invisible_cities.core.random_sampling import normalize_distribution
from invisible_cities.core.random_sampling import general_thresholds
from invisible_cities.core.core_functions  import weighted_mean_and_std



def main():

    """ Plot some pdfs """

    old_pdfs = sys.argv[1]
    new_pdfs = sys.argv[2]

    db1 = DB.DataSiPM(6499)
    db2 = DB.DataSiPM(6562)

    active = (db1.Active & db2.Active) == 1
    changed = (db1.SensorID >= 5000) & (db1.SensorID < 8064)

    sensor_x = db1.X.values
    sensor_y = db1.Y.values

    with tb.open_file(new_pdfs) as newfile, tb.open_file(old_pdfs) as oldfile:
        ## obins  = np.array(oldfile.root.HIST.sipm_mode_bins)
        ## ospecs = np.array(oldfile.root.HIST.sipm_mode).sum(0)
        ## nbins  = np.array(newfile.root.HIST.sipm_mode_bins)
        ## nspecs = np.array(newfile.root.HIST.sipm_mode).sum(0)
        obins  = np.array(oldfile.root.HIST.sipm_adc_bins)
        ospecs = np.array(oldfile.root.HIST.sipm_adc).sum(0)
        nbins  = np.array(newfile.root.HIST.sipm_adc_bins)
        nspecs = np.array(newfile.root.HIST.sipm_adc).sum(0)

        omerms = [weighted_mean_and_std(obins, spec, True) for spec in ospecs]
        nmerms = [weighted_mean_and_std(nbins, spec, True) for spec in nspecs]

        norm_old = np.apply_along_axis(normalize_distribution, 1, ospecs)
        thr_old  = general_thresholds(obins, norm_old, 0.99)

        norm_new = np.apply_along_axis(normalize_distribution, 1, nspecs)
        thr_new  = general_thresholds(nbins, norm_new, 0.99)

        fig, axes = plt.subplots(nrows=3, ncols=3, figsize=(20,6))
        ptitles = ['MEAN new', 'MEAN old', 'Mean new - Mean old',
                   'RMS new' , 'RMS old' , 'RMS new - RMS old'  ,
                   'Thr new' , 'Thr old' , 'Thr new - Thr old'  ]
        vals    = [np.fromiter((v[0] for v in nmerms), np.float),
                   np.fromiter((v[0] for v in omerms), np.float),
                   np.zeros(2),
                   np.fromiter((v[1] for v in nmerms), np.float),
                   np.fromiter((v[1] for v in omerms), np.float),
                   np.zeros(2),
                   thr_new, thr_old, np.zeros(2)]
        for i, (ax, val) in enumerate(zip(axes.flatten(), vals)):
            if not np.any(val):
                diff = vals[i-2] - vals[i-1]
                pinfo = ax.scatter(sensor_x[active & changed],
                                   sensor_y[active & changed],
                                   c=diff[active & changed], s=1.5)
            else:
                pinfo = ax.scatter(sensor_x[active & changed],
                                   sensor_y[active & changed],
                                   c=val[active & changed], s=1.5)
            ax.set_title(ptitles[i])
            ax.set_xlabel('X (mm)')
            ax.set_ylabel('Y (mm)')
            plt.colorbar(pinfo, ax=ax)
        plt.tight_layout()
        fig.show()
        plt.show()

    ## in_base = sys.argv[1]#, sys.argv[2] ]#, sys.argv[3], sys.argv[4] ]
    ## ## assumes order mean data run, mean DC run, mode data run, mode DC run
    ## lbls = [ 'r5906 median ped', 'r5906 mode ped' ]#, 'r4821 mean ped', 'r4821 mode ped' ]

    ## sipmDats = DB.DataSiPM(5906)
    ## noise_sampler = SiPMsNoiseSampler(5906)

    ## binsS = 0
    ## sipmHist = [ 0, 0 ]#, 0, 0 ]
    ## defsDone = [ False, False ]#, False, False ]
    ## #for i in range(4):
    ## #for i in range(2):
    ## for fN in iglob(in_base+'*.h5'):
    ##     with tb.open_file(fN, 'r') as data:
    ##         if not defsDone[i]:
    ##                 #if i == 0:
    ##             binsS = np.array(data.root.HIST.sipm_bins)
    ##             sipmHist[i] = np.zeros((len(data.root.HIST.sipm[0]),
    ##                                     len(data.root.HIST.sipm[0][0])))
    ##             defsDone[i] = True
                
    ##             sipmHist[i] += data.root.HIST.sipm[0]

    ## ## comparisons between plots and saved pdfs
    ## badCh = [ id for id, act in zip(sipmDats.SensorID,sipmDats.Active) if act==0]
    ## channel = input('Which channel do you want to see? ')
    ## if channel in badCh:
    ##     channel = input('Channel bad, pick another ')
    ## while not 'q' in channel:
    ##     inx = sipmDats.SensorID[sipmDats.SensorID==int(channel)].index[0]

    ##     for j, hst in enumerate(sipmHist):
    ##         plt.bar(binsS, hst[inx], width=0.1, log=True, label=lbls[j])
    ##     plt.title('PDF comparison channel '+channel)
    ##     plt.xlabel('pe')
    ##     plt.ylabel('AU')
    ##     plt.legend()
    ##     plt.show(block=False)
    ##     #plt.savefig('SipmSpecCompCh'+channel+'.png')
    ##     channel = input('Move to another channel or quit? ')
    ##     plt.clf()
    ##     plt.close()
    ## exit()


if __name__ == '__main__':
    main()
