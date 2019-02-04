import numpy as np
import tables as tb
import matplotlib.pyplot as plt


import invisible_cities.core.fit_functions as fitf
from   invisible_cities.core.core_functions import weighted_mean_and_std
import invisible_cities.database.load_db   as DB


## def weighted_av_std(values, weights):

##     avg = np.average(values, weights=weights)

##     var = np.average((values-avg)**2, weights=weights)
##     # renormalize
##     var = weights.sum() * var / (weights.sum()-1)

##     return avg, np.sqrt(var)

def check_all_sipms_present(wf_file):


def sipm_connectivity_check(elec_name, dark_name, RUN_NO):
    """ Compare basic parameters of electronics only
    and dark current runs and have the user classify
    channels with little difference """

    sensors = DB.DataSiPM(int(RUN_NO)).SensorID.values

    with tb.open_file(elec_name, 'r') as eFile, \
         tb.open_file(dark_name, 'r') as dFile:

        min_incr = 0.5

        if 'HIST' not in eFile.root:
            print('Error, this check requires spectra')
            exit()

        binsE = np.array(eFile.root.HIST.sipm_dark_bins)
        binsD = np.array(dFile.root.HIST.sipm_dark_bins)

        specsE = np.array(eFile.root.HIST.sipm_dark).sum(axis=0)
        specsD = np.array(dFile.root.HIST.sipm_dark).sum(axis=0)

        ## Bin half widths as x errors
        xerrE = 0.5 * np.diff(binsE)[0]
        xerrD = 0.5 * np.diff(binsD)[0]

        with open('dodgy_channels.txt', 'a') as dChan:
            dChan.write('Run \t Channel \t classification \n')

            for ich, (espec, dspec) in enumerate(zip(specsE, specsD)):

                ## Maybe should be replaced with
                ## database access.
                chan_no = sensors[ich]

                avE, rmsE = weighted_mean_and_std(binsE, espec, unbiased=True)
                avD, rmsD = weighted_mean_and_std(binsD, dspec, unbiased=True)

                mean_low   = avE + min_incr >= avD
                rms_low    = rmsE + min_incr >= rmsD
                #stats_diff = espec.sum() != dspec.sum()

                if mean_low or rms_low:# or stats_diff:
                    print("Possible dodgy channel: "+str(chan_no))
                    print("identified as having: ")
                    if mean_low:   print('low mean: meanE='+str(avE)+', meanD='+str(avD))
                    if rms_low:    print('low rms: rmsE='+str(rmsE)+', rmsD='+str(rmsD))
                    #if stats_diff: print('missing entries: sumE/sumD='+str(espec.sum()/dspec.sum()))

                    plt.yscale('log')
                    plt.errorbar(binsE, espec, xerr=xerrE,
                                 yerr=np.sqrt(espec), fmt='b.',
                                 label='Electronics only')
                    plt.errorbar(binsD, dspec, xerr=xerrD,
                                 yerr=np.sqrt(dspec), fmt='r.', ecolor='r',
                                 label='Dark current')
                    plt.title('Elec/dark current comparison ch'+str(chan_no))
                    plt.xlabel('ADC')
                    plt.ylabel('AU')
                    plt.legend()
                    plt.show(block=False)

                    clif = input("How do you classify this channel? [dead/noisy/suspect/ok]")
                    dChan.write(RUN_NO+' \t '+str(chan_no)+' \t '+clif+' \n');
                    plt.clf()
                    plt.close()

            check_chan = input("Want to check specific channels? [y/n]")
            if check_chan == 'y':
                check_chan = input("Which channel number? [num/stop]")
                while check_chan != 'stop':
                    indx = np.argwhere(sensors==int(check_chan))
                    if len(indx) == 0:
                        print('Channel not found')
                        continue
                    indx = indx[0][0]

                    plt.yscale('log')
                    plt.errorbar(binsE, specsE[indx], xerr=xerrE,
                                 yerr=np.sqrt(specsE[indx]), fmt='b.',
                                 label='Electronics only')
                    plt.errorbar(binsD, specsD[indx], xerr=xerrD,
                                 yerr=np.sqrt(specsD[indx]), fmt='r.',
                                 ecolor='r',
                                 label='Dark current')
                    plt.title('Elec/dark current comparison ch'+str(check_chan))
                    plt.xlabel('ADC')
                    plt.ylabel('AU')
                    plt.legend()
                    plt.show(block=False)

                    clif = input("How do you classify this channel? [dead/noisy/suspect/ok]")
                    dChan.write(RUN_NO+' \t '+str(check_chan)+' \t '+clif+' \n');
                    plt.clf()
                    plt.close()
                    check_chan = input("Another channel? [num/stop]")
    return


def sipm_rms_check(wf_file):
    """
    Simple comparison of raw waveform rms
    to see if we have bad connections
    """

    ## Not great but since big variation in rms between DICES
    dice_rms_mins = [2.2, 2.0, 2.1, 2.0, 3.5, 3.1, 2.7, 2.5, 3.3, 2.8, 2.5, 2.6, 1.9, 2.4, 2.2, 2.1, 3.6, 3.3, 2.6, 2.5, 3.0, 2.7, 2.3, 2.5, 2.0, 2.3, 2.1, 2.1]

    with tb.open_file(wf_file, 'r') as data:

        bad_log = {}

        try:
            sipmrwf = data.root.RD.sipmrwf
        except tb.NoSuchNodeError:
            print('No raw wafeforms in file')
            exit()

        ch_nums = [[si['channel'], si['sensorID']] for si in data.root.Sensors.DataSiPM]

        for evt, tp in enumerate(sipmrwf):
            print('Checking event ', evt)
            for chNo, sipm in zip(ch_nums, tp):
                rms = np.std(sipm, ddof=1)

                if rms < dice_rms_mins[chNo[1] // 1000 -1]:
                    plt.plot(sipm)
                    plt.title('Raw waveform for atca ch '+str(chNo[0])+', sensor id '+str(chNo[1]))
                    plt.xlabel('Sample number')
                    plt.ylabel('ADC')
                    plt.show(block=False)

                    clif = input("How do you classify this channel? [OK/bad] ")
                    plt.clf()
                    plt.close()
                    if 'q' in clif:
                        exit()
                        
                    if 'bad' in clif:
                        if evt not in bad_log.keys():
                            bad_log[evt] = [chNo[1]]
                        else:
                            bad_log[evt].append(chNo[1])

        print('bad channel summary')
        for bkey in bad_log.keys():
            print('Event ', bkey, ' channels ', bad_log[bkey])
                    

def display_pmt_fits(spectra_file):

    with tb.open_file(spectra_file, 'r') as data:

        bins  = np.array(data.root.HIST.pmtspe_bins)
        specs = np.array(data.root.HIST.pmtspe).sum(axis=0)

        ## bin x error
        xerr = 0.5 * np.diff(bins)[0]

        ## Overall fitting function
        respF = fitf.SensorSpeResponse(bins)
        ffuncs = {'ngau':respF.set_gaussians, 'intgau':respF.min_integ_gaussians, 'dfunc': respF.scaled_dark_pedestal, 'conv':respF.dark_convolution}

        for i, spec in enumerate(specs):

            ## Get fit parameters
            vals = [ 22, 23, 24 ]
            plt.errorbar(bins, spec,
                         xerr=xerr, yerr=np.sqrt(spec),
                         fmt='r.', ecolor='r', label='Data')

            ## function name would come from the parameters too
            plt.plot(bins, ffuncs['ngau'](bins, *vals), label='Full function')
