import sys
import numpy as np
import tables as tb
import matplotlib.pyplot as plt

from scipy.signal import find_peaks
from scipy.signal import find_peaks_cwt

from glob import glob

from invisible_cities.database            import               load_db as DB
from invisible_cities.core.core_functions import weighted_mean_and_std

def file_sorter(elem):
    return int(elem[-7:-3])

def plot_pdfs():

    file_name = sys.argv[1]
    if '.h5' in file_name:
        print('Single file mode')
        all_files = [file_name]
        runs = [file_name[-7:-3]]
    else:
        all_files = sorted(glob(file_name + '*.h5'), key=file_sorter)
        runs = [fn[-7:-3] for fn in all_files]

    pdf_type = 'adc'
    if len(sys.argv) > 2:
        pdf_type = sys.argv[2]

    sipmIn = [tb.open_file(fn) for fn in all_files]

    #with tb.open_file(file_name, 'r') as sipmIn:
        
    bins = np.array(sipmIn[0].get_node('/HIST/', pdf_type+'_bins'))
    pdfs = [np.array(siIn.get_node('/HIST/', pdf_type)).sum(axis=0) for siIn in sipmIn]

    if 'Sensors' in sipmIn[0].root:
        atcaNos = np.fromiter((x['channel'] for x in sipmIn[0].root.Sensors.DataSiPM), np.int)
        chNos   = np.fromiter((x['sensorID'] for x in sipmIn[0].root.Sensors.DataSiPM), np.int)
        if np.all(chNos == -1):
            dats = DB.DataSiPM(5166)
            chns = dats.ChannelID
            chNos = np.fromiter((dats.loc[chns == x, 'SensorID'] for x in atcaNos), np.int)
    else:
            ## print('No sensor info in file, assuming DB ordering for run 5166')
            ## atcaNos = DB.DataSiPM(5166).ChannelID.values
            ## chNos = DB.DataSiPM(5166).SensorID.values
        print('No sensor info in file, hardwired')
        dices = [14, 15, 17, 21]
        atcaNos = np.fromiter((d * 1000 + i for d in dices for i in range(64)), np.int)
        chNos = atcaNos
        
    means   = []
    rmss    = []
    peak1   = []
    grads   = []
    cable   = []
    n_cable = [1, 2, 3, 4, 5, 6, 7]
            
    for ich, pdf in enumerate(zip(*pdfs)):
        if chNos[ich] < 21000: continue
        #plt.cla()
        #plt.ion()
        #plt.show()

        for j, pf in enumerate(pdf):
            mean, rms = weighted_mean_and_std(bins, pf, True)
            means.append(mean)
            rmss .append(rms)
            cable.append(n_cable[j])
            peaks = find_peaks(pf, 100, distance=15)
            try:
                peak1.append(bins[peaks[0][1]])
            except IndexError:
                ## Fake to save hassle
                peak1.append(-1)
            #plt.bar(bins, pf, width=1, log=True, label='Run '+runs[j])
        grads.append((peak1[-1] - peak1[-len(pdf)]) / (cable[-1] - cable[-len(pdf)]))
        #plt.title('Zero suppressed PDF for sensor '+str(chNos[ich])+', atca '+str(atcaNos[ich]))
        #plt.xlabel('ADC')
        #plt.ylabel('AU')
        #plt.legend()
        #plt.draw()
        #plt.pause(0.001)
        #plt.show()
        #status = input('next?')
        #if 'q' in status:
        #    exit()

    plt.scatter(cable, means)
    plt.title('Scatter of distribution mean with number of cables')
    plt.xlabel('Number of cables')
    plt.ylabel('PDF mean')
    plt.show()

    plt.scatter(cable, rmss)
    plt.title('Scatter of distribution rms with number of cables')
    plt.xlabel('Number of cables')
    plt.ylabel('PDF rms')
    plt.show()

    plt.scatter(cable, peak1)
    plt.xlabel('Number of cables')
    plt.ylabel('Approximate 1 pe peak position')
    plt.show()

    plt.hist(grads)
    plt.title('histogram of 1 pe peak pos / no. cables gradient')
    plt.xlabel('gradient')
    plt.ylabel('AU')
    plt.show()


if __name__ == '__main__':
    plot_pdfs()
