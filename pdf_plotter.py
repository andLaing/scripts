import sys
import numpy as np
import tables as tb
import matplotlib.pyplot as plt

from invisible_cities.database import load_db as DB


def plot_pdfs():

    file_name = sys.argv[1]

    pdf_type = 'sipm_adc'
    if len(sys.argv) > 2:
        pdf_type = 'sipm_'+sys.argv[2]

    with tb.open_file(file_name, 'r') as sipmIn:
        
        bins = np.array(sipmIn.get_node('/HIST/', pdf_type+'_bins'))
        pdfs = np.array(sipmIn.get_node('/HIST/', pdf_type)).sum(axis=0)

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

        for ich, pdf in enumerate(pdfs):

            plt.cla()
            plt.ion()
            plt.show()

            plt.bar(bins, pdf, width=1, log=True)
            plt.title('Zero suppressed PDF for sensor '+str(chNos[ich])+', atca '+str(atcaNos[ich]))
            plt.xlabel('ADC')
            plt.ylabel('AU')
            plt.draw()
            plt.pause(0.001)
            status = input('next?')


if __name__ == '__main__':
    plot_pdfs()
