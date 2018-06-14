import sys

import tables as tb
import numpy  as np

from glob import iglob as gl


def count_zeros():

    file_base = sys.argv[0]

    with open('zero_channels.txt', 'w') as out_dat:
        out_dat.write('Proportion of samples at zero per event \n \n \n')
        for ifile in gl(file_base+'*.h5'):
            with tb.open_file(ifile) as data_in:
                try:
                    sipmrwf = data_in.root.RD.sipmrwf
                except tb.NoSuchNodeError:
                    print('The files do not contain raw data')
                    exit()
                    
                sens_ids = np.fromiter((ch[1] for ch in data_in.root.Sensors.DataSiPM), np.int)
                evts     = np.fromiter((ev[0] for ev in data_in.root.Run.events), np.int)

                for evt_no, evtrwf in zip(evts, sipmrwf):
                    out = str(evt_no)
                    
                    indx = np.argwhere(sens_ids == 5010)[0][0]
                    prop_zero_5010 = 1 - np.count_nonzero(evtrwf[indx]) / 80
                    out += ' ch 5010 '+str(prop_zero_5010)
                    
                    indx = np.argwhere(sens_ids == 10010)[0][0]
                    prop_zero_10010 = 1 - np.count_nonzero(evtrwf[indx]) / 80
                    out += ' ch 10010 '+str(prop_zero_10010)
                    
                    indx = np.argwhere(sens_ids == 11010)[0][0]
                    prop_zero_11010 = 1 - np.count_nonzero(evtrwf[indx]) / 80
                    out += ' ch 11010 '+str(prop_zero_11010)+'\n'
                    
                    out_dat.write(out)


if __name__ == '__main__':
    count_zeros()
