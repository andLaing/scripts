import sys
import glob as gl

import optPMTcal as opt

start_or_end_fixed = sys.argv[1]
point = sys.argv[2]
run_num = sys.argv[3]

#file_base = '/calibration/pmtCal_2018/spectra/intOpt/pmtCalOptR'+run_num+'_start'
file_base = '../outdats/pmtCalOptR'+run_num+'_start'
if 'start' in start_or_end_fixed:
    file_base += point+'_wid*.h5'
elif 'end' in start_or_end_fixed:
    file_base += '*_wid*_end'+point+'.h5'
else:
    print('invalid fix point')
    exit()

fileList = gl.glob(file_base)
widths = [float(f[f.find('wid')+3:f.find('end')-1]) for f in fileList]

opt.optPMTCal(fileList, widths, 'intgau', 100, 100)
