import sys
import glob as gl

import optPMTcal as opt

start_or_end_fixed = sys.argv[1]
point = sys.argv[2]
run_num = sys.argv[3]

#file_base = '/calibration/sipmCal_2018/spectra/intOpt/sipmCalOptR'+run_num+'_start'
file_base = '../outdats/sipmCalOptR'+run_num+'_start'
if 'start' in start_or_end_fixed:
    file_base += point+'_wid*.h5'
elif 'end' in start_or_end_fixed:
    file_base += '*_wid*_end'+point+'.h5'
else:
    print('invalid fix point')
    exit()

fileList = gl.glob(file_base)
widths = [float(f[f.find('wid')+3:f.find('end')-1]) for f in fileList]

opt.optSiPMCal(fileList, widths, 'dfunc', 0, 10000)
