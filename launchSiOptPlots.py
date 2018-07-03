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
widths = [f[f.find('start')+5:f.find('wid')-1]+'--'+f[f.find('end')+3:-3] for f in fileList]

opt.optSiPMCal(fileList, widths, 'dfunc', 0, 10000)
