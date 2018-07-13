## Variables are city=city executable, conf=base configuration, out_file
## minStart=first start integral point, maxStart=upper limit on start
## minWidth=first integral width, maxWidth=maximum width
import sys
import numpy as np
from subprocess import call


city           = sys.argv[1]
base_conf      = sys.argv[2]
run_num        = sys.argv[3]
out_file_base  = sys.argv[4]
minStart       = sys.argv[5]
maxStart       = sys.argv[6]
minWidth       = sys.argv[7]
maxWidth       = sys.argv[8]

sample_wid = 0.025
cal_folder = 'pmtCal_2018'
if 'sipm' in city:
    sample_wid = 1
    cal_folder = 'sipmCal_2018'

starts = np.arange(float(minStart), float(maxStart), sample_wid)
wids = np.arange(float(minWidth), float(maxWidth), sample_wid)

for st in starts:
    for wid in wids:
        new_conf_base = base_conf[:-5]
        if '../' in new_conf_base:
            new_conf_base = new_conf_base[3:]

        new_conf = '/calibration/'+cal_folder+'/config/'+new_conf_base
        new_conf += 'R'+run_num
        new_conf += '_intStart'+str(st)+'_intWid'+str(wid)+'.conf'
        out_file = 'intOpt/'+out_file_base+'R'+run_num+'_start'+str(st)+'_wid'+str(wid)+'_end'+str(st+wid)+'.h5'
        log = 'opt'+run_num+'_st'+str(st)+'_wid'+str(wid)

        mkConf = 'python makeRunConf.py '+base_conf+' '+new_conf+' '+run_num+' '+out_file+' '+str(st)+' '+str(wid)

        call(mkConf, shell=True, executable='/bin/bash')

        q_sub = 'qsub -N opt'+run_num+' -v city='+city+',conf='+new_conf+',log='+log+' jobTemplateTemp.sh'
        call(q_sub, shell=True, executable='/bin/bash')
