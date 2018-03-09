import sys


base_conf  = sys.argv[1]
new_conf   = sys.argv[2]
run_number = sys.argv[3]

# OUTFILE: pmtSpeSpecs_gainmode_RRUNNO.h5, sipmSpeSpecs_gainmode_RRUNNO.h5 etc
out_file = sys.argv[4]

# Standard integral values
int_start = '5.8'
int_width = '0.5'
if 'sipm' in base_conf:
    int_start = '41'
    int_width = '2'

if len(sys.argv) > 5:
    int_start = sys.argv[5]
    int_width = sys.argv[6]

with open(base_conf, 'r') as masterF, open(new_conf, 'w') as nF:
    for line in masterF:
        if 'RUNNO' in line:
            line = line.replace('RUNNO', run_number)
        if 'OUTFILE' in line:
            line = line.replace('OUTFILE', out_file)
        if 'INTSTART' in line:
            line = line.replace('INTSTART', int_start)
        if 'INTWID' in line:
            line = line.replace('INTWID', int_width)

        nF.write(line)
