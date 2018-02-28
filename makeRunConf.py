import sys


base_conf  = sys.argv[1]
new_conf   = sys.argv[2]
run_number = sys.argv[3]

with open(base_conf, 'r') as masterF, open(new_conf, 'w') as nF:
    for line in masterF:
        if 'RUNNO' in line:
            line = line.replace('RUNNO', run_number)

        nF.write(line)
