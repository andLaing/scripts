import os
import sys
import subprocess

start_f = int(sys.argv[1])
end_f = int(sys.argv[2])

#env1 = os.environ.copy()

#print(env1)

for f1 in range(start_f, end_f):
    config_f = 'wfTest_proc'+str(f1)+'.conf'

    with open('wfTest.conf', 'r') as masterC, open(config_f, 'w') as procC:
        for line in masterC:
            if 'files_in' in line:
                ## line = line.replace('INFILE', 'test'+str(f1)+'pmtpls.h5')
                line = line.replace('INFILE', 'NoPulseTest0.h5')
            elif 'file_out' in line:
                line = line.replace('OUTFILE', 'testSim_nopls_modeHack'+str(f1)+'.h5')
            procC.write(line)
                
