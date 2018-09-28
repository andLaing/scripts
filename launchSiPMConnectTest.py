import sys

from calutils import          sipm_rms_check
from calutils import sipm_connectivity_check

## Launch as python launchSiPMConnectTest <run num> <elec only/wf file>
## opt: <dark current file>

run_num   = sys.argv[1]
elec_file = sys.argv[2]
if len(sys.argv) > 3:
    dark_file = sys.argv[3]

## First yo can do a basic rms check if the files are raw waveform files:
rms_check = input('Do you want to do a check on the rms? [y/n] ')
if 'y' in rms_check:
    sipm_rms_check(elec_file)

## Now you can do a full spectrum comparsison if you have the spectra
spec_check = input('Do you want to do a spectrum comparison? [y/n] ')
if 'y' in spec_check:
    sipm_connectivity_check(elec_file, dark_file, run_num)
