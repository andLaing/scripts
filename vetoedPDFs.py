import sys
import tables as tb
import numpy as np
from glob import glob
from functools import partial

from invisible_cities.core.system_of_units_c import      units
from invisible_cities.io.pmaps_io            import load_pmaps
from invisible_cities.calib                  import calib_functions         as cf
from invisible_cities.reco                   import calib_sensors_functions as csf
from invisible_cities.database               import load_db                 as DB
from invisible_cities.io.hist_io             import hist_writer
from invisible_cities.core.core_functions    import shift_to_bin_centers


def generate_pdfs():
    """
    Generate SiPM PDFs using Kr RAW data.
    Multiple types of PDF are generated:
    Full spectrum PDFs : using full buffer
    Z vetoed PDFs      : Vetoing all sipms for regions in z where there is an identified s1 or s2
    1 ring vetoed PDFs : Vetoing sipms within 1 ring distance of a hit.
    2 ring vetoed PDFs : As above for 2 rings
    ...
    """

    pmap_file_base = sys.argv[1]
    hit_file_base  = sys.argv[2]
    raw_file_base  = sys.argv[3]

    pmap_files = sorted(glob(pmap_file_base+'*.h5'))## Might not be enough
    ##hit_files  = sorted(glob(hit_file_base+'*.h5'))## Might not be enough
    raw_files  = sorted(glob(raw_file_base+'*.h5'))## Might not be enough

    ## Details of raw waveforms to aid vetoing (assumes all same in run)
    with tb.open_file(raw_files[0]) as rwf_in:
        sipmrwf = rwf_in.root.RD.sipmrwf[0][0]
        wf_range = np.arange(len(sipmrwf))

    run_no = int(sys.argv[4])
    ## Gains and sensor positions
    sipm_gains = DB.DataSiPM(run_no).adc_to_pes.values
    sipm_xy    = DB.DataSiPM(run_no)[['X', 'Y']].values

    ## Start assuming KR data and Kdst
    ## For each event [evt_no, list tuples start and end veto areas]
    reduced_pulse_info = []
    for pmf in pmap_files:
        pmap_dict = load_pmaps(pmf)

        for key, pmap in pmap_dict.items():
            mask_list = []
            for s1 in pmap.s1s:
                mask_list.append(wf_range < s1.times[0]  / units.mus - 1)
                mask_list.append(wf_range > s1.times[-1] / units.mus + 1)
            for s2 in pmap.s2s:
                mask_list.append(wf_range < s2.times[0]  / units.mus - 2)
                mask_list.append(wf_range > s2.times[-1] / units.mus + 2)
            reduced_pulse_info.append([key, np.logical_or.reduce(mask_list)])

    ## output
    histbins = np.arange(-10, 1000, 0.1)
    bin_centres = shift_to_bin_centers(histbins)
    pdf_out = tb.open_file('vetoedPDFs_R'+str(run_no)+'.h5', 'w')
    HIST = partial(hist_writer,
                   pdf_out,
                   group_name  = 'HIST',
                   n_sensors   = 1792,
                   n_bins      = len(bin_centres),
                   bin_centres = bin_centres)
    full_spec = HIST(table_name  = 'full_spec')
    z_vetoed  = HIST(table_name  = 'z_vetoed')
    ## empty arrays for histograms
    shape = 1792, len(bin_centres)
    hist_full_spec = np.zeros(shape, dtype=np.int)
    hist_z_vetoed  = np.zeros(shape, dtype=np.int)

    mask_counter = 0
    for rawf in raw_files:
        print(rawf)
        with tb.open_file(rawf) as raw_in:
            revent_nos = np.fromiter((x[0] for x in raw_in.root.Run.events), np.int)

            indx = np.argwhere(revent_nos==reduced_pulse_info[mask_counter][0])
            print(indx, indx[0][0])
            while indx.shape[0] != 0:
                print(indx[0][0])
                rwf = raw_in.root.RD.sipmrwf[indx[0][0]]
                cwf = csf.sipm_processing["subtract_mode_calibrate"](rwf, sipm_gains)

                hist_full_spec += cf.bin_waveforms(cwf, histbins)
                hist_z_vetoed  += cf.bin_waveforms(cwf[:,reduced_pulse_info[mask_counter][1]],
                                                   histbins)

                mask_counter += 1
                indx = np.argwhere(revent_nos==reduced_pulse_info[mask_counter][0])
                #print(indx, indx[0][0])

    full_spec(hist_full_spec)
    z_vetoed(hist_z_vetoed)
    pdf_out.close()

    
if __name__ == '__main__':
    generate_pdfs()
