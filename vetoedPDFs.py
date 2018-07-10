import sys
import tables as tb
import numpy as np
from glob import glob
from functools import partial

from invisible_cities.core.system_of_units_c import      units
from invisible_cities.io.pmaps_io            import load_pmaps_as_df
from invisible_cities.io.dst_io              import load_dsts
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

    pmap_sorter = sorter_func(pmap_file_base)
    pmap_files = sorted(glob(pmap_file_base+'*.h5'), key=pmap_sorter)
    hit_sorter = sorter_func(hit_file_base)
    hit_files  = sorted(glob(hit_file_base+'*.h5'), key=hit_sorter)
    raw_sorter = sorter_func(raw_file_base)
    raw_files  = sorted(glob(raw_file_base+'*.h5'), key=raw_sorter)

    ## Details of raw waveforms to aid vetoing (assumes all same in run)
    with tb.open_file(raw_files[0]) as rwf_in:
        sipmrwf = rwf_in.root.RD.sipmrwf[0][0]
        wf_range = np.arange(len(sipmrwf))

    run_no = int(sys.argv[4])
    ## Gains and sensor positions
    sipm_gains = DB.DataSiPM(run_no).adc_to_pes.values
    sipm_xy    = DB.DataSiPM(run_no)[['X', 'Y']].values

    

    ## output
    histbins = np.arange(-10, 300, 0.1)
    bin_centres = shift_to_bin_centers(histbins)
    with tb.open_file('vetoedPDFs_R'+str(run_no)+'.h5', 'w') as pdf_out:
        HIST = partial(hist_writer,
                       pdf_out,
                       group_name  = 'HIST',
                       n_sensors   = 1792,
                       n_bins      = len(bin_centres),
                       bin_centres = bin_centres)
        full_spec = HIST(table_name  = 'full_spec')
        z_vetoed  = HIST(table_name  = 'z_vetoed')
        one_ring    = HIST(table_name  = 'one_ring_vetoed')
        two_ring    = HIST(table_name  = 'two_ring_vetoed')
        thr_ring    = HIST(table_name  = 'thr_ring_vetoed')

        ## Hit info
        hit_positions = load_dsts(hit_files, 'DST', 'Events')[['event', 'X', 'Y']].values
        ## Start assuming KR data and Kdst
        ## For each event [evt_no, list tuples start and end veto areas]
        reduced_pulse_info = []
        for pmf in pmap_files:
            print(pmf)
            try:
                ## pmap_dict = load_pmaps(pmf)
                s1s, s2s, _, _, _ = load_pmaps_as_df(pmf)
            except (ValueError, tb.exceptions.NoSuchNodeError):
                print("Empty file. Skipping.")
                continue

            ## for key, pmap in pmap_dict.items():
            for evtNo in s1s['event'].unique():
                evtS1 = s1s[s1s['event'] == evtNo]
                evtS2 = s2s[s2s['event'] == evtNo]
                mask_list = []
                ## for s1 in pmap.s1s:
                for is1 in evtS1['peak'].unique():
                    s1 = evtS1[evtS1['peak'] == is1]
                    mask_list.append((wf_range < s1['time'].iloc[0]  / units.mus - 1) |
                                     (wf_range > s1['time'].iloc[-1] / units.mus + 1) )
                for is2 in evtS2['peak'].unique():
                    s2 = evtS2[evtS2['peak'] == is2]
                    mask_list.append((wf_range < s2['time'].iloc[0]  / units.mus - 2) |
                                     (wf_range > s2['time'].iloc[-1] / units.mus + 2) )
                reduced_pulse_info.append([evtNo, np.logical_and.reduce(mask_list)])
        print('masking info stored')
        mask_counter = 0
        pmap_evts = np.fromiter((x[0] for x in reduced_pulse_info), np.int)
        for rawf in raw_files:
            print(rawf)
            if mask_counter >= len(reduced_pulse_info):
                continue
            try:
                ## empty arrays for histograms
                shape = 1792, len(bin_centres)
                hist_full_spec = np.zeros(shape, dtype=np.int)
                hist_z_vetoed  = np.zeros(shape, dtype=np.int)
                hist_1_vetoed  = np.zeros(shape, dtype=np.int)
                hist_2_vetoed  = np.zeros(shape, dtype=np.int)
                hist_3_vetoed  = np.zeros(shape, dtype=np.int)
                with tb.open_file(rawf) as raw_in:
                    revent_nos = np.fromiter((x[0] for x in raw_in.root.Run.events), np.int)
                    
                    #evt_no = reduced_pulse_info[mask_counter][0]
                    #indx = np.argwhere(revent_nos==evt_no)
                    #print(reduced_pulse_info)
                    #pmap_evts = np.array(reduced_pulse_info)[:, 0]
                    pmap_overlap_indx = np.arange(revent_nos.shape[0])[np.in1d(revent_nos, pmap_evts)]
                    hit_overlap_indx  = np.arange(revent_nos.shape[0])[np.in1d(revent_nos,
                                                  hit_positions[:, 0])]
                    hit_indcs = np.arange(hit_positions[:, 0].shape[0])[np.in1d(hit_positions[:, 0]), revent_nos]
                    hindx = 0
                    #print(indx, indx[0][0])
                    #while indx.shape[0] != 0:
                    for indx in pmap_overlap_indx:
                        #print(indx[0][0])
                        #rwf = raw_in.root.RD.sipmrwf[indx[0][0]]
                        ## cwf = csf.sipm_processing["subtract_mode_calibrate"](raw_in.root.RD.sipmrwf[indx[0][0]], sipm_gains)
                        cwf = csf.sipm_processing["subtract_mode_calibrate"](raw_in.root.RD.sipmrwf[indx], sipm_gains)
                        
                        hist_full_spec += cf.bin_waveforms(cwf, histbins)
                        z_veto = reduced_pulse_info[mask_counter][1]
                        hist_z_vetoed  += cf.bin_waveforms(cwf[:, z_veto], histbins)

                        #dst_indx = np.argwhere(hit_positions[:, 0]==evt_no)
                        #if dst_indx.shape[0] != 0:
                        if indx in hit_overlap_indx:
                            ## hit_p = hit_positions[dst_indx[0][0], 1:]
                            hit_p = hit_positions[hit_indcs[hindx], 1:]
                            hindx += 1
                            hist_1_vetoed += cf.bin_waveforms(ring_veto(cwf, 1, z_veto,
                                                                        hit_p,
                                                                        sipm_xy),
                                                            histbins)
                            hist_2_vetoed += cf.bin_waveforms(ring_veto(cwf, 2, z_veto,
                                                                        hit_p,
                                                                        sipm_xy),
                                                            histbins)
                            hist_3_vetoed += cf.bin_waveforms(ring_veto(cwf, 3, z_veto,
                                                                        hit_p,
                                                                        sipm_xy),
                                                            histbins)

                        mask_counter += 1
                        #if mask_counter < len(reduced_pulse_info):
                        #    evt_no = reduced_pulse_info[mask_counter][0]
                        #    indx = np.argwhere(revent_nos==evt_no)
                        #else:
                            ## dummy evt_no to definitely give no info
                        #    indx = np.argwhere(revent_nos==-100)
                        #print(indx, indx[0][0])
                    full_spec(hist_full_spec)
                    z_vetoed(hist_z_vetoed)
                    one_ring(hist_1_vetoed)
                    two_ring(hist_2_vetoed)
                    thr_ring(hist_3_vetoed) 
            except tb.HDF5ExtError:
                print('corrupt file')
                continue

        ## full_spec(hist_full_spec)
        ## z_vetoed(hist_z_vetoed)
        ## one_ring(hist_1_vetoed)
        ## two_ring(hist_2_vetoed)
        ## thr_ring(hist_3_vetoed)
    #pdf_out.close()


def sorter_func(base_name):
    def sorting_func(elem):
        """
        Forces sort according to file number.
        """
        start_num = elem.split(base_name)[1]

        dot_dist = start_num.find('.')
        bar_dist = start_num.find('_')
        num_len = min(dot_dist, bar_dist)
        if num_len < 0:
            num_len = max(dot_dist, bar_dist)      
        return int(start_num[0:num_len])
    return sorting_func
        

def ring_veto(cwf, n_ring, z_veto, hit_pos, xy):

    pitch = 10
    veto_indcs = np.where((xy[:,0] < hit_pos[0] + n_ring * pitch) &
                          (xy[:,0] > hit_pos[0] - n_ring * pitch) &
                          (xy[:,1] < hit_pos[1] + n_ring * pitch) &
                          (xy[:,1] > hit_pos[1] - n_ring * pitch))
    #print('indices = ', veto_indcs)

    ## Forces out of the histo range (or not if you have weird ranges)
    new_section = cwf[veto_indcs[0]]
    new_section[:, np.invert(z_veto)] = -100000
    cwf[veto_indcs[0]] = new_section
    
    return cwf
    
if __name__ == '__main__':
    generate_pdfs()
