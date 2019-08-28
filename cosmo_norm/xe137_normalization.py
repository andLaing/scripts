import  os
import sys

import numpy  as np
import pandas as pd

import matplotlib.pyplot as plt

from invisible_cities.core .configure     import            configure
from invisible_cities.icaro.hst_functions import shift_to_bin_centers


if __name__ == '__main__':

    config = configure(sys.argv).as_namespace

    flux_file = os.path.expandvars(config.flux_file)
    acti_file = os.path.expandvars(config.acti_file)
    out_file  = os.path.expandvars(config.file_out)
    sim_muons = int(config.n_simulated_muons)
    if not hasattr(config, 'log_bins'):
        log_bins = False
    else:
        log_bins = bool(config.log_bins)
    lab_flux  = float(config.lab_flux)
    gen_area  = float(config.gen_area)

    if hasattr(config, 'bin_edges'):
        bins = config.bin_edges
    else:
        bin_range = config.bin_range
        if log_bins:
            bins = np.logspace(np.log10(bin_range[0]),
                               np.log10(bin_range[0]), bin_range[2])
        else:
            bins = np.linspace(*bin_range)
        
    binned_sim_muons, _ = np.histogram(np.random.uniform(bins[0]  ,
                                                         bins[-1] ,
                                                         sim_muons),
                                       bins = bins)

    xe137_df = pd.read_hdf(acti_file)
    xe137_df['GeV'] = xe137_df.Xemunrg * 1e-3

    xe137_count, _ = np.histogram(xe137_df.GeV.values, bins = bins)
    xe137_exp = xe137_count / binned_sim_muons

    ## Get flux in the same bins
    for i in range(10):

        vals = pd.read_hdf(flux_file, 'muon_flux_'+str(i))

        histo, _ = np.histogram(vals.E.values, bins = bins)

        try:
            flux_histo += histo
        except NameError:
            flux_histo  = histo

    norm_flux = flux_histo / flux_histo.sum()
    fluxCMS   = norm_flux * lab_flux
    fluxS     = fluxCMS * gen_area
    xe137S    = xe137_exp * fluxS ## Xe137 per bin per second
    xe137Y    = xe137S * 3.1536e7

    df = pd.DataFrame({"BinMin"         :   bins[:-1],
                       "BinMax"         :    bins[1:],
                       "n_xe137"        : xe137_count,
                       "xe137PerMu"     :   xe137_exp,
                       "NormFlux"       :   norm_flux,
                       "FluxPerCM2PerS" :     fluxCMS,
                       "FluxPerS"       :       fluxS,
                       "xe137PerS"      :      xe137S,
                       "xe137PerY"      :      xe137Y})

    total_xe137PS = xe137S.sum()
    total_xe137PY = xe137Y.sum()
    print('Xe-137 per second = ', total_xe137PS)
    print('Xe-137 per calendar yr = ', total_xe137PY)

    
    plt.plot(shift_to_bin_centers(bins), xe137Y, '^')
    plt.xlabel('Muon energy (GeV)')
    plt.ylabel('Xe-137 expectation per yr per bin')
    plt.show()

    df.to_hdf(out_file, 'xe137')
