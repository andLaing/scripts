import  os
import sys

import numpy  as np
import pandas as pd

import matplotlib.pyplot as plt

from invisible_cities.core .configure     import            configure
from invisible_cities.icaro.hst_functions import shift_to_bin_centers


def xe137_normalisation(conf_list, spec_shift = 0,
                        suppress_df = False):

    config = configure(conf_list).as_namespace

    flux_file  = os.path.expandvars(config.flux_file)
    acti_file  = os.path.expandvars(config.acti_file)
    out_file   = os.path.expandvars(config.file_out)
    sim_muons  = int(config.n_simulated_muons)
    if not hasattr(config, 'log_bins'):
        log_bins = False
    else:
        log_bins = bool(config.log_bins)
    lab_flux   = float(config.lab_flux)
    lab_flux_e = float(config.lab_flux_err)
    gen_area   = float(config.gen_area)

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
    xe137_exp      = xe137_count / binned_sim_muons
    xe137_exp_err  = xe137_exp * np.sqrt(1 / xe137_count + 1 / binned_sim_muons)

    ## Get flux in the same bins
    for i in range(10):

        vals = pd.read_hdf(flux_file, 'muon_flux_'+str(i))

        histo, _ = np.histogram(vals.E.values + spec_shift, bins = bins)

        try:
            flux_histo += histo
        except NameError:
            flux_histo  = histo

    norm_flux = flux_histo / flux_histo.sum()
    flux_err  = norm_flux * np.sqrt(1 / flux_histo + 1 / flux_histo.sum())
    fluxCMS   = norm_flux * lab_flux
    fluxCMS_e = fluxCMS * np.sqrt((flux_err / norm_flux)**2
                                  + (lab_flux_e / lab_flux)**2)
    fluxS     = fluxCMS * gen_area
    fluxS_e   = gen_area * fluxCMS_e
    xe137S    = xe137_exp * fluxS ## Xe137 per bin per second
    xe137S_e  = xe137S * np.sqrt((xe137_exp_err / xe137_exp)**2
                                 + (fluxS_e / fluxS)**2)
    xe137Y    = xe137S * 3.1536e7
    xe137Y_e  = xe137S_e * 3.1536e7

    total_xe137PS = xe137S.sum()
    perSec_err    = np.sqrt(np.sum(xe137S_e**2))
    total_xe137PY = xe137Y.sum()
    perYr_err     = np.sqrt(np.sum(xe137Y_e**2))
    print('Xe-137 per second = ', total_xe137PS, '+/-', perSec_err)
    print('Xe-137 per calendar yr = ', total_xe137PY, '+/-', perYr_err)

    if suppress_df:
        ## Optimisation or error calc,
        ## just return the per yr prediction
        return total_xe137PY, perYr_err


    plt.errorbar(shift_to_bin_centers(bins), xe137Y, fmt='^',
                 xerr = np.diff(bins) / 2, yerr = xe137Y_e)
    plt.xlabel('Muon energy (GeV)')
    plt.ylabel('Xe-137 expectation per yr per bin')
    plt.show()
    
    ## Output the calculation at each stage to file
    df = pd.DataFrame({"BinMin"         :     bins[:-1],
                       "BinMax"         :      bins[1:],
                       "n_xe137"        :   xe137_count,
                       "xe137PerMu"     :     xe137_exp,
                       "Err137PerMu"    : xe137_exp_err,
                       "NormFlux"       :     norm_flux,
                       "FluxErr"        :      flux_err,
                       "FluxPerCM2PerS" :       fluxCMS,
                       "FluxAreaSErr"   :     fluxCMS_e,
                       "FluxPerS"       :         fluxS,
                       "FluxSErr"       :       fluxS_e,
                       "xe137PerS"      :        xe137S,
                       "xe137SErr"      :      xe137S_e,
                       "xe137PerY"      :        xe137Y,
                       "xe137YErr"      :      xe137Y_e})

    df.to_hdf(out_file, 'xe137')


if __name__ == '__main__':
    xe137_normalisation(sys.argv)
