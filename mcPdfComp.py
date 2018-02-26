#!/usr/bin/env python
import numpy as np
import tables as tb
import matplotlib.pyplot as plt


def main():
    """ Compare mc sipm pdfs gneerated in different ways """

    bins = 0
    modeHist = []
    meanMauHist = []
    # Standard sim
    with tb.open_file('outdats/mcpdfTest_stand.h5', 'r') as inF:
        bins = np.array(inF.root.HIST.sipm_bins)
        modeHist.append( np.array(inF.root.HIST.sipm).sum(axis=0) )
        meanMauHist.append( np.array(inF.root.HIST.sipmMAU).sum(axis=0) )
    # +Ped -mean
    with tb.open_file('outdats/mcpdfTest_meanHack.h5', 'r') as inF:
        modeHist.append( np.array(inF.root.HIST.sipm).sum(axis=0) )
        meanMauHist.append( np.array(inF.root.HIST.sipmMAU).sum(axis=0) )
    # +Ped -mode
    with tb.open_file('outdats/mcpdfTest_meanHack.h5', 'r') as inF:
        modeHist.append( np.array(inF.root.HIST.sipm).sum(axis=0) )
        meanMauHist.append( np.array(inF.root.HIST.sipmMAU).sum(axis=0) )
    ##
    for mo1, mo2, mo3, ma1, ma2, ma3 in zip(*modeHist, *meanMauHist):
        fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(20,6))
        axes[0].bar(bins, mo1, width=0.1, log=True, label='standard sim, mode')
        axes[0].bar(bins, mo2, width=0.1, log=True, label='mean sim, mode')
        axes[0].bar(bins, mo3, width=0.1, log=True, label='mode sim, mode')
        hdls1, lbls1 = axes[0].get_legend_handles_labels()
        axes[0].legend(hdls1, lbls1)
        axes[0].set_xlabel('pe')
        axes[0].set_ylabel('entries')
        axes[1].bar(bins, ma1, width=0.1, log=True, label='standard sim, MAU')
        axes[1].bar(bins, ma2, width=0.1, log=True, label='mean sim, MAU')
        axes[1].bar(bins, ma3, width=0.1, log=True, label='mode sim, MAU')
        hdls2, lbls2 = axes[1].get_legend_handles_labels()
        axes[1].legend(hdls2, lbls2)
        axes[1].set_xlabel('pe')
        axes[1].set_ylabel('entries')
        plt.tight_layout()
        fig.show()
        catch = input('test good?')
        if 'q' in catch:
            exit()
        else:
            plt.close(fig)
            plt.clf()
        
        

if __name__ == '__main__':
    main()
