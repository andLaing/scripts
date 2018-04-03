#!/usr/bin/env python
import matplotlib
import numpy as np
from scipy import signal
import matplotlib.pyplot as plt
from IPython.display import Image
#from ROOT import TFile, TH1F
import sys
import tables as tb
from glob import iglob



#pms = [ 0, 1, 4 ,5, 8, 9, 18, 19, 22, 23, 26, 27 ]
#fis_id = [ 1, 2, 5, 6, 9, 'A', 3, 4, 7, 8, 'B', 'C' ]
#pms = [ 1, 4, 5, 8, 9, 18, 19, 22, 23, 30, 31 ]
#fis_id = [ 8, 9, 'A', 'B', 'C', 1, 2, 3, 4, 5, 6 ]
fis_id = {0:1, 1:2, 4:8, 8:5, 9:6, 12:9, 13:'A', 16:3, 17:4, 24:7, 28:'B', 29:'C'}
## grp1 = [ 0, 1, 12, 13, 6, 7, 4, 5 ]
## grp2 = [ 16, 17, 18, 19 ]
grp1 = [ 0, 1, 16, 17, 24, 4 ]
grp2 = [ 12, 13, 28, 29 ]
grp3 = [ 8, 9 ]

#testFREQ = [ 5625., 8437.5, 15625., 16250., 17500., 17812.5, 18125., 18437.5, 18750., 19062.5, 21250., 21875. ]
## testFREQ = [ 5937.5, 7812.5, 15625., 16250., 17500., 17812.5, 18125., 18437.5, 18750., 19062.5, 21250., 21875. ]
##testFREQ = [ 1250., 1562.5, 2500., 2812.5, 3125., 4375., 4687.5, 5625., 5937.5, 6250., 6562.5, 7500., 8125., 8437.5, 8750., 9687.5, 15312.5, 15625., 15937.5, 16250., 16562.5, 17187.5, 17500., 17812.5, 18125., 18437.5, 18750., 19062.5, 19375., 19687.5, 20000., 20312.5, 20625., 20937.5, 21250., 21875., 22187.5, 22500., 22812.5 ]
testFREQ = [ x*312.5 for x in range(1, 80)]

def main():
    """Using Vicente's script with some additions to study
    evolution of oscillation in pmt baselines"""
    #print(testFREQ)
    nevt = 100

    #max_freq = np.zeros( (len(pms), nevt), dtype=np.float64 )
    #max_mag = np.zeros( (len(pms), nevt), dtype=np.float64 )
    #sign_freqs = [ [] for i in pms ]
    #n_frqs = np.zeros( (len(pms), nevt), dtype=np.float64 )
    max_freq = np.zeros( (12, nevt), dtype=np.float64 )
    max_mag = np.zeros( (12, nevt), dtype=np.float64 )
    sign_freqs = [ [] for i in range(12) ]
    n_frqs = np.zeros( (12, nevt), dtype=np.float64 )
    ##
    ## Filter definition
    fSample = 40E6
    freqLPF = 100E3
    freqLPFd = 2*freqLPF / fSample
    b, a = signal.butter(1, freqLPFd, 'low', analog=False)
    ##

    filtering = True
    defsDone = False
    #indxs = np.array()
    pms = []
    bigFreqs = []
    ## MONITORING
    MAGS = { 1:[ np.array([]) for i in range(len(testFREQ)) ], 2:[ np.array([]) for i in range(len(testFREQ)) ], 3:[ np.array([]) for i in range(len(testFREQ)) ] }
    FREQS = { 1:[ np.array([]) for i in range(len(testFREQ)) ], 2:[ np.array([]) for i in range(len(testFREQ)) ], 3:[ np.array([]) for i in range(len(testFREQ)) ] }
    PHAS = { 1:[ np.array([]) for i in range(len(testFREQ)) ], 2:[ np.array([]) for i in range(len(testFREQ)) ], 3:[ np.array([]) for i in range(len(testFREQ)) ] }
    ##
    cEvt = 0
    testWF = np.empty(1)
    testFill = False
    for fN in iglob(sys.argv[1]+'*_waveforms.h5'):
        while cEvt < 200:
            with load_input(fN) as file0:
                if not defsDone:
                    wf_len = get_wf_len(file0)
                #meanAbs = [ np.zeros(int(wf_len/2+1), np.float64) for i in pms ]
                    meanAbs = [ np.zeros(int(wf_len/2+1), np.float64) for i in range(12)]
                    phases = [ [] for i in range(12) ]
                    frq_plot = np.empty(int(wf_len/2+1))
                    pms = np.fromiter((file0.root.Sensors.DataPMT[i][0] for i in range(len(file0.root.Sensors.DataPMT))), np.int)
                    defsDone = True
                for i in range(len(pms)):
                #indxs.append( getPMid(file0, i) )
                    for evt in range(min(100, len(file0.root.RD.pmtrwf))):
                        if i == 0:
                            cEvt = cEvt + 1

                        sg = getWF(file0, i, evt)

                        sg = sg - np.mean(sg)

                        if filtering:
            ## Filtering
                            sg = signal.lfilter(b, a, sg)
            ##
                        sg_rms = np.std(sg)

                        ft = np.fft.rfft(sg)
                        freq = np.fft.rfftfreq(len(sg), d=25E-9)
                        if i == 0 and evt == 0:
                            frq_plot = freq
                        ftAb = np.absolute(ft)
                        ftPha = np.angle(ft)
                    ## Section max. monitoring
                        grp = 1
                        if pms[i] in grp2:
                            grp = 2
                        elif pms[i] in grp3:
                            grp = 3
                        for ishite, fq in enumerate(testFREQ):
                            inx = np.abs(freq-fq).argmin()
                            MAGS[grp][ishite] = np.append(MAGS[grp][ishite], ftAb[inx] )
                            PHAS[grp][ishite] = np.append(PHAS[grp][ishite], ftPha[inx] )

                        if grp == 2 and evt == 0:
                            inx = np.argmax(np.select([(freq>4000)&(freq<6000)], [ftAb]))
                            inx2 = np.argmax(np.select([(freq>7500)&(freq<9500)], [ftAb]))
                        #print(freq[inx], freq[inx2])
                        #print(ftAb[inx+9], freq[inx+9], ftPha[inx+9])
                        ##if pms[indxs[-1]] == 9 and not testFill:
                            if pms[i] == 18 and not testFill:
                                testWF = sg
                                testFill = True
                    ##
                        meanAbs[i] += ftAb#/nevt
                        phases[i].append(ftPha)
                        max_freq[i,evt] = freq[np.argmax(ftAb)]
                        max_mag[i,evt] = np.amax(ftAb)/(sg_rms*len(sg))

            #for j in range(len(imp_freq)):
            #    mag_imp[j,i,evt] = ftAb[np.argwhere(freq==imp_freq[j])]
            #    pha_imp[j,i,evt] = np.angle(ft)[np.argwhere(freq==imp_freq[j])]

                        sign_freqs[i].append( np.fromiter((i for (i,j) in zip(freq, ftAb) if j >= 0.05*len(sg)*sg_rms and i < 50000), np.float64) )
                        n_frqs[i,evt] = len(sign_freqs[i][-1])

    
    #print('Indices: ', indxs)
    fig, axes = plt.subplots(nrows=4,ncols=3, figsize=(20,6))
    for i, ax, row in zip(range(len(pms)), axes.flatten(), max_freq):
        ax.plot(np.arange(1, nevt+1), row)
        #j = indxs[i]#getPMid(file0, i)
        ax.set_title('Max. freq. pmt '+str(pms[i])+' ('+str(fis_id[pms[i]])+')')

    #plt.tight_layout()
    fig.show()

    #fig2, axes2 = plt.subplots(nrows=2,ncols=4)
    #for i, (ax, fr) in enumerate(zip(axes2.flatten(), mag_imp)):
    #    for k, row in enumerate(fr):
    #        j = getPMid(file0, k)
    #        ax.plot(np.arange(1, nevt+1), row/mag_imp[i,0,:])#, label='rel. mag. pmt '+str(pms[j]))
    #    ax.set_xlim(0, 20)
        #ax.set_ylim(-10,10)
    #    ax.set_title('Magnitude of at '+str(imp_freq[i])+' Hz')

    #fig2.show()
    print('nEvt check ', cEvt)
    fig3, axes3 = plt.subplots(nrows=4,ncols=3, figsize=(20,6))
    for i, ax, row in zip(range(len(pms)), axes3.flatten(), meanAbs):
        ax.plot(frq_plot, row/cEvt)
        bigFreqs.append([(mg,fq) for (mg,fq) in zip(row/cEvt,frq_plot) if mg > 1500 and mg < 2000 and fq >4500 and fq<10000])
        ax.set_xlim(500, 25000)
        #j = indxs[i]#getPMid(file0, i)
        #print('pmt ', pms[j], frq_plot[np.argmax(row)])
        ax.set_title('mean fft abs. val. pmt '+str(pms[i])+' ('+str(fis_id[pms[i]])+')')

    plt.tight_layout()
    fig3.show()
    #print(bigFreqs)
    ##fig3.savefig('meanFFTPMTs.pdf')

    #fig4, axes4 = plt.subplots(nrows=2,ncols=4)
    #for k, (ax, row) in enumerate(zip(axes4.flatten(), pha_imp)):
    #    for i, fr in enumerate(row):
    #        j = getPMid(file0, i)
    #        ax.plot(np.arange(1, nevt+1), fr/pha_imp[k,0,:])#, label='rel. phase pmt '+str(pms[j]))
    #    ax.set_xlim(0, 20)
        #ax.set_ylim(-50,50)
        #ax.legend(loc=0)
    #    ax.set_title('Phase '+str(imp_freq[k])+' Hz')

    ## Monitoring
    fig4, axes4 = plt.subplots(nrows=5,ncols=2, figsize=(20,15))
    print('Group 1:')
    # Get the MAGNITUDES and PHASES fro the group of interest.
    indx0 = np.argwhere(pms==0)[0][0]
    #grp = [ 1, 4, 16, 17, 24 ]
    grp = [ 8, 9 ]
    #grp = [ 12, 13, 28, 29 ]
    i#ndxs = [ np.argwhere(pms==1)[0][0], np.argwhere(pms==12)[0][0], np.argwhere(pms==13)[0][0], np.argwhere(pms==16)[0][0], np.argwhere(pms==17)[0][0] ]
    #indxs = [ np.argwhere(pms==7)[0][0], np.argwhere(pms==18)[0][0], np.argwhere(pms==19)[0][0] ]
    indxs = [ np.argwhere(pms==grp[0])[0][0], np.argwhere(pms==grp[1])[0][0] ]
    #for i in range(5):
    for i in range(2):
    #for i in range(3):
        axes4[i][0].plot(meanAbs[indx0]/cEvt, meanAbs[indxs[i]]/cEvt)
        axes4[i][0].set_title('Magnitudes Elecid 0 vs. Magnitudes Elecid '+str(grp[i]))
        axes4[i][0].set_xlabel('Frequency mean magnitudes elecid 0')
        axes4[i][0].set_ylabel('Frequency mean magnitudes elecid '+str(grp[i]))
        axes4[i][1].hist(np.concatenate(np.subtract(phases[indx0], phases[indxs[i]])))
        axes4[i][1].set_title('Phases Elecid 0 Phases Elecid '+str(grp[i]))
        axes4[i][1].set_xlabel('Frequency phases elecid 0 - elcid '+str(grp[i]))
        #axes4[i][1].set_ylabel('Frequency phases elecid '+str(grp[i]))
    plt.tight_layout()
    fig4.show()
    ## for ip1, row in enumerate(axes4):
    ##     ## row[0].hist(MAGS[1][ip1])
    ##     ## row[1].hist(PHAS[1][ip1])
    ##     print(np.mean(MAGS[1][ip1]), ', ', np.std(MAGS[1][ip1], ddof=1))
    ##     row[0].plot(MAGS[1][ip1], PHAS[1][ip1], 'bo')
    ##     if ip1 != len(axes4)-10:
    ##         row[1].plot(range(len(PHAS[1][ip1])), PHAS[1][ip1]/PHAS[1][-10], 'b-')
    ##         row[1].set_ylim(-10, 10)
    ##     else:
    ##         row[1].plot(range(len(PHAS[1][ip1])), PHAS[1][ip1])
    ##freqFile = tb.open_file('osc_magnitudes.h5', 'w')
    ##atca1 = np.array([])
    ##atca2 = np.array([])
    ##t1 = np.arange(0.,3.2e-3,25e-9)
    ##t1 = t1[:-1]
    for ip1 in range(len(MAGS[1])):
        magMean = np.mean(MAGS[1][ip1])
        magRMS = np.std(MAGS[1][ip1], ddof=1)
        ##atca1 = np.append(atca1, round(magMean, 1)/len(t1))
        print(magMean, ', ', magRMS, magRMS/magMean)

    print('correlation atca1: ', np.corrcoef(PHAS[1]), np.max(np.corrcoef(PHAS[1])))
    ## plt.tight_layout()
    ## axes4[0][0].hist(MAGS[1][0], bins=100)
    ## axes4[0][1].hist(FREQS[1][0], bins=100)
    ## axes4[0][2].hist(PHAS[1][0], bins=100)
    ## #axes4[0][2].plot(range(cEvt), PHAS[1][0]/PHAS[1][3])
    ## #axes4[0][2].set_ylim(-5, 5)
    ## axes4[1][0].hist(MAGS[1][1], bins=100)
    ## axes4[1][1].hist(FREQS[1][1], bins=100)
    ## axes4[1][2].hist(PHAS[1][1], bins=100)
    ## #axes4[1][2].plot(range(cEvt), PHAS[1][1]/PHAS[1][3])
    ## #axes4[1][2].set_ylim(-5, 5)
    ## axes4[2][0].hist(MAGS[1][2], bins=100)
    ## axes4[2][1].hist(FREQS[1][2], bins=100)
    ## axes4[2][2].hist(PHAS[1][2], bins=100)
    ## #axes4[2][2].plot(range(cEvt), PHAS[1][2]/PHAS[1][3])
    ## #axes4[2][2].set_ylim(-5, 5)
    ## axes4[3][0].hist(MAGS[1][3], bins=100)
    ## axes4[3][1].hist(FREQS[1][3], bins=100)
    ## axes4[3][2].hist(PHAS[1][3], bins=100)
    #axes4[3][2].plot(range(cEvt), PHAS[1][3])
    ## fig4.show()
    ## #fig4.savefig('corrs_group1.png')
    ## fig5, axes5 = plt.subplots(nrows=len(testFREQ),ncols=2, figsize=(20,15))
    print('Group 2:')
    ## for ip2, row in enumerate(axes5):
    ##     ##row[0].hist(MAGS[2][ip2])
    ##     ##row[1].hist(PHAS[2][ip2])
    ##     print(np.mean(MAGS[2][ip2]), ', ', np.std(MAGS[2][ip2], ddof=1))
    ##     row[0].plot(MAGS[2][ip2], PHAS[2][ip2], 'bo')
    ##     if ip2 != len(axes5)-10:
    ##         row[1].plot(range(len(PHAS[2][ip2])), PHAS[2][ip2]/PHAS[2][-10],'b-')
    ##         row[1].set_ylim(-10, 10)
    ##     else:
    ##         row[1].plot(range(len(PHAS[2][ip2])), PHAS[2][ip2])    
    for ip2 in range(len(MAGS[2])):
        magMean = np.mean(MAGS[2][ip2])
        #atca2 = np.append(atca2, round(magMean, 1)/len(t1))
        magRMS = np.std(MAGS[2][ip2], ddof=1)
        print(magMean, ', ', magRMS, magRMS/magMean)
    print('Group 3:')
    for ip3 in range(len(MAGS[3])):
        magMean = np.mean(MAGS[3][ip3])
        #atca2 = np.append(atca2, round(magMean, 1)/len(t1))
        magRMS = np.std(MAGS[3][ip3], ddof=1)
        print(magMean, ', ', magRMS, magRMS/magMean)
    #freqFile.create_array(freqFile.root, 'atca1_mag', atca1)
    #freqFile.create_array(freqFile.root, 'atca2_mag', atca2)
    #freqFile.close()
    ## plt.tight_layout()
    ## axes5[0][0].hist(MAGS[2][0], bins=100)
    ## axes5[0][1].hist(FREQS[2][0], bins=100)
    ## axes5[0][2].hist(PHAS[2][0], bins=100)
    ## #axes5[0][2].plot(range(cEvt), PHAS[2][0]/PHAS[2][3])
    ## #axes5[0][2].set_ylim(-5, 5)
    ## axes5[1][0].hist(MAGS[2][1], bins=100)
    ## axes5[1][1].hist(FREQS[2][1], bins=100)
    ## axes5[1][2].hist(PHAS[2][1], bins=100)
    ## #axes5[1][2].plot(range(cEvt), PHAS[2][1]/PHAS[2][3])
    ## #axes5[1][2].set_ylim(-5, 5)
    ## axes5[2][0].hist(MAGS[2][2], bins=100)
    ## axes5[2][1].hist(FREQS[2][2], bins=100)
    ## axes5[2][2].hist(PHAS[2][2], bins=100)
    ## #axes5[2][2].plot(range(cEvt), PHAS[2][2]/PHAS[2][3])
    ## #axes5[2][2].set_ylim(-5, 5)
    ## axes5[3][0].hist(MAGS[2][3], bins=100)
    ## axes5[3][1].hist(FREQS[2][3], bins=100)
    ## axes5[3][2].hist(PHAS[2][3], bins=100)
    ## #axes5[3][2].plot(range(cEvt), PHAS[2][3])
    ##fig5.show()
    #fig5.savefig('corrs_group2.png')
    ## print('sec1: ', np.mean(MAGS[1][0]), np.std(MAGS[1][0], ddof=1), np.mean(FREQS[1][0]), np.std(FREQS[1][0], ddof=1))
    ## print('sec2: ', np.mean(MAGS[1][1]), np.std(MAGS[1][1], ddof=1), np.mean(FREQS[1][1]), np.std(FREQS[1][1], ddof=1))
    ## print('sec3: ', np.mean(MAGS[1][2]), np.std(MAGS[1][2], ddof=1), np.mean(FREQS[1][2]), np.std(FREQS[1][2], ddof=1))
    ## print('sec4: ', np.mean(MAGS[1][3]), np.std(MAGS[1][3], ddof=1), np.mean(FREQS[1][3]), np.std(FREQS[1][3], ddof=1))
    
    ## print('sec1: ', np.mean(MAGS[2][0]), np.std(MAGS[2][0], ddof=1), np.mean(FREQS[2][0]), np.std(FREQS[2][0], ddof=1))
    ## print('sec2: ', np.mean(MAGS[2][1]), np.std(MAGS[2][1], ddof=1), np.mean(FREQS[2][1]), np.std(FREQS[2][1], ddof=1))
    ## print('sec3: ', np.mean(MAGS[2][2]), np.std(MAGS[2][2], ddof=1), np.mean(FREQS[2][2]), np.std(FREQS[2][2], ddof=1))
    ## print('sec4: ', np.mean(MAGS[2][3]), np.std(MAGS[2][3], ddof=1), np.mean(FREQS[2][3]), np.std(FREQS[2][3], ddof=1))

    #osc = oscSim(t1, MAGS[2], PHAS[2], testFREQ, 0)
    #print([ x[0] for x in MAGS[1] ])
    #print([ y[0] for y in PHAS[1] ])
    ## osc = 2*MAGS[1][3][0]*np.cos(2*np.pi*FREQS[1][3][0]*t1+PHAS[1][3][0])/len(t1) +2*MAGS[1][2][0]*np.cos(2*np.pi*FREQS[1][2][0]*t1+PHAS[1][2][0])/len(t1) +2*MAGS[1][1][0]*np.cos(2*np.pi*FREQS[1][1][0]*t1+PHAS[1][1][0])/len(t1) +2*MAGS[1][0][0]*np.cos(2*np.pi*FREQS[1][0][0]*t1+PHAS[1][0][0])/len(t1)+2*6416.*np.cos(2*np.pi*17812.5*t1-0.78)/len(t1)+2*6565.*np.cos(2*np.pi*18437.5*t1-1.36)/len(t1)+2*3281.1*np.cos(2*np.pi*21250*t1-2.7)/len(t1)+2*4045*np.cos(2*np.pi*21875*t1-0.48)/len(t1)+2*3121.*np.cos(2*np.pi*15625.0*t1+0.86)/len(t1)+2*4970.2*np.cos(2*np.pi*18750*t1+0.9)/len(t1)+2*4388.3*np.cos(2*np.pi*18125*t1-2.18)/len(t1)+2*3813.1*np.cos(2*np.pi*17500*t1+1.58)/len(t1)##+2*4578.5*np.cos(2*np.pi*17187.5*t1-0.22)/len(t1)+2*3306.9*np.cos(2*np.pi*19687.5*t1+1.58)/len(t1)+2*3591.9*np.cos(2*np.pi*20312.5*t1-1.8)/len(t1)+2*4970.2*np.cos(2*np.pi*18750*t1+0.9)/len(t1)+2*4388.3*np.cos(2*np.pi*18125*t1-2.18)/len(t1)+2*3813.1*np.cos(2*np.pi*17500*t1+1.58)/len(t1)+2*3281.1*np.cos(2*np.pi*21250*t1-2.7)/len(t1)##+2*MAGS[1][1][0]*np.cos(2*np.pi*FREQS[1][1][0]*t1+PHAS[1][1][0])/len(t1) +2*MAGS[1][0][0]*np.cos(2*np.pi*FREQS[1][0][0]*t1+PHAS[1][0][0])/len(t1) ##+2*6416.*np.cos(2*np.pi*17812.5*t1-0.78)/len(t1)+2*6565.*np.cos(2*np.pi*18437.5*t1-1.36)/len(t1)+2*3665.*np.cos(2*np.pi*20312.5*t1-1.6)/len(t1)+2*3351.8*np.cos(2*np.pi*21250.0*t1-2.5)/len(t1)+2*3121.*np.cos(2*np.pi*15625.0*t1+0.86)/len(t1)+2*4578.5*np.cos(2*np.pi*17187.5*t1-0.22)/len(t1)##+2*6416.*np.cos(2*np.pi*17812.5*t1-0.78)/len(t1)
    ##osc = 2*MAGS[1][3][0]*np.cos(2*np.pi*FREQS[1][3][0]*t1+PHAS[1][3][0])/len(t1) +2*MAGS[1][2][0]*np.cos(2*np.pi*FREQS[1][2][0]*t1+PHAS[1][2][0])/len(t1)  +2*6565.*np.cos(2*np.pi*18437.5*t1-1.36)/len(t1)+2*6416.*np.cos(2*np.pi*17812.5*t1-0.78)/len(t1)
    ##6565.71207864 18437.5 -1.36222014797
    ##6416.04072096 17812.5 -0.777775987157
    ##5853.91808889 16250.0 0.658226971655
    ##3665.12974529 20312.5 -1.5941896121
    ##3281.14180306 21250.0 -2.71376460918
    ##3121.00028538 15625.0 0.86317009025
    ##4578.47876823 17187.5 -0.21506481351
    ##3306.89965379 19687.5 1.57968920445
    ##3591.86462791 20312.5 -1.7954724426
    ##4970.24085215 18750.0 0.913235386377
    ##4388.30779441 18125.0 -2.18433202124
    ##3813.10819438 17500.0 1.58188922597
    ##4045.05555794 21875.0 -0.484227927918
    #print('Check: ', FREQS[1][3][0], FREQS[1][2][0], FREQS[1][1][0], FREQS[1][0][0])
    #magFFT = np.absolute(np.fft.rfft(osc))
    #frqs = np.fft.rfftfreq(len(osc), d=25E-9)

    #fig6, axes6 = plt.subplots(nrows=1, ncols=2, figsize=(20,6))
    #axes6[0].plot(t1, testWF)
    #axes6[0].plot(t1, osc)
    #axes6[1].plot(frqs, magFFT)
    #axes6[1].set_xlim(0,25000)
    #fig6.show()
    #fig6.savefig('grp2evt0_simtest.pdf')
        
    #fig4.show()
    #raw_input("ready to move on?")
    #input("ready to move on?")

    #close_file(file0)


def oscSim(tm, magnitude, phase, frequencies, evNo):

    val = 0
    for fq, mg, ph in zip(frequencies, magnitude, phase):
        val += 2*mg[evNo]*np.cos(2*np.pi*fq*tm+ph[evNo])/len(tm)

    return val
    
def load_input(fName):
    """ load and return the input file """

    #if '.h5' in fName:
    return tb.open_file( fName, 'r' )
    #else:
    #    return TFile( sys.argv[1] )


def get_wf_len(inF):
    """ check the length of the wf in samples """

    #if '.h5' in sys.argv[1]:
    return len(inF.root.RD.pmtrwf[0][0])
    #else:
     #   return inF.Get('histP/pmt_0_trace_evt_2').GetNbinsX()


def getWF(inF, ipm, ievt):
    """ Get a specific waveform from file """

    #if '.h5' in sys.argv[1]:
    return inF.root.RD.pmtrwf[ievt][ipm]
    #else:
    #   wf_name = 'histP/pmt_'+str(pms[ipm])+'_trace_evt_'+str(ievt+1) 
    #   histo = inF.Get( wf_name )
    #   max_bin = histo.GetNbinsX()+1
    #   return np.fromiter((histo.GetBinContent(i) for i in range(1,max_bin)), np.int32)


def getPMid(inF, i):
    """ return the index for the pm for identification """

    #if '.h5' in sys.argv[1]:
    return pms.index(inF.root.Sensors.DataPMT[i][0])
    #else:
    #    return i


def close_file(inF):
    #if '.h5' in sys.argv[1]:
    inF.close()
    #else:
    #    inF.Close()

if __name__ == '__main__':
    main()
