# -*- coding: utf-8 -*-
"""
Created on Fri Jun 30 12:16:26 2017

@author: puller
"""

from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import os, sys
import EDI_functions as EDI
import TI_predictor as TI
from datetime import datetime
import seaborn as sns
#sns.set_style('darkgrid')
sns.set(style = 'darkgrid', font = u'Verdana')
#sns.set_style({'axes.grid':True, "axes.facecolor": ".9"})
#sns.despine()

# constants
h = 10**(-8)
fcr = 0.5
Tmin = 0; Tmax = 9
vload_min = None
dilutions_min = None
method = 'LAD'
rframe = 2 #reference frame; set to None to use all sites

fs = 28
fs1 = 42
H = 8

cols = ['b', 'g','r','c', 'm']
marks = ['o', 's', '^', 'd', 'x', 'p', 'v', '<', '>', '1', '2', '3', '4', '*', 'h', '+']
# The genome annotations
datapath = './Frequency_Data/'
head = ['name', 'x1', 'x2', 'width', 'ri']
annot = []
with open(datapath + 'annotations.txt', 'r') as fhandle:
    for line in fhandle:
        l = [x if j ==0 else int(x) for j, x in enumerate(line.split())]
        annot.append({name: l[j] for j, name in enumerate(head)})
coords = {anno['name']: (anno['x1'], anno['x2']) for anno in annot}
feas = ['gag', 'pol', 'env']


#loading frequency data
data = EDI.load_patient_data(patient_names = 'all', filepath = datapath)
Npat = len(data['pat_names'])

def region(j0jL):
    if type(j0jL) is str:
        return coords[j0jL]
    else:
        return j0jL

def translate_measures(measure):
    names = ['polymorphic', 'diversity', 'entropy']
    funcnames = ['ambiguous_above', 'hamming_above', 'entropy_above']
    return funcnames[names.index(measure)]

def ttest_region(func_name, j0jL, cutoff, method, rf = rframe):
    CUT = EDI.window_cutoff(data, func_name, region(j0jL), cutoff, rf = rf)
    print CUT.ttk.shape
    ttk, xxk, jjk = CUT.realdata(Tmin, Tmax,  fcr = fcr, vload_min = vload_min,
                                 dilutions_min = dilutions_min)
    ttk_est = np.zeros(ttk.shape)
    dtdx_t0 = np.zeros((Npat, 2))
    for jpat in xrange(Npat):
        idx_pat = np.where(jjk == jpat)[0]
        idx_data = np.where(jjk != jpat)[0]
        ttk_data, dtdx_t0[jpat,:] = EDI.fitmeth_byname(ttk[idx_data], xxk[idx_data], method = method)
        ttk_est[idx_pat] = dtdx_t0[jpat,0]*xxk[idx_pat] + dtdx_t0[jpat,1]
    return ttk_est, ttk, xxk, jjk, dtdx_t0


def loadK31(reg, filepath, fromHIV = False):
    data = {}
    if fromHIV:
        sys.path.append("/scicore/home/neher/neher/HIV/hivwholeseq")    
        from hivwholeseq.patients.patients import load_patients, Patient
        pats = load_patients(csv=True)
        fmt = "%d/%m/%Y"
        fhandle = open(filepath + 'K31_info_{}.txt'.format(reg), 'w')
        for pcode, pat in pats.iterrows():
            try:
                EDI = datetime.strptime(pat["infect date best"], fmt)
                P = Patient(pat)
                aft = P.get_allele_frequency_trajectories(reg)[0]
                for si, (scode, sample) in enumerate(P.samples.iterrows()):
                    try:
                        date = datetime.strptime(sample["date"], fmt)
                        af = aft[si]
                        TI = date.toordinal() - EDI.toordinal()
                        fhandle.write('{}\t{}\t{}\n'.format(pcode, scode, TI))
                        np.save(filepath + '{}_{}_{}_data.npy'.format(pcode, scode, reg), af.data)
                        np.save(filepath + '{}_{}_{}_mask.npy'.format(pcode, scode, reg), af.mask)
                        data['{}_{}'.format(pcode,scode)] = (date.toordinal() - EDI.toordinal(), af)
                        print(pcode, scode, "WORKED!!!")
                    except:
                            print(scode, "didn't work")
    
            except:
                print("skipping patient ", pcode)
        fhandle.close()
    else:
        with open(filepath + 'K31_info_{}.txt'.format(reg), 'r') as fhandle:
            for line in fhandle:
                words = line.split()
                pat_name = '_'.join(words[:2])
                af_data = np.load(filepath + '{}_{}_data.npy'.format(pat_name, reg))
                af_mask = np.load(filepath + '{}_{}_mask.npy'.format(pat_name, reg))
                af = np.ma.masked_array(af_data, mask = af_mask)
                data[pat_name] = (int(words[2]), af)
    return data
    
    
def K31_diversity(data, cutoff, rf = rframe):
    TT = np.zeros(len(data.keys()))
    DD = np.zeros(len(data.keys()))
    pats = []
    samples = []
    for j, key in enumerate(data.keys()):
        (TT[j], af) = data[key]
        af = af[:4,:]
        if rf is not None:
            af = af[:,rf::3]
        mask = (1.0 - af.max(axis=0))>cutoff
        DD[j] = ((1.0 - (af**2).sum(axis=0))*mask).mean()
        j = key.rfind('_')
        pats.append(key[:j])
        samples.append(key[j+1:])
#        idx = np.where(np.max(af, axis = 0) <= 1. - cutoff)[0]
#        DD[j] = np.mean(1. - np.ma.sum(af[:,idx]**2, axis = 0))
    return TT/365.25, DD, pats, samples
         
if __name__=="__main__":
    '''Predicting infection dates for the 31 additional patients'''
    plt.ioff()
    plt.close('all')

    outdir_name = './plotK31/'
    if not os.path.exists(outdir_name):
        os.makedirs(outdir_name)

    #Creating figures for the manuscript
#    t0 = time.time()
    reg = 'gag'
    j0jL = coords[reg]
    measure = 'diversity'
    meas = translate_measures(measure)
    cutoff1 = 0.01

#    ttk_est, ttk, xxk, jjk, dtdx_t0 = ttest_region(meas, j0jL, cutoff1, method)
    CUT = EDI.window_cutoff(data, meas, region(j0jL), cutoff1, rf = rframe)
    ttk, xxk, jjk = CUT.realdata(Tmin, Tmax,  fcr = fcr, vload_min = vload_min,
                             dilutions_min = dilutions_min)
    ttk_est, dtdx = TI.TI_from_diversity(xxk, j0jL, cutoff1, rf = rframe)
    
    K31data = loadK31(reg, './K31_data/{}/'.format(reg))
    TT, DD, pats, samples = K31_diversity(K31data, cutoff1)
    TTest, dtdx_t0 = TI.TI_from_diversity(DD, j0jL, cutoff1, rf = rframe)
#    print DD.shape, TT.shape
#    print K31data.keys()
#    print data.keys()
    
    Tlim = np.max(TT)
    jj = np.where(ttk <= Tlim)
#    fig, ax = plt.subplots(1,3, figsize = (3*H, H))
#    ax[0].plot(ttk, xxk, '^', label = 'HIVEVO')    
#    ax[0].plot(TT, DD, 'o', label = 'K31')
#    ax[0].plot(dtdx_t0[0]*np.sort(DD) + dtdx_t0[1], np.sort(DD), '--k')
#    ax[1].plot(ttk, ttk_est, '^')
#    ax[1].plot(TT, TTest, 'o')
#    ax[1].plot(np.sort(TT), np.sort(TT), '--k')
##    ax[2].hist(ttk_est-ttk, alpha = 0.5)
#    ax[2].hist(ttk_est[jj]-ttk[jj], alpha = 0.5)
#    ax[2].hist(TTest - TT, alpha = 0.5)
#    ax[0].legend(fontsize = 0.8*fs, loc = 0)
#    ax[0].set_xlabel('TI [years]', fontsize = fs)
#    ax[0].set_ylabel('diversity', fontsize = fs)
#    ax[1].set_xlabel('TI [years]', fontsize = fs)
#    ax[1].set_ylabel('ETI [years]', fontsize = fs)
#    ax[2].set_xlabel('ETI - TI [years]', fontsize = fs)
#    plt.tight_layout()
##    ax[2].legend(('HIVEVO', 'HIVEVO (TI < 5 years)', 'K31'), fontsize = 0.8*fs)
#    fig.savefig(outdir_name + 'K31_{}_rf2.pdf'.format(reg))
#    plt.close()
    
#    fig, ax = plt.subplots(1,1, figsize = (H, H))
#    ax.plot(ttk, xxk, '^', markersize = 12, label = 'training')    
#    ax.plot(TT, DD, 'o', markersize = 12, label = 'validation')
#    ax.plot(dtdx_t0[0]*np.sort(DD) + dtdx_t0[1], np.sort(DD), '--k')
#    ax.legend(fontsize = 0.8*fs, loc = 0)
#    ax.set_xlabel('TI [years]', fontsize = fs)
#    ax.set_ylabel('diversity', fontsize = fs)
#    ax.tick_params(labelsize = .8*fs)
#    plt.savefig(outdir_name + 'K31_{}_diversity_rf2.pdf'.format(reg))
#    plt.close()
    
    fig, ax = plt.subplots(1,1, figsize = (H, H))
    ax.plot(ttk, ttk_est, '^', markersize = 12, label = 'training')
    ax.plot(TT, TTest, 'o', markersize = 12, label = 'validation')
    ax.plot(np.sort(TT), np.sort(TT), '--k')
    ax.legend(fontsize = 0.8*fs, loc = 0)
    ax.set_xlabel('TI [years]', fontsize = fs)
    ax.set_ylabel('ETI [years]', fontsize = fs)
    ax.tick_params(labelsize = .8*fs)
    plt.savefig(outdir_name + 'K31_{}_ETIvsTI_rf2.pdf'.format(reg))
    plt.close()
    
    fig, ax = plt.subplots(1,1, figsize = (H, H))
    ax.hist(ttk_est[jj]-ttk[jj], alpha = 0.5, label = 'training')
    ax.hist(TTest - TT, alpha = 0.5, label = 'validation')
    ax.legend(fontsize = 0.8*fs, loc = 0)
    ax.set_xlabel('ETI - TI [years]', fontsize = fs)
    ax.tick_params(labelsize = .8*fs)
    plt.savefig(outdir_name + 'K31_{}_MAE_rf2.pdf'.format(reg))
    plt.close()
    
    
#    pats_unique = list(set(pats))
#    pairs = [[j for j, p in enumerate(pats) if p ==pat] for pat in pats_unique]
#    fig, ax = plt.subplots(1,1, figsize = (H, H))
#    ax.plot(ttk, xxk, '^', markersize = 12, label = 'training')
#    for jp, p in enumerate(pairs):
#        ax.plot(TT[p], DD[p], '--o', markersize = 12, label = pats_unique[jp])
#    ax.plot(dtdx_t0[0]*np.sort(DD) + dtdx_t0[1], np.sort(DD), '--k')
##    ax.legend(fontsize = 0.8*fs, loc = 0)
#    ax.set_xlabel('TI [years]', fontsize = fs)
#    ax.set_ylabel('diversity', fontsize = fs)
#    ax.tick_params(labelsize = .8*fs)
#    plt.savefig(outdir_name + 'K31_{}_diversity_bypat_rf2.pdf'.format(reg))
#    plt.close()
    
    
    fig, ax = plt.subplots(1, 3, figsize = (3*H, 2*H), sharey = True)
    for j in xrange(3):
        CUT = EDI.window_cutoff(data, meas, region(j0jL), cutoff1, rf = j)
        ttk, xxk, jjk = CUT.realdata(Tmin, Tmax, fcr = fcr)
        ax[j].plot(ttk, xxk, '^', markersize = 12, label = 'training')
            
        T, D, pats, samples = K31_diversity(K31data, cutoff1, rf = j)
        jj = np.argsort(T)
        ax[j].plot(T[jj], D[jj], 'o', markersize = 12, label='validation')
        ax[j].tick_params(labelsize = .8*fs1)
        ax[j].set_xlabel('TI [years]', fontsize = fs1)
        ax[j].set_title('codon pos {}'.format(j+1), fontsize = fs1)
    ax[0].set_ylabel('diversity', fontsize = fs1)
    ax[0].legend(fontsize = 0.8*fs1, loc = 0)
    fig.subplots_adjust(wspace = 0.1)
    plt.savefig(outdir_name + 'K31_{}_diversity.pdf'.format(reg))
    plt.close()
#    return None