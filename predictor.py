# -*- coding: utf-8 -*-
"""
Created on Mon May  8 14:37:05 2017

@author: puller
"""

from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import os, time
import EDI_functions as EDI
import TI_predictor as TI
import seaborn as sns
sns.set_style('darkgrid')

# constants
h = 10**(-8)
fcr = 0.5
Tmin = 0; Tmax = 9
vload_min = None
dilutions_min = None
method = 'LAD'
rframe = 2 #reference frame; set to None to use all sites

fs = 28
H = 8

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

#def TI_from_diversity(DD, func_name, j0jL, cutoff, method, rf = rframe, 
#                      exclude_pat = None, nboot = None):
#    '''estimate the time of infection (TI) from the specified diversity values'''
#    CUT = EDI.window_cutoff(data, func_name, region(j0jL), cutoff, rf = rf)
#    if exclude_pat is not None:
#        CUT.pat_names = list(data['pat_names'])
#        CUT.pat_names.remove(exclude_pat)
#        CUT.assemble_pats()
##        print CUT.pat_names, CUT.ttk.shape
#    ttk, xxk, jjk = CUT.realdata(Tmin, Tmax,  fcr = fcr, vload_min = vload_min,
#                                 dilutions_min = dilutions_min)
#    if nboot is None:
#        ttk_data, dtdx_t0 = EDI.fitmeth_byname(ttk, xxk, method = method)
#        TTest = dtdx_t0[0]*DD + dtdx_t0[1]
#        return TTest, dtdx_t0
#    else:
#        Npat = len(CUT.pat_names)
#        jjboot = np.random.randint(0, high = Npat, size = (nboot, Npat))
#        TTest = np.zeros((nboot, len(DD)))
#        dtdx_t0 = np.zeros((nboot, 2))
#        for jboot, idx_boot in enumerate(jjboot):
#            tk = np.ma.concatenate([ttk[np.where(jjk == j)] for j in idx_boot])
#            xk = np.ma.concatenate([xxk[np.where(jjk == j)] for j in idx_boot])
#            ttk_est, dtdx_t0[jboot,:] = EDI.fitmeth_byname(tk, xk, method = method)
#            TTest[jboot,:] = dtdx_t0[jboot, 0]*DD + dtdx_t0[jboot, 1]
#        return TTest, dtdx_t0

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
    
if __name__=="__main__":
    '''Predicting infection date (for web-application)'''
    plt.ioff()
    plt.close('all')

    outdir_name = './tmp/'
    if not os.path.exists(outdir_name):
        os.makedirs(outdir_name)

    #Creating figures for the manuscript
    t0 = time.time()
    j0jL = coords['pol']
    measure = 'diversity'
    meas = translate_measures(measure)
    cutoff1 = 0.


    ttk_est, ttk, xxk, jjk, dtdx_t0 = ttest_region(meas, j0jL, cutoff1, method)
    CUT = EDI.window_cutoff(data, meas, region(j0jL), cutoff1, rf = rframe)
    ttk, xxk, jjk = CUT.realdata(Tmin, Tmax,  fcr = fcr, vload_min = vload_min,
                             dilutions_min = dilutions_min)
    for j, name in enumerate(data['pat_names']):
        idx_pat = np.where(jjk == j)[0]
        DD = xxk[idx_pat]
        TT = ttk[idx_pat]
        
        TTest, dtdx_t0 = TI.TI_from_diversity(DD, j0jL, cutoff1)

        TTboot, dtdx_boot = TI.TI_from_diversity(DD, j0jL, cutoff1, nboot = 10**3)
        fig, ax = plt.subplots(1, 2, figsize = (2*H, H))
        idx_data = np.where(jjk != j)[0]
        ax[0].scatter(ttk[idx_data], ttk_est[idx_data], color = 'b', marker = 'x', s = 40)
        ax[0].scatter(TT, TTest, color = 'r', marker = 'o', s = 40)
        ax[0].plot(np.sort(ttk), np.sort(ttk), '--k')
        
        ax[0].set_xlabel('TI [years]', fontsize = fs)
        ax[0].set_ylabel('ETI [years]', fontsize = fs)
        ax[0].legend(fontsize = 0.6*fs, loc = 2, ncol = 2)
        ax[0].tick_params(labelsize = .8*fs)
#        ax[0].axis('tight')
        
        for jD, D in enumerate(DD):
            ax[1].hist(TTboot[:,jD])
        ax[1].set_xlabel('ETI [years]', fontsize = fs)
        ax[1].tick_params(labelsize = .8*fs)
        plt.axis('tight')
        plt.savefig(outdir_name + 'check_{}.pdf'.format(name))
        plt.close()
        
        TI.TI_bootstrap_plot(DD, j0jL, cutoff1, outdir_name + 'hist_{}.pdf'.format(name), nboot = 10**3)

    t1 = time.time()
    print t1 - t0
    