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
sys.path.append("/scicore/home/neher/neher/HIV/hivwholeseq")    
from hivwholeseq.patients.patients import load_patients, Patient
#import seaborn as sns
#sns.set_style('darkgrid')

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

def dt(d1,d2):
    	return (d1.toordinal()-d2.toordinal())/365.25

def loadK31(func_name, reg, cutoff, rf = rframe):
    pats = load_patients(csv=True)
    fmt = "%d/%m/%Y"
    res = []
    for pcode, pat in pats.iterrows():
        try:
            EDI = datetime.strptime(pat["infect date best"], fmt)
            P = Patient(pat)
            aft = P.get_allele_frequency_trajectories(reg)[0]
            for si, (scode, sample) in enumerate(P.samples.iterrows()):
                try:
                    date = datetime.strptime(sample["date"], fmt)
                    if rf is None:
                        af = aft[si]
                    else:
                        af = aft[si][:,rf::3]
                    mask = (1.0 - af.max(axis=0))>cutoff
                    div = ((1.0 - (af**2).sum(axis=0))*mask).mean()
                    print(EDI, date, div)
                    res.append((dt(date, EDI), div))
                except:
                        print(scode, "didn't work")

        except:
            print("skipping patient ", pcode)
    res = np.array(res)
    return res[:,0], res[:,1]
    
if __name__=="__main__":
    '''Predicting infection dates for the 31 additional patients'''
    plt.ioff()
    plt.close('all')

    outdir_name = './plotK31/'
    if not os.path.exists(outdir_name):
        os.makedirs(outdir_name)

    #Creating figures for the manuscript
#    t0 = time.time()
    reg = 'pol'
    j0jL = coords['pol']
    measure = 'diversity'
    meas = translate_measures(measure)
    cutoff1 = 0.01

#    ttk_est, ttk, xxk, jjk, dtdx_t0 = ttest_region(meas, j0jL, cutoff1, method)
    CUT = EDI.window_cutoff(data, meas, region(j0jL), cutoff1, rf = rframe)
    ttk, xxk, jjk = CUT.realdata(Tmin, Tmax,  fcr = fcr, vload_min = vload_min,
                             dilutions_min = dilutions_min)
    ttk_est, dtdx = TI.TI_from_diversity(xxk, j0jL, cutoff1, rf = rframe)
    
    TT, DD = loadK31(measure, reg, cutoff1)
    TTest, dtdx_t0 = TI.TI_from_diversity(DD, j0jL, cutoff1, rf = rframe)
    print dtdx
    print dtdx_t0
    
    fig, ax = plt.subplots(1,2, figsize = (2*H, H))
    ax[0].plot(ttk, xxk, 's', label = 'HIVEVO')    
    ax[0].plot(TT, DD, 'o', label = 'K31')
    ax[0].plot(dtdx_t0[0]*np.sort(DD) + dtdx_t0[1], np.sort(DD), '--k')
    ax[1].plot(ttk, ttk_est, 's')
    ax[1].plot(TT, TTest, 'o')
    ax[1].plot(np.sort(TT), np.sort(TT), '--k')
    ax[1].set_xlabel('TI [years]', fontsize = fs)
    ax[1].set_ylabel('ETI [years]', fontsize = fs)
    ax[0].legend(fontsize = 0.8*fs, loc = 0)
    plt.savefig(outdir_name + 'K31_rf2.pdf')
    plt.close()
    