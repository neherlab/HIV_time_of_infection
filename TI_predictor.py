# -*- coding: utf-8 -*-
"""
Created on Tue May  9 11:31:53 2017

@author: puller
"""
from __future__ import division
import numpy as np
import EDI_functions as EDI

#default arguments
measure = 'diversity'
names = ['polymorphic', 'diversity', 'entropy']
funcnames = ['ambiguous_above', 'hamming_above', 'entropy_above']
func_name = funcnames[names.index(measure)]

fcr = 0.5
Tmin = 0; Tmax = 9
vload_min = None
dilutions_min = None
method = 'LAD'
rframe = 2 #reference frame; set to None to use all sites

fs = 28
H = 8


#loading frequency data
datapath = './Frequency_Data/'
data = EDI.load_patient_data(patient_names = 'all', filepath = datapath)
Npat = len(data['pat_names'])


def region(j0jL):
    if type(j0jL) is str:
        # The genome annotations
        head = ['name', 'x1', 'x2', 'width', 'ri']
        annot = []
        with open(datapath + 'annotations.txt', 'r') as fhandle:
            for line in fhandle:
                l = [x if j ==0 else int(x) for j, x in enumerate(line.split())]
                annot.append({name: l[j] for j, name in enumerate(head)})
        coords = {anno['name']: (anno['x1'], anno['x2']) for anno in annot}
        return coords[j0jL]
    else:
        return j0jL

def TI_from_diversity(DD, j0jL, cutoff, nboot = None, rf = rframe):
    '''
    Estimate the time of infection (TI) from the specified diversity values

    Input arguments:
    DD: list/array of diversity values
    j0jL: tuple specifying the genetic region to use
    cutoff: lower cutoff value, xc
    nboot: number of bootstraps over different patients (if None, then no bootstrapping)

    Output arguments:
    TTest: estimated times of infection (with rows corresponding to bootstrap relizations)
    dtdx_t0: slope and intercept values (with rows corresponding to bootstrap relizations)
    '''
    
    CUT = EDI.window_cutoff(data, func_name, region(j0jL), cutoff, rf = rf)
    ttk, xxk, jjk = CUT.realdata(Tmin, Tmax,  fcr = fcr, vload_min = vload_min,
                                 dilutions_min = dilutions_min)
    if nboot is None:
        ttk_data, dtdx_t0 = EDI.fitmeth_byname(ttk, xxk, method = method)
        TTest = dtdx_t0[0]*DD + dtdx_t0[1]
        return TTest, dtdx_t0
    else:
        Npat = len(CUT.pat_names)
        jjboot = np.random.randint(0, high = Npat, size = (nboot, Npat))
        TTest = np.zeros((nboot, len(DD)))
        dtdx_t0 = np.zeros((nboot, 2))
        for jboot, idx_boot in enumerate(jjboot):
            tk = np.ma.concatenate([ttk[np.where(jjk == j)] for j in idx_boot])
            xk = np.ma.concatenate([xxk[np.where(jjk == j)] for j in idx_boot])
            ttk_est, dtdx_t0[jboot,:] = EDI.fitmeth_byname(tk, xk, method = method)
            TTest[jboot,:] = dtdx_t0[jboot, 0]*DD + dtdx_t0[jboot, 1]
        return TTest, dtdx_t0


def TI_bootstrap_plot(DD, j0jL, cutoff, filename, nboot = 10**3):
    '''
    Plot bootstrap histograms for inferred times of infection
    
    Input arguments:
    DD: list/array of diversity values
    j0jL: tuple specifying the genetic region to use
    cutoff: lower cutoff value, xc
    filename: path to file for saving figure
    nboot: number of bootstraps over different patients
    '''
    import matplotlib.pyplot as plt
    plt.ioff()
#    plt.close('all')
    TTboot, dtdx_boot = TI_from_diversity(DD, j0jL, cutoff, nboot = nboot)
    fig, ax = plt.subplots(1, 1, figsize = (2*H, H))
    for jD, D in enumerate(DD):
        ax.hist(TTboot[:,jD])
    ax.set_xlabel('ETI [years]', fontsize = fs)
    ax.tick_params(labelsize = .8*fs)
    plt.axis('tight')
    plt.savefig(filename)
    plt.close()
    return None
    