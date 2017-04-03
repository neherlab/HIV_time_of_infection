# -*- coding: utf-8 -*-
"""
Created on Tue Sep 27 15:50:00 2016

@author: vpuller
"""
from __future__ import division
import numpy as np
import sys

# constants
h = 10**(-8)
err = .002

def load_patient_data(patient_names = 'all', 
                      q = 4, 
                      timescale = 'years', 
                      filepath = None, 
                      fromHIVEVO = False):
    if patient_names == 'all':
        patient_names = ['p{}'.format(j+1) for j in xrange(11)]
    if fromHIVEVO:
        #sys.path.append('/ebio/ag-neher/share/users/vpuller/HIVEVO/HIVEVO_access') 
        sys.path.append('/home/vadim/ebio/users/vpuller/HIVEVO/HIVEVO_access') 
        from hivevo.patients import Patient
        from hivevo.HIVreference import HIVreference
        
        ref = HIVreference(load_alignment=False)
        Lref = len(ref.seq)
        data_all = {}
        for pat_name in patient_names:
            PAT = Patient.load(pat_name)
            tt = PAT.times(unit = timescale)
            vload = PAT.n_templates_viral_load
            dilutions = PAT.n_templates_dilutions
            freqs_raw = PAT.get_allele_frequency_trajectories('genomewide', error_rate=err)[:,:q,:]
            map_to_ref = PAT.map_to_external_reference('genomewide')
            freqs = np.ma.zeros((tt.shape[0], q, Lref)); freqs.mask = True
            freqs[:,:,map_to_ref[:,0]] = freqs_raw[:,:,map_to_ref[:,1]]
            data_all[pat_name] = (tt, freqs, vload, dilutions)
            if filepath is not None:
                np.save(filepath + '{}_data.npy'.format(pat_name), freqs.data)
                np.save(filepath + '{}_mask.npy'.format(pat_name), freqs.mask)
                np.save(filepath + '{}_tt.npy'.format(pat_name), tt) 
                np.save(filepath + '{}_viral_load.npy'.format(pat_name), vload)
                np.save(filepath + '{}_dilutions.npy'.format(pat_name), dilutions.data)
                np.save(filepath + '{}_dilutions_mask.npy'.format(pat_name), dilutions.mask)
        data_all['Lref'] = freqs.shape[2]
        data_all['pat_names'] = patient_names
        
    elif filepath is not None:
        data_all = {}
        for pat_name in patient_names:
            tt = np.load(filepath + '{}_tt.npy'.format(pat_name))
            data = np.load(filepath + '{}_data.npy'.format(pat_name))
            mask = np.load(filepath + '{}_mask.npy'.format(pat_name))
            freqs = np.ma.masked_array(data, mask = mask)
            vload = np.load(filepath + '{}_viral_load.npy'.format(pat_name))
            dilutions = np.load(filepath + '{}_dilutions.npy'.format(pat_name))
            dilutions_mask = np.load(filepath + '{}_dilutions_mask.npy'.format(pat_name))
            dilutions = np.ma.masked_array(dilutions, mask = dilutions_mask)
            data_all[pat_name] = (tt, freqs, vload, dilutions)
        data_all['Lref'] = freqs.shape[2]
        data_all['pat_names'] = patient_names
    else:
        print 'Path to data is not specified'
    return data_all
        
# Diversity measures
def func_byname(funcname):
    funcs = [Ambiguous, Hamming_cut, Entropy_cut]
    heads = ['ambiguous_above', 'hamming_above', 'entropy_above']
    return funcs[heads.index(funcname)]
    

def Entropy_cut(freqs, alpha, q = 4):
    '''
    Calculating entropy with minor allele frequency greater than alpha
    '''
    H = np.ma.zeros(freqs.shape[::2])
    for jt, freq in enumerate(freqs):
        idx = np.where(np.max(freq, axis =0) <= 1. - alpha)[0]
        if len(idx) > 0:
            H[jt,idx] = -np.ma.sum(freq[:,idx]*np.log(freq[:,idx] +h)/np.log(q), axis = 0)            
    return H

def Hamming_cut(freqs, alpha):
    '''
    Calculating Hamming distance with minor allele frequency greater than alpha
    '''
    HM = np.ma.zeros(freqs.shape[::2])
    for jt, freq in enumerate(freqs):
        idx = np.where(np.max(freq, axis = 0) <= 1. - alpha)[0]
        HM[jt, idx] = 1. - np.ma.sum(freq[:,idx]**2, axis = 0)
    return HM

def Ambiguous(freqs, alpha):
    '''
    Calculating the number of sites with minor allele frequency greater than alpha
    '''
    return (np.max(freqs, axis = 1) <= 1. - alpha)
    

# data pre-processing for analysis in a single window
class window_cutoff(object):
    '''
    Evaluating statistic of interest (Hamming distance, entropy, etc.)
    within a genome window with a given cutoff
    '''
    def __init__(self, data, func_name, j0jL, cut, rf = None):
        '''
        data - patient data
        func - function defining statistic of interest (number of ambiguous sites, entropy, etc.)
        j0jL - = (j0, jL) - the nucleotide window to use in averaging
        cut - cutoff
        rf - reference frame, if considering only one position per codon
        '''
        (j0, jL) = j0jL
        func = func_byname(func_name)
        self.pat_names = data['pat_names']
        if rf is not None:
            self.ws = (jL - j0 - rf)//3
        else:
            self.ws = jL - j0 
        
        self.tt_all = {}
        self.vload_all = {}
        self.dilutions_all = {}
        self.xka_all = {}
        self.yka_all = {}
        self.yka_var_all = {}
        self.Nka_all = {}
        for jpat, pat_name in enumerate(self.pat_names):
            self.tt_all[pat_name] = data[pat_name][0]
            self.vload_all[pat_name] = data[pat_name][2]
            self.dilutions_all[pat_name] = data[pat_name][3]
            if rf is None:
                freqs = data[pat_name][1][:,:,j0:jL]
            else:
                freqs = data[pat_name][1][:,:,j0+rf:jL:3]
            msk = freqs.sum(axis = 1) <= 1. - err
            freqs = freqs/np.ma.sum(freqs, axis=1)[:,np.newaxis,:]
            xka = func(freqs, cut)
            xka = np.ma.masked_where(msk, xka)
            if np.sum(xka.mask):
                Nka = xka.shape[1] - np.sum(xka.mask, axis = 1)
            else:
                Nka = xka.shape[1]*np.ones(xka.shape[0])
            yka = np.ma.mean(xka, axis = 1)
            yka_var = np.ma.sum((xka - yka[:,np.newaxis])**2, axis = 1)/(Nka - 1)
            
            self.yka_all[pat_name] = yka
            self.yka_var_all[pat_name] = yka_var
            self.xka_all[pat_name] = xka
            self.Nka_all[pat_name] = Nka
        self.assemble_pats()

    def assemble_pats(self):
        '''
        Collecting patients into a single dataset
        '''
        self.ttk = np.ma.concatenate([self.tt_all[p] for p in self.pat_names])
        self.vload = np.ma.concatenate([self.vload_all[p] for p in self.pat_names])
        self.dilutions = np.ma.concatenate([self.dilutions_all[p] for p in self.pat_names])
        self.xxk = np.ma.concatenate([self.yka_all[p] for p in self.pat_names])
        self.vvark = np.ma.concatenate([self.yka_var_all[p] for p in self.pat_names])
        self.NNk = np.ma.concatenate([self.Nka_all[p] for p in self.pat_names])
        self.jjk = np.ma.concatenate([jp*np.ones(self.tt_all[p].shape[0], dtype = 'int') 
        for jp, p in enumerate(self.pat_names)])
        return None
    
    def realdata(self, Tmin, Tmax, fcr = 0., vload_min = None, dilutions_min = None):
        if dilutions_min is not None:
            idx = np.where((1-self.xxk.mask)*(self.ttk > Tmin)*(self.ttk < Tmax)*\
            (self.NNk/self.ws >= fcr)*(self.dilutions > dilutions_min))[0]
        elif vload_min is not None:
            idx = np.where((1-self.xxk.mask)*(self.ttk > Tmin)*(self.ttk < Tmax)*\
            (self.NNk/self.ws >= fcr)*(self.vload > vload_min))[0]
        else:
            idx = np.where((1-self.xxk.mask)*(self.ttk > Tmin)*(self.ttk < Tmax)*\
            (self.NNk/self.ws >= fcr))[0]
        return self.ttk[idx], self.xxk[idx], self.jjk[idx]


class sliding_window(object):
    '''
    Infection date inference within a genome window 
    '''
    def __init__(self, data, func_name, ws, cut):
        '''
        pat_class - initialized class with patient data
        func - function defining statistic of interest (number of ambiguous sites, entropy, etc.)
        ws - window size (# nucleotides)
        args - additional arguments of func
        '''
        func = func_byname(func_name)
        self.pat_names = data['pat_names']
        Npat = len(self.pat_names)
        Lref = data['Lref']
        self.tt_all = {}
        self.xka_all = {}
        self.yka_all = {}
        self.yka_var_all = {}
        self.Nka_all = {}
        self.ws = ws
        self.slopes = np.ma.zeros((Npat, Lref))
        self.tt0 = np.ma.zeros((Npat, Lref))
        self.var = np.ma.zeros((Npat, Lref))
        for jpat, pat_name in enumerate(self.pat_names):
            self.tt_all[pat_name] = data[pat_name][0]
            freqs = data[pat_name][1]
            msk = freqs.sum(axis = 1) <= 1.-err
            freqs = freqs/np.ma.sum(freqs, axis=1)[:,np.newaxis,:]
            xka = func(freqs, cut)
            xka = np.ma.masked_where(msk, xka)
            yka = np.ma.zeros((xka.shape[0], Lref - ws + 1))
            yka_var = np.ma.zeros((xka.shape[0], Lref - ws + 1))
            Nka = np.ma.zeros((xka.shape[0], Lref - ws + 1))
            for jt, xk in enumerate(xka):
                num = np.convolve(np.ones(ws, dtype=float), xk.filled(0.), mode='valid')
                if np.sum(xk.mask):
                    den = np.convolve(np.ones(ws, dtype=float), (1. - xk.mask), mode='valid')
                else:
                    den = np.convolve(np.ones(ws, dtype=float), np.ones(xk.shape), mode='valid')
                yk = num/np.ma.array(den)
                
                num_var = np.convolve(np.ones(ws, dtype=float), xk.filled(0.)**2, mode='valid')
                vark = (num_var - den*yk**2)/np.ma.array(den - 1)
                
                yka[jt,:] = yk #np.ma.masked_invalid(yk)
                yka_var[jt,:] = vark
                Nka[jt,:] = den
            self.yka_all[pat_name] = yka
            self.yka_var_all[pat_name] = yka_var
            self.xka_all[pat_name] = xka
            self.Nka_all[pat_name] = Nka
        self.assemble_pats()
        
    def assemble_pats(self):
        '''
        Collecting patients into a single dataset
        '''
        self.ttk = np.ma.concatenate([self.tt_all[p] for p in self.pat_names])
        self.xxk = np.ma.concatenate([self.yka_all[p] for p in self.pat_names], axis = 0)
        self.vvark = np.ma.concatenate([self.yka_var_all[p] for p in self.pat_names])
        self.NNk = np.ma.concatenate([self.Nka_all[p] for p in self.pat_names])
        self.jjk = np.ma.concatenate([jp*np.ones(self.tt_all[p].shape[0], dtype = 'int') 
        for jp, p in enumerate(self.pat_names)])
        return None
    

#fitting procedures
def fitmeth_byname(ttk, xxk, jjk = None, method = 'LAD', dtdx = True):
    '''
    TODO
    '''
    if method == 'LAD':
        ttk_est, ab = EDI_LAD(ttk, xxk)
    elif method == 'LAD_slope':
        ttk_est, ab = EDI_LAD(ttk, xxk, onlyslope = True)
    else:
        print 'WARNING from EDI_functions: unknown regression method'
    return ttk_est, ab

def EDI_LAD(ttk, xxk, jjk = None, bypat = False, onlyslope = False):
    '''
    Least absolute deviations (LAD) fitting by checking all point pairs
    
    Input arguments: 
    ttk - time points
    xxk - corresponding data points
    jjk - patient numbers (necessary, if bypat == True)
    
    Output arguments:
    inferred times for every time point, pair slope-intercept
    '''
    if onlyslope:
        slope = np.array([ti/xxk[i] for i, ti in enumerate(ttk) if xxk[i] != 0])
        dtdx_t0 = np.zeros((slope.shape[0],2))
        dtdx_t0[:,0] = slope
        error = np.array([np.sum(np.abs(ttk - beta*xxk - t0)) for beta, t0 in dtdx_t0])
    elif bypat:
        dtdx_t0 = []
        for jpat in xrange(np.max(jjk) + 1):
            idx = np.where(jjk == jpat)[0]
            dtdx_t0.extend([[(ttk[i] - ttk[j])/(xxk[i] - xxk[j]),
                             (xxk[i]*ttk[j] - xxk[j]*ttk[i])/(xxk[i] - xxk[j])] 
            for i in idx for j in idx if (i > j and xxk[i] != xxk[j])])
        dtdx_t0 = np.array(dtdx_t0)
        error = np.array([np.sum(np.abs(ttk - beta*xxk - t0)) for beta, t0 in dtdx_t0])
    else:
        N = ttk.shape[0]
        Mtt = ttk[:,np.newaxis] - ttk[np.newaxis,:]
        Mxx = xxk[:,np.newaxis] - xxk[np.newaxis,:]
        Mxt = xxk[:,np.newaxis]*ttk[np.newaxis,:] - xxk[np.newaxis,:]*ttk[:,np.newaxis]
        Mdtdx = np.ma.masked_where(Mxx == 0, Mtt/Mxx)
        Mt0 = np.ma.masked_where(Mxx == 0, Mxt/Mxx)
        dtdx_t0 = np.ma.zeros((N*(N-1)//2,2))
        dtdx_t0[:,0] = Mdtdx[np.triu_indices(N, 1)]
        dtdx_t0[:,1] = Mt0[np.triu_indices(N, 1)]            
        A = ttk[np.newaxis,:] - dtdx_t0[:,:1]*xxk[np.newaxis,:] - dtdx_t0[:,1:]
        error = np.sum(np.abs(A), axis = 1)
    if error.shape[0] > 0:
        tk_est = np.polyval(dtdx_t0[np.argmin(error),:], xxk)
        return tk_est, dtdx_t0[np.argmin(error),:]
    else:
        return np.array([None]*xxk.shape[0]), np.array([None, None])

    
def EDI_LAD_multisite(ttk, xxk):
    TTk = np.tile(ttk,(xxk.shape[1],1)).T
    def error(dtdx_t0):
        return np.ma.sum(np.ma.abs(dtdx_t0[0,:]*xxk + dtdx_t0[1,:] - TTk), axis = 0)
    (Nk, Na) = xxk.shape
    slopes = np.ma.zeros((Nk*(Nk-1)//2, 2, Na))
    idx = np.tril_indices(Nk, -1)
    tij = ttk[:, np.newaxis, np.newaxis] - ttk[np.newaxis,:, np.newaxis]
    xxk_ij = xxk[:,np.newaxis,:] - xxk[np.newaxis,:,:]
    dtdx = tij/xxk_ij
    slopes[:,0,:] = dtdx[idx[0], idx[1],:]
    xt_ij = xxk[:,np.newaxis,:]*ttk[np.newaxis,:, np.newaxis] - xxk[np.newaxis,:,:]*ttk[:, np.newaxis, np.newaxis]
    t0 = xt_ij/xxk_ij    
    slopes[:,1,:] = t0[idx[0], idx[1],:]

    err = np.ma.array([error(dtdx_t0) for dtdx_t0 in slopes])

    jjmin = np.argmin(err, axis = 0)
    dtdx_t0 = slopes[jjmin,:,range(jjmin.shape[0])].T
    tk_est = dtdx_t0[0,:]*xxk + dtdx_t0[1,:]
    return tk_est, dtdx_t0
