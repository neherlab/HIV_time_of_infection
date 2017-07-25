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
import pandas as pd
sns.set(style = 'darkgrid', font = u'Verdana')


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

cols = ['b', 'g','r','c', 'y', 'm']
marks = ['o', 's', '^', 'd', '*', 'p', 'v', '<', '>', 'h', '1', '2', '3', '4', '*', '+', 'x']
ms = 10
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

#def ttest_region(func_name, j0jL, cutoff, method, rf = rframe):
#    CUT = EDI.window_cutoff(data, func_name, region(j0jL), cutoff, rf = rf)
#    print CUT.ttk.shape
#    ttk, xxk, jjk = CUT.realdata(Tmin, Tmax,  fcr = fcr, vload_min = vload_min,
#                                 dilutions_min = dilutions_min)
#    ttk_est = np.zeros(ttk.shape)
#    dtdx_t0 = np.zeros((Npat, 2))
#    for jpat in xrange(Npat):
#        idx_pat = np.where(jjk == jpat)[0]
#        idx_data = np.where(jjk != jpat)[0]
#        ttk_data, dtdx_t0[jpat,:] = EDI.fitmeth_byname(ttk[idx_data], xxk[idx_data], method = method)
#        ttk_est[idx_pat] = dtdx_t0[jpat,0]*xxk[idx_pat] + dtdx_t0[jpat,1]
#    return ttk_est, ttk, xxk, jjk, dtdx_t0


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
                aft = P.get_allele_frequency_trajectories(reg, cov_min=500)[0]
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
    
    
def K31_diversity(data, cutoff, rf = rframe, fcr = 0.5, verbose = False):
    TT = []
    DD = []
    pats = []
    samples = []
    for j, key in enumerate(data.keys()):
        (Tloc, af) = data[key]
#        if np.sum(af.mask) > 0:
#            print key, np.mean(np.sum(af.mask, axis =0) >0)
        if np.mean(np.sum(af.mask, axis=0)>0) > fcr:
            if verbose:
                print 'Omitting {}'.format(key)
            continue
        af = af[:4,:]
        if rf is not None:
            af = af[:,rf::3]
        mask = (1.0 - af.max(axis=0))>cutoff
        Dloc = ((1.0 - (af**2).sum(axis=0))*mask).mean()
        TT.append(Tloc)
        DD.append(Dloc)
        j = key.rfind('_')
        pats.append(key[:j])
        samples.append(key[j+1:])
#        idx = np.where(np.max(af, axis = 0) <= 1. - cutoff)[0]
#        DD[j] = np.mean(1. - np.ma.sum(af[:,idx]**2, axis = 0))
    return np.array(TT)/365.25, np.array(DD), pats, samples

def prob_bins(bins, lamb):
    Pbins = .5*np.abs(np.exp(-np.abs(bins[:-1])/lamb) - np.exp(-np.abs(bins[1:])/lamb))
    j = np.where((bins[:-1] <0)*(bins[1:]>0))[0]
    if j:
        Pbins[j] = .5*(np.exp(bins[:-1][j]/lamb) + np.exp(-bins[1:][j]/lamb))
    return Pbins
    
         
if __name__=="__main__":
    '''Predicting infection dates for the 31 additional patients'''
    plt.ioff()
    plt.close('all')

    outdir_name = './plotK31/'
    if not os.path.exists(outdir_name):
        os.makedirs(outdir_name)

    #Creating figures for the manuscript
    measure = 'diversity'
    meas = translate_measures(measure)
    cutoff1 = 0.002
    bypairs = False
    pairs_legend = True
    label_type = 'templates' #'subtype' #'dilutions' #'templates' #'route' #'id'
#    reg = 'pol'
    for reg in ['gag', 'pol']:
        print reg
        j0jL = coords[reg]
    
        #Loading and processing the training set data (11 patients)
        CUT = EDI.window_cutoff(data, meas, region(j0jL), cutoff1, rf = rframe)
        ttk, xxk, jjk = CUT.realdata(Tmin, Tmax,  fcr = fcr, vload_min = vload_min,
                                 dilutions_min = dilutions_min)
        ttk_est, dtdx = TI.TI_from_diversity(xxk, j0jL, cutoff1, rf = rframe)
        
        
        #Loading and processing the validation dataset data (31 patients)
        K31data = loadK31(reg, './K31_data/{}/'.format(reg), fromHIV = False)
        TT, DD, pats, samples = K31_diversity(K31data, cutoff1, verbose = True)
        TTest, dtdx_t0 = TI.TI_from_diversity(DD, j0jL, cutoff1, rf = rframe)
        
        meta = pd.read_csv('./K31_data/patients.csv')
        meta_sample = pd.read_csv('./K31_data/sample_timeline_sequenced.csv')
        meta_table = pd.read_excel('./K31_data/Suppl table 31 patients true TI 170629.xlsx', skiprows = 2)
        
        TTmax = np.max(TT)
        jj = np.where(ttk <= TTmax)        
        if label_type in ['id', 'route', 'subtype']:     
            pats_unique = list(set(pats))
            f = (len(pats_unique)//len(cols)+1)
            cc = cols*f
            mm = [s for s in marks[:f] for j in xrange(len(cols))]
            pairs = [[j for j, p in enumerate(pats) if p ==pat] for pat in pats_unique]
        elif label_type in ['samples', 'templates', 'dilutions']:
            f = (len(pats)//len(cols)+1)
            cc = cols*f
            mm = [s for s in marks[:f] for j in xrange(len(cols))]
            

            
        
        
        fig, ax = plt.subplots(1,1, figsize = (1.3*H, 1.3*H))
        ax.plot(ttk, ttk_est, linestyle = '', marker = '^', markersize = 0.8*ms, markerfacecolor = 'gray', label = 'training data')
        if label_type in ['id', 'route', 'subtype']:
#            ax.plot(TT, TTest, 'ob', markersize = 12, label = 'validation data')
            for jp, p in enumerate(pairs):
                if label_type == 'id':
                    lab = pats_unique[jp]
                elif label_type == 'route':
                    lab = meta.at[np.where(meta['id'] == pats_unique[jp])[0][0], 'transmission route']
                elif label_type == 'subtype':
                    lab = meta_table.at[np.where(meta['id'] == pats_unique[jp])[0][0], 'HIV-1 subtyp']
                else:
                    lab = None
                    
                ax.plot(TT[p], TTest[p], ':' + mm[jp] + cc[jp], markersize = ms, label = lab)
            if pairs_legend:
                ax.legend(fontsize = 0.4*fs, loc = 0)
                ax.set_title(label_type, fontsize = fs)
        elif label_type in ['samples', 'templates', 'dilutions']:
            for j, sample in enumerate(samples):
                if label_type == 'samples':
                    lab = sample
                elif label_type == 'templates':
                    lab = meta_sample.at[np.where(meta_sample['id'] == sample)[0][0], 'templates approx']
                elif label_type == 'dilutions':
                    lab = meta_sample.at[np.where(meta_sample['id'] == sample)[0][0], 'dilutions']
                ax.plot(TT[j], TTest[j], mm[j] + cc[j], markersize = ms, label = lab)
#            ax.plot(TT, TTest, 'ob', markersize = 12, label = 'validation data')
            ax.legend(fontsize = 0.3*fs, loc = 0, ncol = 3)
        else:
            ax.plot(TT, TTest, 'ob', markersize = 12, label = 'validation data')
        ax.plot(np.sort(ttk), np.sort(ttk), linestyle = '--', color = 'gray')
        ax.plot(np.sort(ttk), 2.+ np.sort(ttk), linestyle = '--', color = 'gray')
        ax.plot(np.sort(ttk), -2.+ np.sort(ttk), linestyle = '--', color = 'gray')
        ax.set_xlabel('TI [years]', fontsize = fs)
        ax.set_ylabel('ETI [years]', fontsize = fs)
        ax.set_xlim([0, 10.])
        ax.set_ylim([0, 10.])
        ax.tick_params(labelsize = .8*fs)
        fig.tight_layout()
        plt.savefig(outdir_name + 'K31_{}_ETIvsTI_rf2_{}.pdf'.format(reg, label_type))
        plt.close()


#        templates = [meta_sample.at[np.where(meta_sample['id'] == sample)[0][0], 'templates approx'] for j, sample in enumerate(samples)]
#        templates = []
#        for j, sample in enumerate(samples):
#            temp = meta_sample.at[np.where(meta_sample['id'] == sample)[0][0], 'templates approx']
#            try:
#                templates.append(float(temp))
#            except:
#                templates.append(None)
#        
        fig, ax = plt.subplots(1,2, figsize = (2.6*H, 1.3*H))
        for j, sample in enumerate(samples):
            temp = meta_sample.at[np.where(meta_sample['id'] == sample)[0][0], 'templates approx']
            try:
                ax[0].semilogx(float(temp), TTest[j] - TT[j], mm[j] + cc[j], markersize = ms, label = sample)
                ax[1].semilogx(float(temp), np.abs(TTest[j] - TT[j]), mm[j] + cc[j], markersize = ms, label = sample)
            except:
                None
        ax[0].legend(fontsize = 0.3*fs, loc = 0)
        ax[0].set_ylabel('ETI - TI [years]', fontsize = fs)
        ax[1].set_ylabel('|ETI - TI| [years]', fontsize = fs)
        for j in xrange(2):        
            ax[j].set_xlabel('# templates', fontsize = fs)
            ax[j].tick_params(labelsize = .8*fs)
        fig.tight_layout()
        plt.savefig(outdir_name + 'K31_{}_error_vs_templates.pdf'.format(reg))
        plt.close()


        dTT = TTest - TT
        dTTmin = np.min(dTT)//1
        dTTmax = np.max(dTT)//1 + 1.
        TTbins = np.linspace(dTTmin, dTTmax, 2*int(dTTmax - dTTmin) +1)
        fig, ax = plt.subplots(1,1, figsize = (H, H))
#        ax.hist(ttk_est[jj]-ttk[jj], alpha = 0.5, color = 'gray', label = 'training data')
#        ax.hist(TTest - TT, alpha = 0.5, color = 'b', label = 'validation data')
        ax.hist(ttk_est-ttk, bins = TTbins, color = 'gray', label = 'training data')
        ax.hist(TTest - TT, bins = TTbins, alpha = 0.5, color = 'b', label = 'validation data')
        ax.legend(fontsize = 0.8*fs, loc = 0)
        ax.set_xlabel('ETI - TI [years]', fontsize = fs)
        ax.tick_params(labelsize = .8*fs)
        fig.tight_layout()
        plt.savefig(outdir_name + 'K31_{}_hist_rf2.pdf'.format(reg))
        plt.close()
        
        
        fig, ax = plt.subplots(1,1, figsize = (H, H))
        nn, bins = np.histogram(ttk_est-ttk, bins = TTbins)
        cbins = (bins[:-1] + bins[1:])/2.
        lamb = np.mean(np.abs(ttk_est-ttk))
        Pbins = prob_bins(bins, lamb)
        ax.plot((bins[:-1] + bins[1:])/2., nn, marker = 'd', color = 'gray', label = 'training data')
        ax.plot(cbins, np.sum(nn)*Pbins, linestyle = '--', color = 'gray')
#        ax.hist(ttk_est[jj]-ttk[jj], alpha = 0.5, color = 'gray', label = 'training data')
        nn, bins = np.histogram(TTest-TT, bins = TTbins)
        cbins = (bins[:-1] + bins[1:])/2.
        lamb = np.mean(np.abs(TTest-TT))
        Pbins = prob_bins(bins, lamb)
        ax.plot(cbins, nn, marker = 'o', color = 'b', label = 'validation data')
        ax.plot(cbins, np.sum(nn)*Pbins, '--b')
        ax.legend(fontsize = 0.8*fs, loc = 0)
        ax.set_xlabel('ETI - TI [years]', fontsize = fs)
        ax.tick_params(labelsize = .8*fs)
        fig.tight_layout()
        plt.savefig(outdir_name + 'K31_{}_hist_rf2_tmp.pdf'.format(reg))
        plt.close()
        
        
        fig, ax = plt.subplots(1,1, figsize = (2*H, 2*H))
        ax.plot(ttk, xxk, linestyle = '', marker = '^', markersize = 0.8*ms, markerfacecolor = 'gray', label = 'training data')
        if bypairs:
            for jp, p in enumerate(pairs):
                ax.plot(TT[p], DD[p], ':' + mm[jp] + cc[jp], markersize = 12, label = pats_unique[jp])
            if pairs_legend:
                ax.legend(fontsize = 0.4*fs, loc = 0)
        else:
            ax.plot(TT, DD, 'ob', markersize = 12, label = 'validation data')
            ax.legend(loc = 0, fontsize = 0.8*fs)
        ax.plot(dtdx_t0[0]*np.sort(xxk) + dtdx_t0[1], np.sort(xxk), linestyle = '--', color = 'gray')
        ax.set_xlabel('TI [years]', fontsize = fs)
        ax.set_ylabel('diversity', fontsize = fs)
        ax.tick_params(labelsize = .8*fs)
        fig.tight_layout()
        plt.savefig(outdir_name + 'K31_{}_diversity_bypat_rf2.pdf'.format(reg))
#        plt.savefig(outdir_name + 'K31_{}_diversity_bypat_rf2_N{}.pdf'.format(reg, len(pats_unique)))
        plt.close()
        
        
        # Outliers        
        fig, ax = plt.subplots(1,1, figsize = (H,H), sharey = True)
        dTT = np.abs(TTest - TT)
        dttk = np.abs(ttk_est-ttk)
        dTmax = np.max(np.concatenate((dTT, dttk)))//1 + 1
        Tbins = np.linspace(0., dTmax, 2*int(dTmax) + 1)
        nn, bins = np.histogram(dTT, Tbins)
        ax.plot(bins[1:], 1. - np.cumsum(nn)/TT.shape[0])
        nn, bins = np.histogram(dttk, Tbins)
        ax.plot(bins[1:], 1. - np.cumsum(nn)/ttk.shape[0])
               
        dTT = TTest - TT
        dTT = dTT[np.where(dTT >=0)]
        nn, bins = np.histogram(dTT, Tbins)
        ax.plot(bins[1:], (np.sum(nn) - np.cumsum(nn))/TT.shape[0], '--b')
        dttk = ttk_est-ttk
        dttk = dttk[np.where(dttk >=0)]
        nn, bins = np.histogram(dttk, Tbins)
        ax.plot(bins[1:],  (np.sum(nn) - np.cumsum(nn))/ttk.shape[0], '--g')
        
        dTT = TT - TTest
        dTT = dTT[np.where(dTT >0)]
        nn, bins = np.histogram(dTT, Tbins)
        ax.plot(bins[1:], (np.sum(nn) - np.cumsum(nn))/TT.shape[0], ':b')
        dttk = ttk-ttk_est
        dttk = dttk[np.where(dttk >0)]
        nn, bins = np.histogram(dttk, Tbins)
        ax.plot(bins[1:], (np.sum(nn) - np.cumsum(nn))/ttk.shape[0], ':g')
        
        ax.set_xlabel('|ETI - TI| [years]', fontsize = fs)
        ax.set_ylabel('fraction')
#        ax.set_ylim((0.,1.))
        plt.savefig(outdir_name + 'K31_{}_Nout.pdf'.format(reg))
        
        # diversity by codon position
        fig, ax = plt.subplots(1, 3, figsize = (3*H, 2*H), sharey = True)
        for j in xrange(3):
            CUT = EDI.window_cutoff(data, meas, region(j0jL), cutoff1, rf = j)
            ttk, xxk, jjk = CUT.realdata(Tmin, Tmax, fcr = fcr)
            ax[j].plot(ttk, xxk, marker = '^', linestyle = '', markerfacecolor = 'gray', markersize = 0.8*ms, label = 'training data')
                
            T, D, pats, samples = K31_diversity(K31data, cutoff1, rf = j)
            jj = np.argsort(T)
            ax[j].plot(T[jj], D[jj], 'ob', markersize = ms, label='validation data')
            ax[j].tick_params(labelsize = .8*fs1)
            ax[j].set_xlabel('TI [years]', fontsize = fs1)
            ax[j].set_title('codon pos {}'.format(j+1), fontsize = fs1)
        ax[0].set_ylabel('diversity', fontsize = fs1)
        ax[0].legend(fontsize = 0.8*fs1, loc = 0)
        fig.subplots_adjust(wspace = 0.1)
        plt.savefig(outdir_name + 'K31_{}_diversity.pdf'.format(reg))
        plt.close()
        

        
