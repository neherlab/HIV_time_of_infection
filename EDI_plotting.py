# -*- coding: utf-8 -*-
"""
Created on Mon Mar 20 13:39:45 2017

@author: vadim
"""
from __future__ import division
import numpy as np
import scipy.stats.mstats as mstats
import matplotlib.pyplot as plt
import seaborn as sns
import EDI_functions as EDI
sns.set_style('darkgrid')


# constants
h = 10**(-8)
fcr = 0.5
Tmin = 0; Tmax = 9
method = 'LAD'

cols = sns.color_palette(n_colors=6)*2
#marks = ['o', 'v', '^', '<', '>', '1', '2', '3', '4', 's', 'p', '*', 'h', '+', 'x', 'd']
marks = ['o', 's', '^', 'd', 'x', 'p', 'v', '<', '>', '1', '2', '3', '4', '*', 'h', '+']
marks1 = ['o']*6 + ['^']*6
styles = ['-']*6 + ['--']*6
fs = 28
H = 8

# The genome annotations
datapath = '/home/vadim/ebio_vadim/HIV_time_of_infection/Frequency_Data/'  
head = ['name', 'x1', 'x2', 'width', 'ri']
annot = []
with open(datapath + 'annotations.txt', 'r') as fhandle:
    for line in fhandle:
        l = [x if j ==0 else int(x) for j, x in enumerate(line.split())]
        annot.append({name: l[j] for j, name in enumerate(head)})
coords = {anno['name']: (anno['x1'], anno['x2']) for anno in annot}
feas = ['gag', 'pol', 'env'] 


#loading frequency data
data = EDI.load_patient_data(filepath = datapath)
        
def leg_byname(funcname):
    legs = ['ambiguous sites', 'diversity', 'site entropy']
    heads = ['ambiguous_above', 'hamming_above', 'entropy_above']
    return legs[heads.index(funcname)]

def region(j0jL):
    if type(j0jL) is str:
        return coords[j0jL]
    else:
        j0jL
    
def plot_traj_xt(j0jL, measure, cutoff, filename):
    '''
    plot diversity time trajectories by codon position
    Input arguments
    j0jL: tuple of initial and final positions of the genome window
    measure: diversity measure
    cutoff: low frequency cutoff value
    filename: file path to save the figure
    '''
#    j0, jL = region(j0jL)
    def ax_traj_xt(ax, rf = None): 
        CUT = EDI.window_cutoff(data, measure, region(j0jL), cutoff, rf = rf)
        ttk, xxk, jjk, Npat = CUT.realdata(Tmin, Tmax, fcr)
        for jpat in xrange(Npat):
            jj = np.where(jjk == jpat)
            ax.plot(ttk[jj], xxk[jj], '--' + marks1[jpat], c=cols[jpat], markersize = 12)
        return ax
        
    fig, ax = plt.subplots(1, 3, figsize = (3*H, 2*H), sharey = True)
    for j in xrange(3):
        ax_traj_xt(ax[j], rf = j)
        ax[j].tick_params(labelsize = .8*fs)
        ax[j].set_xlabel('TI [years]', fontsize = fs)
        ax[j].set_title('codon pos {}'.format(j+1), fontsize = fs)
    ax[0].set_ylabel(leg_byname(measure), fontsize = fs)
    ax[0].legend(data['pat_names'], fontsize = 0.8*fs, loc = 0)
    fig.subplots_adjust(wspace = 0.1)
    plt.savefig(filename)
    plt.close()
    return None

def ttest_region(func_name, j0jL, cutoff, method,\
                 return_slope = False, return_all = False):
    CUT = EDI.window_cutoff(data, func_name, region(j0jL), cutoff)
    ttk, xxk, jjk, Npat = CUT.realdata(Tmin, Tmax, fcr)
    
    ttk_est = np.zeros(ttk.shape)
    dtdx_t0 = np.zeros((Npat, 2))
    for jpat in xrange(Npat):
        idx_pat = np.where(jjk == jpat)[0]
        idx_data = np.where(jjk != jpat)[0]
        ttk_data, dtdx_t0[jpat,:] = EDI.fitmeth_byname(ttk[idx_data], xxk[idx_data], method = method)
        ttk_est[idx_pat] = dtdx_t0[jpat,0]*xxk[idx_pat] + dtdx_t0[jpat,1]
    if return_all:
        return ttk_est, ttk, xxk, jjk, dtdx_t0
    elif return_slope:
        return ttk_est, ttk, dtdx_t0
    else:
        return ttk_est, ttk


def plot_median_new(j0jL, func_names, cutoffs, filehead):
    '''
    plotting absolute error as as a function of the low frequency cutoff
    for several diversity measures
    
    Input arguments:
    j0jL: tuple of initial and final positions of the genome window
    func_names: list of diversity measures
    cutoffs: array of the cutoff values
    filehead: common path for the output files
    '''
    err = []
    if cutoffs.ndim == 1:
        cutoffs = np.tile(cutoffs, (len(func_names),1))
    Ncut = cutoffs.shape[1]
    dtdx_t0 = np.zeros((len(data['pat_names']),2, len(func_names), Ncut))
    for jf, name in enumerate(func_names):
        ttk_abserr = np.zeros((2, Ncut))
        for jcut, cut in enumerate(cutoffs[jf,:]):
            ttk_est, ttk, dtdx_t0[:,:,jf,jcut] = ttest_region(name, j0jL,
                                        cut, method, return_slope = True) 
            dttk = ttk_est - ttk
            ttk_abserr[:, jcut] = np.array([np.median(dttk), np.mean(np.abs(dttk))])    
        err.append(ttk_abserr[1,:])

    fig, ax = plt.subplots(1, 1,figsize = (H, H))
    for jf, name in enumerate(func_names):
        ax.plot(cutoffs[jf,:], err[jf], '-', color = cols[jf])
    ax.set_ylabel('ETI - TI, mean abs. error, [years]', fontsize = fs)
    ax.legend([leg_byname(name) for name in func_names], loc = 0, fontsize = 0.8*fs)
    ax.set_xlabel(r'$x_{c}$', fontsize = fs)
    ax.tick_params(labelsize = .8*fs)
    ax.set_xticks(np.arange(0.,.5,.1))
    plt.savefig(filehead + 'cut_abserr.pdf')
    plt.close()
    
    fig, ax = plt.subplots(1, 1,figsize = (1.2*H, H))
    dtdx_med = np.median(dtdx_t0, axis = 0)
    for jf, name in enumerate(func_names):
        ax.plot(cutoffs[jf,:], dtdx_med[0,jf,:])
    ax.legend([leg_byname(name) for name in func_names], loc = 0, fontsize = 0.8*fs)
    ax.set_ylabel('slope [years/diversity]', fontsize = fs)
    ax.set_xlabel(r'$x_{c}$', fontsize = fs)
    ax.tick_params(labelsize = .8*fs)
    ax.set_xticks(np.arange(0.,.5,.1))
    plt.savefig(filehead + 'cut_s.pdf')
    plt.close()
    
    fig, ax = plt.subplots(1, 1,figsize = (1.2*H, H))
    dtdx_med = np.median(dtdx_t0, axis = 0)
    for jf, name in enumerate(func_names):
        ax.plot(cutoffs[jf,:], dtdx_med[1,jf,:])
    ax.legend([leg_byname(name) for name in func_names], loc = 0, fontsize = 0.8*fs)
    ax.set_ylabel('intercept [years]', fontsize = fs)
    ax.set_xlabel(r'$x_{c}$', fontsize = fs)
    ax.tick_params(labelsize = .8*fs)
    ax.set_xticks(np.arange(0.,.5,.1))
    plt.savefig(filehead + 'cut_t0.pdf')
    plt.close()
    return None

def plot_tEDI_vs_tDI_bypat(j0jL, func_name, cutoff, filehead):
    '''
    Plotting estimated versus actual date of infection
    
    Input arguments:
    j0jL: tuple of initial and final positions of the genome window
    func_name: diversity measure
    cutoff: low frequency cutoff value
    filename: file path to save the figures
    '''
    ttk_est, ttk, xxk, jjk, dtdx_t0 =\
    ttest_region(func_name, j0jL, cutoff, method, return_all = True) 
        
    # scatter-plot EDI vs DI
    fig, ax = plt.subplots(1, 1, figsize = (H, H)) 
    for j in xrange(np.max(jjk) + 1):
        jj = np.where(jjk == j)
        ax.scatter(ttk[jj], ttk_est[jj], color = cols[j], marker = marks1[j],\
        s = 40, label = data['pat_names'][j])
    ax.plot(np.sort(ttk), np.sort(ttk), '--k')
    
    ax = draw_ellipse(ax, xy = (5.3, 1.6), ab = (2.7,.7), psi = .25*np.pi,\
    lw = 1, c = 'r')
    ax = draw_ellipse(ax, xy = (1., 3.35), ab = (1.,.4), psi = .25*np.pi,\
    lw = 1, c = 'b')
    ax.set_xlabel('TI [years]', fontsize = fs)
    ax.set_ylabel('ETI [years]', fontsize = fs)
    ax.legend(fontsize = 0.6*fs, loc = 2, ncol = 2)
    ax.tick_params(labelsize = .8*fs)
    ax.axis('tight')
    plt.savefig(filehead + 'ETI_vs_TI.pdf')
    plt.close()
    
    fig, ax = plt.subplots(1, 1, figsize = (H, H)) 
    dttk = np.abs(ttk_est - ttk)            
    dtkave, dtkvar, tk = moving_average(ttk, dttk)
    ax.plot(tk, dtkave)
    ax.set_xlabel('TI [years]', fontsize = fs)
    ax.set_ylabel(r'$|$ETI - TI$|$, [years]', fontsize = fs)
    ax.tick_params(labelsize = .8*fs)
    ax.axis('tight')
    plt.savefig(filehead + 'error_vs_TI.pdf')
    plt.close()
    
    fig, ax = plt.subplots(1, 1, figsize = (H, H))
    ax.hist(ttk_est - ttk, alpha = 0.5)
    ax.set_xlabel(r'ETI - TI, [years]', fontsize = fs)
    ax.tick_params(labelsize = .8*fs)
    plt.savefig(filehead + 'hist.pdf')
    plt.close()
    return ttk_est, ttk

def draw_ellipse(ax, xy = (0.,0), ab = (1.,1.), psi = 0., n = 50, lw = 1., c = 'b'):
    pphi = np.linspace(0., 2*np.pi, num = n)
    xxyy = np.array([ab[0]*np.cos(pphi), ab[1]*np.sin(pphi)])
    Spsi = np.array([[np.cos(psi), -np.sin(psi)], [np.sin(psi), np.cos(psi)]])
    xxyy1 = Spsi.dot(xxyy)
    ax.plot(xy[0] + xxyy1[0,:],xy[1] + xxyy1[1,:], lw = lw, c = c)
    return ax

def moving_average(ttk, dttk, n = 25):
    idx_sort = np.argsort(ttk)
    tk = ttk[idx_sort].astype('int')
    dtk = dttk[idx_sort]
    dtk_ave = np.array([np.mean(dtk[j:j+n]) for j, t in enumerate(tk[:-n])])
    dtk_var = np.array([np.mean(dtk[j:j+n]**2) for j, t in enumerate(tk[:-n])])-\
    dtk_ave**2
    tk_ave = np.array([np.mean(tk[j:j+n]) for j, t in enumerate(tk[:-n])])
    return dtk_ave, dtk_var/n, tk_ave

def ROC_curves(func_name, j0jL, cutoff, Trecent, filename):
    '''recent infection
    
    Plotting receiver operating characteristic (ROC curve)
    Input arguments:
    j0jL: tuple of initial and final positions of the genome window
    func_name: diversity measure
    cutoff: low frequency cutoff value
    Trecent: time threshold definig 
    filename: file path to save the figure
    '''
    fig, ax = plt.subplots(1, 1, figsize = (H, H))
    legs = []
    for cut in cutoff:
        MM = ROC_curve(func_name, j0jL, cut, Trecent, ax)
        AUC = np.sum(np.diff(MM[:,1,0])*MM[1:,0,0])
        legs.append('xc = {}, AUC = {:.2g}'.format(cut, AUC))
    ax.set_xlabel('1-specificity', fontsize = fs)
    ax.set_ylabel('sensitivity', fontsize = fs)
    ax.tick_params(labelsize = 0.8*fs)
    ax.legend(legs, fontsize = 0.8*fs, loc = 0)
    plt.savefig(filename)
    plt.close()
    return None

def ROC_curve(func_name, j0jL, cutoff, tcr, ax = None):
    def contable(ttk, xxk, tcr, xcr):
        TP = np.count_nonzero((ttk < tcr)*(xxk < xcr))
        FP = np.count_nonzero((ttk >= tcr)*(xxk < xcr))
        FN = np.count_nonzero((ttk < tcr)*(xxk >= xcr))
        TN = np.count_nonzero((ttk >= tcr)*(xxk >= xcr))
        return np.array([[TP, FN], [FP, TN]])    
    CUT = EDI.window_cutoff(data, func_name, region(j0jL), cutoff)
    ttk, xxk, jjk, Npat = CUT.realdata(Tmin, Tmax, fcr)
    
    xxcr = np.sort(np.unique(xxk))
    M = np.array([contable(ttk, xxk, tcr, xcr) for xcr in xxcr])
    MM = M/np.sum(M, axis=2, keepdims = True)
    if ax is not None:
        ax.plot(MM[:,1,0], MM[:,0,0])
    return MM

def maketable_slopes(j0jL, func_name, cutoffs, methods, filehead):
    '''
    Save table of slopes, intercepts and errors for different values of 
    low frequency cutoff
    
    Input arguments:
    j0jL: tuple of initial and final positions of the genome window
    func_names: list of diversity measures
    cutoffs: array of the cutoff values
    methods: fitting methods
    filehead: common path for the output files
    '''
    err = []
    Ncut = cutoffs.shape[0]
    dtdx_t0 = np.zeros((len(data['pat_names']),2, len(methods), Ncut))
    for jmeth, meth in enumerate(methods):
        ttk_median = np.zeros((3, Ncut))
        ttk_abserr = np.zeros((2, Ncut))
        for jcut, cut in enumerate(cutoffs):
            ttk_est, ttk, dtdx_t0[:,:,jmeth,jcut] = ttest_region(func_name,\
            j0jL, cut, meth, return_slope = True) 
            dttk = ttk_est - ttk
            ttk_median[:,jcut] = np.array([np.percentile(dttk, 25),\
            np.percentile(dttk, 50), np.percentile(dttk, 75)])
            ttk_abserr[:, jcut] = np.array([np.median(dttk), np.mean(np.abs(dttk))])
            
        err.append(ttk_abserr[1,:])
    
    dtdx_med = np.median(dtdx_t0, axis = 0)    
            
    with open(filehead + 'cut_st0.txt', 'w') as filehandle:
        filehandle.write('\t' + '\t\t\t\t\t\t'.join(methods) + '\n')
        head = 'x_c'
        for jmeth, meth in enumerate(methods):
            if meth != 'LAD_slope':
                head += '\ts[years/D]\tt_0[years]\tMAE [years]'
            else:
                head += '\ts[years/D]\tMAE [years]'
        filehandle.write(head + '\n')
        for jcut, cut in enumerate(cutoffs):
            line = '{:.2f}'.format(cut)
            for jmeth, meth in enumerate(methods):
                if meth != 'LAD_slope':
                    line += ' &\t{:.2f} &\t{:.2f} &\t\t{:.2f}\t'.format(\
                    dtdx_med[0,jmeth,jcut], dtdx_med[1,jmeth,jcut], err[jmeth][jcut])
                else:
                    line += ' &\t{:.2f} &\t\t{:.2f}\t'.format(\
                    dtdx_med[0,jmeth,jcut], err[jmeth][jcut])
            filehandle.write(line + '\\\\\n')
    return None


def plot_slope_bootstrap(j0jL, func_name, cutoff, filename, nboot = 10**3):
    '''
    Bootstrap plot of slope and intercept values
    
    Input arguments:
    j0jL: tuple of initial and final positions of the genome window
    func_name: diversity measure
    cutoff: low frequency cutoff
    filename: path to the file for saving the figure
    nboot: number of bootstrap realizations
    '''
    CUT = EDI.window_cutoff(data, func_name, region(j0jL), cutoff)
    ttk, xxk, jjk, Npat = CUT.realdata(Tmin, Tmax, fcr)
    dtdx_t0 = np.zeros((nboot, 2))

    jjboot = np.random.randint(0, high = Npat, size = (nboot, Npat))
    for jboot, idx_boot in enumerate(jjboot):
        tk = np.ma.concatenate([ttk[np.where(jjk == j)] for j in idx_boot])
        xk = np.ma.concatenate([xxk[np.where(jjk == j)] for j in idx_boot])
        ttk_est, dtdx_t0[jboot,:] = EDI.fitmeth_byname(tk, xk, method = method)

    fig, ax = plt.subplots(1, 2, figsize = (2*H, H), sharey = True)
    ax[0].hist(dtdx_t0[:,0], alpha = 0.5)
    ax[1].hist(dtdx_t0[:,1], alpha = 0.5)
        
    ax[0].set_xlabel('slope [years]', fontsize = fs)
    ax[0].set_ylabel(method, fontsize = fs)
    ax[0].tick_params(labelsize = .8*fs)
    
    ax[1].set_xlabel('intercept', fontsize = fs)
    ax[1].tick_params(labelsize = .8*fs)
    plt.savefig(filename)
    plt.close()
    
    fig, ax = plt.subplots(1, 1, figsize = (1.2*H, H))
    Hist, xedges, yedges, cax = ax.hist2d(dtdx_t0[:,0], dtdx_t0[:,1],\
    cmap = plt.cm.Blues)
    ax.set_xlabel('slope [years]', fontsize = fs)
    ax.set_ylabel('intercept', fontsize = fs)
    ax.tick_params(labelsize = .8*fs)
    cbar = fig.colorbar(cax)
    cbar.ax.tick_params(labelsize = .8*fs)
    plt.savefig(filename[:-4] + '_2d.pdf')
    plt.close()
    return None

def plot_corrcoeff0(j0jL, measures, cutoffs, filename):
    '''
    Plot pearson correlation coefficients between times 
    and corresponding diversity values
    
    Input arguments:
    j0jL: tuple of initial and final positions of the genome window
    measures: diversity measures
    cutoffs: low frequency cutoffs
    filename: path to the file for saving the figure
    nboot: number of bootstrap realizations
    '''
    fig, ax = plt.subplots(1, len(measures), figsize = (H*len(measures), 2*H),\
    sharey = True)
    titls = [leg_byname(name) for name in measures]
        
    for j, measure in enumerate(measures):
        rxt = np.zeros((cutoffs.shape[0], len(data['pat_names'])))
        for jcut, cut in enumerate(cutoffs):
            CUT = EDI.window_cutoff(data, measure, region(j0jL), cut)
            ttk_all, xxk_all, jjk, Npat = CUT.realdata(Tmin, Tmax, fcr)
            for jpat in xrange(Npat):
                idx = np.where(jjk == jpat)
                ttk = ttk_all[idx]
                xxk = xxk_all[idx]
                rxt[jcut, jpat] = np.corrcoef(ttk, xxk)[0,1]

        for jr, r in enumerate(rxt.T):
            ax[j].plot(cutoffs, r, styles[jr])
        ax[j].set_title(titls[j], fontsize = fs)
        ax[j].tick_params(labelsize = .8*fs)
        ax[j].set_xlabel(r'$x_c$', fontsize = fs)
        ax[j].set_xticks(np.arange(0.,.5,.1))
    ax[0].legend(data['pat_names'], fontsize = 0.8*fs, loc = 0)
    ax[0].set_ylabel(r'$r_{xt}$', fontsize = fs)
    fig.subplots_adjust(hspace = 0.1)
    plt.savefig(filename)
    plt.close()
    return None
        
#Sliding window plot functions
def plot_sliding_ws(func_name, cutoff, wws, filename, lstep = 10):
    '''
    Plot absolute arror as a function of position in the genome 
    (also plots coverage of valid time points)
    
    func_name: diversity measure
    cutoff: low-frequency cutoff
    wws: genome widnow sizes to use
    filename: file path to save the main figure
    lstep: step along the genome (# nucleotides)
    '''
    dtdx_all = []
    cov_sites_all = []
    cov_points_all = []
    ll_all = []
    
    f = 4; hsp = 0.
    plt.figure(10, figsize = (2*H, f/(f-1)*H + hsp))
    ax0 = plt.subplot2grid((f, 1), (0,0), rowspan = f-1)
    ax1 = plt.subplot2grid((f, 1), (f-1,0), rowspan = 1, sharex = ax0)
    ax = [ax0, ax1]
    for jws, ws in enumerate(wws): 
        ttk_est, dtdx_t0, ttk, xxk, NNk, jjk =\
        ttest_sliding_window(func_name, cutoff, ws, lstep = lstep)
        
        L = ttk_est.shape[1]*lstep
        ll = range(ws//2, L + (ws+1)//2, lstep)
        tabserr = np.mean(np.abs(ttk_est - ttk[:, np.newaxis]), axis = 0)
        ax[0].plot(ll, tabserr, color = cols[jws], label = 'ws = {}'.format(ws))
        ax[0].set_ylabel('ETI - TI, mean abs. error [years]', fontsize = fs)
            
        dtdx_all.append(dtdx_t0)
        cov_sites_all.append(NNk.mean(axis = 0)/ws)
        cov_points_all.append(np.mean(1 - xxk.mask, axis = 0))
        ll_all.append(ll)

    handles, labels = ax[0].get_legend_handles_labels()
    ax[0].legend(fontsize = 0.8*fs, loc = 0)
    ax[0].tick_params(labelsize = 0.8*fs)
    
    draw_genome(annot, ax[1], fs = 18, pad = 1)
    ax[1].set_xlim((0, L + ws))
    ax[1].set_axis_off()
    
##    coord = np.array(anno.loc[anno['name'].isin(feas), ['x1', 'x2']])
#    coord = np.array([[anno['x1'], anno['x2']] for anno in annot 
#    if anno['name'] in feas])
    
    for jgene, gene in enumerate(feas):
        tk_est, tk = ttest_region(func_name, gene, cutoff, 'LAD')
        dtk = tk_est - tk
        ttk_range = np.mean(np.abs(dtk))
        (xL, xR) = coords[gene]
        ax[0].axhline(y = ttk_range, xmin = xL/(L + ws), xmax = xR/(L + ws),\
        color = 'k')
        ax[0].text((xL + xR)/2, ttk_range + .05,
                gene,
                color='k',
                fontsize = .8*fs,
                ha='center')
    
    plt.subplots_adjust(hspace = hsp)
    plt.savefig(filename)
    plt.close()
    
    #plotting coverage
    plt.figure(10, figsize = (2*H, f/(f-1)*H + hsp))
    ax0 = plt.subplot2grid((f, 1), (0,0), rowspan = f-1)
    ax1 = plt.subplot2grid((f, 1), (f-1,0), rowspan = 1, sharex = ax0)
    for jws, ws in enumerate(wws):
        ax0.plot(ll_all[jws], cov_sites_all[jws], label = 'ws = {}'.format(ws))
    ax0.set_ylabel('coverage, sites %', fontsize = fs)
    ax0.legend(fontsize = 0.8*fs, loc = 0)
    ax0.tick_params(labelsize = .8*fs)
    
    draw_genome(annot, ax1, fs = 18, pad = 1)
    ax1.set_xlim((0, L + ws))
    ax1.set_axis_off()
    plt.subplots_adjust(hspace = hsp)
    plt.savefig(filename[:-4] + '_cov_sites.pdf')
    plt.close(10)
    
    plt.figure(10, figsize = (2*H, f/(f-1)*H + hsp))
    ax0 = plt.subplot2grid((f, 1), (0,0), rowspan = f-1)
    ax1 = plt.subplot2grid((f, 1), (f-1,0), rowspan = 1, sharex = ax0)
    for jws, ws in enumerate(wws):
        ax0.plot(ll_all[jws], cov_points_all[jws], label = 'ws = {}'.format(ws))
    ax0.set_ylabel('coverage, sdatapoints %', fontsize = fs)
    ax0.legend(fontsize = 0.8*fs, loc = 0)
    ax0.tick_params(labelsize = .8*fs)
    
    draw_genome(annot, ax1, fs = 18, pad = 1)
    ax1.set_xlim((0, L + ws))
    ax1.set_axis_off()
    plt.subplots_adjust(hspace = hsp)
    plt.savefig(filename[:-4] + '_cov_points.pdf')
    plt.close(10)
        
    return None

def ttest_sliding_window(func_name, cutoff, ws, lstep = 10, Ncr = 5):
    SW= EDI.sliding_window(data, func_name, ws, cutoff)
    idx = np.where((SW.ttk > Tmin)*(SW.ttk < Tmax))[0]
    ttk = SW.ttk[idx]
    jjk = SW.jjk[idx]; Npat = np.max(SW.jjk) + 1
    xxk = SW.xxk[idx,:][:,::lstep]
    NNk = SW.NNk[idx,:][:,::lstep]
    xxk = np.ma.masked_where(NNk/ws < fcr, xxk)

    ttk_est = np.ma.zeros(xxk.shape)
    for jpat in xrange(Npat):
        idx_pat = np.where(jjk == jpat)[0]
        idx_data = np.where(jjk != jpat)[0]
        ttk_data, dtdx_t0 = EDI.EDI_LAD_multisite(ttk[idx_data], xxk[idx_data,:])
        ttk_est[idx_pat,:] = dtdx_t0[0,:]*xxk[idx_pat,:] + dtdx_t0[1,:]
    msk = np.zeros_like(xxk)
    msk[:,np.where(np.sum(1-xxk.mask, axis=0) < Ncr)[0]] = 1
    ttk_est = np.ma.masked_where(msk, ttk_est)
    return ttk_est, dtdx_t0, ttk, xxk, NNk, jjk

def ma_quantiles(ttk, prob = None):
    if prob is None:
        tquant = np.ma.zeros((3, ttk.shape[1]))
        tquant[1,:] = np.ma.median(ttk, axis = 0)
        ttk1 = np.ma.masked_where(ttk > tquant[1,:], ttk)
        tquant[0,:] = np.ma.median(ttk1, axis = 0)
        ttk1 = np.ma.masked_where(ttk < tquant[1,:], ttk)
        tquant[2,:] = np.ma.median(ttk1, axis = 0)
    else:
        tquant = mstats.mquantiles(ttk, prob, axis = 0)
    return tquant
    
def draw_genome(anno_elements,
                ax=None,
                rows=4,
                readingframe=True, fs=9,
                y1=0,
                height=1,
                pad=0.2):
    '''Draw genome boxes'''
    from matplotlib.patches import Rectangle
#    from Bio.SeqFeature import CompoundLocation
#    import pandas as pd

    if ax is None:
        fig, ax = plt.subplots(1, 1)

    ax.set_ylim([-pad,rows*(height+pad)])
#    anno_elements = []
#    for name, feature in annotations.iteritems():
#        if type(feature.location) is CompoundLocation:
#            locs = feature.location.parts
#        else:
#            locs = [feature.location]
#        for li,loc in enumerate(locs):
#            x = [loc.nofuzzy_start, loc.nofuzzy_end]
#            anno_elements.append({'name': name,
#                                  'x1': x[0],
#                                  'x2': x[1],
#                                  'width': x[1] - x[0]})
#            if name[0]=='V':
#                anno_elements[-1]['ri']=3
#            elif li:
#                anno_elements[-1]['ri']=(anno_elements[-2]['ri'] + ((x[0] - anno_elements[-2]['x2'])%3))%3
#            else:
#                anno_elements[-1]['ri']=x[0]%3

    anno_elements.sort(key = lambda x:x['x1'])
    for ai, anno in enumerate(anno_elements):
        if readingframe:
            anno['y1'] = y1 + (height + pad) * anno['ri']
        else:
            anno['y1'] = y1 + (height + pad) * (ai%rows)
        anno['y2'] = anno['y1'] + height
        anno['height'] = height

    for anno in anno_elements:
        r = Rectangle((anno['x1'], anno['y1']),
                      anno['width'],
                      anno['height'],
                      facecolor=[0.8] * 3,
                      edgecolor='k',
                      label=anno['name'])

        xt = anno['x1'] + 0.5 * anno['width']
        yt = anno['y1'] + 0.2 * height + height * (anno['width']<500)
        anno['x_text'] = xt
        anno['y_text'] = yt

        ax.add_patch(r)
        ax.text(xt, yt,
                anno['name'],
                color='k',
                fontsize=fs,
                ha='center')

    return None