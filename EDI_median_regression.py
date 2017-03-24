# -*- coding: utf-8 -*-
"""
Created on Fri Oct  7 11:34:24 2016

@author: vpuller
"""

from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import os, time
import EDI_plotting as EDI_plot

if __name__=="__main__":
    '''Infection date cutoff dependence analysis'''
    plt.ioff()
    plt.close('all')

    outdir_name = '/home/vadim/ebio_vadim/HIV_time_of_infection/article_figures/'
    if not os.path.exists(outdir_name):
        os.makedirs(outdir_name)

    #Creating figures for the manuscript
    t0 = time.time()
    j0jL = 'pol'
    measure = 'hamming_above' 
    cutoff1 = 0.
    filehead = outdir_name + 'pol_' + measure + '_LAD_cut05_'

    #Trajectories
    for gene in EDI_plot.feas:
        print gene, EDI_plot.coords[gene]
        EDI_plot.plot_traj_xt(gene, measure, cutoff1,\
        outdir_name + '{}_hamming_trajxt.pdf'.format(gene))
        

    #Sliding window: absolute error for different ws
    wws = [1000, 3000, 5000, 7000]
    EDI_plot.plot_sliding_ws(measure, cutoff1, wws,\
    outdir_name + '{}_LAD_slidingwindow_abserr.pdf'.format(measure[:-6]))
    t1 = time.time()
    print t1 - t0
    
    #Dependence on cutoff
    measures = ['ambiguous_above', 'hamming_above', 'entropy_above']
    cc = np.array([np.linspace(.01, .45), np.linspace(0., .45), np.linspace(0., .45)])    
    for gene in EDI_plot.feas:
        EDI_plot.plot_median_new(gene, measures,\
        cc, outdir_name + '{}_LAD_'.format(gene))

    # accuracy
    EDI_plot.plot_tEDI_vs_tDI_bypat(j0jL, measure, cutoff1, \
    outdir_name + 'pol_{}_LAD_bypat_'.format(measure[:-6]))
    

    #ROC curve
    EDI_plot.ROC_curves(measure, j0jL, [0., .25], .5, 
                        outdir_name + 'pol_ROC_cut.pdf')
    
               
    #Table of slopes and intercepts
    cc = np.linspace(.0, .45, 10)
    for meas in measures:
        EDI_plot.maketable_slopes(j0jL, meas, cc,\
        ['LAD', 'LAD_slope'], outdir_name + 'table_pol_{}_'.format(meas))
    

    #Additional supplementary figure
    #distribution of slopes and intercepts
    EDI_plot.plot_slope_bootstrap(j0jL, measure, cutoff1,\
    outdir_name + 'slope_bootstrap.pdf')

    # Pearson correlation vs. cutoff
    cc = np.linspace(.01, .4)
    EDI_plot.plot_corrcoeff0(j0jL, measures, cc,\
    outdir_name + 'pol_corrcoeff.pdf')
    t1 = time.time()
    print t1 - t0
