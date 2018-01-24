#! /usr/bin/env python
# -*- coding: utf-8 -*

##########################################################################
# @author: yann.leguen@cea.fr
# Copyright (C) CEA, 2017
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html for details.
##########################################################################

import os
import json
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from sklearn import mixture
import nibabel.gifti.giftiio as gio
from pylab import exp

# Define constants for the plots
LABEL_SIZE = 16
mpl.rcParams['xtick.labelsize'] = LABEL_SIZE
mpl.rcParams['ytick.labelsize'] = LABEL_SIZE 
TEXT_SIZE = 16

def gauss(x,mu,sigma,A):
    return A*exp(-(x-mu)**2/2/sigma**2)


def DPF_distrib_histogram(ax, x):
    """
    Parameters
    ax: contains the reference to the plot
    x: array containing all the pits DPF values in the considered areal
    """
    
    X  = x.reshape(-1, 1)
    lowest_bic = np.inf
    bic = []
    # Find the mixture with lowest BIC (Best information criterion)
    for n_components in range(1,3):
        # Fit a mixture of Gaussians with EM
        gmm = mixture.GMM(n_components=n_components)
        gmm.fit(X) # train it!
        bic.append(gmm.bic(X))
        if bic[-1] < lowest_bic:
            lowest_bic = bic[-1]
            best_gmm = gmm
            print ("Best gmm with components "+str(n_components)+" and BIC "+
                   str(bic[-1]))
    gmm = best_gmm

    COLOR, label, alpha, NORMED = 'blue', 'pits', 1, False
    number_bins = int(x.shape[0]/16)
    if number_bins <= 0:
        number_bins = 1
    print number_bins
    a,b = min(x), max(x)
    ax2 = ax.twinx()
    ax2.spines['right'].set_visible(False)
    ax2.spines['top'].set_visible(False)
    ax2.xaxis.set_ticks_position('none')
    ax2.yaxis.set_ticks_position('none')
    ax2.spines['left'].set_visible(False)
    for ylabel_i in ax2.get_yticklabels():
        ylabel_i.set_visible(False)
        ylabel_i.set_fontsize(0.0)

    n, bins, patches = ax.hist(X, number_bins, facecolor=COLOR, alpha=alpha,
                               range=(a,b), label=label+": "+str(x.shape[0]),
                               normed=NORMED)
    x_bins=(bins[1:]+bins[:-1])/2
    
    linspace = np.linspace(-2, 2, 1000).reshape(-1, 1)  
    ax2.plot(linspace, np.exp(gmm.score_samples(linspace)[0]), 'r')
    mu1 = gmm.means_[0]
    std1 = np.sqrt(gmm.covars_[0])

    if gmm.n_components == 1:
        threshold = mu1-2*std1
    elif gmm.n_components == 2:
        mu2 = gmm.means_[1]
        std2 = np.sqrt(gmm.covars_[1])
        A2 = gmm.weights_[1]
        A1 = gmm.weights_[0]
        x_samples = [mu2+(mu1-mu2)/250.0*p for p in range(250)]
        gauss1 = gauss(x_samples, mu1, std1, 1)
        gauss2 = gauss(x_samples, mu2, std2, 1)
        threshold = x_samples[np.argmin(np.abs(gauss1-gauss2))]

    if gmm.n_components == 1:
        ax.plot(np.repeat(threshold,200), np.linspace(min(n), max(n), num=200),
                color='magenta', lw=3, label = 'threshold')
    elif gmm.n_components == 2:
        ax.plot(np.repeat(threshold,200), np.linspace(min(n), max(n), num=200),
                color='green', lw=3, label = 'threshold')
    ax.set_xlim([-2,2])
    return threshold



if __name__ == '__main__':
    # In the following, lh refers to the symmetric template side

    ROOT_DIR = ""
    file_areals = ('clusters_total_average_pits_smoothed0.7_60_sym_lh_dist15.0'
                   '_area100.0_ridge2.0.gii')
    areals = os.path.join(ROOT_DIR, '2016_HCP_pits_cleaned',
                          'pits_density_sym_lh', file_areals)

    array_areals = gio.read(areals).darrays[0].data
    # Obtain a list of areal numbers
    areals_list =  np.unique(array_areals)
    # Exclude the first areal corresponding to corpus callosum and fornix
    areals_list = areals_list[1:]

    # input directory
    indir = os.path.join(ROOT_DIR, '2016_HCP_pits_cleaned',
                         'pits_analysis_lh', 'pits_DPF_distribution')

    for side in ['R', 'L']:
        INPUT = os.path.join(indir, side)
        count  = 0
        thresholds = []

        # Loop over all areals, to get the "best" threshold in each areal
        # We plot one histogram for each areal pits DPF distribution
        # We arbitrarily consider figures with grid (5, 6)
        for k, areal in enumerate(areals_list):
            if count%30 == 0:
                count = 0
                gs = gridspec.GridSpec(5, 6, width_ratios=[8,8,8,8,8,8],
                                       height_ratios=[6,6,6,6,6])
                fig = plt.figure()
                fig.subplots_adjust(hspace=0.2, wspace=0.2)
                ax = []
            ax.append(fig.add_subplot(gs[count]))
            file_DPF_distribution = os.path.join(INPUT,
                                                 'Areal_'+str(areal)+'.txt')
            X  = np.loadtxt(file_DPF_distribution)
            thresholds.append(DPF_distrib_histogram(ax[count], X))

            # If you know the areal names you can display them
            # if int(areal) in df_labels.index:
            #    ax[count].set_title(df_labels.loc[int(parcel)]['Name'],
            #    fontsize =TEXT_SIZE-4, fontweight = 'bold',
            #    verticalalignment="bottom")

            ax[count].locator_params(axis='x', nbins=4)
            ax[count].locator_params(axis='y', nbins=5)

            ax[count].spines['right'].set_visible(False)
            ax[count].spines['top'].set_visible(False)
            # Remove the ticks for both axis
            ax[count].xaxis.set_ticks_position('none')
            ax[count].yaxis.set_ticks_position('none')
            # Set the ticks on the left for y-axis
            ax[count].yaxis.set_ticks_position('left')

            # Only display ylabel for the plot on the left of the grid (6,5)
            if count%6==0:
                ax[count].set_ylabel('Number of pits', fontsize=TEXT_SIZE,
                                     fontweight = 'bold', labelpad=0)

            # Remove x_tick for plots which are not at the bottom
            if count < 24:
                ax[count].spines['bottom'].set_visible(False)
                for xlabel_i in ax[count].get_xticklabels():
                    xlabel_i.set_visible(False)
                    xlabel_i.set_fontsize(0.0)
            # Set the xlabel for plots at the bottom
            else:
                ax[count].xaxis.set_ticks_position('bottom')
                ax[count].set_xlabel('DPF', fontsize=TEXT_SIZE,
                                     fontweight='bold', labelpad=0)

            #fig.subplots_adjust(hspace=0.2,wspace=0.2)
            count +=1
        # Save the list of thresholds to file
        output_file = os.path.join(INPUT,'thresholds'+side+'.txt')
        np.savetxt(output_file, thresholds)
        
    # Close all the figures or show them with plt.show()
    plt.show()
    # plt.close('all')
