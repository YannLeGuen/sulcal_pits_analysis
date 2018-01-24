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
import pandas as pd
import nibabel.gifti.giftiio as gio

# Boolean decide wether or not threshold,
# before selecting the deepest pit in the areal
THRESHOLDING = False


def create_phenotype(database, s_ids, array_areals, areal_list, side,
                     indir, outdir, sd_template='lh'):
    """ 
    Parameters

    database: path of Morphologist database
    array_areals: array containing the areals information for each vertex
    areal_list: contain the areal numbers
    s_ids: list of s_ids to consider
    side: either R or L
    indir: directory containing the thresholds table
    outdir: output directory
    sd_template: side of the template
    """
    if not os.path.isdir(outdir):
        os.makedirs(outdir)

    outdir_cc = os.path.join(outdir, 'case_control')
    if not os.path.isdir(outdir_cc):
        os.makedirs(outdir_cc)

    INPUT = os.path.join(indir, side)
    if THRESHOLDING:
        thr_file = os.path.join(INPUT, 'Areal_'+str(areal)+'.txt')
        thresholds = np.loadtxt(thr_file)

    NB_AREALS = len(areal_list)
    # Matrix containing the retained pit (if any) DPF value for each subject 
    DATA_DPF = np.zeros((len(s_ids), NB_AREALS))*np.nan
    
    for j, s_id in enumerate(s_ids):
        print "Currently processing s_id "+str(s_id)
        file_pits = os.path.join(database, s_id, "t1mri", "BL",
                                 "default_analysis", "segmentation", "mesh",
                                 "surface_analysis_"+sd_template,
                                 s_id+"_"+side+"white_pits_on_atlas.gii")
        file_DPF = os.path.join(database, s_id, "t1mri", "BL",
                                "default_analysis", "segmentation", "mesh",
                                "surface_analysis_"+sd_template,
                                s_id+"_"+side+"white_DPF_on_atlas.gii")
        array_pits = gio.read(file_pits).darrays[0].data
        array_DPF = gio.read(file_DPF).darrays[0].data


        # Remove the pits with a DPF below threshold
        if THRESHOLDING:
            for k, areal in enumerate(areal_list):
                ind = np.where(array_areals == areal)
                array_pits[ind] = (array_pits[ind]
                                    *(array_DPF[ind]>thresholds[k]))

        if False:
            # Locate the remaining pits
            index_pits = np.nonzero(array_pits)[0]
            # And their corresponding areals
            areals = array_areals[index_pits]
            for k, areal in enumerate(areal_list):
                ind = np.where(areals == areal)[0]
                # If the subject has pit in this areal we consider the deepest one
                if ind.size:
                    index_max_DPF = np.argmax(array_DPF[index_pits[ind]])
                    # MULTIPLY BY 20 BECAUSE
                    # SOLAR REQUIRES THIS TO ALLOW LARGE ENOUGH STD FOR ITS MODEL
                    # 20 is chosen to have an std sufficient in each areal
                    # It's somewhat arbitrary but doesn't affect the heritability
                    DATA_DPF[j,k] = array_DPF[index_pits[ind]][index_max_DPF]
        elif True:
            for k, areal in enumerate(areal_list):
                ind = np.where(array_areals == areal)[0]
                if ind.size:
                    index_max_DPF = np.argmax(array_DPF[ind])
                else:
                    DATA_DPF[j,k] = array_DPF[ind][index_max_DPF]
    """
    Process the DATA_DPF matrix 3 steps:
    1st filter out areals with less than 50% subjects having a pit
    2nd create DPF quantitative phenotype file for each areal kept
    3rd create a case control phenotype stating if a subject has pit or not
    """
    # We do not consider subject with pit DPF = 0,
    # because the pit DPF must be > 0.
    # Else use find zeros of numpy and replace them with almost 0
    DATA_DPF = np.nan_to_num(DATA_DPF)
    index_columns_kept = []
    # Identify the columns with at least 50% subjects having a sulcal pit
    for j in range(DATA_DPF.shape[1]):
        if np.count_nonzero(DATA_DPF[:,j]) > DATA_DPF.shape[0]*0.5:
            index_columns_kept.append(j)

    # For the columns kept create a phenotype file containing subjects
    # with at least one pit
    for index in index_columns_kept:
        num = str(int(areal_list[index]))
        df3 = pd.DataFrame()
        df3['FID'] = np.asarray(s_ids)[np.nonzero(DATA_DPF[:,index])].tolist()
        df3['IID'] = df3['FID']
        df3['Areal_'+num] = DATA_DPF[:,index][np.nonzero(DATA_DPF[:,index])]
        output = os.path.join(outdir,'DPF_pit'+num+'side'+side+'.csv')
        df3.to_csv(output, sep=' ',  header=True, index=False)


    # Create a case control phenotype file, stating if a subject has a pit
    df = pd.DataFrame()
    for j in range(DATA_DPF.shape[1]):
        df['Areal_'+str(int(areal_list[j]))] = DATA_DPF[:,j]
    df[df != 0] = 2
    df[df == 0] = 1
    df['IID'] = np.asarray(s_ids)
    output_cc = os.path.join(outdir_cc, 'all_pits_side'+side+'.csv')
    df.to_csv(output_cc, sep= ',',  header=True, index=False)

            
if __name__ == '__main__':    
    # In the following, lh refers to the symmetric template side
    ROOT_DIR = ""
    file_areals = ('clusters_total_average_pits_smoothed0.7_60_sym_lh_dist15.0'
                   '_area100.0_ridge2.0.gii')
    areals = os.path.join(ROOT_DIR,  '2016_HCP_pits_cleaned',
                          'pits_density_sym_lh', file_areals)
    
    array_areals = gio.read(areals).darrays[0].data
    # Obtain a list of areal numbers
    areal_list =  np.unique(array_areals)
    # Exclude the first areal corresponding to corpus callosum and fornix
    areal_list = areal_list[1:]

    # Path to the Morphologist database
    database = ""
    s_ids = os.listdir(database)
    s_ids = sorted(s_ids)

    # Directory containing the thresholds if needed
    indir = os.path.join(ROOT_DIR,  '2016_HCP_pits_cleaned',
                          'pits_analysis_lh', 'pits_DPF_distribution')
    # Output directory for the phenotypes
    outdir = os.path.join(ROOT_DIR,  '2016_HCP_pits_cleaned',
                          'pits_analysis_lh', 'phenotype_pits_max_DPF')

    # Loop over side and create phenotypes files
    # Can be parallelized
    # Note to avoid excessive gio.read contrary to pits_DPF_distributrion.py
    # We do not pass the considered areal as parameter, it's questionable
    # This also avoids repeating the thresholding when necessary
    for side in ['R', 'L']:
        create_phenotype(database, s_ids, array_areals, areal_list, side,
                         indir, outdir, sd_template='lh')
