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
import nibabel.gifti.giftiio as gio
from multiprocessing import cpu_count
from multiprocessing import Pool

def pits_DPF_distribution(parameters):

    """pits_DPF_distribution(database, array_areals, s_ids,
    areal, side,  sd_template="lh",
    outdir='/tmp/'):
    """
    database, array_areals, s_ids, areal, side,  sd_template, outdir = parameters
    """ 
    Parameters

    database: path of Morphologist database
    array_areals: array containing the areals information for each vertex
    s_ids: list of s_ids to consider
    areal: considered areal number in the areal file
    side: either R or L
    sd_template: side of the template
    outdir: output directory
    """
    if not (side == 'R' or side == 'L'):
        raise ValueError("argument side must be either 'R' or 'L'")
        
    OUTPUT = os.path.join(outdir, side)
    if not os.path.isdir(OUTPUT):
        os.makedirs(OUTPUT)

    # Array containing the pits DPF values in the areal for all subjects
    X = np.array([])    
    for s_id in s_ids:
        file_pits = os.path.join(database, s_id, "t1mri", "BL",
                                 "default_analysis", "segmentation", "mesh",
                                 "surface_analysis_"+sd_template,
                                 s_id+"_"+side+"white_pits_on_atlas.gii")
        file_DPF = os.path.join(database, s_id, "t1mri", "BL",
                                "default_analysis", "segmentation", "mesh",
                                "surface_analysis_"+sd_template,
                                s_id+"_"+side+"white_DPF_on_atlas.gii")
        #print file_DPF
        if os.path.isfile(file_pits) and os.path.isfile(file_DPF):
            array_pits = gio.read(file_pits).darrays[0].data
            array_DPF = gio.read(file_DPF).darrays[0].data
            # Find the index of pits
            index_pits = np.nonzero(array_pits)[0]
            # Identify their corresponding areals
            areals = array_areals[index_pits]
            # Find index of pits in the considered areal
            ind = np.where(areals ==  areal)[0]
            # test if the subject has pit in this parcel
            if ind.size:
                X = np.concatenate((X,array_DPF[index_pits[ind]]))
        else:
            print str(s_id)+" is missing either pits or DPF files"
            
    # Save the array X to output
    output = os.path.join(OUTPUT, 'Areal_'+str(areal)+'.txt')
    print "saving "+output
    np.savetxt(output, X)

    
if __name__ == '__main__':

    ROOT_DIR = ""
    # In the following, lh refers to the symmetric template side
    file_areals = ('clusters_total_average_pits_smoothed0.7_60_sym_lh_dist15.0'
                   '_area100.0_ridge2.0.gii')
    areals = os.path.join(ROOT_DIR, '2016_HCP_pits_cleaned',
                          'pits_density_sym_lh', file_areals)

    array_areals = gio.read(areals).darrays[0].data
    # Obtain a list of areal numbers
    areals_list =  np.unique(array_areals)
    # Exclude the first areal corresponding to corpus callosum and fornix
    areals_list = areals_list[1:]


    # Path to the Morphologist database
    database = ''
    s_ids = os.listdir(database)
    s_ids = sorted(s_ids)
        
    # output directory
    outdir = os.path.join(ROOT_DIR, '2016_HCP_pits_cleaned',
                          'pits_analysis_lh', 'pits_DPF_distribution')

    parameters = []
    # Compute the pits DPF distribution for each side and areal
    # Note: It's possible to run the code in parallel
    for side in ['R', 'L']:
        for areal in areals_list:
            parameters.append([database, array_areals, s_ids,
                               areal, side, 'lh', outdir])
    number_CPU = cpu_count()-1
    pool = Pool(processes = number_CPU)
    pool.map(pits_DPF_distribution, parameters)
    pool.close()
    pool.join()
