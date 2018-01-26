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
import argparse

from multiprocessing import Pool

def texture_to_atlas(parameters):
    
    ROOT_DIR, s_id = parameters
    subj_dir = os.path.join(ROOT_DIR, "UKB", s_id)
    features = ['white_DPF', 'white_pits_smoothed0.7_60',
                'white_pits', 'geodesic_depth']
    # Native space hemisphere side
    for sd in ['L', 'R']:
        # Template symmetric (fsaverage_sym) hemisphere side
	for sd2 in ['lh', 'rh']:
 	    for feat in features:
                # File texture (input)
                tex = os.path.join(subj_dir, "t1mri", "BL",
                                   "default_analysis", "segmentation",
                                   "mesh", "surface_analysis",
                                   s_id+ "_"+sd+feat+".gii")
                # File texture on atlas (output)
                tex_atlas = os.path.join(subj_dir, "t1mri", "BL",
                                         "default_analysis", "segmentation",
                                         "mesh",
                                         "surface_analysis_sym_"+sd2,
                                         s_id+ "_"+sd+feat+"_on_atlas.gii")

                # File white mesh atlas spherical (input)
                atlas_spher = os.path.join(ROOT_DIR,
                                           "folder_gii", "sym",
                                           sd2+"h.sphere.reg.gii")
                # File white mesh sphere (input)
                white_spher = os.path.join(subj_dir, "t1mri", "BL",
                                           "default_analysis",
                                           "segmentation",
                                           "mesh", "surface_analysis",
                                           (s_id+"."+sd+".sphere.reg."+
                                            sd2+".fsaverage_sym.surf.gii"))
                                                   
                cmd = " ".join(['python -m brainvisa.axon.runprocess',
                               'ProjectTextureOntoAtlas',
                               tex,
                               white_spher,
                               atlas_spher,
                               tex_atlas])
                                        
                print cmd
                os.system(cmd)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=__doc__)
    parser.add_argument('-a', '--array_id', type=int,
                        help="Array ID from PBS")
    parser.add_argument('-j', '--job', type=int,
                        help="Number of job per array id")
    parser.add_argument('--ncore',
                        help='Nb of cpu cores to use (default %i)'
                        % 1, type=int)
    options = parser.parse_args()
    a = int(options.array_id)
    j = int(options.job)
    if options.ncore is not None:
        number_CPU = options.ncore
    else:
        number_CPU = 1
 
    ROOT_DIR = ""
    s_ids = sorted(os.listdir(ROOT_DIR+"/UKB"))
    parameters = []
    for s_id in s_ids[a*j:j+a*j]:
        parameters.append([ROOT_DIR, s_id])
    
    pool = Pool(processes = number_CPU)
    pool.map(texture_to_atlas, parameters)
    pool.close()
    pool.join()

