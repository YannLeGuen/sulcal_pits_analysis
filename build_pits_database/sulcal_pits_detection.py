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


def pits_extraction(parameters):
    DATABASE, s_id = parameters
    subj_dir = os.path.join(DATABASE, s_id)

    # Default parameters
    DPF_alpha, thresh_ridge, thresh_dist, thresh_area = 0.03, 1.5, 20.0, 50.0
    # side, group average Fiedler length, group average surface area
    # Ideally these two group parameters should be recomputed on new cohorts
    param = {'R': ['right', 91950.37, 237.54],
             'L': ['left', 91312.20, 236.18]}
    
    for sd in ['L', 'R']:
        pole_cingular = os.path.join(subj_dir, "t1mri", "BL",
                                     "default_analysis", "segmentation",
                                     "mesh", "surface_analysis",
                                     s_id+ "_"+sd+"white_pole_cingular.gii")
        white_mesh = os.path.join(subj_dir, "t1mri", "BL",
                                  "default_analysis", "segmentation", "mesh",
                                  s_id+ "_"+sd+"white.gii")
        DPF = os.path.join(subj_dir, "t1mri", "BL",
                           "default_analysis", "segmentation", "mesh",
                           "surface_analysis", s_id+ "_"+sd+"white_DPF.gii")
        pits = os.path.join(subj_dir, "t1mri", "BL",
                            "default_analysis", "segmentation", "mesh",
                            "surface_analysis", s_id+ "_"+sd+"white_pits.gii")
        noisy_pits = os.path.join(subj_dir, "t1mri", "BL",
                                  "default_analysis", "segmentation", "mesh",
                                  "surface_analysis",
                                  s_id+ "_"+sd+"white_noisy_pits.gii")
        ridges = os.path.join(subj_dir, "t1mri", "BL",
                              "default_analysis", "segmentation", "mesh",
                              "surface_analysis",
                              s_id+ "_"+sd+"white_ridges.gii")
        basins = os.path.join(subj_dir, "t1mri", "BL",
                              "default_analysis", "segmentation", "mesh",
                              "surface_analysis",
                              s_id+ "_"+sd+"white_basins.gii")

        cmd = " ".join(['python -m brainvisa.axon.runprocess',
                        'SulcalPitsExtraction'
                        white_mesh, "%s" %  param[sd][0],
                        pole_cingular,
                        '%s %s %s' % (DPF_alpha, thresh_ridge, thresh_dist),
                        '%s %s %s' % (param[sd][2], thresh_area, param[sd][1]),
                        DPF,
                        pits,
                        noisy_pits,
                        ridges,
                        basins])
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
                        help='Nb of cpu cores to use (default %i)'%1, type=int)
    options = parser.parse_args()
    a = int(options.array_id)
    j = int(options.job)
    if options.ncore is not None:
        number_CPU = options.ncore
    else:
        number_CPU = 1

    DATABASE = ""
    s_ids = sorted(os.listdir(DATABASE))
    parameters = []
    for s_id in s_ids[a*j:j+a*j]:
        parameters.append([DATABASE, s_id])

    pool = Pool(processes = number_CPU)
    pool.map(pits_extraction, parameters)
    pool.close()
    pool.join()

