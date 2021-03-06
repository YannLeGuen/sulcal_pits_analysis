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

def smooth_pits(parameters):
    DATABASE, s_id = parameters
    subj_dir = os.path.join(DATABASE, s_id)

    # Default values
    dt, nb_iterations = 0.7, 60
    for sd in ['L', 'R']:
        white_mesh = os.path.join(subj_dir, "t1mri", "BL",
                                  "default_analysis", "segmentation", "mesh",
                                  s_id+ "_"+sd+"white.gii")
        pits = os.path.join(subj_dir, "t1mri", "BL",
                            "default_analysis", "segmentation", "mesh",
                            "surface_analysis", s_id+"_"+sd+"white_pits.gii")
        pits_smooth = os.path.join(subj_dir, "t1mri", "BL",
                                   "default_analysis", "segmentation", "mesh",
                                   "surface_analysis",
                                   s_id+ "_"+sd+"white_pits_smoothed0.7_60.gii")
        cmd = " ".join(['python -m brainvisa.axon.runprocess',
                        'PitsTextureSmoothing',
                        pits,
                        white_mesh,
                        '%s %s' % (dt, nb_iterations),
                        pits_smooth])
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
    pool.map(smooth_pits, parameters)
    pool.close()
    pool.join()

