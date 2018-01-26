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

def watershed_parcels(parameters):
    # f_density should contain the type of fsaverage template
    # either "_asym_"(fsaverage normal template),
    # "_sym_lh" (fsaverage_sym left side),
    # "_sym_rh" (fsaverage_sym right side)
    # this will be stored in tp_sd
    sd, ROOT_DIR, f_density= parameters

    # Default threshold parameters
    thr_dist, thr_area, thr_ridge = 15.0, 100.0, 2.0
    thr_dist, thr_area, thr_ridge = str(thr_dist), str(thr_area), str(thr_ridge)

    if "_sym_" in f_density:
        tp_sd = f_density[-7:]
        template = os.path.join(ROOT_DIR, "folder_gii", "sym",
                                sd.lower()+'h.inflated4.smooth10.white.gii')
        f_pole_cingular = os.path.join(ROOT_DIR, "folder_gii", "sym",
                                        sd+'white_pole_cingular.gii')
    else:
        tp_sd = "_asym"
        template =  os.path.join(ROOT_DIR, "folder_gii",
                                 sd.lower()+'h.inflated.white.gii')
        f_pole_cingular = os.path.join(ROOT_DIR, "folder_gii",
                                       sd+'white_pole_cingular.gii')

    output = os.path.join(ROOT_DIR, "pits_density"+tp_sd,
                          ('clusters_'+f_density+'_dist'+thr_dist+
                           '_area'+thr_area+'_ridge'+thr_ridge+'.gii'))
    f_density = os.path.join(ROOT_DIR, "pits_density"+tp_sd, f_density+'.gii')


    cmd = " ".join(['python -m brainvisa.axon.runprocess Group_watershed'
                    template,
                    f_density;
                    f_pole_cingular,
                    '%s %s %s' % (thr_dist, thr_area, thr_ridge),
                    output])
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

    ROOT_DIR = ""
    parameters = []
    for sd in ['L', 'R']:
        #parameters.append([sd, ROOT_DIR, sd+'_average_pits_smoothed0.7_60_asym'])
        for sd2 in ['lh', 'rh']:
            parameters.append([sd, ROOT_DIR,
                               sd+'_average_pits_smoothed0.7_60_sym_'+sd2])

            
    for sd2 in ['lh', 'rh']:
        parameters.append([sd, ROOT_DIR,
                           'total_average_pits_smoothed0.7_60_sym_'+sd2])

    pool.map(watershed_parcels, parameters[a*j:j+a*j])
    pool.close()
    pool.join()
    print "Elapsed time: " + str(time.time()-t)
