#! /usr/bin/env python
# -*- coding: utf-8 -*

##########################################################################
# @author: yann.leguen@cea.fr
# Copyright (C) CEA, 2017
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html for details.
##########################################################################

from soma import aims
from brainvisa.cortical_surface.surface_tools import texture_tools as textureTls
import os
import argparse
import numpy as np
from multiprocessing import Pool

def cingular_project(parameters):
    DATABASE, s_id = parameters
    
    subj_dir = os.path.join(DATABASE, s_id)

    for sd in ['L', 'R']:
        file_white_mesh = os.path.join(subj_dir, "t1mri", "BL",
                                       "default_analysis", "segmentation",
                                       "mesh", s_id+ "_"+sd+"white.gii")
        pole_cingular = os.path.join(subj_dir, "t1mri", "BL",
                                     "default_analysis", "segmentation",
                                     "mesh", "surface_analysis",
                                     s_id+ "_"+sd+"white_pole_cingular.gii")
        #aparc = os.path.join(subj_dir, "label", 
        #                     s_id+ "."+sd+".aparc.native.label.gii")
        fs_cortex_label = os.path.join(subj_dir, "label", 
                                       sd.lower()+'h.cortex.label')
        re = aims.Reader()
        ws = aims.Writer()
        mesh = re.read(file_white_mesh)
        nbv = len(mesh.vertex())
        data = np.ones(nbv)
        
        f = open(fs_cortex_label,'r')
        lines = f.readlines()
        for a in lines[2:]:
            a_spl = a.split()
            data[int(a_spl[0])] = a_spl[4]
        f.close()
        #data_in = gio.read(file_aparc).darrays[0].data
        #data = np.zeros_like(data_in)
        #data[data_in==-1]=1
        a, b = textureTls.textureTopologicalCorrection(mesh, data, 1)
        cingular_tex_clean, cingular_tex_boundary = a, b
        tex_out = aims.TimeTexture_S16()
        tex_out[0].assign(cingular_tex_clean)
        ws.write(tex_out, pole_cingular)
        


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
    pool.map(cingular_project, parameters)
    pool.close()
    pool.join()



    # Below convert cingulate template to giifti
    """ ROOT_DIR = ""
    sds = ['R', 'L']
    path = os.path.join(ROOT_DIR, 'folder_gii/sym')
    for sd in ['L', 'R']:
        white_mesh = os.path.join(path, sd.lower()+"h.white.gii")
        pole_cingular = os.path.join(path, sd+"white_pole_cingular.gii")
        re = aims.Reader()
        ws = aims.Writer()
        mesh = re.read(white_mesh)
        nbv = len(mesh.vertex())
        data = np.ones(nbv)
        fs_cortex_label = os.path.join(ROOT_DIR, 'fsaverage_sym/label',
                                       sd.lower()+'h.cortex.label')
        f = open(fs_cortex_label,'r')
        lines = f.readlines()
        for a in lines[2:]:
            a_spl = a.split()
            data[int(a_spl[0])] = a_spl[4]
        f.close()

        a, b = textureTls.textureTopologicalCorrection(mesh, data, 1)
        cingular_tex_clean, cingular_tex_boundary = a, b
        tex_out = aims.TimeTexture_S16()
        tex_out[0].assign(cingular_tex_clean)
        ws.write(tex_out, pole_cingular)"""
