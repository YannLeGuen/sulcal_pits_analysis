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
import nibabel.gifti.giftiio as gio

if __name__ == '__main__':

    ROOT_DIR = ""
    s_ids = sorted(os.listdir(ROOT_DIR))

    for sd in ['rh', 'lh']:
        OUTDIR = os.path.join(ROOT_DIR, 'pits_density_sym_'+sd)
        if not os.path.exists(OUTPUT):
            os.makedirs(OUTPUT)
        fname_R = os.path.JOIN(OUTDIR,
                               'R_average_pits_smoothed0.7_60_sym_'+sd+'.gii')
        fname_L = os.path.JOIN(OUTDIR,
                               'L_average_pits_smoothed0.7_60_sym_'+sd+'.gii')
        fname_LR = os.path.JOIN(OUTDIR,
                                'total_average_pits_smoothed0.7_60_sym_'+sd+'.gii')

        pits_data_R = np.array([])
        pits_data_L = np.array([])
        count_R = 0
        count_L = 0
        total_density = np.array([])
        for k, s_id in enumerate(s_ids):
            path_s =  os.path.join(ROOT_DIR, s_id)
            if "imagen" in path and "Freesurfer" not in path:
                s_id = s_id[len(s_id)-12:]
            pits_R = os.path.join(path_s, 't1mri', 'BL', 'default_analysis'
                                  'segmentation', 'mesh',
                                  'surface_analysis_sym_'+sd,
                                  s_id+'_Rwhite_pits_smoothed0.7_60_on_atlas.gii')
            pits_L = os.path.join(path_s, 't1mri', 'BL', 'default_analysis'
                                  'segmentation', 'mesh',
                                  'surface_analysis_sym_'+sd,
                                  s_id+'_Lwhite_pits_smoothed0.7_60_on_atlas.gii')
            if os.path.isfile(pits_R):
                count_R +=1
                if pits_data_R.size == 0:
                    pits_data_R = gio.read(pits_R).darrays[0].data
                else:
                    pits_data_R += gio.read(pits_R).darrays[0].data
            else:
                print k
                print s_id
            if os.path.isfile(pits_L):
                count_L +=1       
                if pits_data_L.size == 0:
                    pits_data_L = gio.read(pits_L).darrays[0].data
                else:
                    pits_data_L += gio.read(pits_L).darrays[0].data

        pits_data_R_temp = pits_data_R/count_R
        pits_data_L_temp = pits_data_L/count_L

        total_density = (pits_data_R_temp+pits_data_L_temp)/2
        g_L_R = gio.read(pits_L)
        g_L_R.darrays[0].data = total_density
        gio.write(g_L_R, fname_LR)

        g_R = gio.read(pits_L)
        g_L = gio.read(pits_L)
        g_R.darrays[0].data = pits_data_R_temp
        g_L.darrays[0].data = pits_data_L_temp
        gio.write(g_R, fname_R)
        gio.write(g_L, fname_L)
