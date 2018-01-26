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
# Use freesurfer wrapper from A. Grigis
from freesurferW.exceptions import FreeSurferRuntimeError
from freesurferW.wrappers import FSWrapper
from soma import aims

def mri_convert( input_volume_file, output_volume_file,
                  fsconfig="/i2bm/local/freesurfer/SetUpFreeSurfer.sh"):
    cmd = ["mri_convert", input_volume_file,
           output_volume_file]

    # Execute the FS command
    recon = FSWrapper(cmd, shfile=fsconfig)
    recon()
    if recon.exitcode != 0:
        raise FreeSurferRuntimeError(
            recon.cmd[0], " ".join(recon.cmd[1:]), recon.stderr + recon.stdout)
    

def compute_geodesic_depth_map(parameters):
    DIR_F, DIR_BV, s_id = parameters
    path_nu = os.path.join(DIR_F, s_id, 'T1w', s_id, 'mri', 'nu.mgz')
    path_nat = os.path.join(DIR_F, s_id, 'T1w', s_id, 'mri', 'orig', '001.mgz')
    nu_nii = os.path.join(DIR_BV, s_id, 'mri', 'nu.nii.gz')
    native_nii = os.path.join(DIR_BV, s_id, 'mri', 'native.nii.gz')
    outdir = os.path.dirname(nu_nii)
    if not os.path.isdir(outdir):
        os.makedirs(outdir)
    #mri_convert(path_nu, nu_nii)
    #mri_convert(path_nat, native_nii)
    path_ribbon = os.path.join(DIR_F, s_id, 'T1w', s_id, 'mri', 'ribbon.mgz')
    ribbon_nii = os.path.join(DIR_BV, s_id, 'mri', 'ribbon.nii.gz')
    #mri_convert(path_ribbon, ribbon_nii)
    for vol in [nu_nii, native_nii, ribbon_nii]:
        # Convert file to Aims format
        cmd = " ".join(['AimsFileConvert',
                        '-i %s' % vol,
                        '-o %s' % vol,
                        '-t S16'])
        print cmd
        os.system(cmd)
        
    # Replace the label in the ribbon with Morphologist expected label
    ribbon_bv = os.path.join(DIR_BV, s_id, 'mri', 'ribbon_bv.nii.gz')
    cmd = " ".join(['AimsReplaceLevel',
                    '-i %s' % ribbon_nii,
                    '-o %s' % ribbon_bv,
                    '-g 42 41 2 3',
                    '-n 100 200 200 100'])
    print cmd
    os.system(cmd)
    # Split the hemisphere and replace labels
    l_gw = os.path.join(DIR_BV, s_id, 'mri', 'Lgrey_white'+s_id+'.nii.gz')
    r_gw = os.path.join(DIR_BV, s_id, 'mri', 'Rgrey_white'+s_id+'.nii.gz')
    cmd = " ".join(['AimsReplaceLevel',
                    '-i %s' % ribbon_nii,
                    '-o %s' % l_gw,
                    '-g 42 41 2 3',
                    '-n 0 0 200 100'])
    print cmd
    os.system(cmd)
    cmd = " ".join(['AimsReplaceLevel',
                    '-i %s' % ribbon_nii,
                    '-o %s' % r_gw,
                    '-g 42 41 2 3',
                    '-n 100 200 0 0'])
    print cmd
    os.system(cmd)

    # Compute the grey/white histogram
    histo = os.path.join(DIR_BV, s_id, 'mri', 'nu_'+s_id+'.his')
    cmd = " ".join(['VipGreyStatFromClassif',
                    '-i %s' % nu_nii,
                    '-c %s' % ribbon_bv,
                    '-a %s' % histo,
                    '-g 100 -w 200'])
    print cmd
    os.system(cmd)
    # Create the hemi cortex required for the geodesic depth process
    l_hemcor = os.path.join(DIR_BV, s_id, 'mri', 'Lcortex_'+s_id+'.nii.gz')
    r_hemcor = os.path.join(DIR_BV, s_id, 'mri', 'Rcortex_'+s_id+'.nii.gz')
    for gw, hemcor in zip([l_gw, r_gw], [l_hemcor, r_hemcor]):
        cmd = " ".join(['VipHomotopic',
                        '-i %s' % nu_nii,
                        '-cl %s' % gw,
                        '-h %s' % histo+'.han',
                        '-o %s' % hemcor,
                        '-m C -w t'])
        print cmd
        os.system(cmd)

    # Freesurfer surface are in native space, but not the ribbon
    # Thus, we tranform the obtained hemicortex to the native space
    ribbon = aims.read(ribbon_nii) 
    native = aims.read(native_nii)
    s2sb = aims.AffineTransformation3d(ribbon.header()['transformations'][0])
    n2sb = aims.AffineTransformation3d(native.header()['transformations'][0])
    s2n = n2sb.inverse() * s2sb
    transfo_fs2nat = os.path.join(DIR_BV, s_id, 'mri', 'fs2nat.trm')
    aims.write(s2n, transfo_fs2nat)
    for hemcor in [l_hemcor, r_hemcor]:
        cmd = " ".join(['AimsResample',
                        '-i %s' % hemcor,
                        '-o %s' % hemcor,
                        '-m %s' % transfo_fs2nat,
                        '-t  0',
                        '-r %s' % native_nii])
        print cmd
        os.system(cmd)
    
    l_white = os.path.join(DIR_BV, s_id, 't1mri/BL/default_analysis',
                           'segmentation/mesh', s_id+'_Lwhite.gii')
    l_white_aims = os.path.join(DIR_BV, s_id, 't1mri/BL/default_analysis',
                                'segmentation/mesh', s_id+'_Lwhite_aims.gii')
    l_depth = os.path.join(DIR_BV, s_id, 't1mri/BL/default_analysis',
                           'segmentation/mesh/surface_analysis',
                           s_id+'_Lgeodesic_depth.gii')
    r_white = os.path.join(DIR_BV, s_id, 't1mri/BL/default_analysis',
                           'segmentation/mesh', s_id+'_Rwhite.gii')
    r_white_aims = os.path.join(DIR_BV, s_id, 't1mri/BL/default_analysis',
                            'segmentation/mesh', s_id+'_Rwhite_aims.gii')
    r_depth = os.path.join(DIR_BV, s_id, 't1mri/BL/default_analysis',
                           'segmentation/mesh/surface_analysis',
                           s_id+'_Rgeodesic_depth.gii')

    # Convert mesh to Aims Referential
    from freesurfer.freesurferMeshToAimsMesh import freesurferMeshToAimsMesh
    for white, white_aims in zip([l_white, r_white],
                                 [l_white_aims, r_white_aims]):
        freesurferMeshToAimsMesh(white, native_nii, white_aims)

        

    # Finally compute the geodesic depth of interest 
    for white, hemcor, depth in zip([l_white_aims, r_white_aims],
                                    [l_hemcor, r_hemcor], [l_depth, r_depth]):
            
        cmd = " ".join(['python -m brainvisa.axon.runprocess',
                        'whitemeshdepthmap',
                        white,
                        hemcor,
                        depth,
                        '10.0']) # closing_size
        print cmd
        os.system(cmd)


        # Another option is to compute the geodesic depth map in the volume
        # The code is commented below, however it performs sighlty worse
        """
        cmd = " ".join(['AimsThreshold',
                        hemcor,
                        white,
                        '-m eq -t 0 -b'])
        print cmd
        os.system(cmd)

        closedwhite = os.path.join('/tmp', os.path.basename(white[:-4]+
                                                            '_closed.gii'))
        cmd = " ".join(['AimsMorphoMath -m clo',
                        '-i %s' % white,
                        '-o %s' % closedwhite,
                        '-r 10']) # closing_size
        print cmd
        os.system(cmd)
        cmd = " ".join(['AimsMorphoMath -m ero',
                        '-i %s' % closedwhite,
                        '-o %s' % closedwhite,
                        '-r 2.']) # erosion_size?
        print cmd
        os.system(cmd)
        
        hulltex = os.path.join('/tmp', os.path.basename(white[:-4]+
                                                         '_hull.gii'))
        # How to call context from Python?
        fusion = context.runProcess('fusion3Dmesh', closedwhite, white,
                                    hulltex, 0)
                                 
        cmd = " ".join(['cartoLinearComb.py',
                        '-i %s' % hulltex,
                        '-o %s' % hulltex,
                        '-f -I1 + 32767'])
        print cmd
        os.system(cmd)
        # Convert file to Aims format
        cmd = " ".join(['AimsFileConvert',
                        '-i %s' % hulltex,
                        '-o %s' % hulltex,
                        '-t S16'])
        print cmd
        os.system(cmd)

        cmd = " ".join(['AimsMeshDistance',
                        '-i %s' % white,
                        '-o %s' % depth,
                        '-t %s' % hulltex])
        print cmd
        os.system(cmd)

        """    




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

    # Directory Freesurfer database
    DIR_F = ''
    # Directory Brainvisa/Morphologist database
    DIR_BV = ''
    s_ids = sorted(os.listdir(DIR_F))
    parameters = []
    for s_id in s_ids[a*j:j+a*j]:
        parameters.append([DIR_F, DIR_BV, s_id])

    pool = Pool(processes = number_CPU)
    pool.map(pits_extraction, parameters)
    pool.close()
    pool.join()










        
