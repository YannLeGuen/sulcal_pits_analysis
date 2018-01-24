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
import pandas as pd
import numpy as np
import nibabel.gifti.giftiio as gio
# Import to create a poster figure with the snapshots
import PIL.Image
# Import related to Anatomist object manipulation
from soma import aims
import anatomist.api 
ana = anatomist.api.Anatomist()
import paletteViewer
# Define the referential and transformation necessary,
# because Freesurfer meshes are indexed in the reverse order.
ref = ana.createReferential()
tr = ana.createTransformation([0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, -1],
                              ref, ana.centralReferential())

# Method to notify the Anatomist window of changes and remove cursor
def updateWindow(window, obj):
    window.addObjects(obj)
    obj.setChanged()
    obj.notifyObservers()
    ana.execute('WindowConfig', windows=[window], cursor_visibility=0)

def compute_pits_frequencies(indir, threshold_freq):
    """
    Parameters
    indir: directory containing the case/control pits phenotype
    threshold_freq: threshold frequency
    """
    dict_freq = {}
    for side in ['R', 'L']:
        dict_freq[side] = {}
        file_in = os.path.join(indir, 'case_control',
                               'all_pits_side'+side+'.csv')
        df = pd.read_csv(file_in)
        nb_subj = df.shape[0]
        # Loop over the areals
        for col in df.columns:
            # Ignore the IID column
            if col != 'IID':
                df[col] = df[col]-1 # Code case/control is 2/1 so substract 1
                nb_pits = np.count_nonzero(df[col]) # Pits correspond to 1
                freq = float(nb_pits)/nb_subj
                # Make all areal with freq < threshold to 0
                dict_freq[side][col] = freq*(freq>threshold_freq)

    return dict_freq


def map_pits_heritability(template, side, areals, indir, cb, interval_cb,
                          dict_freq, pval_thresh, outdir):
    """
    Parameters
    template: path to template mesh file
    areals: path to gii file containing the texture of areals
    side: hemisphere side from which SOLAR estimates are displayed
    indir: directory containing h2 and pval dictionaries
    cb: colorbar names must be an array with [cb_h2, cb_pval]
    interval_cb: number boundaries for the colobar display
                 must be a double array [inter_h2[0,1], inter_pval[0,1]]
    dict_freq: dictionary containing pits frequency per areal
    pval_thresh: p-value threshold
    outdir: output directory for the snapshots
    """
    if not os.path.isdir(outdir):
        os.makedirs(outdir)
    # Define the quaternions that are gonna be used for the snapshots
    if (side == "L"):
        view_quaternions = {'intern' : [0.5, 0.5, 0.5, 0.5],
                            'extern' : [0.5, -0.5, -0.5, 0.5],
                            'bottom' : [1, 0, 0, 0],
                            'top': [0, 0, 1, 0]}
    else:
        view_quaternions = {'intern' : [0.5, 0.5, 0.5, 0.5],
                            'extern' : [0.5, -0.5, -0.5, 0.5],
                            'bottom' : [1, 0, 0, 0],
                            'top': [0, 0, 1, 0]}
    
    # Load h2 and pval dictionaries
    path_h2 = os.path.join(indir, 'h2_dict.json')
    with open(path_h2, 'r') as f:
        data = json.load(f)
    dict_h2 = json.loads(data)
    path_pval = os.path.join(indir, 'pval_dict.json')
    with open(path_pval, 'r') as f:
        data = json.load(f)
    dict_pval = json.loads(data)

    # Areal value in each vertex
    array_areals = gio.read(areals).darrays[0].data
    # Create arrays that will contain h2 and pval estimates on each vertex
    array_h2 = np.zeros(len(array_areals))
    array_pval = np.zeros(len(array_areals))
    list_areals =  np.unique(array_areals)
    
    # Load the mesh in anatomist
    mesh = ana.loadObject(template)

    for areal in list_areals:
        # make sure areal has the right format
        areal = str(int(areal))
        if (dict_h2[side].has_key(areal) and
            # Check if areal passed the threshold frequency in both hemispheres
            dict_freq['L']["Areal_"+areal]*dict_freq['R']["Areal_"+areal] > 0):
            # Finally, checl if the p-value pass the pval threshold
            if dict_pval[side][areal] < pval_thresh :
                # Select the index of vertices belonging to the areal
                ind = np.where(array_areals == float(areal))[0]
                # Give them the h2 estimate of this areal
                array_h2[ind] = dict_h2[side][areal]
                # And the log10 pvalue
                array_pval[ind] = -np.log10(dict_pval[side][areal])

    ### Block performing the display of h2 estimates ###
    # Create anatomist window
    window = ana.createWindow('3D', geometry=[0, 0, 584, 584])
    tex = aims.TimeTexture(dtype='FLOAT')
    tex[0].assign(array_h2)                  
    pits_tex = ana.toAObject(tex)
    tex_mesh = ana.fusionObjects([mesh, pits_tex], method='FusionTexSurfMethod')
    tex_mesh.assignReferential(ref)
    tex_mesh.setMaterial(front_face='counterclockwise')
    pits_tex.setPalette(cb[0], minVal=interval_cb[0][0],
                        maxVal=interval_cb[0][1], absoluteMode=True)
    updateWindow(window, tex_mesh)

    # Loop through the quaternions and do the snapshot
    for sd in view_quaternions.keys():
            q = aims.Quaternion(view_quaternions[sd])
            window.camera(view_quaternion=view_quaternions[sd],
                                 zoom=0.65)
            # Snapshot file output
            output = os.path.join(outdir, 'snapshot_h2_'+side+'_'+sd+'.png')
            ana.execute('WindowConfig', windows=[window], snapshot=output)


    #### Block performing the display of pvalues estimates ###
    # Create anatomist window, second time (else zoom problem for snapshots)
    window = ana.createWindow('3D', geometry=[0, 0, 584, 584])
    tex = aims.TimeTexture(dtype='FLOAT')
    tex[0].assign(array_pval)                  
    pits_tex = ana.toAObject(tex)
    tex_mesh = ana.fusionObjects([mesh, pits_tex],
                                 method='FusionTexSurfMethod')
    tex_mesh.assignReferential(ref)
    tex_mesh.setMaterial(front_face='counterclockwise')
    pits_tex.setPalette(cb[1], minVal=interval_cb[1][0],
                        maxVal=interval_cb[1][1], absoluteMode=True)
    updateWindow(window, tex_mesh)
    
    # Loop through the quaternions and do the snapshot
    for sd in view_quaternions.keys():
            q = aims.Quaternion(view_quaternions[sd])
            window.camera(view_quaternion=view_quaternions[sd],
                                 zoom=0.65)
            # Snapshot file output
            output = os.path.join(outdir, 'snapshot_pval_'+side+'_'+sd+'.png')
            ana.execute('WindowConfig', windows=[window], snapshot=output)



def create_poster(images, output_image, tile=4, padding=4, bgcolor=(255,255,255),
                  align_grid=False):
    # print 'calculating poster size...'
    rowsizes = []
    maxsize = [0, 0]
    n = 0
    for imagename in images:
        image = PIL.Image.open( imagename )
        if n % tile == 0:
            rowsizes.append( ( 0, 0 ) )
        rowsizes[-1] = ( rowsizes[-1][0] + image.size[0] + padding,
            max( rowsizes[-1][1], image.size[1] + padding ) )
        maxsize = (max((maxsize[0], image.size[0])),
                   max((maxsize[1], image.size[1])))
        n += 1
    print 'sizes:', rowsizes
    if align_grid:
        size = ((maxsize[0] + padding) * min(tile, len(images)) - padding,
                (maxsize[1] + padding) * len(rowsizes) - padding)
    else:
        size = ( max( [ rs[0] for rs in rowsizes ] ) +padding,
                 sum( [ rs[1] for rs in rowsizes ] ) +padding)

    print 'size:', size
    outimage = PIL.Image.new( 'RGB', size, bgcolor )
    print 'creating image...'
    n = 0
    line = 0
    xpos = 0
    ypos = 0
    #rowsizes.insert( 0, ( 0, 0 ) )
    for imagename in images:
        image = PIL.Image.open( imagename )
        if align_grid:
            bbox = ((maxsize[0] + padding) * n, (maxsize[1] + padding) * line,
                    0, 0)
            bbox = (bbox[0], bbox[1],
                    bbox[0] + image.size[0], bbox[1] + image.size[1])
        else:
            y = ypos + ( rowsizes[line][1] - image.size[1] ) / 2
            bbox = ( xpos, y, xpos+image.size[0], y + image.size[1] )
        outimage.paste( image, bbox )
        xpos = bbox[2] + padding
        n += 1
        if n == tile:
            n = 0
            line += 1
            xpos = 0
            ypos = bbox[3] + padding
    path = os.path.dirname(output_image)
    outimage.save( open( output_image, 'w' ) )
    # print 'done.'

def poster_h2_pval(indir, output):
    """
    Parameters
    indir: directory containing h2 and pval snapshots
    output: filename of the poster output
    """
    outdir = os.path.dirname(output)
    if not os.path.isdir(outdir):
        os.makedirs(outdir)
    padding = 6
    # List containing the snapshots filename ordered for the poster display
    images = []
    view = ['extern', 'intern', 'bottom', 'top']
    for ft in ['h2', 'pval']:
        for sd in ['L', 'R']:
            for v in view:
                image = os.path.join(indir, 'snapshot_'+ft+'_'+sd+'_'+v+'.png')
                images.append(image)
    # Create a poster with the snapshots
    create_poster(images,  output, padding=padding, tile=8, align_grid=False)

    # Code section to add the colorbars to the image
    # The colorbars need to have previously been saved in a file
    # Unfortunately Anatomist does not enable to do it through line of code
    # So it needs to be done manually
    image = PIL.Image.open(output)
    fts = ['herit45', 'pval4']
    colorbar = {}
    max_large = 0
    for ft in fts:
        colorbar[ft] = PIL.Image.open('/home/yl247234/Images/PALETTES/'
                                           'palette_'+ft+'.png')
        if int(colorbar[ft].size[0]) > max_large:
            max_large = int(colorbar[ft].size[0])
        colorbar[ft] = colorbar[ft].resize((int(colorbar[ft].size[0]),
                                            int((image.size[1]-
                                                 (len(fts)-1)*padding)
                                                /len(fts))))
    bgcolor = (255,255,255)
    size =( image.size[0]+padding+max_large, 
            image.size[1])
    outimage = PIL.Image.new( 'RGB', size, bgcolor)
    bbox = (0, 0, image.size[0], image.size[1])
    outimage.paste( image, bbox )
    for k,ft in enumerate(fts):
        if ft != None:
            bbox = (image.size[0]+padding,
                    int(image.size[1]*k/len(fts))+k*padding, size[0],
                    int(image.size[1]*k/len(fts))+k*padding+
                    int((image.size[1]-(len(fts)-1)*padding)/len(fts)))
            outimage.paste(colorbar[ft], bbox)

    path = os.path.dirname(output)
    outimage.save(open(output, 'w'))

     
if __name__ == '__main__':
    # Careful when running this script with colorbar interval between 0 and 1
    # You must use brainvisa bug fix and not the current release as of 06/2017
    ROOT_DIR = ""
    
    indir  = os.path.join(ROOT_DIR, '2016_HCP_pits_cleaned',
                          'pits_analysis_lh', 'phenotype_pits_DPF')
    # Minimum pits frequency in the areal, in both hemispheres, for the h2 of
    # this areal to be displayed on the template
    THRESHOLD_FREQ = 0.5
    # Dictoniary containing the pits frequency per areal
    dict_freq = compute_pits_frequencies(indir, THRESHOLD_FREQ)

    template  = os.path.join(ROOT_DIR, 'folder_gii', 'lh.inflated.white.gii')
    path_areals = os.path.join(ROOT_DIR, '2016_HCP_pits_cleaned',
                               'pits_density_sym_lh')
    areals = os.path.join(path_areals,
                          ('clusters_total_average_pits_smoothed0.7_60_sym_lh'
                           '_dist15.0_area100.0_ridge2.0.gii'))
    indir = os.path.join(ROOT_DIR, '2016_HCP_pits_cleaned',
                          'pits_analysis_lh', 'dictionaries_heritability')
    outdir = os.path.join(ROOT_DIR, '2016_HCP_pits_cleaned',
                          'pits_analysis_lh', 'snapshots')
    # Define colorbar names for h2 and pval respectively
    cb = ['Yellow-red-fusion', 'green_yellow_red']
    interval_cb = [[0,0.45], [1, 8]]
    # Divide by nb areals (58)* nb hemispheres (2) to get Bonferroni threshold
    PVAL_THRESH = 5e-2#/(58*2)
    for side in ['R', 'L']:
        # Map heritability and obtained single snapshots
        map_pits_heritability(template, side, areals, indir, cb, interval_cb,
                              dict_freq, PVAL_THRESH, outdir)


    indir = outdir
    outdir = os.path.join(ROOT_DIR, '2016_HCP_pits_cleaned',
                          'pits_analysis_lh', 'poster_results')
    output = os.path.join(outdir, 'figure4_lh_HCP_poster.png')
    # Nicely crop the previous png snapshots
    os.system('ls '+indir+'/snapshot*.png|while read input;'
              'do convert  $input -trim /tmp/toto.png;'
              'convert  /tmp/toto.png -transparent white $input;'
              'done')
    # Remove intermediate file
    os.system('rm /tmp/toto.png')
    # Create the poster with snapshots and colorbars
    # Need to have saved the colorbar snapshots manually before see code above
    poster_h2_pval(indir, output)
