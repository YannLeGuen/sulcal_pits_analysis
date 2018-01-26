

#0 Organize the file from Freesurfer output and convert them to gii

#1 Identified cingular pole on each individual cingular_projection.py
Using Freesurfer annotation in ?h.cortex.label

#2 Sulcal Pits Extraction sulcal_pits_detection.py
Using Brainvisa process SulcalPitsExtraction

#3 Smooth Sulcal Pits texture tex_pits_smoothing.py
Using Brainvisa process PitsTextureSmoothing

# Compute geodesic depth map on subject white mesh geodesic_depth_map.py
NB: This process is not required in the sulcal pits workflow

#4 Projecture textures onto the template project_to_template.py
Using Brainvisa process ProjectTextureOntoAtlas
Project the pits, DPF, geodesic depth and pits smoothed textures

#5 Compute Group pits density pits_density.py

#6 Group watershed on the group pits density group_watershed.py
to obtain the group-clusters of pits (i.e. areals)

# job_submission.pbs to submit jobs on a cluster