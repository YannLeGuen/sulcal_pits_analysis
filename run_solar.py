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
import glob
import re
import numpy as np
import pandas as pd

# verbose for details of race and ethnicity
verbose = False

def create_solar_pedigree(org_ped_restrict, org_ped_unrestrict, output):
    """
    Parameters
    org_ped_restricted: file containing the restricted pedigree from HCP
    org_ped_unrestricted: file containing the unrestricted pedigree from HCP
    output: output file containing SOLAR formatted pedigree
    """
    outdir = os.path.dirname(output)
    if not os.path.isdir(outdir):
        os.makedirs(outdir)
    # SOLAR pedigree required columns
    columns = ['ID', 'FA', 'MO', 'Age', 'SEX']
    # Load the full restricted HCP data
    df_src = pd.read_csv(org_ped_restrict)
    # Select columns of interest
    df = df_src[['Subject', 'Father_ID', 'Mother_ID', 'Age_in_Yrs']]
    df.index = df['Subject']
    # Gender is in unrestricted data
    df_src2 = pd.read_csv(org_ped_unrestrict)
    df2 = df_src2[['Subject', 'Gender']]
    df2.index = df2['Subject']
    # Merge the gender in first dataframe
    df['Gender']  = df2['Gender']
    # Assign the SOLAR expected column names
    df.columns = columns

    # SOLAR expect a column with the twin pair number for MZ twins
    # The code below prepares this column (not optimized)
    # We also prepare household ID column for the pair of twin only if needed
    df_short = df_src[['Subject', 'Zygosity_corr', 'Mother_ID', 'Father_ID']]
    df_short.index = df_short['Subject']
    mz_zyg = np.zeros(len(df_short.index))
    twin_zyg_hid = np.zeros(len(df_short.index))
    mztwins_done = []
    dztwins_done = []
    count = 1
    count_hid = 1
    for i,zyg in enumerate(df_short['Zygosity_corr']):
        if zyg == 'MZ' and df_short.index[i] not in mztwins_done:
            MO = df_short.loc[df_short.index[i]]['Mother_ID']
            FA = df_short.loc[df_short.index[i]]['Father_ID']
            df_sub = df_short.loc[(df_short['Mother_ID'] == MO) &
                                  (df_short['Father_ID'] == FA) ]

            df_sub = df_sub.loc[df_sub['Zygosity_corr'] == 'MZ']
            # Check if at least 2 twins in the pair,
            # because sometimes HCP has singleton pair (at least in S900)..
            if len(df_sub.index) >= 2:
                for subj in df_sub.index:
                    mztwins_done.append(subj)
                    ind = list(df_short.index).index(subj)
                    mz_zyg[ind] = count
                    twin_zyg_hid[ind] = count_hid
                count += 1
                count_hid += 1
        elif zyg == 'NotMZ' and df_short.index[i] not in dztwins_done:
            MO = df_short.loc[df_short.index[i]]['Mother_ID']
            df_sub = df_short.loc[df_short['Mother_ID'] == MO ]
            df_sub = df_sub.loc[df_sub['Zygosity_corr'] == 'NotMZ']
            if len(df_sub.index) >= 2:
                for subj in df_sub.index:
                    mztwins_done.append(subj)
                    ind = list(df_short.index).index(subj)
                    twin_zyg_hid[ind] = count_hid
                count_hid += 1
        else:
            pass
    df['MZTWIN'] = mz_zyg
    df['MZTWIN'] = [ int(i) for i in df['MZTWIN'] ]
    # Format the columns ID, FA and MO
    for column in ['ID', 'FA', 'MO']:
        df[column] = ['%06d' % int(i) if not np.isnan(i)
                      else 0  for i in df[column]]

    # Create required family ID column by concatenating MO and FA
    l  = list(df['MO'].astype(str)+df['FA'].astype(str))
    df['FAMID'] = l
    for column in ['FAMID']:
        df[column] = ['%012d' % i for i in df[column].astype(int)]

    # Save SOLAR formatted pedigree (we did not include HID yet)
    df.to_csv(output, header=True, index=False)


def convert_pheno_to_solar(pheno, org_ped_restricted, restrict_to_white=True):
    """
    Parameters
    pheno: file path containing the phenotype
    restrict_to_white: boolean to restrict the analysis to white
    org_ped_restricted: file containing the restricted pedigree from HCP
    """
    cov = org_ped_restricted

    df_cov = pd.read_csv(cov)
    df_cov  = df_cov[['Subject', 'Age_in_Yrs', 'Race', 'Ethnicity']]
    df_cov.index = df_cov['Subject']
    # Identify all the available Race and Ethinicity labels
    race_types = list(set(df_cov['Race']))
    ethnicity_types = list(set(df_cov['Ethnicity']))
    
    # Print some statistics for the whole HCP dataset
    if verbose:
        print "\nTotal individual "+str(df_cov.shape[0])
        print " RACE"
        for race in race_types:
            print race +" "+str(df_cov.loc[df_cov['Race']== race].shape[0])
        print "\n ETHNI"
        for eth in ethnicity_types:
            print eth+" "+str(df_cov.loc[df_cov['Ethnicity']== ethni].shape[0])


    df_phen = pd.read_csv(pheno)
    df_phen = df_phen.dropna()
    df_phen.index = df_phen['IID'].astype('int')
    # Add the Age column as covariate
    df_phen['Age'] = df_cov['Age_in_Yrs']
    # Find the phenotype column name
    columns = list(set(df_phen.columns)-set(['IID']))
    df_phen = df_phen[[u'IID', u'Age']+columns]
                      
    # Exclude individuals with race and ethnicity not included
    df_cov = df_cov.drop(df_cov.loc[df_cov['Race']==
                                    'Unknown or Not Reported'].index)
    df_cov = df_cov.drop(df_cov.loc[df_cov['Ethnicity']==
                                    'Unknown or Not Reported'].index)
    df_cov = df_cov.drop(df_cov.loc[df_cov['Race']==
                                    'More than one'].index)
    df_cov = df_cov.drop(df_cov.loc[df_cov['Race']==
                                    'Am. Indian/Alaskan Nat.'].index)
    df_phen = df_phen.loc[df_cov.index]

    # Bit of code to restrict the analysis to white not hispanic
    if restrict_to_white:
        df_phen = df_phen.loc[df_cov.loc[df_cov['Race']=='White'].index]
        df_cov = df_cov.loc[df_phen.index]


    # Update the list of races and ethinicities
    race_types = list(set(df_cov['Race']))
    ethnicity_types = list(set(df_cov['Ethnicity']))

    if verbose:
        # Display statistic for this phenotype
        print "\nTotal individual "+str(df_cov.shape[0])
        print " RACE"
        for race in race_types:
            print race+" "+str(df_cov.loc[df_cov['Race']== race].shape[0])
        print " ETHNI"
        for eth in ethnicity_types:
            print eth+" "+str(df_cov.loc[df_cov['Ethnicity']== ethni].shape[0])
        
    count = 0
    count_max = None
    max_val = 0
    print "\n"
    for race in race_types:
        for ethni in ethnicity_types:
            if df_cov.loc[ (df_cov['Ethnicity']== ethni) &
                           (df_cov['Race']== race) ].shape[0] != 0:
                df_phen['Group'+str(count)] = np.asarray(
                    ((df_cov['Ethnicity']== ethni) &
                     (df_cov['Race']== race))).astype('int')

                if max_val < sum(df_phen['Group'+str(count)]):
                    max_val = sum(df_phen['Group'+str(count)])
                    count_max = count
                count += 1
    # We withdraw the main group
    df_phen = df_phen.drop('Group'+str(count_max),1)
    # If any drop Nan values
    df_phen = df_phen.dropna()
    df_phen['IID'] = df_phen['IID'].astype('int')
    # Do not overwrite the previous pheno file
    dirout = os.path.join(os.path.dirname(pheno), "solar_ready")
    if not os.path.isdir(dirout):
        os.makedirs(dirout)    
    pheno = os.path.join(dirout, os.path.basename(pheno))
    df_phen.to_csv(pheno, header=True, index=False)

    
def run_solar(areal, side, pheno, work_dir, org_ped_restricted):
    """
    Parameters
    areal: considered areal number
    side: hemisphere side
    pheno: file path containing the phenotype
    work_dir: root folder in which solar pedigree and outputs are written
    org_ped_restricted: file containing the restricted pedigree from HCP
    """
    trait = "Areal_"+str(areal)
    outdir = os.path.join(work_dir, side, trait)
    if not os.path.isdir(outdir):
        os.makedirs(outdir)

    # Convert phenotype to solar
    convert_pheno_to_solar(pheno, org_ped_restricted,
                           restrict_to_white=True)
    pheno = os.path.join(os.path.dirname(pheno), "solar_ready",
                         os.path.basename(pheno))
    # Run solar
    cmd = "solar pheno_analysis "+work_dir+" "+outdir+" "+trait+" "+pheno
    print cmd
    # Note: The code for Group covariates is not self adaptive yet
    # i.e. need to check pheno_analysis.tcl file to adjust accordingly
    os.system(cmd)

def parse_solar_out(work_dir, dict_areals, outdir='/tmp/'):
    """
    Parameters
    work_dir: root folder in which solar pedigree and outputs are written
    dict_areals: list of areals retained for each side
    outdir: output directory for h2 and pval dictionaries
    """
    if not os.path.isdir(outdir):
        os.makedirs(outdir)
    
    # Dictionaries containing h2 and associated p-values estimates
    dict_h2 = {}
    dict_pval = {}
    # Reference list of numbers for parsing
    nb = ['.', '-', 'e', '0', '1', '2', '3', '4', '5', '6', '7', '8', '9', '10']
    
    for side in ['R', 'L']:
        dict_h2[side] = {}
        dict_pval[side] = {}
        for areal in dict_areals[side]:
            solar_output = os.path.join(work_dir, side, 'Areal_'+areal,
                                        'polygenic.out')
            for line in open(solar_output, 'r'):
                # Find the heritability estimate and pval line
                if 'H2r is' in line and '(Significant)' in line:
                    print line[4:len(line)-15]
                    h2 = line[11:len(line)-30] 
                    p = line[26:len(line)-15]                
                    for k,l in enumerate(h2):
                        if not (l  in nb):
                            break
                    h2 = float(h2[:k])
                    p = float(p)
                    print "We extracted h2: "+str(h2)+" pval: "+str(p)
                    dict_h2[side][areal] = h2
                    dict_pval[side][areal] = p
                    """
                    Code below parse output in case of common environment model
                    It needs to be adapted to new variables

                    if 'C2 is' in line and '(Significant)' in line:
                        #print line[4:len(line)-15]
                        c2 = line[12:len(line)-30] 
                        p = line[26:len(line)-15]                
                        for k,l in enumerate(c2):
                            if not (l  in nb):
                                break
                        c2 = float(c2[:k])
                        p = float(p)
                        #if p<5e-2/10.0:
                        #print "We extracted c2: "+str(c2)+" pval: "+str(p)
                        #dict_h2[pheno][areal][i] = h2
                        #dict_pval[pheno][areal][i] = p
                        dict_h2[pheno][trait] = c2
                        dict_pval[pheno][trait] = p

                    elif 'C2 is' in line:# and '(Significant)' in line:
                        #print line
                        if len(line) < 22:
                            c2 = line[11:len(line)-1]
                            p = 1
                            #print "We extracted c2: "+str(c2)+" pval: "+str(p)
                        else:
                            c2 = line[12:len(line)-30] 
                            if 'Not Significant' in line:
                                p = float(line[26:len(line)-19])
                            else:
                                p = float(line[26:len(line)-15])
                            for k,l in enumerate(c2):
                                if not (l  in nb):
                                    break
                            c2 = float(c2[:k])
                            p = float(p)
                            #print "We extracted c2: "+str(c2)+" pval: "+str(p)
                        dict_h2[pheno][trait] = c2
                        dict_pval[pheno][trait] = p
                    """
    encoded = json.dumps(dict_h2)
    output_h2 = os.path.join(outdir, 'h2_dict.json')
    with open(output_h2, 'w') as f:
        json.dump(encoded, f)

    encoded = json.dumps(dict_pval)
    output_pval = os.path.join(outdir, 'pval_dict.json')
    with open(output_pval, 'w') as f:
        json.dump(encoded, f)

        
if __name__== '__main__':

    ROOT_DIR = ""
    
    # Original HCP files restricted and unrestricted data
    restrict = os.path.join(ROOT_DIR, '2016_HCP_pits_cleaned',
                         'pedigree', 'RESTRICTED_HCP_file.csv')
    unrestrict = os.path.join(ROOT_DIR, '2016_HCP_pits_cleaned',
                              'pedigree', 'unrestricted_HCP_file.csv')
    # Pedigree output formatted for SOLAR
    pedigree = os.path.join(ROOT_DIR, '2016_HCP_pits_cleaned/pedigree/',
                            'HCP_S1200_pedigree_no_hid.csv')
    # Create the pedigree
    create_solar_pedigree(restrict, unrestrict, pedigree)
    
    # Folder in which solar pedigree and outputs are written.
    # cf scripts makeped.tcl & pheno_analysis.tcl
    work_dir = os.path.join(ROOT_DIR, '2016_HCP_pits_cleaned',
                            'pits_analysis_lh/solar_work_dir')
    if not os.path.isdir(work_dir):
        os.makedirs(work_dir)
    
    # Contains the phenotype saved by pits_extraction.py
    pheno_dir  = os.path.join(ROOT_DIR, '2016_HCP_pits_cleaned',
                              'pits_analysis_lh', 'phenotype_pits_DPF')

    os.system("solar makeped "+work_dir+" "+pedigree)
    # Dictionary containing list of retained areals
    dict_areals = {}
    # There is likely a more efficient solution than to use glob and re.
    # For example save an intermediary dictionary in pits_extraction.py
    for side in ['R', 'L']:
        dict_areals[side] = []
        for filename in glob.glob(os.path.join(pheno_dir,'*.csv')):
            if 'side'+side in filename:
                m = re.search(pheno_dir+'/DPF_pit(.+?)side'+side, filename)
                if m:
                    areal = m.group(1)
                    dict_areals[side].append(str(areal))
                    # NB: SOLAR encounter issues when running
                    # several threads in parallel
                    run_solar(areal, side, filename, work_dir, restrict)
    
    outdir = os.path.join(ROOT_DIR, '2016_HCP_pits_cleaned',
                          'pits_analysis_lh', 'dictionaries_heritability')
    parse_solar_out(work_dir, dict_areals, outdir=outdir)
