"""
solar_heritability.py

Functions to compute heritability.

@author Nicholas P. Tatonetti
"""

import os
import sys
import csv
import time
import gzip
import numpy
import random
import shutil
from collections import defaultdict

from tqdm import tqdm

from familial_recurrence import assign_family_ethnicities
from familial_recurrence import load_demographics
from familial_recurrence import load_family_ids

from solar_setup_h2_analysis import fixed_effects_averages
from solar_setup_h2_analysis import load_generic_pedigree

common_data_path = './common_data/'
_ccs_map_file_path = os.path.join(common_data_path, 'ccs_icd9_map.txt')


tcl_load_string = """
proc loadped {} {
    load pedigree pedigree.ped
    load phenotypes phenotypes.phen
}
"""

tcl_analysis_string_b = """
proc runanalysis {} {
    trait pheno
    covariates sex age
    polygenic -screen
}
"""

# We assume that all of our quantitative variables will not strictly follow
# a normal distribution. Most are skewed one way or another. SOLAR provides
# some functions for helping to correct for biases introduced by non normal
# traits. We use tdist as it is the more computationally efficeint.
tcl_analysis_string_q = """
proc runanalysis {} {
    trait pheno
    covariates sex age
    tdist
    polygenic -screen
}
"""

TRAIT_TYPE_QUANTITATIVE = 0
TRAIT_TYPE_BINARY = 1

def build_solar_directories(h2_path, iid2ped, empi2trait, ethnicities, min_ascertained, fam2empi, fam2eth, fam2count, fam2proband, use_proband, verbose = True, family_ids_only = None, trait_type = TRAIT_TYPE_BINARY):
    
    run_list = dict()
    
    if not os.path.exists(h2_path):
        os.mkdir(h2_path)
    
    for ethnicity in ethnicities:
        
        trait_ped = list()
        
        for fid in family_ids_only:
            if fam2count[fid] < min_ascertained:
                continue
            if ethnicity != 'ALL' and not fam2eth[row[0]] == ethnicity:
                continue
            
            for iid in fam2empi[fid]:
                trait_value = empi2trait.get(iid, None)
                trait_ped.append( iid2ped[iid] + [trait_value] )
        
        if len(trait_ped) > 32000:
            print >> sys.stderr, "ERROR: Too many individuals."
            return False
        
        run_list[ethnicity] = [1,]
        
        for rid in run_list[ethnicity]:
            solar_ped = list()
            solar_phen = list()
            for row in trait_ped:
                famid, iid, fid, mid, sex, age, trait = row
                
                solar_ped.append( [famid, iid, fid, mid, sex] )
                if trait is None:
                    trait = ''
                    age = ''
                if use_proband:
                    solar_phen.append([iid, trait, age, 1 if iid == fam2proband[famid] else 0])
                else:
                    solar_phen.append([iid, trait, age])
            
            ethnicity_directory = os.path.join(h2_path, ethnicity)
            if not os.path.exists(ethnicity_directory):
                os.mkdir(ethnicity_directory)
            
            solar_working_path = os.path.join(ethnicity_directory, 'wd%d' % rid)
            
            if not os.path.exists(solar_working_path):
                os.mkdir(solar_working_path)
            
            ped_fh = open(os.path.join(solar_working_path, 'pedigree.ped'), 'w')
            writer = csv.writer(ped_fh, delimiter=',')
            writer.writerow(['famid', 'id', 'fa', 'mo', 'sex'])
            writer.writerows(solar_ped)
            ped_fh.close()
            
            phen_fh = open(os.path.join(solar_working_path, 'phenotypes.phen'), 'w')
            writer = csv.writer(phen_fh, delimiter=',', quoting=csv.QUOTE_NONE)
            if use_proband:
                writer.writerow(['id', 'pheno', 'age', 'proband'])
            else:
                writer.writerow(['id', 'pheno', 'age'])
            writer.writerows(solar_phen)
            phen_fh.close()
            
            # load_pedigree.tcl
            tcl_fh = open(os.path.join(solar_working_path, 'load_pedigree.tcl'), 'w')
            tcl_fh.write(tcl_load_string)
            tcl_fh.close()
            
            # run_analysis.tcl
            tcl_fh = open(os.path.join(solar_working_path, 'run_analysis.tcl'), 'w')
            if trait_type == TRAIT_TYPE_BINARY:
                tcl_fh.write(tcl_analysis_string_b)
            else:
                tcl_fh.writer(tcl_analysis_string_q)
            tcl_fh.close()
    
    return run_list

def single_solar_run(h2_path, run_list, verbose = True):
    
    ethnicity = 'ALL'
    for rid in run_list[ethnicity]:
        solar_working_path = os.path.join(h2_path, ethnicity, 'wd%d' % rid)
        
        if os.path.exists(os.path.join(solar_working_path, 'pheno', 'polygenic.out')):
            print >> sys.stderr, "Solar appears to already have been run. File exists at %s. Skipping." % os.path.join(solar_working_path, 'pheno', 'polygenic.out')
            continue
        
        if verbose:
            print >> sys.stderr, "%s..." % solar_working_path
        
        #print "cd %s && solar loadped > /dev/null && cd -" % solar_working_path
        os.system("cd %s > /dev/null && solar loadped > /dev/null && cd - > /dev/null " % solar_working_path)
        os.system("cd %s > /dev/null && solar runanalysis > /dev/null && cd - > /dev/null " % solar_working_path)
    
    if verbose:
        print >> sys.stderr, "ok."
    
    if verbose:
        print >> sys.stderr, "Parsing polygenic output to get h2 values..."
    
    h2s = list()
    h2errs = list()
    pvals = list()
    
    for rid in run_list[ethnicity]:
        
        polygenic_out_fn = os.path.join(h2_path, ethnicity, 'wd%d' % rid, 'pheno', 'polygenic.out')
        if not os.path.exists(polygenic_out_fn):
            print >> sys.stderr, "WARNING: No polygenic outfile found at: %s" % polygenic_out_fn
            continue
        
        polygenic_out_fh = open(polygenic_out_fn)

        results = polygenic_out_fh.read().split('\n')
        h2r_raw = [row.strip() for row in results if row.find('H2r') != -1]
        
        try:
            h2r = float(h2r_raw[0].split()[2])
            p = float(h2r_raw[0].split()[5])
        except IndexError:
            if verbose:
                print >> sys.stderr, "FAILED. SOLAR failed to run. Could be a convergence error."
            h2r = None
            p = None
        
        h2r_err = None
        
        try:
            h2r_err = float(h2r_raw[1].split()[3])
        except IndexError:
            h2r_err = None
        
        h2s.append(h2r)
        h2errs.append(h2r_err)
        pvals.append(p)
        
    return h2s, h2errs, pvals

if __name__ == '__main__':
    print >> sys.stderr, "This script is not intended to be run directly."
    sys.exit(100)