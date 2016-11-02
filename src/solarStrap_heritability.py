#!/usr/bin/env python

"""
solarStrap_heritability.py

Compute heritability for one or more traits.

@author Nicholas P. Tatonetti, Rami Vanguri

USAGE EXAMPLE
-------------
python solarStrap_heritability.py trait=path/to/traitfile.txt.gz type=D demog=demogtable.txt.gz fam=family_ids.txt.gz ped=generic_pedigree.txt.gz [sd=path/to/working_directory] [ace=yes] [verbose=yes]

"""

import os
import sys
import csv
import time
import gzip
import numpy
import string
import random
import shutil
from collections import defaultdict

from tqdm import tqdm

from h2o_utility import assign_family_ethnicities
from h2o_utility import load_demographics
from h2o_utility import load_family_ids
from h2o_utility import load_relationships
from h2o_utility import load_generic_pedigree

from h2o_utility import build_solar_directories
from h2o_utility import single_solar_run
from h2o_utility import solar
from h2o_utility import random_string
from h2o_utility import solar_strap
from h2o_utility import prevelance
from h2o_utility import estimate_h2o

from h2o_utility import TRAIT_TYPE_BINARY
from h2o_utility import TRAIT_TYPE_QUANTITATIVE

common_data_path = ''

def main(demographic_file, family_file, pedigree_file, trait_path, solar_dir, trait_type, diag_slice=None, ethnicities=None, verbose=False, house=False):
    
    if trait_type == 'D':
        trait_type_code = TRAIT_TYPE_BINARY
        use_proband = True
        min_ascertained = 1
        _DT_ = True
        _QT_ = False
    elif trait_type == 'Q':
        trait_type_code = TRAIT_TYPE_QUANTITATIVE
        use_proband = False
        min_ascertained = 2
        _DT_ = False
        _QT_ = True
    else:
        raise Exception("Unknown trait type provided (%s). Must be either D or Q." % trait_type)
    
    if not os.path.exists(solar_dir):
        raise Exception("A directory must exist at: %s" % solar_dir)
    
    # load family ids
    empi2fam, fam2empi = load_family_ids(os.path.join(common_data_path, family_file))
    print >> sys.stderr, "(%d families loaded)" % len(fam2empi)
    
    # load demographic data (sex, birth_decade, race, ethnicity, age)
    demog_file_path = os.path.join(common_data_path, demographic_file)
    empi2demog = load_demographics(demog_file_path)
    
    # set up some convenience dictionaries
    empi2sex = dict()
    empi2age = dict()
    for empi, data in empi2demog.items():
        empi2sex[empi] = data['sex']
        empi2age[empi] = data['age'] if not data['age'] == 'NULL' else ''
    
    print >> sys.stderr, "Loading phenotype (trait) data..."
    all_traits = defaultdict(dict)
    
    # Load the trait data
    fh = gzip.open(trait_path)
    reader = csv.reader(fh, delimiter='\t')
    reader.next()
    
    for empi, icd9, value in reader:
        if empi in empi2fam:
            all_traits[icd9][empi] = float(value)
    
    available_diagnoses = sorted(all_traits.keys())
    diags_to_process = available_diagnoses
    
    print >> sys.stderr, "Loaded data for %d phenotypes." % len(all_traits.keys())
    
    if _DT_:
        print >> sys.stderr, "%10s %13s %13s" % ("Trait", "N Cases", "N Controls")
        for trait in all_traits.keys():
            ncases = sum(all_traits[trait].values())
            nctrls = len(all_traits[trait])-ncases
            print >> sys.stderr, "%10s %13d %13d" % (trait, ncases, nctrls)
            if nctrls == 0:
                print >> sys.stderr, "Error: There are no controls available for the dichotomous trait."
    if _QT_:
        print >> sys.stderr, "%10s %13s" % ("Trait", "N Samples")
        for trait in all_traits.keys():
            ncases = len(all_traits[trait])
            print >> sys.stderr, "%10s %13d" % (trait, ncases)
    
    if not diag_slice is None:
        if len(diag_slice) == 1:
            start = diag_slice[0]
            stop = diag_slice[0] + 1
        elif len(diag_slice) == 2:
            start = diag_slice[0]
            stop = diag_slice[1]
        else:
            raise Exception("ERROR: Slice provided in unexpected format: %s" % diag_slice)
        
        diags_to_process = available_diagnoses[start:stop]
    
    print >> sys.stderr, "Loading ascertainment counts..."
    
    all_fam2count = defaultdict(dict)
    all_fam2case = defaultdict(dict)
    families_with_case = defaultdict(set)
    
    for fam_id, members in fam2empi.items():
        for icd9 in diags_to_process:
            all_fam2count[icd9][fam_id] = sum([1 for e in members if e in all_traits[icd9]])
            all_fam2case[icd9][fam_id] = sum([all_traits[icd9][e] for e in members if e in all_traits[icd9]])
            if _DT_ and all_fam2case[icd9][fam_id] >= 1:
                families_with_case[icd9].add(fam_id)
            if _QT_ and all_fam2count[icd9][fam_id] >= min_ascertained:
                families_with_case[icd9].add(fam_id)
    
    all_fam2proband = defaultdict(dict)
    if _DT_:
        print >> sys.stderr, "Assigning a proband for each family for each trait..."
        for icd9 in diags_to_process:
            for iid, trait in all_traits[icd9].items():
                if trait == 1.0:
                    if not empi2fam[iid] in all_fam2proband[icd9]:
                        all_fam2proband[icd9][empi2fam[iid]] = iid
        
        print >> sys.stderr, "Checking proband assignments..."
        for icd9 in diags_to_process:
            for famid in fam2empi.keys():
                if all_fam2case[icd9][famid] > 0 and not famid in all_fam2proband[icd9]:
                    raise Exception("Family: %s did not have a proband!" % famid)
    
    generic_ped_path = os.path.join(common_data_path, pedigree_file)
    iid2ped = load_generic_pedigree(generic_ped_path, empi2sex, empi2age)
    eth2fam, fam2eth = assign_family_ethnicities(fam2empi, empi2demog)
    
    if ethnicities is None:
        ethnicities = ['ALL']
    elif ethnicities == 'each':
        ethnicities = ['ALL'] + eth2fam.keys()
    
    print >> sys.stderr, "Evaluating %d ethnicities: %s" % (len(ethnicities), ', '.join(ethnicities))
    
    if _DT_:
        min_ascertained = 1
    if _QT_:
        min_ascertained = 2
    
    num_attempts = 200
    
    #step_factor = 1 # every 100, 1000, 10000
    #step_factor = 2 # every 200, 2000, 20000
    step_factor = 3 # every 300, 3000, 30000
    
    num_families_range = (i*10**exp for exp in range(2, 4) for i in range(1, 11, step_factor))
    
    results_file = open(os.path.join(solar_dir, 'solar_strap_results.csv'), 'w')
    results_writer = csv.writer(results_file)

    results_writer.writerow(['trait', 'ethnicity', 'num_families', 'model', 'h2o', 'h2o_lower', 'h2o_upper', 'solarerr', 'solarpval', 'num_attempts', 'num_converged', 'num_significant', 'posa'])
    
    for num_families in num_families_range:
        for icd9 in diags_to_process:
            
            print >> sys.stderr, "Running solarStrap analysis for %s, num_fam = %d" % (icd9, num_families)
            print >> sys.stderr, " AE: yes, ACE: %s" % ('yes' if house else 'no')
            
            icd9_path = os.path.join(solar_dir, icd9)
            if not os.path.exists(icd9_path):
                os.mkdir(icd9_path)
            
            h2_path = os.path.join(solar_dir, icd9, 'h2')
            
            print >> sys.stderr, "Number of families with case: %d" % (len(families_with_case[icd9]))
            
            if num_families > len(families_with_case[icd9]):
                print >> sys.stderr, "Not enough families available, skipping."
                continue
            
            for eth in ethnicities:
                
                ae_h2r_results, ace_h2r_results = solar_strap(num_families,
                                                              families_with_case,
                                                              icd9,
                                                              trait_type_code,
                                                              num_attempts,
                                                              solar_dir,
                                                              iid2ped,
                                                              all_traits,
                                                              eth,
                                                              min_ascertained,
                                                              fam2empi,
                                                              fam2eth,
                                                              all_fam2count,
                                                              all_fam2proband,
                                                              use_proband,
                                                              house,
                                                              verbose)
                estimates = estimate_h2o(ae_h2r_results)
                if estimates:
                    h2o, h2olo, h2ohi, solarerr, solarpval, num_converged, num_significant, posa = estimates
                    print "%10s %10s %5d %4s %7.2f %7.2f %7.2e %7d %7d %7d %7.2f" % (icd9, eth, num_families, 'AE', h2o, solarerr, solarpval, num_attempts, num_converged, num_significant, posa)
                    results_writer.writerow([icd9, eth, num_families, 'AE', h2o, h2olo, h2ohi, solarerr, solarpval, num_attempts, num_converged, num_significant, posa])
                
                if house:
                    estimates = estimate_h2o(ace_h2r_results)
                    if estimates:
                        h2o, h2olo, h2ohi, solarerr, solarpval, num_converged, num_significant, posa = estimates
                        print "%10s %10s %5d %4s %7.2f %7.2f %7.2e %7d %7d %7d %7.2f" % (icd9, eth, num_families, 'ACE', h2o, solarerr, solarpval, num_attempts, num_converged, num_significant, posa)
                        results_writer.writerow([icd9, eth, num_families, 'AE', h2o, h2olo, h2ohi, solarerr, solarpval, num_attempts, num_converged, num_significant, posa])

    results_file.close()

if __name__ == '__main__':
    args = dict([x.split('=') for x in sys.argv[1:]])
    
    if not 'sd' in args:
        args['sd'] = './working'
    
    main(demographic_file = args['demog'],
        family_file = args['fam'],
        pedigree_file = args['ped'],
        trait_path = args['trait'],
        solar_dir = args['sd'],
        trait_type = args['type'],
        diag_slice = None if not 'slice' in args else map(int, args['slice'].split('-')),
        ethnicities = None if not 'eth' in args else args['eth'],
        verbose = False if not 'verbose' in args else args['verbose'].lower() == 'yes',
        house = False if not 'ace' in args else args['ace'].lower() == 'yes')
