#!/usr/bin/env python

"""
solarStrap_heritability.py

Compute heritability for one or more traits.

@author Nicholas P. Tatonetti, Rami Vanguri

USAGE EXAMPLE
-------------
python solarStrap_heritability.py trait=path/to/traitfile.txt.gz type=D demog=demogtable.txt.gz fam=family_ids.txt.gz ped=generic_pedigree.txt.gz [sd=path/to/working_directory]

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

from familial_recurrence import assign_family_ethnicities
from familial_recurrence import load_demographics
from familial_recurrence import load_family_ids
from familial_recurrence import load_relationships

from solar_heritability import build_solar_directories
from solar_heritability import single_solar_run

from solar_setup_h2_analysis import load_generic_pedigree

common_data_path = ''

def solar(input_args):
    h2_path, families_with_case, icd9, \
        num_families, iid2ped, all_traits, ethnicity, \
        min_ascertained, fam2empi, fam2eth, all_fam2count, \
        all_fam2proband, use_proband, binarytrait = input_args
    
    if os.path.exists(h2_path):
        shutil.rmtree(h2_path)
    
    if binarytrait:
        trait_type = 1
    else:
        trait_type = 0
    
    chosen_families = random.sample(families_with_case[icd9], num_families)
    run_list = build_solar_directories(h2_path,
                                       iid2ped,
                                       all_traits[icd9],
                                       [ethnicity],
                                       min_ascertained,
                                       fam2empi,
                                       fam2eth,
                                       all_fam2count[icd9],
                                       all_fam2proband[icd9],
                                       use_proband,
                                       verbose=False,
                                       family_ids_only = chosen_families)
    
    results = single_solar_run(h2_path, run_list, ethnicity, verbose=False)
    if len(results[0]) != 0:
        h2, h2err, pval = zip(*results)[0]
    else:
        h2, h2err, pval = None, None, None

    return (h2, h2err, pval)

def random_string(N):
    """http://stackoverflow.com/questions/2257441/random-string-generation-with-upper-case-letters-and-digits-in-python/23728630#23728630"""
    return ''.join(random.choice(string.ascii_uppercase + string.digits) for _ in range(N))

def solar_strap(input_args):
    """ 
    Runs the solar strap procedure. The arguments are accepted like this to make it easier for parallel processing.
    """
    
    num_families, families_with_case, icd9, _max_significant_, _max_attempts_, \
        solar_dir, iid2ped, all_traits, ethnicity, min_ascertained, fam2empi, \
        fam2eth, all_fam2count, all_fam2proband, use_proband, verbose, binarytrait = input_args
    
    if num_families > len(families_with_case[icd9]):
        print >> sys.stderr, "Number of families large than what is available."
        return num_families, list(), 0, 0
    
    h2s = list()
    num_attempts = 0
    num_successes = 0
    num_significant = 0
    
    while num_significant < _max_significant_ and num_attempts < _max_attempts_:
        
        num_attempts += 1
        
        start_time = time.time()
        
        h2_path = os.path.join(solar_dir, icd9, 'h2_%d_%s' % (num_families, random_string(5)))
        h2, h2err, pval = solar((h2_path,
                        families_with_case,
                        icd9,
                        num_families,
                        iid2ped,
                        all_traits,
                        ethnicity,
                        min_ascertained,
                        fam2empi,
                        fam2eth,
                        all_fam2count,
                        all_fam2proband,
                        use_proband,
                        binarytrait))
        
        shutil.rmtree(h2_path)
        run_time = time.time()-start_time
        
        ####
        # NOTE: These hyperparameters are set by "Heritability - SOLARStrap - 84phenos - Hyperparamters"
        pcutoff = 0.05
        edge_eps = 1e-9
        denoise_eps = 0.05
        
        if h2 is None or h2err is None or float(h2) < edge_eps or float(h2) > (1-edge_eps):
            print "%20s %7d %10s %10s %10s %5d %5d %7.2f %7.2fs" % (ethnicity, num_families, h2, h2err, pval, num_significant, num_successes, num_significant/max(1.0,float(num_successes)), run_time)
            continue
        
        if float(h2err) < denoise_eps*float(h2):
            continue
        
        num_successes += 1
        
        if pval < pcutoff:
            num_significant += 1
            h2s.append(h2)
            if verbose:
                print "%20s %7d %10s %10s %10s %5d %5d %7.2f %7.2fs" % (ethnicity, num_families, h2, h2err, pval, num_significant, num_successes, num_significant/float(num_successes), run_time)
    
    return num_families, h2s, num_significant, num_successes

def prevelance(trait):
    n_affected = 0
    n_unaffected = 0
    for iid, value in trait.items():
        if value == 1:
            n_affected += 1
        elif value == 0:
            n_unaffected += 1
    prevelance = n_affected/float(n_unaffected+n_affected)
    return n_affected, n_unaffected, prevelance

def main(demographic_file, family_file, pedigree_file, trait_path, solar_dir, trait_type, diag_slice=None, ethnicities=None):
    
    if trait_type == 'D':
        use_proband = True
        min_ascertained = 1
        _DT_ = True
        _QT_ = False
    elif trait_type == 'Q':
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
        print >> sys.stderr, "%10s %13s %13s" % ("Trait", "N Controls", "N Cases")
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
        start = diag_slice[0]
        stop = diag_slice[1]
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
    generic_pedigree = load_generic_pedigree(generic_ped_path, empi2sex, empi2age)
    
    iid2ped = dict()
    for row in generic_pedigree:
        iid2ped[row[1]] = row
    del(generic_pedigree)
    
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
    
    _max_significant_ = 200
    _max_attempts_ = 200
    
    #step_factor = 1 # every 100, 1000, 10000
    #step_factor = 2 # every 200, 2000, 20000
    step_factor = 3 # every 300, 3000, 30000
    
    num_families_range = (i*10**exp for exp in range(2, 4) for i in range(1, 11, step_factor))
    
    for num_families in num_families_range:
        for icd9 in diags_to_process:
            
            print >> sys.stderr, "Running solarStrap analysis for %s" % icd9
        
            icd9_path = os.path.join(solar_dir, icd9)
            if not os.path.exists(icd9_path):
                os.mkdir(icd9_path)
            
            h2_path = os.path.join(solar_dir, icd9, 'h2')
            
            print >> sys.stderr, "Number of families with case: %d" % (len(families_with_case[icd9]))
            
            for eth in ethnicities:
                
                nf, h2values, num_significant, num_successes = solar_strap((num_families,
                                                                            families_with_case,
                                                                            icd9,
                                                                            _max_significant_,
                                                                            _max_attempts_,
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
                                                                            True,
                                                                            _DT_))
                if len(h2values) >= 3:
                    medi = int(numpy.floor(len(h2values)/2))
                    print >> sys.stdout, "SUMMARY: ", icd9, eth, nf, num_significant, num_successes, h2values[medi], numpy.percentile(h2values, 2.5), numpy.percentile(h2values, 97.5)

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
        ethnicities = None if not 'eth' in args else args['eth'])
