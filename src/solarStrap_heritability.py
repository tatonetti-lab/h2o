#!/usr/bin/env python

"""
solarStrap_heritability.py

Compute heritability for one or more traits.

@author Nicholas P. Tatonetti, Rami Vanguri

USAGE EXAMPLE
-------------
python solarStrap_heritability.py trait=path/to/traitfile.txt.gz type=D demog=demogtable.txt.gz fam=family_ids.txt.gz ped=generic_pedigree.txt.gz [nfam=0.15] [sd=path/to/working_directory] [ace=no] [verbose=no] [samples=200] [buildonly=no] [proband=yes]

"""
__version__ = 0.9


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

def main(demographic_file, family_file, pedigree_file, trait_path, solar_dir, trait_type, num_families_range, diag_slice=None, ethnicities=None, verbose=False, house=False, prefix='', nprocs=1, num_attempts=200, buildonly=False, use_proband=True):

    if trait_type == 'D':
        trait_type_code = TRAIT_TYPE_BINARY
        _DT_ = True
        _QT_ = False
    elif trait_type == 'Q':
        trait_type_code = TRAIT_TYPE_QUANTITATIVE
        if use_proband:
            print >> sys.stderr, "Setting use_proband=False for quantitative trait."
        use_proband = False
        _DT_ = False
        _QT_ = True
    else:
        raise Exception("Unknown trait type provided (%s). Must be either D or Q." % trait_type)

    print >> sys.stderr, "Use proband? %s " % ('yes' if use_proband else 'no')

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
    diag2idx = dict([(diag, idx) for idx, diag in enumerate(diags_to_process)])
    print >> sys.stderr, "Loaded data for %d phenotypes." % len(all_traits.keys())

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
            # There are two scenarios for deciding which families to sample from:
            # 1) The preferred one for observational heritability is that only families
            # with at least one individual with the trait (disease) are included. We do
            # this because we have more confidence in the ascertainment of cases than
            # controls. We then control for having only families with traits  by setting
            # a "proband" in each family. Solar then will use the proband to adjust the
            # heritability estimate.
            #
            # 2) The second scenario is that we include any family that has at least two
            # members with ascertained traits (either case or control). You want to use this
            # if you are very confidence in your ascertainment of controls. This is usually
            # not the case if your data are coming from observational sources. In this
            # case we leave out the proband assigment.

            if not use_proband and all_fam2count[icd9][fam_id] >= 2:
                # quantitative traits must have least two individuals with measured traits
                # or dichotomous traits with high confidence controls
                families_with_case[icd9].add(fam_id)
            elif _DT_ and all_fam2case[icd9][fam_id] >= 1 and all_fam2count[icd9][fam_id] >= 2:
                # dichotomous traits with noisy/uncertain controls
                families_with_case[icd9].add(fam_id)

    if _DT_:
        fmt_string_h = "%10s %" + str(max(map(len, diags_to_process))+1) + "s %10s %10s %10s %10s %10s"
        fmt_string_r = "%10s %" + str(max(map(len, diags_to_process))+1) + "s %10d %10d %10d %10.4f %10.4f"

        print >> sys.stderr, fmt_string_h % ('Trait Idx', 'Trait', 'N Ascert', 'N Cases', 'N Families', 'Avg Ctl/Fam', 'Avg Aff/Fam', )
        summary_stats = list()
        for idx, icd9 in enumerate(diags_to_process):
            apf = numpy.mean([all_fam2case[icd9][famid] for famid in families_with_case[icd9]])
            cpf = numpy.mean([(all_fam2count[icd9][famid]-all_fam2case[icd9][famid]) for famid in families_with_case[icd9]])

            summary_stats.append( [apf, diag2idx[icd9], icd9, len(all_traits[icd9]), sum(all_traits[icd9].values()), len(families_with_case[icd9]), cpf] )

        summary_stats = sorted(summary_stats)
        for row in summary_stats:
            args = row[1:]
            args.append(row[0])
            print >> sys.stderr, fmt_string_r % tuple(args)

    if _QT_:
        print >> sys.stderr, "%10s %13s %10s" % ("Trait", "N Samples", "Avg Aff/Fam")
        for trait in all_traits.keys():
            apf = numpy.mean([all_fam2count[icd9][famid] for famid in families_with_case[icd9]])
            ncases = len(all_traits[trait])
            print >> sys.stderr, "%10s %13d %10.4f" % (trait, ncases, apf)


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

    results_file = open(os.path.join(solar_dir, '%s_solar_strap_results.csv' % prefix), 'w')
    results_writer = csv.writer(results_file)
    results_writer.writerow(['trait', 'ethnicity', 'num_families', 'model', 'h2o', 'h2o_lower', 'h2o_upper', 'solarerr', 'solarpval', 'num_attempts', 'num_converged', 'num_significant', 'posa'])

    runs_file = open(os.path.join(solar_dir, '%s_solar_strap_allruns.csv' % prefix), 'w')
    runs_writer = csv.writer(runs_file)
    runs_writer.writerow(['trait', 'ethnicity', 'num_families', 'model', 'h2o', 'solarerr', 'pvalue'])

    if num_families_range is None:
        num_families_range = [0.15,]

    for num_families in num_families_range:
        for icd9 in diags_to_process:

            if type(num_families) == float and num_families < 1:
                num_families = max(500, int(num_families*len(families_with_case[icd9])))
            else:
                num_familes = int(num_families)

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
                                                              fam2empi,
                                                              fam2eth,
                                                              all_fam2count,
                                                              all_fam2proband,
                                                              use_proband,
                                                              house,
                                                              nprocs,
                                                              verbose,
                                                              buildonly)
                for h2o, err, pval in ae_h2r_results:
                    runs_writer.writerow([icd9, eth, num_families, 'AE', h2o, err, pval])

                estimates = estimate_h2o(ae_h2r_results)
                if estimates:
                    h2o, h2olo, h2ohi, solarerr, solarpval, num_converged, num_significant, posa = estimates
                    print "%10s %10s %5d %4s %7.2f %7.2f %7.2e %7d %7d %7d %7.2f" % (icd9, eth, num_families, 'AE', h2o, solarerr, solarpval, num_attempts, num_converged, num_significant, posa)
                    results_writer.writerow([icd9, eth, num_families, 'AE', h2o, h2olo, h2ohi, solarerr, solarpval, num_attempts, num_converged, num_significant, posa])

                if house:
                    for h2o, err, pval in ace_h2r_results:
                        runs_writer.writerow([icd9, eth, num_families, 'ACE', h2o, err, pval])

                    estimates = estimate_h2o(ace_h2r_results)
                    if estimates:
                        h2o, h2olo, h2ohi, solarerr, solarpval, num_converged, num_significant, posa = estimates
                        print "%10s %10s %5d %4s %7.2f %7.2f %7.2e %7d %7d %7d %7.2f" % (icd9, eth, num_families, 'ACE', h2o, solarerr, solarpval, num_attempts, num_converged, num_significant, posa)
                        results_writer.writerow([icd9, eth, num_families, 'ACE', h2o, h2olo, h2ohi, solarerr, solarpval, num_attempts, num_converged, num_significant, posa])

    results_file.close()
    runs_file.close()

if __name__ == '__main__':
    args = dict([x.split('=') for x in sys.argv[1:]])

    if not 'sd' in args:
        args['sd'] = './working'

    if not 'name' in args:
        args['name'] = random_string(5)

    print >> sys.stderr, ""
    print >> sys.stderr, "SolarStrap v%4.1f - Estimate heritability of disease using observational data." % __version__
    print >> sys.stderr, "-----------------------------------------------------------------------------"
    print >> sys.stderr, "Summary results will be saved in %(sd)s/%(name)s_solar_strap_results.csv" % args
    print >> sys.stderr, "Results from each bootstrap will be saved at %(sd)s/%(name)s_solar_strap_allruns.csv" % args
    print >> sys.stderr, ""

    main(demographic_file = args['demog'],
        family_file = args['fam'],
        pedigree_file = args['ped'],
        trait_path = args['trait'],
        solar_dir = args['sd'],
        trait_type = args['type'],
        num_families_range = None if not 'nfam' in args else map(float, args['nfam'].split(',')),
        diag_slice = None if not 'slice' in args else map(int, args['slice'].split('-')),
        ethnicities = None if not 'eth' in args else args['eth'],
        verbose = False if not 'verbose' in args else args['verbose'].lower() == 'yes',
        house = False if not 'ace' in args else args['ace'].lower() == 'yes',
        prefix = args['name'],
        nprocs = 1 if not 'nprocs' in args else int(args['nprocs']),
        num_attempts = 200 if not 'samples' in args else int(args['samples']),
        buildonly = False if not 'buildonly' in args else args['buildonly'].lower() == 'yes',
        use_proband = True if not 'proband' in args else args['proband'].lower() == 'yes')
