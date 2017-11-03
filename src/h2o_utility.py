"""
Utility functions for loading data, setting up directories, and running solar.

"""

import os
import sys
import csv
import gzip
import numpy
import shutil
import string
import random

from tqdm import tqdm
from collections import defaultdict

import multiprocessing as mp
import time

# Define some constants
TRAIT_TYPE_QUANTITATIVE = 0
TRAIT_TYPE_BINARY = 1

# the typical id used for household when it is not directly available is
# the mother id (mo)
# we also explore the use of the family id
tcl_load_string = """
proc loadped {} {
    field hhid mo
    load pedigree pedigree.ped
    load phenotypes phenotypes.phen
}
"""

tcl_analysis_string_b = """
proc runanalysis {} {
    model new
    trait pheno
    covariates sex age
    polygenic -screen
}
proc runanalysishouse {} {
    model new
    trait pheno
    covariates sex age
    house
    polygenic -screen
}
"""

tcl_analysis_string_q = """
proc runanalysis {} {
    trait pheno
    covariates sex age
    tdist
    polygenic -screen
}
proc runanalysishouse {} {
    model new
    trait pheno
    covariates sex age
    house
    tdist
    polygenic -screen
}
"""

def blood_proportion(rel):
    return float(blood_fraction(rel))/100.

def blood_fraction(r):
    if r == 'Aunt/Uncle':
        return 25
    elif r == 'Child':
        return 50
    elif r == 'Cousin':
        return 12.5
    elif r == 'First cousin once removed':
        return 6.25
    elif r == 'Grandaunt/Granduncle':
        return 12.5
    elif r == 'Grandchild':
        return 25
    elif r == 'Grandnephew/Grandniece':
        return 12.5
    elif r == 'Grandparent':
        return 25
    elif r == 'Great-grandaunt/Great-granduncle':
        return 6.25
    elif r == 'Great-grandchild':
        return 12.5
    elif r == 'Great-grandnephew/Great-grandniece':
        return 6.25
    elif r == 'Great-grandparent':
        return 12.5
    elif r == 'Great-great-grandchild':
        return 6.25
    elif r == 'Great-great-grandparent':
        return 6.25
    elif r == 'Nephew/Niece':
        return 25
    elif r == 'Parent':
        return 50
    elif r == 'Sibling':
        return 50
    elif r == 'Spouse':
        return 0
    elif r == 'Child/Nephew/Niece':
        return 25
    elif r == 'Parent/Aunt/Uncle':
        return 25
    elif r == 'Sibling/Cousin':
        return 12.5
    else:
        #rint >> sys.stderr, "Encountered relationshp not found: %s" % r
        return 0

def random_string(N):
    return ''.join(random.choice(string.ascii_uppercase + string.digits) for _ in range(N))

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

def assign_family_ethnicities(fam2empi, empi2demog, print_breakdown=True):
    """
    Assign each family a single ethnicity. This makes the (incorrect) assumption that everyone in the same
    family will have the same ethnicity. However, because there are so much missing data, this is better
    than nothing at all.

    fam2empi        dict, map from family id to a set of patients identifiers that belong in that family
    empi2demog      dict, map from patient identifiers to the demographic data (including self-reported race/ethnicity data)
    """

    print >> sys.stderr, "Identifying most commonly reported ethnicity for each family..."

    eth2fam = defaultdict(set)
    fam2eth = dict()
    for fid, members in tqdm(fam2empi.items()):
        family_ethnicities = [empi2demog[e]['race'] for e in members if empi2demog[e]['race'] != 'Unknown']
        if len(family_ethnicities) == 0:
            mode = 'Unknown'
        else:
            mode = max(set(family_ethnicities), key=family_ethnicities.count)

        eth2fam[mode].add(fid)
        fam2eth[fid] = mode

    if print_breakdown:
        for eth, families in eth2fam.items():
            print >> sys.stderr, "\t%s: %d" % (eth, len(families))

    return eth2fam, fam2eth

def load_demographics(demographic_file_path):
    """
    Load the demographic patient data from file.

    """

    print >> sys.stderr, "Loading patient demographic data..."

    # load demographic data (sex, birth_decade, race, ethnicity, age)
    fh = gzip.open(demographic_file_path)
    reader = csv.reader(fh, delimiter='\t')
    reader.next()
    demog_header = ['sex', 'birthdec', 'race', 'ethnicity', 'age']

    empi2demog = dict()
    for row in tqdm(reader):

        empi = row[0]
        empi2demog[empi] = dict(zip(demog_header, row[1:]))

        # re-map the race codes
        if empi2demog[empi]['race'] == 'NA':
            if empi2demog[empi]['ethnicity'] == 'W':
                empi2demog[empi]['race'] = 'White'
            elif empi2demog[empi]['ethnicity'] == 'B':
                empi2demog[empi]['race'] = 'Black'
            elif empi2demog[empi]['ethnicity'] in ('H', 'S', 'O'):
                empi2demog[empi]['race'] = 'Hispanic'
            elif empi2demog[empi]['ethnicity'] in ('', 'U', 'D', '2'):
                empi2demog[empi]['race'] = 'Unknown'
            else:
                empi2demog[empi]['race'] = 'Other'
        else:
            if empi2demog[empi]['race'] == 'W':
                if empi2demog[empi]['ethnicity'] == 'H':
                    empi2demog[empi]['race'] = 'Hispanic'
                else:
                    empi2demog[empi]['race'] = 'White'
            elif empi2demog[empi]['race'] == 'B':
                if empi2demog[empi]['ethnicity'] == 'H':
                    empi2demog[empi]['race'] = 'Hispanic'
                else:
                    empi2demog[empi]['race'] = 'Black'
            elif empi2demog[empi]['race'] == 'O':
                empi2demog[empi]['race'] = 'Hispanic'
            elif empi2demog[empi]['race'] in ('U','D'):
                if empi2demog[empi]['ethnicity'] == 'H':
                    empi2demog[empi]['race'] = 'Hispanic'
                else:
                    empi2demog[empi]['race'] = 'Unknown'
            else:
                empi2demog[empi]['race'] = 'Other'

        # recast age as float
        empi2demog[empi]['age'] = float(empi2demog[empi]['age']) if empi2demog[empi]['age'] != 'NULL' else numpy.nan

    print >> sys.stderr, "Loaded demographic data for %d patients." % len(empi2demog)

    return empi2demog

def load_family_ids(family_file_path):
    """
    Load family identifiers form file.
    """

    print >> sys.stderr, "Loading family identifiers..."

    fh = gzip.open(family_file_path)
    reader = csv.reader(fh, delimiter='\t')
    reader.next()
    empi2fam = dict()
    fam2empi = defaultdict(set)
    for fid, iid in reader:
        empi2fam[iid] = fid
        fam2empi[fid].add( iid )
    fam2empi = dict(fam2empi)

    return empi2fam, fam2empi

def load_relationships_by_pairs(relationships_file):
    # load the relationships
    fh = gzip.open(os.path.join(relationships_file))
    rel_data = csv.reader(fh, delimiter='\t')
    rel_data.next()
    
    relationships = defaultdict(dict)
    
    for empi1, rel, empi2, provided in rel_data:
        rel = rel.strip()
        relationships[empi1][empi2] = blood_proportion(rel)
    
    return relationships

def load_relationships(relationship_file_path, print_breakdown=True):
    """
    Load the patient relationships from file.
    """

    # load the relationships
    print >> sys.stderr, "Loading the relationships..."
    fh = gzip.open(relationship_file_path)
    rel_data = csv.reader(fh, delimiter='\t')
    rel_data.next()

    relationships = defaultdict(lambda: defaultdict(set))

    for empi1, rel, empi2, provied in tqdm(rel_data):
        rel = rel.strip().lower()
        if rel.find('in-law') != -1 or rel.find('spouse') != -1:
            # we explicitly exclude any in-law relationships
            continue

        relationships[rel][empi1].add(empi2)
        relationships['All'][empi1].add(empi2)

    if print_breakdown:
        print >> sys.stderr, "Relationship breakdown: "
        for rel, rels in relationships.items():
            nrels = sum([len(relatives) for relatives in rels.values()])
            print >> sys.stderr, "%50s %10d" % (rel, nrels)

    return relationships

def load_generic_pedigree(generic_ped_path, empi2sex, empi2age):
    """
    Load the pedigree file from the path.
    generic_ped_path    string, path to pedigree file
    empi2sex            dict, map from patient identifier to sex (M|F)
    empi2age            dict, map from patient identifier to age (int)
    """

    print >> sys.stderr, "Loading generic pedigree file..."

    # load the generic pedigree file
    fh = gzip.open(generic_ped_path)
    reader = csv.reader(fh, delimiter='\t')
    header = reader.next()
    ped_data = list()
    own_ancestor = set()
    conflicted_fams = set()
    for fam_id, ind_id, fat_id, mot_id, own_a in tqdm(reader):
        if own_a == '1':
            own_ancestor.add(ind_id)
            conflicted_fams.add(fam_id)

        if empi2sex.get(fat_id, None) is None:
            empi2sex[fat_id] = 'M'
        if empi2sex.get(mot_id, None) is None:
            empi2sex[mot_id] = 'F'

        ped_data.append((fam_id, ind_id, fat_id, mot_id))

    print >> sys.stderr, "ok. (%d individuals in pedigree)" % len(ped_data)

    print >> sys.stderr, "Building pedigree. Found %d ancestor conflicts..." % len(own_ancestor)

    generic_pedigree = list()
    next_id = 1
    for fam_id, ind_id, fat_id, mot_id in tqdm(ped_data):
        #if fam_id in conflicted_fams:
        #    continue

        if ind_id in own_ancestor:
            continue

        if fat_id in own_ancestor:
            fat_id = 0
        if mot_id in own_ancestor:
            mot_id = 0

        if (fat_id == 0 or mot_id == 0) and not (fat_id == mot_id == 0):
            fat_id = 'S%d' % next_id
            next_id += 1
            mot_id = 'S%d' % next_id
            next_id += 1
            generic_pedigree.append([fam_id, fat_id, 0, 0, 1, ''])
            generic_pedigree.append([fam_id, mot_id, 0, 0, 2, ''])

        sexcode = empi2sex.get(ind_id, None)
        if sexcode is None:
            sex = 0
        elif sexcode == 'M':
            sex = 1
        elif sexcode == 'F':
            sex = 2
        else:
            sex = 0

        age = empi2age.get(ind_id, '')

        generic_pedigree.append([fam_id, ind_id, fat_id, mot_id, sex, age])

    iid2ped = dict()
    for row in generic_pedigree:
        iid2ped[row[1]] = row
    del(generic_pedigree)

    return iid2ped

class SolarException(Exception):
    pass

def build_gcta_directories(h2_path, empi2demog, empi2trait, fam2empi, fam2count, relationships, trait_type, verbose = True, family_ids_only = None):
    """
    Build the directoies with the files need to run GCTA a single time.
    
    GCTA requires four files:
    1. *.phen file of the phenotyes, columns are family ID, individual ID and phenotypes
    2. *.covar file of the covariates, columns are family ID, individual ID and discrete covariates
    3. *.qcovar file of the quantitative covariates, columns are family ID, individual ID and quantitative covariates
    4. *.grm.gz file of the genetic relatedness between all pairs of indivdiuals, 
                columns are indices of pairs of individuals (row numbers of the test.grm.id), number of non-missing SNPs and the estimate of genetic relatedness
    
    """
    
    if not os.path.exists(h2_path):
        os.mkdir(h2_path)
    
    if verbose:
        print >> sys.stderr, "Building files at path: %s" % h2_path
    
    if os.path.exists(os.path.join(h2_path, 'wd')):
        raise Exception("Directory already exists at path: %s" % os.path.join(h2_path, 'wd'))
    
    if verbose:
        print >> sys.stderr, "Building trait-specfic pedigree file...",
    
    patient_data = dict()
    patients = set()
    
    for famid in family_ids_only:
        if fam2count[famid] < 2:
            continue
        
        for pid in fam2empi[famid]:
            if trait_type == TRAIT_TYPE_BINARY:
                trait_value = empi2trait.get(pid, None)
                trait_value = int(trait_value) if not trait_value is None else trait_value
            else:
                trait_value = empi2trait.get(pid, None)
            #if trait_value is None:
            #    continue
            
            patient_data[pid] = [famid, pid, empi2demog[pid]['sex'], empi2demog[pid]['age'], trait_value]
            patients.add( pid )
    
    patients = sorted(patients)
    pat2index = dict([(pid,i+1) for i,pid in enumerate(patients)])
    pat2patid = dict([(pid,i+10) for i, pid in enumerate(patients)])
    pat2famid = dict([(pid,i+100) for i, pid in enumerate(patients)])
    
    if verbose:
        print >> sys.stderr, "ok."
        print >> sys.stderr, "Found %d individuals that are members of families with at least 2 acertained individual(s)." % len(patient_data)
    
    if verbose:
        print >> sys.stderr, "Building the phenotype and covariate files for GCTA..."
    
    gcta_phen = list()
    gcta_covar = list()
    gcta_qcovar = list()
    
    for pid in tqdm(patients):
        famid, iid, sex, age, trait = patient_data[pid]
        
        if trait is None:
            trait = -9
        
        gcta_phen.append([pat2index[iid], pat2patid[iid], trait])
        gcta_covar.append([pat2index[iid], pat2patid[iid], sex, empi2demog[iid]['race']])
        gcta_qcovar.append([pat2index[iid], pat2patid[iid], age])
    
    if verbose:
        print >> sys.stderr, "Building the GRM file for covar..."
    
    gcta_grm = list()
    gcta_grm_id = list()
    
    for i in tqdm(range(len(patients))):
        pid1 = patients[i]
        gcta_grm_id.append([i+1, pat2patid[pid1]])
        for j, pid2 in enumerate(patients):
            if i == j:
                gcta_grm.append( [i+1, j+1, 1000, 1.0])
            elif i > j:
                gcta_grm.append( [i+1, j+1, 1000, relationships[pid1].get(pid2, 0.0)] )
    
    gcta_working_path = os.path.join(h2_path, 'working')
    
    if not os.path.exists(gcta_working_path):
        os.mkdir(gcta_working_path)
    
    if verbose:
        print >> sys.stderr, "All files loaded. Writing to disk at %s..." % gcta_working_path
    
    fh = open(os.path.join(gcta_working_path, 'trait.phen'), 'w')
    writer = csv.writer(fh, delimiter=' ', quoting=csv.QUOTE_NONE)
    writer.writerows(gcta_phen)
    fh.close()
    
    fh = open(os.path.join(gcta_working_path, 'trait.covar'), 'w')
    writer = csv.writer(fh, delimiter=' ', quoting=csv.QUOTE_NONE)
    writer.writerows(gcta_covar)
    fh.close()
    
    fh = open(os.path.join(gcta_working_path, 'trait.qcovar'), 'w')
    writer = csv.writer(fh, delimiter=' ', quoting=csv.QUOTE_NONE)
    writer.writerows(gcta_qcovar)
    fh.close()
    
    fh = gzip.open(os.path.join(gcta_working_path, 'trait.grm.gz'), 'w')
    writer = csv.writer(fh, delimiter=' ', quoting=csv.QUOTE_NONE)
    writer.writerows(gcta_grm)
    fh.close()
    
    fh = open(os.path.join(gcta_working_path, 'trait.grm.id'), 'w')
    writer = csv.writer(fh, delimiter=' ', quoting=csv.QUOTE_NONE)
    writer.writerows(gcta_grm_id)
    fh.close()
    
    if verbose:
        print >> sys.stderr, "Finished building working directories."
    
    sys.exit(101)
    return
    

def build_solar_directories(h2_path, iid2ped, empi2trait, fam2empi, fam2count, fam2proband, use_proband, trait_type, verbose = True, family_ids_only = None):
    """
    Build the directories with files needed to run SOLAR a single time.
    """


    if not os.path.exists(h2_path):
        os.mkdir(h2_path)

    if verbose:
        print >> sys.stderr, "Building files at path: %s" % h2_path

    if os.path.exists(os.path.join(h2_path, 'wd')):
        raise Exception("Directory already exists at path: %s" % os.path.join(h2_path, 'wd'))

    if verbose:
        print >> sys.stderr, "Building trait-specfic pedigree file...",

    trait_ped = list()
    for famid in family_ids_only:
        if fam2count[famid] < 2:
            continue

        for pid in fam2empi[famid]:
            row = iid2ped[pid]
            trait_value = empi2trait.get(pid, None)
            trait_ped.append( row + [trait_value] )

    if verbose:
        print >> sys.stderr, "ok."
        print >> sys.stderr, "Found %d individuals that are members of families with at least 2 acertained individual(s)." % len(trait_ped)

    # solar can't handle more than 32K individuals in a pedigree file, if you are using solar you should turn this back on
    if len(trait_ped) > 32000:
        raise SolarException("Too many individuals in a single pedigree. SOLAR max is 32,000.")

    num_runs = 1

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

    if verbose:
        print >> sys.stderr, "ok. (%d individuals in pedigree.)" % len(solar_ped)

    solar_working_path = os.path.join(h2_path, 'working')

    if verbose:
        print >> sys.stderr, "Creating solar working directory at %s..." % solar_working_path

    if not os.path.exists(solar_working_path):
        os.mkdir(solar_working_path)

    if verbose:
        print >> sys.stderr, "Writing .ped and .phen files for solar...",

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

    if verbose:
        print >> sys.stderr, "ok."
        print >> sys.stderr, "Writing out tcl scripts to run solar...",

    # load_pedigree.tcl
    tcl_fh = open(os.path.join(solar_working_path, 'load_pedigree.tcl'), 'w')
    tcl_fh.write(tcl_load_string)
    tcl_fh.close()

    # run_analysis.tcl
    tcl_fh = open(os.path.join(solar_working_path, 'run_analysis.tcl'), 'w')
    if trait_type == TRAIT_TYPE_QUANTITATIVE:
        tcl_fh.write(tcl_analysis_string_q)
    elif trait_type == TRAIT_TYPE_BINARY:
        tcl_fh.write(tcl_analysis_string_b)
    else:
        raise Exception("Error: Impossible trait type encountered: %s" % trait_type)
    tcl_fh.close()

    if verbose:
        print >> sys.stderr, "ok."

    if verbose:
        print >> sys.stderr, "Finished building working directories."

    return

def extract_convered_estimates(h2r_results, edge_eps, denoise_eps):
    converged = list()
    for h2, h2err, pval in h2r_results:
        if h2 == None or h2err == None or h2 < edge_eps or h2 > (1-edge_eps):
            continue

        if h2err < denoise_eps*h2:
            continue

        converged.append( (h2, h2err, pval) )
    return converged

def estimate_h2o(h2r_results, ci = 95., show_warnings=True, show_errors=True):

    num_converged = 0
    num_significant = 0
    sig_h2s = list()

    pcutoff = 0.05
    edge_eps = 1e-9
    denoise_eps = 0.05

    converged = extract_convered_estimates(h2r_results, edge_eps, denoise_eps)
    num_converged = len(converged)

    for h2, h2err, pval in converged:
        if pval < pcutoff:
            num_significant += 1
            sig_h2s.append((h2, h2err, pval))

    if num_significant == 0:
        if show_errors:
            print >> sys.stderr, "ERROR: There are no significant and converged estimates available."
        return False

    if num_significant < 30:
        if show_warnings:
            print >> sys.stderr, "WARNING: There are fewer than 30 (%d) significant and converged estimates." % num_significant

    h2o, solarerr, solarpval = sorted(sig_h2s)[len(sig_h2s)/2]
    h2o_estimates = zip(*sig_h2s)[0]

    cidiff = (100.-ci)/2.
    h2olo = numpy.percentile(h2o_estimates, cidiff)
    h2ohi = numpy.percentile(h2o_estimates, 100.-cidiff)
    posa = num_significant/float(num_converged)

    return h2o, h2olo, h2ohi, solarerr, solarpval, num_converged, num_significant, posa

def solar_strap(num_families, families_with_case, icd9, trait_type, num_attempts, solar_dir, iid2ped, all_traits, eth, fam2empi, fam2eth, all_fam2count, all_fam2proband, use_proband, house=False, nprocs=1, verbose=False, buildonly=False):
    """
    Run the bootstrapping algorithm to estimate the observational heritability for both
    the AE and ACE (if house=True) models of heritability. Results of this funciton can be parsed with
    estimate_h2o().
    """

    ae_h2r_results = list()
    ace_h2r_results = list()
    def log_solar_results(single_results):
        ae_h2r_results.append((single_results['AE']['h2r'], single_results['AE']['err'], single_results['AE']['pvalue']))
        ace_h2r_results.append((single_results['ACE']['h2r'], single_results['ACE']['err'], single_results['ACE']['pvalue']))
        if verbose and not buildonly:
                aeh2 = single_results['AE']['h2r']
                aeer = single_results['AE']['err']
                aepv = single_results['AE']['pvalue']
                if aeh2 is None:
                    aeh2 = numpy.nan
                if aeer is None:
                    aeer = numpy.nan
                if aepv is None:
                    aepv = numpy.nan

                aceh2 = single_results['ACE']['h2r']
                aceer = single_results['ACE']['err']
                acepv = single_results['ACE']['pvalue']
                if aceh2 is None:
                    aceh2 = numpy.nan
                if aceer is None:
                    aceer = numpy.nan
                if acepv is None:
                    acepv = numpy.nan

                print >> sys.stderr, "%10s %15s %5d %5d %7.2f %7.2f %10.2e %7.2f %7.2f %10.2e %10.4f" % (icd9, eth, num_families, len(ae_h2r_results), aeh2, aeer, aepv, aceh2, aceer, acepv, single_results['APF'])

    if num_families > len(families_with_case[icd9]):
        print >> sys.stderr, "Number of families to be sampled (%d) is larger than what is available (%d)." % (num_families, len(families_with_case[icd9]))
    else:
        if verbose and not buildonly:
            print >> sys.stderr, "%10s %15s %5s %5s %7s %7s %10s %7s %7s %10s %10s" % ('Trait', 'Ethnicity', 'NFam', 'Samp', 'AE h2', 'err', 'pval', 'ACE h2', 'err', 'pval', 'Sample AFP')

        pool = mp.Pool(nprocs)

        for i in range(num_attempts):
            h2_path = os.path.join(solar_dir, icd9, 'h2_%d_%s' % (num_families, random_string(5)))
            solar_args = (h2_path,
                            families_with_case,
                            icd9,
                            trait_type,
                            num_families,
                            iid2ped,
                            all_traits,
                            eth,
                            fam2empi,
                            fam2eth,
                            all_fam2count,
                            all_fam2proband,
                            use_proband,
                            house,
                            verbose,
                            buildonly)
            if nprocs == 1:
                log_solar_results(solar(*solar_args))
            else:
                pool.apply_async(solar, args = solar_args, callback = log_solar_results)

        pool.close()
        pool.join()

    return ae_h2r_results, ace_h2r_results

def gcta_strap(num_families, families_with_case, icd9, trait_type, num_attempts, solar_dir, iid2ped, all_traits, eth, empi2demog, fam2empi, fam2eth, all_fam2count, relationships, house=False, nprocs=1, verbose=False, buildonly=False):
    """
    Run the bootstrapping algorithm to estimate the observational heritability for both
    the AE and ACE (if house=True) models of heritability. Results of this funciton can be parsed with
    estimate_h2o().
    """
    
    h2_path = os.path.join(solar_dir, icd9, 'h2_%d_%s' % (num_families, random_string(5)))
    
    if num_families > len(families_with_case[icd9]):
        print >> sys.stderr, "Number of families to be sampled (%d) is larger than what is available (%d)." % (num_families, len(families_with_case[icd9]))
    else:
        gcta_args = (h2_path,
                        families_with_case,
                        icd9,
                        trait_type,
                        num_families,
                        all_traits,
                        eth,
                        empi2demog,
                        fam2empi,
                        fam2eth,
                        all_fam2count,
                        relationships,
                        verbose,
                        buildonly)
        gcta(*gcta_args)
    
    ae_h2r_results = list()
    ace_h2r_results = list()
    
    return ae_h2r_results, ace_h2r_results

def gcta(h2_path, families_with_case, icd9, trait_type, num_families, all_traits, eth, empi2demog, fam2empi, fam2eth, all_fam2count, relationships, verbose, buildonly):
    """
    Setup the data for GCTA and run GCTA for the given condition.
    """
    
    # if the h2_path exists, we remove it
    if os.path.exists(h2_path):
        shutil.rmtree(h2_path)
    
    available_families = [fid for fid in families_with_case[icd9] if eth == 'ALL' or fam2eth[fid] == eth]
    chosen_families = random.sample(available_families, num_families)
    if trait_type == TRAIT_TYPE_BINARY:
        apf = numpy.mean([numpy.sum([all_traits[icd9].get(iid, 0) for iid in fam2empi[famid]]) for famid in chosen_families])
    elif trait_type == TRAIT_TYPE_QUANTITATIVE:
        apf = numpy.mean([all_fam2count[icd9][famid] for famid in chosen_families])
    
    build_gcta_directories(h2_path,
                           empi2demog,
                           all_traits[icd9],
                           fam2empi,
                           all_fam2count[icd9],
                           relationships,
                           trait_type,
                           verbose=True, # this one is kindof annoying if it is on
                           family_ids_only=chosen_families)

    print >> sys.stderr, "BUILDONLY! We haven't implemented GCTA yet."
    

def solar(h2_path, families_with_case, icd9, trait_type, num_families, iid2ped, all_traits, eth, fam2empi, fam2eth, all_fam2count, all_fam2proband, use_proband, house, verbose, buildonly):
    """
    Setup the data for solar and run solar for the given condition.

    """

    # if the h2_path exists, we remove it
    if os.path.exists(h2_path):
        shutil.rmtree(h2_path)

    available_families = [fid for fid in families_with_case[icd9] if eth == 'ALL' or fam2eth[fid] == eth]
    chosen_families = random.sample(available_families, num_families)
    if trait_type == TRAIT_TYPE_BINARY:
        apf = numpy.mean([numpy.sum([all_traits[icd9].get(iid, 0) for iid in fam2empi[famid]]) for famid in chosen_families])
    elif trait_type == TRAIT_TYPE_QUANTITATIVE:
        apf = numpy.mean([all_fam2count[icd9][famid] for famid in chosen_families])

    build_solar_directories(h2_path,
                           iid2ped,
                           all_traits[icd9],
                           fam2empi,
                           all_fam2count[icd9],
                           all_fam2proband[icd9],
                           use_proband,
                           trait_type,
                           verbose=False, # this one is kindof annoying if it is on
                           family_ids_only=chosen_families)

    if not buildonly:
        #print >> sys.stderr, h2_path
        results = single_solar_run(h2_path, house, verbose)
        results['APF'] = apf
        shutil.rmtree(h2_path)
    else:
        results = {'AE':{'h2r':None, 'err':None, 'pvalue':None}, 'ACE':{'h2r':None, 'err':None, 'pvalue':None}, 'APF': apf}

    return results

def parse_polygenic_out(polygenic_out_fn, verbose):

    if verbose:
        print >> sys.stderr, "Parsing polygenic output to get h2 values..."

    if not os.path.exists(polygenic_out_fn):
        print >> sys.stderr, "WARNING: No polygenic outfile found at: %s" % polygenic_out_fn
        h2r, h2r_err, p = None, None, None
    else:
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

    return {'h2r':h2r, 'err':h2r_err, 'pvalue':p}

def single_solar_run(h2_path, house=False, verbose=False, really_verbose=False):
    """
    Runs solar for the AE and ACE (if house=True) models, returns parsed results.

    """

    solar_working_path = os.path.join(h2_path, 'working')

    if really_verbose:
        print >> sys.stderr, "Loading data into solar in %s with:" % solar_working_path
        print >> sys.stderr, "cd %s > /dev/null && solar loadped > /dev/null && cd - > /dev/null " % solar_working_path
    os.system("cd %s > /dev/null && solar loadped > /dev/null && cd - > /dev/null " % solar_working_path)

    if really_verbose:
        print >> sys.stderr, "Running analysis without house effects..."
        print >> sys.stderr, "cd %s > /dev/null && solar runanalysis > /dev/null && cd - > /dev/null " % solar_working_path
    os.system("cd %s > /dev/null && solar runanalysis > /dev/null && cd - > /dev/null " % solar_working_path)

    polygenic_out_fn = os.path.join(solar_working_path, 'pheno', 'polygenic.out')

    ae_results = parse_polygenic_out(polygenic_out_fn, verbose=really_verbose)

    if house:
        if really_verbose:
            print >> sys.stderr, "Running analysis with house effects..."
            print >> sys.stderr, "cd %s > /dev/null && solar runanalysishouse > /dev/null && cd - > /dev/null " % solar_working_path
        os.system("cd %s > /dev/null && solar runanalysishouse > /dev/null && cd - > /dev/null " % solar_working_path)
        ace_results = parse_polygenic_out(polygenic_out_fn, verbose=really_verbose)
    else:
        ace_results = {'h2r':None, 'err':None, 'pvalue':None}

    if really_verbose:
        print >> sys.stderr, 'AE', "%(h2r)s %(err)s %(pvalue)s" % ae_results
        if house:
            print >> sys.stderr, 'ACE', "%(h2r)s %(err)s %(pvalue)s" % ace_results

    return {'AE': ae_results, 'ACE': ace_results}
