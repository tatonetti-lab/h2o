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

        family_ethnicities = [empi2demog[e]['race'] for e in members if empi2demog.get(e,0) != 0 if empi2demog[e]['race'] != 'Unknown' ]
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


def load_covariates(cov_file_path):
    """
    Load covariates data from file.

    """

    print >> sys.stderr, "Loading patient covariates..."

    # load covariates data (empi, diabetes status etc.)
    fh = gzip.open(cov_file_path)
    reader = csv.reader(fh, delimiter='\t')
    header_row = reader.next()
    #get list of covariates, ignore first column (empi)
    cov_header = header_row[1:]

    empi2cov = dict()
    for row in tqdm(reader):

        empi = row[0]
        empi2cov[empi] = dict(zip(cov_header, row[1:]))

    return empi2cov,cov_header


#KLB hhid_4
def load_hhid(cov_file_path):
    """
    Load hhid data from file.

    """

    print >> sys.stderr, "Loading patient household IDs..."
    fh = gzip.open(cov_file_path)
    reader = csv.reader(fh, delimiter='\t')
    header_row = reader.next()

    empi2hhid = dict()
    for row in tqdm(reader):
        empi = row[0]
        empi2hhid[empi] = row[1]

    return empi2hhid


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
        if empi2demog[empi]['race'] in ('NA', 'NULL', '') :
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
        empi2demog[empi]['age'] = float(empi2demog[empi]['age']) if empi2demog[empi]['age'] != 'NULL' else None

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

def build_solar_directories(h2_path, iid2ped, empi2trait,
fam2empi, fam2count, fam2proband, use_proband, trait_type, empi2cov, cov_list, empi2hhid,
verbose = True, family_ids_only = None):
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
        #Household ID default is mother ID
        if empi2hhid is None:
            solar_ped.append( [famid, iid, fid, mid, sex, mid] )

        else:
            #If missing household ID when user-defined leave blank, blank or 0 makes SOLAR assign a new unique hhid for that patient
            solar_ped.append( [famid, iid, fid, mid, sex, empi2hhid.get(iid,'')] )

        if trait is None:
            trait = ''
            age = ''

        if empi2cov == None:
            if use_proband:
                solar_phen.append([iid, trait, age, 1 if iid == fam2proband[famid] else 0])
            else:
                solar_phen.append([iid, trait, age])
        else:
            iidcov_dict = empi2cov.get(iid,'')
            if use_proband:
                holder = [iid, trait, age, 1 if iid == fam2proband[famid] else 0]
            else:
                holder = [iid, trait, age]

            for each in cov_list:
                if type(iidcov_dict) == dict:
                    holder.append(iidcov_dict.get(each,''))
                else:
                    holder.append('')

            solar_phen.append(holder)

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

    writer.writerow(['famid', 'id', 'fa', 'mo', 'sex', 'hhid'])
    writer.writerows(solar_ped)
    ped_fh.close()

    phen_fh = open(os.path.join(solar_working_path, 'phenotypes.phen'), 'w')
    writer = csv.writer(phen_fh, delimiter=',', quoting=csv.QUOTE_NONE)


    write_holder = ['id', 'pheno', 'age']
    if use_proband:
        write_holder.append('proband')
    if cov_list != None:
        for each in cov_list:
            write_holder.append(each)

    writer.writerow(write_holder)
    writer.writerows(solar_phen)
    phen_fh.close()

    if verbose:
        print >> sys.stderr, "ok."
        print >> sys.stderr, "Writing out tcl scripts to run solar...",

    # the typical id used for household when it is not directly available is
    # the mother id (mo)
    # other options are household ID (bulidingID+SteetID+ZipCode+familyID) or familyID

    #make cov_list a list if None for tcl

    #normalize quant trait with inorm instead of tdist (more robust)

    tcl_load_string = """
    proc loadped {} {
        load pedigree pedigree.ped
        load phenotypes phenotypes.phen
    }
    """
    if len(cov_list) == 0:
        tcl_analysis_string_b = """
        proc runanalysis {} {
            model new
            trait pheno
            polygenic -screen
        }

        proc runanalysishouse {} {
            model new
            trait pheno
            house
            polygenic -screen
        }
        """

        tcl_analysis_string_q = """
        proc runanalysis {} {
            define ipheno = inorm_pheno
            trait ipheno
            polygenic -screen
        }
        proc runanalysishouse {} {
            model new
            define ipheno = inorm_pheno
            trait ipheno
            house
            polygenic -screen
        }
        """
    else:
        tcl_analysis_string_b = """
        proc runanalysis {} {
            model new
            trait pheno
            covariates sex age %s
            polygenic -screen
        }

        proc runanalysishouse {} {
            model new
            trait pheno
            covariates sex age %s
            house
            polygenic -screen
        }
        """ % (" ".join(cov_list)," ".join(cov_list))

        tcl_analysis_string_q = """
        proc runanalysis {} {
            define ipheno = inorm_pheno
            trait ipheno
            covariates sex age %s
            polygenic -screen
        }
        proc runanalysishouse {} {
            model new
            define ipheno = inorm_pheno
            trait ipheno
            covariates sex age %s
            house
            polygenic -screen
        }
        """ % (" ".join(cov_list)," ".join(cov_list))

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

def extract_converged_sig_estimates(results, params, process_flag, pcutoff):
    """
    extracts all runs that converge and that converge and are significant
    """
    converged = list()
    sig_converged = list()
    num_significant = 0

    if process_flag == "h2":
        edge_eps, denoise_eps = params

        for h2, h2_err, pval in results:
            if h2 == None or h2_err == None or h2 < edge_eps or h2 > (1-edge_eps):
                continue

            if h2_err < denoise_eps*h2:
                continue

            converged.append( (h2, h2_err, pval) )

            if pval < pcutoff:
                num_significant += 1
                sig_converged.append( (h2, h2_err, pval) )



    elif process_flag == "c2":
        edge_eps_c2, denoise_eps_c2 = params
        #Note, edge eps and denoise eps are not estimated yet.  Placeholder = 0

        for c2, c2_err, c2_pval in results:
            if c2 == None or c2_err == None or c2 <= edge_eps_c2 or c2 >= (1-edge_eps_c2):
                #Update <= to < and >= to > to mirror h2 formatting when c2 params estimated. Param currently 0
                continue

            if c2_err <= denoise_eps_c2*c2:
                #Update <= to < to mirror h2 formatting when c2 params estimated.  Param currently 0
                continue

            converged.append( (c2, c2_err, c2_pval) )

            if c2_pval < pcutoff:
                num_significant += 1
                sig_converged.append( (c2, c2_err, c2_pval) )



    elif process_flag == "h2c2":

        edge_eps, denoise_eps, edge_eps_c2, denoise_eps_c2 = params

        for h2, h2_err, pval, c2, c2_err, c2_pval in results:
            if h2 == None or h2_err == None or h2 < edge_eps or h2 > (1-edge_eps) or\
            c2 == None or c2_err == None or c2 <= edge_eps_c2 or c2 >= (1-edge_eps_c2):
            #For c2 part, update <= to < and >= to > to mirror h2 formatting when c2 params estimated. Param currently 0
                continue

            if h2_err < denoise_eps*h2 or c2_err <= denoise_eps_c2*c2:
            #For c2 part, update <= to < to mirror h2 formatting when c2 params estimated.  Param currently 0
                continue

            converged.append( (h2, h2_err, pval, c2, c2_err, c2_pval) )

            if pval < pcutoff and c2_pval < pcutoff:
                num_significant += 1
                sig_converged.append((h2, h2_err, pval, c2, c2_err, c2_pval))



    return converged, sig_converged,  len(converged), num_significant



def estimate_h2o(h2r_results, process_flag = "h2", ci = 95., show_warnings=True, show_errors=True):

    num_converged = 0
    num_significant = 0
    sig_h2s = list()

    pcutoff = 0.05
    cidiff = (100.-ci)/2.
    edge_eps = 1e-9
    denoise_eps = 0.05

    edge_eps_c2 = 0
    denoise_eps_c2 = 0
    #These two parameters are not estimated yet, using 0 as a placeholder

    if process_flag == "h2":
        params = [edge_eps, denoise_eps]
        converged, sig_converged, num_converged , num_significant = extract_converged_sig_estimates(h2r_results, params, "h2", pcutoff)

    elif process_flag == "c2":
        params = [edge_eps_c2,denoise_eps_c2]
        converged, sig_converged, num_converged, num_significant = extract_converged_sig_estimates(h2r_results, params, "c2", pcutoff)

    elif process_flag == "h2c2":
        params = [edge_eps, denoise_eps, edge_eps_c2, denoise_eps_c2]
        converged, sig_converged, num_converged, num_significant = extract_converged_sig_estimates(h2r_results, params, "h2c2", pcutoff)
# Don't have to return length of converged and significant, exclude in function

    num_converged = len(converged)
    num_significant = len(sig_converged)


    if num_significant == 0:
        if show_errors:
            print >> sys.stderr, "ERROR: There are no significant and converged estimates available."
        return False

    if num_significant < 30:
        if show_warnings:
            print >> sys.stderr, "WARNING: There are fewer than 30 (%d) significant and converged estimates." % num_significant

    posa = num_significant/float(num_converged)

    if process_flag == "h2" or process_flag == "c2":
        estimate, solarerr, solarpval = sorted(sig_converged)[len(sig_converged)/2]

        estimates = zip(*sig_converged)[0]

        estlo = numpy.percentile(estimates, cidiff)
        esthi = numpy.percentile(estimates, 100.-cidiff)

        return estimate, estlo, esthi, solarerr, solarpval, num_converged, num_significant, posa
        #Update naming here

    else:
        h2o, h2_err, h2_pval, c2o, c2_err, c2_pval = sorted(sig_converged)[len(sig_converged)/2]

        h2o_estimates = zip(*sig_converged)[0]
        h2olo = numpy.percentile(h2o_estimates, cidiff)
        h2ohi = numpy.percentile(h2o_estimates, 100.-cidiff)

        c2o_estimates = zip(*sig_converged)[3]
        c2olo = numpy.percentile(c2o_estimates, cidiff)
        c2ohi = numpy.percentile(c2o_estimates, 100.-cidiff)


        return h2o, h2olo, h2ohi, h2_err, h2_pval, num_converged, num_significant, posa, c2o, c2olo, c2ohi, c2_err, c2_pval




def solar_strap(num_families, families_with_case, icd9, trait_type, num_attempts, solar_dir,
iid2ped, all_traits, eth, fam2empi, fam2eth,all_fam2count, all_fam2proband,
use_proband, house=False, nprocs=1, verbose=False, buildonly=False, empi2cov=None, cov_list=None, empi2hhid=None, output_fams=False):

    """
    Run the bootstrapping algorithm to estimate the observational heritability for both
    the AE and ACE (if house=True) models of heritability. Results of this funciton can be parsed with
    estimate_h2o().
    """

    ae_h2r_results = list()
    ace_h2r_results = list()
    h2r_families = list()

    def log_solar_results(single_results):
        ae_h2r_results.append((single_results['AE']['h2r'], single_results['AE']['err'], single_results['AE']['pvalue']))

        ace_h2r_results.append((single_results['ACE']['h2r'], single_results['ACE']['err'], single_results['ACE']['pvalue'],single_results['ACE']['c2'], single_results['ACE']['c2_err'], single_results['ACE']['c2_pvalue']  ))

        if output_fams:
            h2r_families.append(single_results['chosen_families'])
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

                acec2 = single_results['ACE']['c2']
                acec2er = single_results['ACE']['c2_err']
                acec2pv = single_results['ACE']['c2_pvalue']
                if acec2 is None:
                    acec2 = numpy.nan
                if acec2er is None:
                    acec2er = numpy.nan
                if acec2pv is None:
                    acec2pv = numpy.nan

                print >> sys.stderr, "%10s %15s %5d %5d %7.2f %7.2f %10.2e %7.2f %7.2f %10.2e %7.2f %7.2f %10.2e %10.4f" % (icd9, eth, num_families, len(ae_h2r_results), aeh2, aeer, aepv, aceh2, aceer, acepv, acec2, acec2er, acec2pv, single_results['APF'])

    if num_families > len(families_with_case[icd9]):
        print >> sys.stderr, "Number of families to be sampled (%d) is larger than what is available (%d)." % (num_families, len(families_with_case[icd9]))
    else:

        if verbose and not buildonly:
            print >> sys.stderr, "%10s %15s %5s %5s %4s %5s %7s %5s %3s %4s %15s %6s %7s %10s" % ('Trait', 'Ethnicity', 'NFam', 'Samp', 'AE h2', 'err', 'pval', 'ACE h2', 'err', 'pval','Shared Env c2', 'c2 err', 'c2 pval', 'Sample AFP')

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
                            buildonly,
                            empi2cov,
                            cov_list,
                            empi2hhid,
                            output_fams)
            if nprocs == 1:
                log_solar_results(solar(*solar_args))
            else:
                pool.apply_async(solar, args = solar_args, callback = log_solar_results)

        pool.close()
        pool.join()


    return ae_h2r_results, ace_h2r_results, h2r_families

def solar(h2_path, families_with_case, icd9, trait_type, num_families,
iid2ped, all_traits, eth, fam2empi, fam2eth,all_fam2count, all_fam2proband,
use_proband, house, verbose, buildonly, empi2cov, cov_list, empi2hhid, output_fams):
    """
    Setup the data for solar and run solar for the given condition.

    """

    #if the h2_path exists, we remove it
    if not buildonly and os.path.exists(h2_path):
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
                           empi2cov,
                           cov_list,
                           empi2hhid,
                           verbose=False, # this one is kindof annoying if it is on
                           family_ids_only=chosen_families)

    if not buildonly:
        #print >> sys.stderr, h2_path

        results = single_solar_run(h2_path, trait_type, house, verbose)
        results['APF'] = apf

        if output_fams:
            results["chosen_families"] = chosen_families
        #Deletes directory
        shutil.rmtree(h2_path)
    else:

        results = {'AE':{'h2r':None, 'err':None, 'pvalue':None}, 'ACE':{'h2r':None, 'err':None, 'pvalue':None, 'c2':None, 'c2_err':None, 'c2_pvalue':None}, 'APF': apf}

        if output_fams:
            results["chosen_families"] = chosen_families

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

        c2_raw = [row.strip() for row in results if row.find('C2') != -1]

        try:
            c2 = float(c2_raw[0].split()[2])
            c2_pval = float(c2_raw[0].split()[5])

        #NOTE There is sometimes no standard error in c2 output!
        except IndexError:
            if verbose:
                print >> sys.stderr, "SOLAR failed to converge on a shared environment estimate. Could be a convergence error."
            c2 = None
            c2_pval = None
        c2_err = None

        try:
            c2_err = float(c2_raw[1].split()[3])
        except (ValueError, IndexError):
            c2_raw = None


    return {'h2r':h2r, 'err':h2r_err, 'pvalue':p,'c2':c2, 'c2_err':c2_err, 'c2_pvalue':c2_pval }

def single_solar_run(h2_path, trait_type, house=False, verbose=False, really_verbose=False):
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


    if trait_type == TRAIT_TYPE_QUANTITATIVE:
        polygenic_out_fn = os.path.join(solar_working_path, 'ipheno', 'polygenic.out')

    ae_results = parse_polygenic_out(polygenic_out_fn, verbose=really_verbose)

    if house:
        if really_verbose:
            print >> sys.stderr, "Running analysis with house effects..."
            print >> sys.stderr, "cd %s > /dev/null && solar runanalysishouse > /dev/null && cd - > /dev/null " % solar_working_path
        os.system("cd %s > /dev/null && solar runanalysishouse > /dev/null && cd - > /dev/null " % solar_working_path)
        ace_results = parse_polygenic_out(polygenic_out_fn, verbose=really_verbose)

    else:

        ace_results = {'h2r':None, 'err':None, 'pvalue':None, 'c2':None,  'c2_err':None, 'c2_pvalue':None }

    if really_verbose:
        print >> sys.stderr, 'AE', "%(h2r)s %(err)s %(pvalue)s" % ae_results
        if house:
            print >> sys.stderr, 'ACE', "%(h2r)s %(err)s %(pvalue)s %(c2)s %(c2_err)s %(c2_pvalue)s" % ace_results

    return {'AE': ae_results, 'ACE': ace_results}
