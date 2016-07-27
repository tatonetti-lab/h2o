"""
solar_setup_h2_analysis.py

@author Nicholas P. Tatonett, 2016
"""

import os
import sys
import csv
import gzip
import math
import numpy
import pystan
import random

from tqdm import tqdm
from scipy import stats
from collections import defaultdict

from familial_recurrence import assign_family_ethnicities
from familial_recurrence import load_demographics
from familial_recurrence import load_family_ids


USAGE_PROMPT = """
setup_h2_analysis.py v0.1

USAGE:
python setup_h2_analysis.py path_to_trait_folder q|quant|b|bin

"""

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

def bayesian_averages(h2s, h2errs, iter=1000):
    h2_code = """
    data {
        int<lower=0> J; // number of relationships
        real y[J]; // h2 estimates
        real<lower=0> sigma[J]; // s.e. of h2 estimates
    }
    parameters {
        real mu;
        real<lower=0> tau;
        real eta[J];
    }
    transformed parameters {
        real theta[J];
        for (j in 1:J)
        theta[j] <- mu + tau * eta[j];
    }
    model {
        eta ~ normal(0, 1);
        y ~ normal(theta, sigma);
    }
    """
    
    h2_dat = {'J': len(h2s),
              'y': h2s,
              'sigma': h2errs
    }
    
    fit = pystan.stan(model_code=h2_code, data=h2_dat, iter=iter, chains=4)
    estimates = fit.extract()['mu']
    results = {
        'mean': estimates.mean(),
        'lower': numpy.percentile(estimates, 2.5),
        'upper': numpy.percentile(estimates, 97.5)
    }
        
    print results['mean']
    
    return results

def fixed_effects_averages(h2s, h2errs):
    # taken from https://www.meta-analysis.com/downloads/Meta-analysis%20fixed%20effect%20vs%20random%20effects.pdf
    
    ws = [1.0/e**2 for e in h2errs]
    hbar = sum([h*w for h, w in zip(h2s, ws)])/sum(ws)
    vbar = 1.0/sum(ws)
    
    results = {
        'mean': hbar,
        'lower': hbar-1.96*numpy.sqrt(vbar),
        'upper': hbar+1.96*numpy.sqrt(vbar)
    }
    
    return results

def load_generic_pedigree(generic_ped_path, empi2sex, empi2age):
    
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
    
    return generic_pedigree

def build_solar_directories(h2_path, generic_pedigree, empi2trait, ethnicities, min_ascertained, fam2empi, fam2eth, fam2count, fam2proband, use_proband, trait_type, verbose = True, family_ids_only = None):
    
    run_list = dict()
    
    if not os.path.exists(h2_path):
        os.mkdir(h2_path)
    
    for ethnicity in ethnicities:
        if verbose:
            print >> sys.stderr, "Building files for: %s" % ethnicity
        
        if os.path.exists(os.path.join(h2_path, ethnicity, 'wd1')):
            if verbose:
                print >> sys.stderr, "Solar directories have already been built. Will not regenerate."
                print >> sys.stderr, "Building run list from existing directory structure...",
            solar_wds = [dn for dn in os.listdir(os.path.join(h2_path, ethnicity)) if dn.startswith('wd')]
            run_list[ethnicity] = range(1,len(solar_wds)+1)
            if verbose:
                print >> sys.stderr, "ok. Found %d directories." % len(solar_wds)
            continue
        
        if verbose:
            print >> sys.stderr, "Building trait-specfic pedigree file...",
        trait_ped = list()
        for row in generic_pedigree:
            # uncomment to exclude any families that do not have ascertained individuals
            if fam2count[row[0]] < min_ascertained:
                continue
            
            if not family_ids_only is None and not row[0] in family_ids_only:
                continue
            
            # uncomment to exclude any families with too many unknown family members
            #if fam2propunk[row[0]] > 0.4:
            #    continue
            
            if ethnicity != 'ALL' and not fam2eth[row[0]] == ethnicity:
                continue
            
            trait_value = empi2trait.get(row[1], None)
            trait_ped.append( row + [trait_value] )
        
        if verbose:
            print >> sys.stderr, "ok."
            print >> sys.stderr, "Found %d individuals (ethnicty = %s) that are members of families with at least %d acertained individual(s)." % (len(trait_ped), ethnicity, min_ascertained)
        
        # solar can't handle more than 32K individuals in a pedigree file, if you are using solar you should turn this back on
        if len(trait_ped) > 32000:
            if verbose:
                print >> sys.stderr, "*********"
                print >> sys.stderr, "There are more than 32,000 individuals in this pedigree. We will need to split them and " + \
                    "then reaggregate the results togehter using a meta-analysis."
                
            num_runs = math.ceil(len(trait_ped)/32000.)
            run2fid = defaultdict(set)
            
            for fid in fam2empi.keys():
                run_id = random.randint(1,num_runs)
                run2fid[run_id].add(fid)
                
            run2fid = dict(run2fid)
            
            if verbose:
                print >> sys.stderr, "Divided families into %d runs." % num_runs
                print >> sys.stderr, "*********"
        else:
            num_runs = 1
            run2fid = dict()
            run2fid[1] = fam2empi.keys()
        
        fid2run = dict()
        for rid, families in run2fid.items():
            for fid in families:
                fid2run[fid] = rid
        
        run_list[ethnicity] = sorted(run2fid.keys())
        
        for rid in run_list[ethnicity]:
            families = run2fid[rid]
            if verbose:
                print >> sys.stderr, "\nBuilding run-specific (run_id = %d) pedigree file..." % rid,
            
            solar_ped = list()
            solar_phen = list()
            for row in trait_ped:
                famid, iid, fid, mid, sex, age, trait = row
                
                if fid2run[famid] != rid:
                    continue
                
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
            
            ethnicity_directory = os.path.join(h2_path, ethnicity)
            if not os.path.exists(ethnicity_directory):
                os.mkdir(ethnicity_directory)
            
            solar_working_path = os.path.join(ethnicity_directory, 'wd%d' % rid)
            
            if verbose:
                print >> sys.stderr, "Creating solar working directory at %s..." % solar_working_path
            
            while os.path.exists(solar_working_path):
                if verbose:
                    print >> sys.stderr, "FAILED.\n\nYou already have a solar working directory. " + \
                        "You'll need to remove the directory. Press any key once you have removed the directory to continue.." + \
                        "You can remove this directory with `rm -rf %s` or all solar directories " % solar_working_path + \
                        "with `rm -rf %s`" % os.path.join(h2_path, 'wd*')
                raw_input()
            
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
                raise Exception("Error: Impossible trait type encountered.")
            tcl_fh.close()
            
            if verbose:
                print >> sys.stderr, "ok."
            
        if verbose:
            print >> sys.stderr, "Finished building working directories."
    
    return run_list

def run_solar_directories(h2_path, run_list, ethnicities, verbose = True):
    
    if verbose:
        print >> sys.stderr, "Running solar in 2 steps."
    
    h2_results = dict()
    
    for ethnicity in ethnicities:
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
                print >> sys.stderr, "FAILED. SOLAR failed to run. Could be a convergence error."
                h2r = None
                p = None
            
            h2r_err = None
            
            if p < 0.05:
                try:
                    h2r_err = float(h2r_raw[1].split()[3])
                except IndexError:
                    h2r_err = None
                h2s.append(h2r)
                h2errs.append(h2r_err)
                pvals.append(p)
            
            print >> sys.stderr, "\t", h2r, h2r_err, p
        
        if verbose:
            print >> sys.stderr, "ok. (found %d h2 estimates of %d expected)" % (len(h2s), len(run_list[ethnicity]))
        
        failed = False
        if len(h2s) == 0:
            if verbose:
                print >> sys.stderr, "WARNING. There are no available estimates of heritability. This run failed."
            failed = True
            
        elif len(h2s) > 1:
            if verbose:
                print >> sys.stderr, "More than one estimate received. We will aggregate using a meta analysis..."
            #print >> sys.stderr, "Aggregating using STAN..."
            #results = bayesian_averages(h2s, h2errs)
            if verbose:
                print >> sys.stderr, "Aggregating using Fixed Effects Model..."
            results = fixed_effects_averages(h2s, h2errs)
            
            h2r = results['mean']
            h2r_lower = results['lower']
            h2r_upper = results['upper']
        else:
            h2r = h2s[0]
            h2r_stderr = h2errs[0]
            h2r_lower = h2r-h2r_stderr
            h2r_upper = h2r+h2r_stderr
        
        if not failed:
            if verbose:
                print >> sys.stderr, "Heritabiity estimates completed: h2r = %.3f (%.3f - %.3f)" % (h2r, h2r_lower, h2r_upper)
        
        results_fh = open(os.path.join(h2_path, '%s_heritability_results.txt' % ethnicity), 'w')
        writer = csv.writer(results_fh)
        writer.writerow(['run', 'h2r', 'h2r_lower', 'h2r_upper'])
        
        if not failed:
            writer.writerow(['AVG', h2r, h2r_lower, h2r_upper])
            h2_results[ethnicity] = [h2r, h2r_lower, h2r_upper]
        
        for hi, (h2r, h2r_stderr) in enumerate(zip(h2s, h2errs)):
            writer.writerow([hi, h2r, h2r-1.96*h2r_stderr, h2r+1.96*h2r_stderr])
        
        results_fh.close()
    
    return h2_results

if __name__ == '__main__':
    print >> sys.stderr, "This script is not intended to be run directly."
    sys.exit(100)