"""
familial_recurrence.py

Estimate the rate of familial recurrence for a given trait. This script can estimate recurrence rates
based on anyone in the same family (familial) or based on a particular relationships type (e.g. sibling).

@author Nicholas P. Tatonetti and Rami Vanguri

USAGE
-----
python src/familial_recurrence.py trait=traits/binary/t2dm

python src/familial_recurrence.py trait=traits/binary/t2dm

"""

import os
import sys
import csv
import gzip
import numpy
import random
from collections import defaultdict

from tqdm import tqdm

def assign_family_ethnicities(fam2empi, empi2demog, print_breakdown=True):
    
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

def compute_recurrence(siblings, trait):
    
    # determine if this trait has negatives and positives, or just positives
    nvals = len(set(trait.values()))
    if nvals == 0 or nvals > 2:
        raise Exception("Traits have too few or too many values: %d" % nvals)        

    N = 0
    N_affected = 0
    N_sibs = 0
    N_sibs_affected = 0
    
    concordant = 0
    discordant = 0
    
    random.seed(0)
    
    for iid, sibs in siblings.items():
        
        if nvals == 2 and not iid in trait:
            continue
        
        sib = list(sibs)[random.randint(0,len(sibs)-1)]
        
        if nvals == 2 and not sib in trait:
            continue
        
        if random.randint(0,1) == 1:
            tmp = iid
            iid = sib
            sib = tmp
        
        if trait.get(iid, False):
            N_affected += 1
        if trait.get(sib, False):
            N_affected += 1
        
        N += 2
        
        if trait.get(iid, False):
            N_sibs += 1
            if trait.get(sib, False):
                N_sibs_affected += 1
        
        if trait.get(iid, False) or trait.get(sib, False):
            if trait.get(iid, False) == trait.get(sib, False):
                concordant += 1
            else:
                discordant += 1
    
    # sib_recurrence = (N_sibs_affected/float(N_sibs))
    # sib_recurrence_err = (sib_recurrence*(1-sib_recurrence)/float(N_sibs))**0.5
    # prevelance = (N_affected/float(N))
    # prevelance_err = (prevelance*(1-prevelance)/float(N))**0.5
    
    return N_sibs_affected, N_sibs, N_affected, N, concordant, discordant

def main(demographic_file, family_file, relationships_file, recurrence_type, trait_path, data_directory='./common_data/'):
    
    # load demographic data
    empi2demog = load_demographics(os.path.join(data_directory, demographic_file))
    
    # load family ids
    empi2fam, fam2empi = load_family_ids(os.path.join(data_directory, family_file))
    
    # load relationships
    relationships = load_relationships(os.path.join(data_directory, relationships_file), True)
    available_relationships = sorted(relationships.keys())
    
    print >> sys.stderr, "Working with trait at path: %s" % trait_path
    
    if not os.path.exists(trait_path):
        print >> sys.stderr, "\n\nNo directory at given path. You must first create a directory with a phenotyping sql file."
        sys.exit(10)
    
    print >> sys.stderr, "Checking for trait file..."
    trait_file_path = os.path.join(trait_path, 'trait_file.txt.gz')
    
    run_mysql = True
    
    trait_query_file_path = os.path.join(trait_path, 'query_trait.sql')
    
    if os.path.exists(trait_file_path):
        print >> sys.stderr, "Found trait file. Checking if the query changed..."
        if os.path.exists(os.path.join(trait_path, '.query_trait.sql')):
            cached_query = open(os.path.join(trait_path, '.query_trait.sql')).read()
            current_query = open(trait_query_file_path).read()
            if cached_query == current_query:
                print >> sys.stderr, "\tCached query matches current query. Will skip querying the database again."
                run_mysql = False
            else:
                print >> sys.stderr, "\tNew query found. Will freshly query the database.."
    
    if run_mysql:
        print >> sys.stderr, "Running MySQL query against database..."
        cmd = "mysql < %s | gzip > %s" % (trait_query_file_path, trait_file_path)
        print >> sys.stderr, "\trunning command: %s" % cmd
        res = os.system(cmd)
        os.system("cp %s %s" % (trait_query_file_path, os.path.join(trait_path, '.query_trait.sql')))
    
    if not os.path.exists(trait_file_path) or os.stat(trait_file_path).st_size == 0:
        print >> sys.stderr, "FAILED\n. MySQL query failed. Check the syntax in your query_trait.sql file."
        sys.exit(13)
    
    print >> sys.stderr, "Loading phenotype (trait) data from %s" % trait_file_path
    
    fh = gzip.open(trait_file_path)
    reader = csv.reader(fh, delimiter='\t')
    reader.next()
    
    empi2trait = dict()
    
    for empi, value in tqdm(reader):
        if empi in empi2fam:
            empi2trait[empi] = True if int(value) == 1 else False
    
    print >> sys.stderr, "Trait loaded, founded %d cases of %d total." % (sum(empi2trait.values()), len(empi2trait))
    
    print >> sys.stderr, "Recurrence calculations:"
    
    print >> sys.stdout, ",".join(map(str, ['Relationship', 'N Rel Affected', 'N Relationships', 'N Affected', 'N Total', 'Concordant', 'Discordant']))
    
    for rel in available_relationships:
        if not (recurrence_type == 'each' or recurrence_type == rel):
            continue
        
        N_sibs_affected, N_sibs, N_affected, N, concordant, discordant = compute_recurrence(relationships[rel], empi2trait)
        
        print >> sys.stdout, ",".join(map(str, [rel, N_sibs_affected, N_sibs, N_affected, N, concordant, discordant]))

if __name__ == '__main__':
    
    args = dict([x.split('=') for x in sys.argv[1:]])
    
    if not 'demog' in args:
        args['demog'] = 'patient_demog_data_with_age.txt.gz'
    if not 'fam' in args:
        args['fam'] = 'family_ids.txt.gz'
    if not 'rel' in args:
        args['rel'] = 'actual_and_inf_rel_part1_unique_clean.txt.gz'
    if not 'type' in args:
        args['type'] = 'each'
    
    main(demographic_file = args['demog'],
        family_file = args['fam'],
        relationships_file = args['rel'],
        recurrence_type = args['type'],
        trait_path = args['trait'])
