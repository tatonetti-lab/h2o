"""
rel2ped.py 
Convert relationships to pedigree format. This script that inferences has been run on the relationships
data so that all possible pairs of relationships are available.

@author Nicholas P. Tatonetti and Fernanda Polubriaginof

USAGE
-----
# west site pedigree file
python src/rel2ped.py demog=patient_demog_data_with_age.txt.gz fam=family_ids.txt.gz rel=actual_and_inf_rel_part1_unique_clean.txt.gz out=west_generic_pedigree_file.txt.gz

# east site pedigree file
python src/rel2ped.py demog=east_patient_demog_data_with_age.txt.gz fam=east_family_ids.txt.gz rel=east_actual_and_inf_rel_part1_unique_clean.txt.gz out=east_generic_pedigree_file.txt.gz

NOTES
-----
To create the relationships file use 
bash src/build_relationship_file.sh
"""

import os
import csv
import sys
import gzip
from copy import deepcopy
from collections import defaultdict

from tqdm import tqdm

def add_new_relations_parent(relationships, e, mid, sex):
    for rel_type, relatives in relationships[e].items():
        mother2rel = None
        rel2mother = None
        if rel_type == 'Aunt/Uncle':
            continue
            #mother2rel = 'Sibling/Brother-in-law/Sister-in-law'
            #rel2mother = mother2rel
        elif rel_type == 'Child':
            mother2rel = 'Grandparent'
            rel2mother = 'Grandchild'
        elif rel_type == 'Child/Child-in-law' or rel_type == 'Child/Nephew/Niece':
            continue
            #mother2rel = 'Grandparent/Grandparent-in-law'
            #rel2mother = 'Grandchild/Grandchild-in-law'
        elif rel_type == 'Cousin':
            continue
            #mother2rel = 'Aunt/Uncle/Aunt-in-law/Uncle-in-law'
            #rel2mother = 'Nephew/Niece/Newphew-in-law/Niece-in-law'
        elif rel_type == 'First cousin once removed':
            continue            
        elif rel_type == 'Grandaunt/Granduncle':
            continue
        elif rel_type == 'Grandchild':
            mother2rel = 'Great-grandparent'
            rel2mother = 'Great-grandchild'
        elif rel_type == 'Grandnephew/Grandniece':
            continue
        elif rel_type == 'Grandparent':
            # This is ambiguous because we don't know if it's on the mother's or 
            # father's side. If we assume mother, then we can do the following:
            # if sex == 'F':
            #     # the relative is the parent of the mother
            #     rel2mother = 'Parent'
            #     # the mother is the child of the relative
            #     mother2rel = 'Child'
            continue
        elif rel_type == 'Great-grandaunt/Great-granduncle':
            continue
        elif rel_type == 'Great-grandchild':
            mother2rel = 'Great-great-grandparent'
            rel2mother = 'Great-great-grandchild'
        elif rel_type == 'Great-grandnephew/Great-grandniece':
            continue
        elif rel_type == 'Great-grandparent':
            continue
        elif rel_type == 'Great-great-grandchild':
            continue
        elif rel_type == 'Great-great-grandparent':
            continue
        elif rel_type == 'Nephew/Niece':
            continue
        elif rel_type == 'Parent':
            continue
        elif rel_type == 'Parent/Aunt/Uncle':
            continue
        elif rel_type == 'Parent-in-law':
            continue
        elif rel_type == 'Parent/Parent-in-law':
            continue
        elif rel_type == 'Sibling/Sibling-in-law':
            continue
        elif rel_type == 'Sibling':
            mother2rel = 'Parent'
            rel2mother = 'Child'
        elif rel_type == 'Sibling/Cousin':
            continue
        elif rel_type == 'Spouse':
            continue
        elif rel_type == 'Nephew/Niece/Nephew-in-law/Niece-in-law':
            continue
        elif rel_type == 'Grandchild/Grandchild-in-law':
            continue
        elif rel_type == 'Grandaunt/Granduncle/Grandaunt-in-law/Granduncle-in-law':
            continue
        elif rel_type == 'Grandnephew/Grandniece/Grandnephew-in-law/Grandniece-in-law':
            continue
        elif rel_type == 'Great-grandchild/Great-grandchild-in-law':
            continue
        elif rel_type == 'Great-grandparent/Great-grandparent-in-law':
            continue
        elif rel_type == 'Grandparent/Grandparent-in-law':
            continue
        elif rel_type == 'Aunt/Uncle/Aunt-in-law/Uncle-in-law':
            continue
        elif rel_type == 'Child-in-law':
            continue
        else:
            raise Exception("Unexpected relationship type: %s" % rel_type)
        
        if not mother2rel is None:
            for rel, provided in relatives:
                # check if this relationship has already been added
                if not (mid, False) in relationships[rel][mother2rel]:
                    relationships[rel][mother2rel].add((mid, False))
                
                if not (rel, False) in relationships[mid][rel2mother]:
                    relationships[mid][rel2mother].add((rel, False))
    
    relationships[e]['Parent'].add((mid, False))
    relationships[mid]['Child'].add((e, False))

def main(demographic_file, family_file, relationships_file, pedigree_outfile, data_directory='./common_data/'):
    
    print >> sys.stderr, "Loading patient demographic data..."
    
    # load demographic data (sex, birth_decade, race, ethnicity, age)
    fh = gzip.open(os.path.join(data_directory, demographic_file))
    reader = csv.reader(fh, delimiter='\t')
    reader.next()
    demog_header = ['sex', 'birthdec', 'race', 'ethnicity', 'age']
    
    demog_groups = defaultdict(set)

    empi2demog = dict()
    for row in tqdm(reader):
        
        empi = row[0]
        empi2demog[empi] = dict(zip(demog_header, row[1:]))
        
        # re-map the race codes
        if empi2demog[empi]['race'] == 'W':
            empi2demog[empi]['race'] = 'White'
        elif empi2demog[empi]['race'] == 'B':
            empi2demog[empi]['race'] = 'Black'
        elif empi2demog[empi]['race'] == 'O':
            empi2demog[empi]['race'] = 'Hispanic'
        elif empi2demog[empi]['race'] == 'U':
            empi2demog[empi]['race'] = 'Unknown'
        else:
            empi2demog[empi]['race'] = 'Other'
        
        for key, val in empi2demog[empi].items():
            demog_groups[key].add(val)
    
    print >> sys.stderr, "Loaded demographic data for %d patients." % len(empi2demog)
    
    # load family ids
    print >> sys.stderr, "Loading family identifiers and relationships..."
    
    fh = gzip.open(os.path.join(data_directory, family_file))
    reader = csv.reader(fh, delimiter='\t')
    reader.next()
    empi2fam = dict()
    fam2empi = defaultdict(set)
    for fid, iid in reader:
        empi2fam[iid] = fid
        fam2empi[fid].add( iid )
    fam2empi = dict(fam2empi)
    len(fam2empi), len(empi2fam)
    
    # load the relationships
    fh = gzip.open(os.path.join(data_directory, relationships_file))
    rel_data = csv.reader(fh, delimiter='\t')
    rel_data.next()

    relationships = defaultdict(lambda: defaultdict(set))

    for empi1, rel, empi2, provided in rel_data:
        rel = rel.strip()
        relationships[empi1][rel].add((empi2, True if provided == '1' else False))

    len(relationships)
    
    
    print >> sys.stderr, "Identifying missing fathers and mothers..."
    
    new_counter = 1
    for fam_id in tqdm(fam2empi.keys()):

        members = deepcopy(fam2empi[fam_id])

        for e in members:
            parents = relationships[e].get('Parent', None)
            mother_ids = set()
            father_ids = set()
        
            if parents is not None:
                for pid, provided in parents:
                    if empi2demog[pid]['sex'] == 'F':
                        mother_ids.add(pid)
                    elif empi2demog[pid]['sex'] == 'M':
                        father_ids.add(pid)
                    else:
                        #print "Unknown sex: %s, assuming mother" % empi2demog[pid]['sex']
                        mother_ids.add(pid)

            if len(mother_ids) == 0:
                #print "Could not find mother, creating new."
                mid = 'N%d' % new_counter
                new_counter += 1
            
                empi2demog[mid] = {'sex': 'F'}
                add_new_relations_parent(relationships, e, mid, 'F')
                fam2empi[fam_id].add(mid)
                empi2fam[mid] = fam_id
            
            if len(father_ids) == 0:
                #print "Could not find father, creating new."
                fid = 'N%d' % new_counter
                new_counter += 1
            
                empi2demog[fid] = {'sex': 'M'}
                add_new_relations_parent(relationships, e, fid, 'M')
                fam2empi[fam_id].add(fid)
                empi2fam[fid] = fam_id
    
    print >> sys.stderr, "Inferring %d missing parents." % (new_counter-1)
    
    print >> sys.stderr, "Putting the data in raw pedigree format..."
    
    ped_data = list()
    for fam_id in tqdm(fam2empi.keys()):

        members = deepcopy(fam2empi[fam_id])

        for e in members:
            parents = relationships[e].get('Parent', None)
            mother_ids = set()
            father_ids = set()
        
            if parents is not None:
                for pid, provided in parents:
                    if empi2demog[pid]['sex'] == 'F':
                        mother_ids.add((pid, provided))
                    elif empi2demog[pid]['sex'] == 'M':
                        father_ids.add((pid, provided))
                    else:
                        #print "Unknown sex: %s, assuming mother" % empi2demog[pid]['sex']
                        mother_ids.add((pid, provided))
        
            if len(mother_ids) == 0:
                mother_ids.add((0, False))
            if len(father_ids) == 0:
                father_ids.add((0, False))
        
            for mid, mp in mother_ids:
                for fid, fp in father_ids:
                    ped_data.append( (fam_id, e, mid, 'P' if mp else 'I', fid, 'P' if fp else 'I') )
    
    print >> sys.stderr, "Built pedigree for %d individuals." % len(ped_data)
    
    print >> sys.stderr, "Identifying multiple mother/multiple father conflicts..."
    
    mother2count = defaultdict(int)
    father2count = defaultdict(int)
    individual2mothers = defaultdict(set)
    individual2fathers = defaultdict(set)
    
    for fam_id, e, mid, mp, fid, fp in tqdm(ped_data):
        mother2count[mid] += 1
        father2count[fid] += 1
        individual2mothers[e].add(mid)
        individual2fathers[e].add(fid)
        
    print len([e for e, moms in individual2mothers.items() if len(moms) > 1])
    print len([e for e, dads in individual2fathers.items() if len(dads) > 1])
    
    pedigree_outfile_path = os.path.join(data_directory, pedigree_outfile)
    
    print >> sys.stderr, "Writing pedigree data to file: %s" % pedigree_outfile_path
    
    if pedigree_outfile_path.endswith('gz'):
        fh = gzip.open(pedigree_outfile_path, 'w')
    else:
        fh = open(pedigree_outfile_path, 'w')
    
    writer = csv.writer(fh, delimiter='\t')
    #writer.writerow(['family_id', 'individual_id', 'mother_id', 'mother_source', 'father_id', 'father_source'])
    writer.writerow(['family_id', 'individual_id', 'father_id', 'mother_id', 'own_ancestor'])
    
    for fam_id, e, mid, mp, fid, fp in tqdm(ped_data):
    
        if len(individual2mothers[e]) == 1:
            best_mid = mid
        else:
            best_mid = sorted([(mother2count[m],m) for m in individual2mothers[e]])[-1][1]
    
        if len(individual2fathers[e]) == 1:
            best_fid = fid
        else:
            best_fid = sorted([(father2count[f],f) for f in individual2fathers[e]])[-1][1]
    
        if mid == best_mid and fid == best_fid:
            #writer.writerow([fam_id, e, mid, mp, fid, fp])
            writer.writerow([fam_id, e, fid, mid, 0])
    fh.close()
    
    print >> sys.stderr, "Finished successfully."
    

if __name__ == '__main__':
    
    args = dict([x.split('=') for x in sys.argv[1:]])
    
    main(demographic_file = args['demog'],
        family_file = args['fam'],
        relationships_file = args['rel'],
        pedigree_outfile = args['out'])
    
    
