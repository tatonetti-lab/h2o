"""
rel2grm.py
Convert relationships to a sparse genetic relationship matrix. This matrix can then be loaded and used to run GCTA.

@author Nicholas P. Tatonetti

USAGE
-----
python src/rel2grm.py demog=patient_demog_data_with_age.txt.gz fam=family_ids.txt.gz rel=actual_and_inf_rel_part1_unique_clean.txt.gz out=west_generic_grm_file.txt.gz

"""

import os
import csv
import sys
import gzip
import scipy.sparse
from tqdm import tqdm

from h2o_utility import load_family_ids, load_relationships_by_pairs

def main(demographic_file, family_file, relationships_file, grm_outfile, data_directory='./common_data/'):
    
    if not grm_outfile.endswith('.txt.gz'):
        print >> sys.stderr, "ERROR: GRM outfile must end with .txt.gz"
        sys.exit(10)
    
    empi2fam, fam2empi = load_family_ids(os.path.join(data_directory, family_file))
    relationships = load_relationships_by_pairs(os.path.join(data_directory, relationships_file))
    
    # the grm must be sorted by the patient ids
    patients = sorted(empi2fam.keys())
    
    gcta_grm = list()
    gcta_grm_id = list()
    print >> sys.stderr, "Building the GRM..."
    
    for i in tqdm(range(len(patients))):
        pid1 = patients[i]
        gcta_grm_id.append([i, pid1])
        for j, pid2 in enumerate(patients):
            if i == j:
                gcta_grm.append( [j, i, 1.0] )
            elif i < j:
                if empi2fam[pid1] == empi2fam[pid2]:
                    rel = relationships[pid1].get(pid2, 0.0)
                    if rel != 0.0:
                        gcta_grm.append( [j, i, rel] )
        
    fn = os.path.join(data_directory, grm_outfile)
    print >> sys.stderr, "Saving sparse matrix to file: %s..." % fn
    fh = gzip.open(fn, 'w')
    writer = csv.writer(fh, delimiter='\t', quoting=csv.QUOTE_NONE)
    writer.writerows(gcta_grm)
    fh.close()
    
    fn = os.path.join(data_directory, grm_outfile.replace('.txt.gz', '.id.txt'))
    print >> sys.stderr, "Saving sparse matrix index file: %s..." % fn
    fh = open(fn, 'w')
    writer = csv.writer(fh, delimiter='\t', quoting=csv.QUOTE_NONE)
    writer.writerows(gcta_grm_id)
    fh.close()
    
    print >> sys.stderr, "Finished."
            
if __name__ == "__main__":
    args = dict([x.split('=') for x in sys.argv[1:]])

    main(demographic_file = args['demog'],
        family_file = args['fam'],
        relationships_file = args['rel'],
        grm_outfile = args['out'])
