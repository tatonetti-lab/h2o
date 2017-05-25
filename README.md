# h2o
Observational Heritability

Code and notebooks for running observational heritability analysis and generating figures for the manuscript. The data required to run the notebooks and the examples are available at the following URL.

http://riftehr.tatonettilab.org/

## Install SOLAR

SOLARStrap requires that SOLAR is installed. You can do so from the following link:

https://www.nitrc.org/frs/?group_id=558

## Estimating Observational Heritability using SOLARStrap

Estimating observational heritability using SOLARStrap requires four data files. The files for `Rhinitis` are provided as an example at the supporting website. Retrieve them with `wget http://riftehr.tatonettilab.org/h2o_data.tar.gz`.

1. Trait file
2. Patient demographics
3. Family identifiers
4. Pedigree file

### Trait File

The trait file is a gzipped tab delimited file containing three columns: `ptid`, `pheno`, and `value`. The `ptid` is a local identifier for the patient. The `pheno` is a unique identifier for the trait (e.g. an ICD9/10 code, a LOINC code, or a local identifier). The `value` is either `1`, `0`, or `NA` if the trait is dichotomous or the quantitative value if the trait is quantiative (again with `NA` used if the value is unknown). Note that all individuals in the pedigree that are not in the trait file will be assigned `NA` as the trait value.

More than one trait can be coded in a single trait file. In this case SOLARStrap will run heritability estimates for all traits in the file. 

`gzcat example/rhinitis/67_trait_data.txt.gz | head`

```
ptid	pheno	value
3000213365	67	1
3000267263	67	1
3000127670	67	1
3000084624	67	1
3000347920	67	1
3000179426	67	1
3000197653	67	1
3000272716	67	1
3000360729	67	1
```

### Patient Demographics

The patient demographics file is a gzipped tab delimited file containing 6 columns: `ptid`, `sex`, `birth_decade`, `race_code`, `ethnic_code`, and `age`. If any of the information is unknown it is coded as `NULL`. Birth decade is not used at this time. `race_code` may be one of `NA` or `U` (unknown), `W` (white), or `B` (black). `ethnic_code` represents self-reporting of hispanic or not hispanic -- `H`, `S`, `O` will all be mapped to hispanic, other codes to non-hispanic. The `race_code` and the `ethnic_code` are used together to group patients and families into self-reported race/ethnicity groups. An `ethnic_code` of hispanic overrides the reported race and a `race_code` will be mapped to hispanic. Note thate this logic is built for the demographic diversity of New York City. Uses in other areas may require different logic. In this case you may want to edit the `load_demographics` method in `h2o_utility.py`. 

```
ptid	sex	birth_decade	race_code	ethnic_code	age
3000149927      F       NULL    W       N       68
3000191095      M       NULL    B       N       75
3000276396      M       NULL    W       N       63
3000181834      F       NULL    W       H       55
3000130794      F       NULL    W       H       60
3000341624      F       NULL    O       D       13
3000267077      M       NULL    W       H       56
3000089601      F       NULL    W       H       52
3000234386      F       NULL    B       N       45
3000176550      M       NULL    B       N       70
```

### Family Identifiers

The family identifiers file is a gzipped tab delimited file containing two columns: `famid` and `ptid`. 

```
famid   ptid
30186335        3000313768
30186335        3000313771
30186335        3000313769
30186335        3000313770
30186335        3000313767
30186334        3000313763
30186334        3000313762
30186334        3000313765
30186334        3000313764
30186334        3000313766
30186337        3000313773
```

### Pedigree file

The pedigree file is a gzipped tab delimited file containing X columns: `family_id`, `individual_id`, `father_id`, `mother_id`, and `own_ancestor`. The `own_ancestor` column is now deprecated but must, for now, still be present. This column can be safely set to `0` for every row.

```
family_id	individual_id	father_id	mother_id	own_ancestor
30000851	3000068718	3000370355	3000068717	0
30000851	3000068717	3000370356	3000370357	0
30000851	3000068716	3000370355	3000068717	0
30000854	3000068719	3000370358	3000370359	0
30000854	3000068721	3000370360	3000068720	0
30000854	3000068720	3000370361	3000370362	0
30000854	3000370358	0	0	0
30000854	3000370359	0	0	0
30000854	3000370361	0	0	0
```

### Running SOLARStrap

There are five required arguments when running SOLARStrap, `trait`, `demog`, `fam`, `ped`, and `type`. The first four are paths to the trait, demographics, family ids, and pedigree files (as described above). `type` is either `D` for dichotomous or `Q` for quantitative. Here is an example using rhinitis. You will be prompted to make a `./working` directory, do so with `mkdir working`.

```
python src/solarStrap_heritability.py trait=example/rhinitis/67_trait_data.txt.gz demog=example/rhinitis/patient_demog_data_with_age.txt.gz fam=example/rhinitis/family_ids.txt.gz ped=example/rhinitis/west_generic_pedigree_file.txt.gz type=D
```

The remaining arguments are optional. They are 
```
- `ace`       Set to `yes` to run estimates modeling the household effect. The mother id will be used as the household id. Default is `no`.
- `verbose`   Set to `yes` to see a lot more output. Default is `no`.
- 'nfam`      The number of families to be sampled in each iteration, can be a float (proportion of total families) or an integer. Default is `0.15` (15%).
- `samples`   The number of iterations to run. Default is `200`. 
- `buildonly` Set to yes to only build the directories, do not run SOLAR. Default is `no`.
- `proband`   Set to yes to use a proband. Default is `yes`.
- `sd`        The path to the working directory. Default is `./working`.
```

### Output files

In addition to the error and log messages that are printed to the screen, SOLARStrap will produce two results files: (1) a file of aggregated results (these are the h2o estimates) and (2) a log of all of the heritability estimates for each iteration. The former is aggregated version of the latter. 

These files will be placed in the working directory provided to SOLARStrap when run (default is `./working`). In addition, when SOLARStrap is run it will print out the path where these files will be saved and their names. For example:

`python src/solarStrap_heritability.py trait=example/rhinitis/67_trait_data.txt.gz demog=example/rhinitis/patient_demog_data_with_age.txt.gz fam=example/rhinitis/family_ids.txt.gz ped=example/rhinitis/west_generic_pedigree_file.txt.gz type=D verbose=yes ace=yes`

```
SolarStrap v 0.9 - Estimate heritability of disease using observational data.
-----------------------------------------------------------------------------
Summary results will be saved in ./working/F6ZF3_solar_strap_results.csv
Results from each bootstrap will be saved at ./working/F6ZF3_solar_strap_allruns.csv
```

If you run solar in verbose mode you will see a log of the h2 estimats from SOLAR if SOLARStrap is running successfully.

```
Number of families with case: 3686
     Trait       Ethnicity  NFam  Samp   AE h2     err       pval  ACE h2     err       pval Sample AFP
        67             ALL   552     1    1.00     nan   2.90e-06    0.86    0.09   1.06e-04     1.3207
        67             ALL   552     2    0.17    0.09   1.73e-01    0.10    0.31   2.96e-01     1.2609
        67             ALL   552     3    1.00     nan   2.90e-06    0.87    0.09   3.58e-05     1.2373
        67             ALL   552     4    1.00     nan   1.50e-06    0.83    0.09   7.70e-05     1.2790
        67             ALL   552     5    1.00     nan   3.20e-06    0.77    0.10   9.62e-05     1.2627
        67             ALL   552     6    1.00     nan   1.04e-08    0.86    0.06   4.60e-06     1.2880
        67             ALL   552     7    1.00     nan   1.34e-09    0.87    0.10   5.00e-07     1.2754
        67             ALL   552     8    1.00     nan   1.30e-06    0.77    0.07   8.95e-05     1.2681
...
```
