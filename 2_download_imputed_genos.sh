#!/bin/bash

# these are examples

# imputed genotypes
scp v1mstoff@eddie3.ecdf.ed.ac.uk:/exports/csce/eddie/biology/groups/pemberton/martin/sheep_imputation/cv_run_3_hap_2019/extract_results/geno_imp/* results/cv_3_hap_2019/

# original genotypes (unimputed)
scp v1mstoff@eddie3.ecdf.ed.ac.uk:/exports/eddie/scratch/v1mstoff/all_genos_org/* data/genos_unimputed_2018/

# download full run (not cv)
scp v1mstoff@eddie3.ecdf.ed.ac.uk:/exports/csce/eddie/biology/groups/pemberton/martin/sheep_imputation/full_1_5_2019/extract_results/genos/* results/full_1_5_2019/

# download 2020 run
scp mstoffel@eddie3.ecdf.ed.ac.uk:/exports/csce/eddie/biology/groups/pemberton/martin/sheep_imputation/full_1_5_2020/extract_results/genos/* results/full_1_5_2020/

scp mstoffel@eddie.ecdf.ed.ac.uk:/exports/csce/eddie/biology/groups/pemberton/martin/sheep_imputation/cv_run_full_2020/extract_results/geno_imp/*.txt results/full_1_5_2020/