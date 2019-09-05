#!/bin/bash

# Merge plink HD and LD sheep datasets, keep individual order as in original datasets
# merge-mode is default so mismatches are set to missing
# if two variants have the same position, try to merge them (fails if allele names arent matching)


plink --bfile ../sheep/data/SNP_chip/Plates_1-2_HD_QC2 --bmerge ../sheep/data/SNP_chip/Plates_1to87_QC3.bed ../sheep/data/SNP_chip/Plates_1to87_QC3.bim ../sheep/data/SNP_chip/Plates_1to87_QC3.fam --make-bed --indiv-sort 0 --sheep --merge-equal-pos --out ../sheep/data/SNP_chip/merged_sheep_geno
