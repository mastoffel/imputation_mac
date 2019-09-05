#!/bin/bash

# remove 3 individuals from HD chip (they are duplicates or have bad quality)

plink --bfile ../sheep/data/SNP_chip/Plates_1-2_HD_QC1 --remove ../sheep/data/SNP_chip/ids_to_remove --sheep --make-bed --out ../sheep/data/SNP_chip/Plates_1-2_HD_QC2
