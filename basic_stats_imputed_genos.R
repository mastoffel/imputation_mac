# full run with heuristic method and 1 to 5% 
library(data.table)
library(purrr)
library(tidyverse)
library(dplyr)
library(snpStats)

# plink name
sheep_plink_name <- "../sheep/data/SNP_chip/ramb_mapping/sheep_geno_imputed_ram_27092019"
# read merged plink data
sheep_bed <- paste0(sheep_plink_name, ".bed")
sheep_bim <- paste0(sheep_plink_name, ".bim")
sheep_fam <- paste0(sheep_plink_name, ".fam")
full_sample <- read.plink(sheep_bed, sheep_bim, sheep_fam)

summary(full_sample$genotypes)

sheep_geno <- as(full_sample$genotypes, Class = "numeric")
sheep_geno_log <- is.na(sheep_geno)
rm(sheep_geno)

missing_snps <- colSums(sheep_geno_log)
missing_snps_per_ind <- rowSums(sheep_geno_log)
hist(missing_snps_per_ind, breaks = 1000)
hist(missing_snps, breaks = 100)
rownames(sheep_geno_log)[missing_snps_per_ind > 100000]
