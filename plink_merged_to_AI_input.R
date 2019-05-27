# Creating the genotype files and other files for AlphaImpute
# this script outputs three files:
# (1) Genotypes.txt: complete genotype data, can be filter for chromosomes
# (2) to_be_imputed_index.txt: index of SNPs which have to be masked and imputed
# (3) to_be_imputed.txt: names of SNPs
# (4) AlphaImputeLinux
# (5) Pedigree.txt

library(snpStats)
library(tidyverse)
library(data.table)
source("create_spec_file_merged.R")
#library(gdata)
#library("WGCNA")
# browseVignettes("snpStats")

### INPUT FOLDER ###
### Contains PLINK FILES, AlphaImputeLinux and Pedigree.txt ###

# on mac
plink_geno_path <- "../sheep/data/SNP_chip/"

### OUTPUT FOLDER ###

# on mac
#output_path_chr <- paste0("all_chr_cv/chr_", chr_num)
#output_path_main_files <- paste0(output_path_chr, "/AI_main_files/")

# plink name
sheep_plink_name <- "merged_sheep_geno"
# read merged plink data
sheep_bed <- paste0(plink_geno_path, sheep_plink_name, ".bed")
sheep_bim <- paste0(plink_geno_path, sheep_plink_name, ".bim")
sheep_fam <- paste0(plink_geno_path, sheep_plink_name, ".fam")
full_sample <- read.plink(sheep_bed, sheep_bim, sheep_fam)

# filter names of snps on one chromosome
all_chr_snps <- full_sample$map %>% filter(chromosome %in% 1:26) %>% .$snp.name

# filter those snps from full dataset and coerce from raw to numeric
sheep_geno_org <- as(full_sample$genotypes[, all_chr_snps], Class = "numeric")
dim(sheep_geno_org)
ids_org <- rownames(sheep_geno_org)
snp_names_org <- colnames(sheep_geno_org)

# plink puts double ids when merging, extract unique ids here
sheep_ids <- unlist(lapply(str_split(rownames(full_sample$genotypes), "\\.", 2), function(x) x[[2]]))
rownames(sheep_geno_org) <- sheep_ids

# clear some space
#rm(full_sample)

# make tibble and put rownames as ID column
sheep_geno <- as_tibble(sheep_geno_org, rownames = "ID")

###### merge individuals on both ld and hd chip #####
#dup_ids <-  which(duplicated(sheep_ids))
setDT(sheep_geno)


# no merging of individuals here which are on HD and LD chip
# only simply deleting the LD individuals which are on the hd chip too

sheep_geno <- sheep_geno[!duplicated(ID)]

# function to merge SNP data from the same individual, when it is both 
# in the HD and LD chip
# if genotypoe is missing on one chip, take the existing genotype
# if genotypes differ between chips, set NA

# merge_geno <- function(vec) {
#    # vec <- as.numeric(vec)
#     if (length(vec) == 1) return(as.numeric(vec))
#     if (sum(is.na(vec)) == 2) return(as.numeric(NA))
#     if (sum(is.na(vec)) == 1) return(vec[!is.na(vec)])
#     if (sum(is.na(vec)) == 0) {
#         if (vec[1] == vec[2]){
#             return(vec[1]) 
#         } else {
#             print(vec)
#             return(as.numeric(NA))
#         }
#     }
# }
# sheep_geno_merged <- sheep_geno[, lapply(.SD, merge_geno), by=ID]

# which SNPs are present in the LD but not the HD SNP chip and have to be imputed?
geno_missing <- colSums(is.na(sheep_geno))
# which SNPs are missing in more than 50% individuals (LD chip SNPs)
to_be_imputed <- names(geno_missing[geno_missing > 0.5 * nrow(sheep_geno)])
write_lines(to_be_imputed, path = "data/to_be_imputed_7386_400638.txt")

# filter individuals which are not in pedigree due to some ID error
not_in_ped <- as.character(c(39,4302,9240,10446,10448,10449,10450,
                             10451,11076,11077,11079,11388))

sheep_geno <- sheep_geno[!(ID %chin% not_in_ped)]

# write to file with col names for masking script
fwrite(sheep_geno, "genotypes_full_734_400638.txt", 
       sep = " ", col.names = TRUE)



