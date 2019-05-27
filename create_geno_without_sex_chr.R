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
source("create_spec_file_AI.R")

### INPUT FOLDER ###
# for mac
# plink_geno_path <- "../sheep/data/SNP_chip/"
# for server, also contains the AlphaImputeLinux file and the Pedigree.txt
plink_geno_path <- "../sheep/data/SNP_chip/"
##############
output_path_main_files <- "genotypes_no_sex_chr/"

# get low density SNP chip data and take only SNPs from chromosome 1
sheep_bed <- paste0(plink_geno_path, "20180209_SoayPlates1-83.bed")
sheep_bim <- paste0(plink_geno_path, "20180209_SoayPlates1-83.bim")
sheep_fam <- paste0(plink_geno_path, "20180209_SoayPlates1-83.fam")
ld_sample <- read.plink(sheep_bed, sheep_bim, sheep_fam)

# coerce from raw to numeric
ld_sheep_geno <- as(ld_sample$genotypes, Class = "numeric")
ld_sheep_geno[is.na(ld_sheep_geno)] <- 9
ld_sheep_geno <- as_tibble(ld_sheep_geno, rownames = "ID")

# get high density SNP chip data 
hd_sheep_bed <- paste0(plink_geno_path, "20140214_SheepHD_QC1_Polym.bed")
hd_sheep_bim <- paste0(plink_geno_path, "20140214_SheepHD_QC1_Polym.bim")
hd_sheep_fam <- paste0(plink_geno_path, "20140214_SheepHD_QC1_Polym.fam")
hd_sample <- read.plink(hd_sheep_bed, hd_sheep_bim, hd_sheep_fam)

# coerce from raw to numeric
hd_sheep_geno <- as(hd_sample$genotypes, Class = "numeric")
hd_sheep_geno[is.na(hd_sheep_geno)] <- 9
hd_sheep_geno <- as_tibble(hd_sheep_geno, rownames = "ID")

# filter all but sex chromosome
all_chr_snps_hd <- hd_sample$map %>% filter((chromosome != 27) & (!is.na(chromosome))) %>% .$snp.name
all_chr_snps_ld <- ld_sample$map %>% filter((chromosome != 27) & (!is.na(chromosome))) %>% .$snp.name

#all_chr_snps_hd <- hd_sample$map %>% filter((chromosome == 18) & (!is.na(chromosome))) %>% .$snp.name
#all_chr_snps_ld <- ld_sample$map %>% filter((chromosome == 18) & (!is.na(chromosome))) %>% .$snp.name
                                            
# filter genotypes
hd_sheep_geno <- hd_sheep_geno[, c("ID", all_chr_snps_hd)]
ld_sheep_geno <- ld_sheep_geno[, c("ID", all_chr_snps_ld)]

# filter individuals from ld which are in hd
ld_sheep_geno_red <- ld_sheep_geno %>% filter(!(ID %in% hd_sheep_geno$ID))

# which loci have to be imputed (which are present in HD but NOT in LD data)
to_be_imputed <- names(hd_sheep_geno)[!(names(hd_sheep_geno) %in% names(ld_sheep_geno_red))]
to_be_imputed_index <- which(!(names(hd_sheep_geno) %in% names(ld_sheep_geno_red)))

# write to files
write_lines(to_be_imputed_index, path = paste0(output_path_main_files, "to_be_imputed_index.txt"))
#write_lines(to_be_imputed, path = paste0(output_path_main_files,"to_be_imputed.txt"))

rm(hd_sample)
rm(ld_sample)

# create table with both hd and ld data
#sheep_geno <- dplyr::bind_rows(hd_sheep_geno, ld_sheep_geno_red)
sheep_geno <- data.table::rbindlist(list(hd_sheep_geno, ld_sheep_geno_red), 
                                    use.names = TRUE, fill = TRUE)

# replace NAs with 9
repl_na <- function(DT) {
    for (j in seq_len(ncol(DT)))
        set(DT,which(is.na(DT[[j]])),j,9)
}
repl_na(sheep_geno)

# filter individuals which are not in pedigree due to some ID error
not_in_ped <- as.character(c(39,4302,9240,10446,10448,10449,10450,
                             10451,11076,11077,11079,11388))

sheep_geno_filt <- sheep_geno[!(ID %chin% not_in_ped)]

fwrite(sheep_geno_filt, paste0(output_path_main_files, "Genotypes.txt"), 
       sep = " ", col.names = TRUE)

#fwrite(hd_sheep_geno_t, "main_files_cv/Genotypes_snps_in_row.txt", sep = " ", 
#col.names = TRUE)



# check 
ind_1097 <- sheep_geno[ID==1097]
ind_1097[, to_be_imputed] <- 9
fwrite(ind_1097, file="ind_1097.txt", col.names = FALSE, sep = " ")

