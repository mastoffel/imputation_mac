# file to run lots of alpha impute runs
library(dplyr)
library(readr)
library(data.table)

# folder with all the files
main_files <- "main_files_cv/"
# load genotypes
geno_org1 <- read_delim(paste0(main_files, 'Genotypes_snps_in_row.txt'), delim = " ")
ids <- colnames(geno_org1)[-1]
geno_org <- data.table::transpose(geno_org1)
colnames(geno_org) <- geno_org[1, ]
geno_org <- geno_org[-1, ]
geno_org <- cbind(data.frame("ID" = ids), geno_org)

#geno_org <- data.table::fread(paste0(main_files, 'full_hd_genotypes.txt'), sep = " ")

# geno_org <- read_delim(paste0(main_files, 'full_hd_genotypes.txt'), delim = " ")
# load SNPs to mask
to_be_imputed <- read_lines(paste0(main_files, 'to_be_imputed.txt'))


create_cv_folder <- function(run_num) {
    
    dir.create(paste0("run_", run_num))
    #file.copy(from = paste0(main_files, "AlphaImputeLinux"), to = paste0("run_", run_num))
    file.copy(from = paste0(main_files, "AlphaImputeSpec.txt"), to = paste0("run_", run_num))
    
    # mask genotypes of 18 random individuals
    # set seed to know which individuals
    geno <- geno_org
    set.seed(run_num)
    geno[sample(1:nrow(geno), 18), to_be_imputed] <- 9
    write_delim(geno, paste0("run_", run_num, "/Genotypes.txt"), delim = " ", col_names = FALSE)
}

sapply(1:2, create_cv_folder)


