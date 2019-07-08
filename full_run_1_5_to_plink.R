# Gather imputed data from full run and make the plink file

# full run with heuristic method and 1 to 5% 
library(data.table)
library(purrr)
library(tidyverse)
library(dplyr)
library(snpStats)

#### imputed genotypes ####
read_imputed <- function(chr) {
    geno_imp <- fread(paste0("results/full_1_5/genos/ImputeGenotypes_", chr, ".txt"), header = FALSE)
    if ((chr != chrs[1])) geno_imp <- geno_imp[, -1]
    setDF(geno_imp)
}

# Read all imputed chromosomes
chrs <- 1:26
geno_imp_list <- map(chrs, read_imputed)

# chromosomes
chr_lengths <- map_chr(geno_imp_list, function(x) ncol(x))
chr_lengths <- as.numeric(chr_lengths)
chr_lengths[1] <- chr_lengths[1] - 1 # minus ID column

# colbind all and 
geno_imp <- setDT(unlist(geno_imp_list, recursive = FALSE), check.names = TRUE)[]
ind_ids <- geno_imp[[1]]
geno_imp[, V1:=NULL]

# set all genotypes which weren't imputed 100% as NA
# not necessary anymore as genotypes and not dosages are used
for(j in seq_along(geno_imp)){
    set(geno_imp, i = which(!(geno_imp[[j]] %in% c(0,1,2))), j=j, value = NA_real_)
}

# remove list and data.table
rm(geno_imp_list)
rm(geno_imp)

# take SNP names and individuals from other data
# load filtered genotypes
sheep_geno <- fread("data/genotypes_full_734_400638.txt", sep = " ", header = TRUE, fill = TRUE)
#ind_ids_2 <- sheep_geno[[1]]
snp_names <- names(sheep_geno)[-1]
rm(sheep_geno)

# transform to matrix
geno_imp_mat <- as.matrix(geno_imp)
rownames(geno_imp_mat) <- ind_ids
colnames(geno_imp_mat ) <- snp_names

# snpmatrix coding is 1,2,3
Sys.setenv('R_MAX_VSIZE'=52000000000)
geno_imp_mat  <- geno_imp_mat + 1
geno_imp_mat_snpstat <- new("SnpMatrix", geno_imp_mat)

# on mac
plink_geno_path <- "../sheep/data/SNP_chip/"
# plink name
sheep_plink_name <- "merged_sheep_geno"
# read merged plink data
sheep_bed <- paste0(plink_geno_path, sheep_plink_name, ".bed")
sheep_bim <- paste0(plink_geno_path, sheep_plink_name, ".bim")
sheep_fam <- paste0(plink_geno_path, sheep_plink_name, ".fam")
full_sample <- read.plink(sheep_bed, sheep_bim, sheep_fam)
full_sample$fam
full_sample$map

# filter individuals which are not in pedigree due to some ID error
not_in_ped <- as.character(c(39,4302,9240,10446,10448,10449,10450,
                             10451,11076,11077,11079,11388))

# plink puts double ids when merging, extract unique ids here
sheep_ids <- unlist(lapply(str_split(rownames(full_sample$fam), "\\.", 2), function(x) x[[2]]))
#
ind_data <- full_sample$fam %>% 
                mutate(id = sheep_ids) %>% 
                filter(!duplicated(id)) %>% 
                filter(!(id %in% not_in_ped))
snp_map <- full_sample$map %>% 
                filter(snp.name %in% snp_names)

full_sample$fam

# which individuals
write.plink(file.base = "data/sheep_geno_imputed_new",
            snps = geno_imp_mat_snpstat,
            #subject.data = ind_data,
            sex = ind_data$sex,
            id = ind_data$member,
            #snp.data = snp_map,
            chromosome = snp_map$chromosome,
            genetic.distance = snp_map$cM,
            position = snp_map$position,
            allele.1 = snp_map$allele.1,
            allele.2 = snp_map$allele.2,
            human.genome = FALSE)


?write.plink


# load plink full genos
# on mac
plink_geno_path <- "data/"
# plink name
sheep_plink_name <- "sheep_geno_imputed"
# read merged plink data
sheep_bed <- paste0(plink_geno_path, sheep_plink_name, ".bed")
sheep_bim <- paste0(plink_geno_path, sheep_plink_name, ".bim")
sheep_fam <- paste0(plink_geno_path, sheep_plink_name, ".fam")
full_sample <- read.plink(sheep_bed, sheep_bim, sheep_fam)
summary(full_sample$genotypes)




# some diagnostics
geno_imp_logical <- !is.na(geno_imp)

# how many SNPs per individual?
all_snps <- rowSums(!is.na(geno_imp))
all_snps_df <- data.frame(inds = ind_ids, imputed_snps = all_snps)

# how many individuals per SNP?
ind_per_snps <- colSums(!is.na(geno_imp))
sum(ind_per_snps < 2088)
all_ind_per_snps_df <- data.frame(inds_per_snp = ind_per_snps)

library(ggplot2)
p1 <- ggplot(all_snps_df, aes(imputed_snps)) + 
    geom_histogram(bins = 1000) +
    theme_minimal() +
    xlab("Imputed SNPs") +
    ylab("Individual count")
p1


p2 <- ggplot(all_ind_per_snps_df, aes(inds_per_snp)) + 
    geom_histogram(bins = 1000) +
    theme_minimal() +
    xlab("Imputed individuals") +
    ylab("SNP count")

p2

p_out <- p1/p2

# ggsave(p_out, filename = "Full_imputed_data_overview.jpg")

