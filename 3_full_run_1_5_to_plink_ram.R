# Gather imputed data from full run and make the plink file

# full run with heuristic method and 1 to 5% 
library(data.table)
library(purrr)
library(tidyverse)
library(dplyr)
library(snpStats)

#### imputed genotypes ####
read_imputed <- function(chr) {
    geno_imp <- fread(paste0("results/full_1_5_2019/ImputeGenotypes_", chr, ".txt"), header = FALSE) # genos/
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

# take SNP names and individuals from other data
# load filtered genotypes
snp_names <- read_delim( "data/snp_names_merged_geno_2019.txt", delim=' ')
# sheep_geno <- fread("data/genotypes_2018.txt", sep = " ", header = TRUE, fill = TRUE)

# transform to matrix
geno_imp_mat <- as.matrix(geno_imp)
rownames(geno_imp_mat) <- ind_ids
colnames(geno_imp_mat ) <- snp_names$snps

# make space
rm(geno_imp)

# snpmatrix coding is 1,2,3
Sys.setenv('R_MAX_VSIZE'=52000000000)
geno_imp_mat  <- geno_imp_mat + 1
geno_imp_mat_snpstat <- new("SnpMatrix", geno_imp_mat)

# plink name
sheep_plink_name <- "../sheep/data/SNP_chip/ramb_mapping/merged_sheep_geno_ram"
# read merged plink data
sheep_bed <- paste0(sheep_plink_name, ".bed")
sheep_bim <- paste0(sheep_plink_name, ".bim")
sheep_fam <- paste0(sheep_plink_name, ".fam")
full_sample <- read.plink(sheep_bed, sheep_bim, sheep_fam)

not_in_ped <- as.character(c(7658, 7628, 7217, 5371, -112, -6, 1791, 5986, 7717))

subject_data <- full_sample$fam %>% 
                    rename(id = member,
                           phenotype = affected) %>% 
                    filter(!(id %in% not_in_ped))
rownames(subject_data) <- subject_data$id
                    

snp_data <- full_sample$map %>% 
                filter(chromosome %in% c(1:26)) %>% 
                dplyr::select(-snp.name) %>% 
                dplyr::rename(genetic.distance = cM)
rownames(snp_data) <- full_sample$map %>% 
                      filter(chromosome %in% c(1:26)) %>%
                      .[["snp.name"]]
#rm(geno_imp_mat)
#rm(full_sample)
# which individuals
# colnames(geno_imp_mat_snpstat) <- snp_names$snp.name
write.plink(file.base = "data/sheep_geno_imputed_ram_27092019",
            snps = geno_imp_mat_snpstat,
            pedigree = subject_data$pedigree,
            id = subject_data$id,
            father = subject_data$father,
            mother = subject_data$mother,
            sex = subject_data$sex,
            phenotype = subject_data$phenotype,
            chromosome = snp_data$chromosome,
            genetic.distance = snp_data$genetic.distance,
            position = snp_data$position,
            allele.1 = snp_data$allele.1,
            allele.2 = snp_data$allele.2,
            human.genome = FALSE)

write.plink(file.base = "data/sheep_geno_imputed_ram_27092019",
            snps = geno_imp_mat_snpstat,
            subject.data = subject_data,
            snp.data = snp_data,
            human.genome = FALSE)


# load plink full genos
# on mac
# plink name
sheep_plink_name <- "data/sheep_geno_imputed_04092019"
# read merged plink data
sheep_bed <- paste0(sheep_plink_name, ".bed")
sheep_bim <- paste0(sheep_plink_name, ".bim")
sheep_fam <- paste0(sheep_plink_name, ".fam")
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

