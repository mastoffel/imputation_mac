# Gather imputed data from full run and make the plink file

# full run with heuristic method and 1 to 5% 
library(data.table)
library(purrr)
library(tidyverse)
library(dplyr)
library(snpStats)

#### imputed genotypes ####
read_imputed <- function(chr) {
    geno_imp <- fread(paste0("results/full_1_5_2020/ImputeGenotypes_", chr, ".txt"), header = FALSE) # genos/
    #if ((chr != chrs[1])) geno_imp <- geno_imp[, -1]
    geno_imp <- geno_imp[, -1]
    geno_imp <- as.matrix(geno_imp) + 1
    geno_imp <- as.data.table(geno_imp)
    #as.matrix(geno_imp)
}

# get individual ids from one of the genotype files
ind_ids <- fread(paste0("results/full_1_5_2020/ImputeGenotypes_1.txt"), header = FALSE, select = 1) %>% 
             unlist()

# Read all imputed chromosomes
chrs <- 1:26
geno_imp_list <- map(chrs, read_imputed)

# chromosomes
chr_lengths <- map_chr(geno_imp_list, function(x) ncol(x))
chr_lengths <- as.numeric(chr_lengths)
#chr_lengths[1] <- chr_lengths[1] - 1 # minus ID column

# colbind all and 
geno_imp <- setDT(unlist(geno_imp_list, recursive = FALSE), check.names = TRUE)[]
#ind_ids <- geno_imp[[1]]
#geno_imp[, V1:=NULL]
rm(geno_imp_list)
# set all genotypes which weren't imputed 100% as NA
# not necessary anymore as genotypes and not dosages are used
#for(j in seq_along(geno_imp)){
#    set(geno_imp, i = which(!(geno_imp[[j]] %in% c(0,1,2))), j=j, value = NA_real_)
#}

# take SNP names and individuals from other data
# load filtered genotypes
# plink name
snp_names <- read_delim("../sheep/data/SNP_chip/ramb_mapping/merged_sheep_geno_ram_snpfilt.bim", 
           delim = "\t", col_names = FALSE) %>% 
            filter(X1 %in% c(1:26)) %>% 
            .$X2

#snp_names <- read_delim( "data/snp_names_merged_geno_2019.txt", delim=' ')
# sheep_geno <- fread("data/genotypes_2018.txt", sep = " ", header = TRUE, fill = TRUE)
ncols_dt <- 1:ncol(geno_imp)
fwrite(geno_imp, "output/geno_imp.txt", row.names = FALSE, col.names = FALSE)

rm(geno_imp)

ncols_dt <- 1:397867
col_inds <- split(ncols_dt, ceiling(seq_along(ncols_dt)/10000))
read_mats <- function(col_ind) {
    tbl <- fread("output/geno_imp.txt", select = col_ind, na.strings = "10")
    as.matrix(tbl)
}

geno_imp <- do.call(cbind, lapply(col_inds, read_mats))
# transform to matrix
#geno_imp <- as.matrix(geno_imp)
rownames(geno_imp) <- ind_ids
colnames(geno_imp) <- snp_names

# make space
#rm(geno_imp)

# snpmatrix coding is 1,2,3
#Sys.setenv('R_MAX_VSIZE'=52000000000)
#geno_imp  <- geno_imp + 1
geno_imp_mat_snpstat <- new("SnpMatrix", geno_imp)
rm(geno_imp)

# plink name
sheep_plink_name <- "../sheep/data/SNP_chip/ramb_mapping/merged_sheep_geno_ram_snpfilt"
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

# after saving, check that MAP has saved correctly
# for a SNP with position 172000000 it saved a scientific number once
options(scipen=999)
write.plink(file.base = "data/sheep_geno_imputed_ram_01052020",
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

# check that no scientific number has been saved and correct
library(data.table)
options(scipen=999)
bim <- fread("data/sheep_geno_imputed_ram_01052020.bim")
bim[which(bim$V2 == "oar3_OAR3_160215438"), ]
#bim[which(bim$V2 == "oar3_OAR3_160215438"), "V2"] <- 172000000
fwrite(bim, file = "data/sheep_geno_imputed_ram_01052020.bim", col.names = FALSE,
        sep = "\t")

# write.plink(file.base = "data/sheep_geno_imputed_ram_27092019",
#             snps = geno_imp_mat_snpstat,
#             subject.data = subject_data,
#             snp.data = snp_data,
#             human.genome = FALSE)


# load plink full genos
# on mac
# plink name
sheep_plink_name <- "data/sheep_geno_imputed_ram_01052020"
# read merged plink data
sheep_bed <- paste0(sheep_plink_name, ".bed")
sheep_bim <- paste0(sheep_plink_name, ".bim")
sheep_fam <- paste0(sheep_plink_name, ".fam")
full_sample <- read.plink(sheep_bed, sheep_bim, sheep_fam)
full_sample$map %>% as_tibble()
summary(full_sample$genotypes)


sheep_plink_name <- "../sheep_ID/data/sheep_geno_imputed_ram_400k_filt"
# read merged plink data
sheep_bed <- paste0(sheep_plink_name, ".bed")
sheep_bim <- paste0(sheep_plink_name, ".bim")
sheep_fam <- paste0(sheep_plink_name, ".fam")
full_sample_old <- read.plink(sheep_bed, sheep_bim, sheep_fam)
full_sample$map %>% as_tibble()
summary(full_sample$genotypes)

geno_old <- full_sample_old$genotypes[full_sample_old$fam$member %in% full_sample$fam$member, ]
geno_new <- full_sample$genotypes[full_sample$fam$member %in% full_sample_old$fam$member, ]
geno_old
geno_new

geno_old_re <- geno_old[, match(colnames(geno_new), colnames(geno_old))]

out <- purrr::map_int(1:ncol(geno_old_re), function(x) sum(geno_old_re[, x] == geno_new[, x], na.rm = TRUE))  
plot(out)
out_sub <- which(out<3500)

geno_old_re[, out_sub[1000]]
geno_new[, out_sub[1000]]

snp1 <- as(geno_old_re[, out_sub[1000]], "numeric")
snp2 <- as(geno_new[, out_sub[1000]], "numeric")

sum(snp1 == snp2, na.rm = TRUE)

df <- cbind(snp1, snp2)
head(df)

table(snp1)
table(snp2)

match(colnames(geno_old), colnames(geno_new))
sum(colnames(geno_old) == colnames(geno_new))

x <- 100000
sum(geno_old_re[, x] == geno_new[, x], na.rm = TRUE)

geno_imp <- full_sample$genotypes

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

