# Gather imputed data from full run and make the plink file

# full run with heuristic method and 1 to 5% 
library(data.table)
library(purrr)
library(tidyverse)
library(dplyr)
library(snpStats)

#### imputed genotypes ####
read_imputed <- function(chr) {
    geno_imp <- fread(paste0("results/full_1_5_2018/ImputeGenotypes_", chr, ".txt"), header = FALSE) # genos/
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
# rm(geno_imp)

#fwrite(geno_imp, file = "data/genotypes_2018.txt", sep = " ")

# take SNP names and individuals from other data
# load filtered genotypes
snp_names <- read_delim( "data/snp_names_merged_geno_2018.txt", delim=' ')
# sheep_geno <- fread("data/genotypes_2018.txt", sep = " ", header = TRUE, fill = TRUE)

# transform to matrix
geno_imp_mat <- as.matrix(geno_imp)
rownames(geno_imp_mat) <- ind_ids
colnames(geno_imp_mat ) <- snp_names$snps

rm(geno_imp)

# snpmatrix coding is 1,2,3
Sys.setenv('R_MAX_VSIZE'=52000000000)
geno_imp_mat  <- geno_imp_mat + 1
geno_imp_mat_snpstat <- new("SnpMatrix", geno_imp_mat)


# plink name
sheep_plink_name <- "../sheep/data/SNP_chip/merged_sheep_geno"
# read merged plink data
sheep_bed <- paste0(sheep_plink_name, ".bed")
sheep_bim <- paste0(sheep_plink_name, ".bim")
sheep_fam <- paste0(sheep_plink_name, ".fam")
full_sample <- read.plink(sheep_bed, sheep_bim, sheep_fam)

# # filter names of snps on one chromosome
# #all_chr_snps <- full_sample$map %>% filter(chromosome == chr_num) %>% .$snp.name
# 
# # filter those snps from full dataset and coerce from raw to numeric
# sheep_geno <- as(full_sample$genotypes[, all_chr_snps], Class = "numeric")
# 
# # plink puts double ids when merging, extract unique ids here
sheep_ids <- unlist(lapply(str_split(rownames(full_sample$genotypes), "\\.", 2), function(x) x[[2]]))
#rownames(full_sample$fam) <- sheep_ids
full_sample$fam$id <- sheep_ids

# make a df with ID and sex
full_sample$fam %>% 
    group_by(id) %>% 
    mean(sex, na.rm=TRUE)

# add sex and get the sequence of ids right
ind_df <- data.frame("id" = rownames(geno_imp_mat)) %>% 
    dplyr::left_join(full_sample$fam[!duplicated(full_sample$fam$id), ], by = "id")

# 
# # clear some space
# rm(full_sample)
# 
# # make tibble and put rownames as ID column
# sheep_geno <- as_tibble(sheep_geno, rownames = "ID")
# 
# ###### merge individuals on both ld and hd chip #####
# #dup_ids <-  which(duplicated(sheep_ids))
# setDT(sheep_geno)
# 
# # function to merge SNP data from the same individual, when it is both 
# # in the HD and LD chip
# # if genotypoe is missing on one chip, take the existing genotype
# # if genotypes differ between chips, set NA
# 
# merge_geno <- function(vec) {
#     # vec <- as.numeric(vec)
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
# 
# sheep_geno_merged <- sheep_geno[, lapply(.SD, merge_geno), by=ID]
# 
# # for the new pedigree (2018)
# not_in_ped <- as.character(c(7658, 7628, 7217, 5371, -112, -6, 1791, 5986, 7717))
# sheep_geno_filt <- sheep_geno_merged[!(ID %chin% not_in_ped)]
# 

# get SNP infos from plink merged file
snp_map <- full_sample$map %>%
                filter(snp.name %in% snp_names$snps)

# check whether SNPs names have the correct order
sum(snp_names$snps == snp_map$snp.name)

#full_sample$fam

#length(row.names(geno_imp_mat))
# which individuals
write.plink(file.base = "data/sheep_geno_imputed_04092019",
            snps = geno_imp_mat_snpstat,
            #subject.data = ind_data,
            sex = ind_df$sex,
            id = row.names(geno_imp_mat), 
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

