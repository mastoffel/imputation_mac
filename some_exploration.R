#### exploration
library(snpStats)
library(tidyverse)
library(sheepDB)
# LD inds
sheep_bed <- "../sheep/data/SNP_chip/20180209_SoayPlates1-83.bed"
sheep_bim <- "../sheep/data/SNP_chip/20180209_SoayPlates1-83.bim"
sheep_fam <- "../sheep/data/SNP_chip/20180209_SoayPlates1-83.fam"
# 
sample <- read.plink(sheep_bed, sheep_bim, sheep_fam)
ld_sheep_names <- rownames(sample$genotypes)
ld_sheep_loci <- colnames(sample$genotypes)

# HD inds
# get high density SNP chip data 
hd_sheep_bed <- "../sheep/data/SNP_chip/20140214_SheepHD_QC1_Polym.bed"
hd_sheep_bim <- "../sheep/data/SNP_chip/20140214_SheepHD_QC1_Polym.bim"
hd_sheep_fam <- "../sheep/data/SNP_chip/20140214_SheepHD_QC1_Polym.fam"

hd_sample <- read.plink(hd_sheep_bed, hd_sheep_bim, hd_sheep_fam)
hd_sheep_names <- rownames(hd_sample$genotypes)
hd_sheep_loci <- colnames(hd_sample$genotypes)

# how many ld loci are on the hd chip? ~ 80%
sum(hd_sheep_loci %in% ld_sheep_loci) / length(ld_sheep_loci)

# how many hd sheep were previously typed on the ld chip ? seems like all of them.
sum(hd_sheep_names %in% ld_sheep_names)

# are the SNP genotypes from ld and hd the same for the shared loci?
# to be checked soon

# pedigree
sheep_pedigree <- read_delim("../sheep/data/SNP_chip/20190208_Full_Pedigree.txt", delim = "\t")

# how many individuals can be perfectly imputed?
sub_ped_perfect <- pedigree %>%
    filter((MOTHER %in% hd_sheep_names) & (FATHER) %in% hd_sheep_names)

perfect_geno <- c(hd_sheep_names, sub_ped_perfect$ID)

sub_ped_perfect2 <- pedigree %>% 
    filter((MOTHER %in% perfect_geno) & (FATHER %in% perfect_geno))

perfect_geno2 <- c(perfect_geno, sub_ped_perfect2$ID)

sub_ped_perfect3 <- pedigree %>% 
    filter((MOTHER %in% perfect_geno2) & (FATHER %in% perfect_geno2))

perfect_geno3 <- c(perfect_geno2, sub_ped_perfect2$ID)

sub_ped_perfect4 <- pedigree %>% 
    filter((MOTHER %in% perfect_geno3) & (FATHER %in% perfect_geno3))


full_genotypes <- length(hd_sheep_names)


# start with the intial 188
hd_sheep <- hd_sheep_names
increase <- 1

while (increase != 0) {
    hd_sheep_before <- length(hd_sheep)
    # check individuals which have both parents from HD sheep
    hd_sheep_new <- sheep_pedigree %>%
        filter((MOTHER %in% hd_sheep) | (FATHER %in% hd_sheep))
    
    hd_sheep <- unique(c(hd_sheep_new$ID, hd_sheep))
    
    increase <- length(hd_sheep) - hd_sheep_before
    
}

# get sheep data
sheep_data <- get_sheep_data(dbpath = "../sheep/data/db/", datasets = c("CensusData", "Sheep"))
names(sheep_data)
# check pedigree
sheep_ped <- left_join(sheep_pedigree, select(sheep_data, ID, Birth , by = "ID"))
# shped <- kinship2::pedigree(id = sheep_pedigree$ID, dadid = sheep_pedigree$FATHER, 
#                             momid = sheep_pedigree$MOTHER, 
#                             sex = rep(3, length(sheep_pedigree$ID)) )


library(ggpedigree)
sum(duplicated(sheep_pedigree$ID))

ggpedigree(as.data.frame(sheep_pedigree), line.alpha = 0.1, point.alpha = 0.1)
