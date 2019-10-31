# X imputation didn't work properly. This is because the mapping to 
# the rambouillet genome does not seem very good. 

library(snpStats)
library(tidyverse)
library(data.table)
source("create_spec_file_merged_sexchr.R")
source("../sheep_id/theme_clean.R")

chr_num <- 27

# on mac
plink_geno_path <- "../sheep/data/SNP_chip/ramb_mapping/"
#####################

############################ CHECK HD chip #####################################

# plink name
#sheep_plink_name <- "merged_sheep_geno_ram"
sheep_plink_name <- "../sheep/data/SNP_chip/ramb_mapping/Plates_1-2_HD_QC3_ram"
#sheep_plink_name <- "../sheep/data/SNP_chip/Plates_1-2_HD_QC2"

# read merged plink data
sheep_bed <- paste0(sheep_plink_name, ".bed")
sheep_bim <- paste0(sheep_plink_name, ".bim")
sheep_fam <- paste0(sheep_plink_name, ".fam")
full_sample_hd <- read.plink(sheep_bed, sheep_bim, sheep_fam)
inds_hd <- full_sample_hd$fam
snps_on_x <- full_sample_hd$map %>% filter(chromosome == 27) 

geno_x <- as_tibble(as(full_sample_hd$genotypes[, snps_on_x$snp.name], Class = "numeric"),
                      rownames = "id") %>% 
          # filter only males
          filter(inds_hd$sex == 1)

geno_x_mat <- geno_x %>% 
    select(-id) %>% 
    as.matrix() 

geno_x_mat[!(geno_x_mat == 1)] <- NA
all_het <- colSums(geno_x_mat, na.rm = TRUE)

df <- tibble(snp.name = names(all_het), num_het = all_het) %>% 
      left_join(snps_on_x, by = "snp.name")# %>% 
      #left_join(par_snps) %>% 
      #mutate(par_snp = replace_na(as.character(par_snp), 0))

# whats marked as par_snp already and where should the line be?  
df %>% ggplot(aes(position, num_het)) + geom_point() +
        xlim(c(0, 10000000)) +
        geom_vline(xintercept = 7700000) # probably at 7700000


# Here we do several filtering steps to not include PAR SNPs in the imputation
# even when they are wrongly mapped 

# save new PAR snps
df %>% filter(position < 7700000) %>% select(snp.name) -> par_snps_hd
par_snps_hd %>% write_delim("data/Ram_PAR_SNPs_HD.txt", col_names = FALSE)

# find SNPs that were previously unmapped and which are now wrongly mapped
unmapped_snps <- read_delim("data/old_unmmapped_X_hd.txt", " ") %>% 
                  rename(snp.name = snp)

df %>% left_join(unmapped_snps, by = "snp.name") %>% 
  ggplot(aes(position, num_het, color = factor(pos_old))) + geom_point()

# they are potentially par snps and need to be removed 
par_snps_hd_pot1 <- unmapped_snps %>% select(snp.name)
par_snps_hd_pot1 %>%  write_delim("data/Ram_PAR_SNPs_HD_pot1.txt", col_names = FALSE)

# find SNPs with lots of heterozygotes 
df %>% filter(!(snp.name %in% par_snps_hd$snp.name)) %>% 
  filter(!(snp.name %in% unmapped_snps$snp.name)) %>% 
  filter(num_het > 0)  %>% select(snp.name) -> par_snps_hd_pot2

par_snps_hd_pot2 %>% write_delim("data/Ram_PAR_SNPs_HD_pot2.txt", col_names = FALSE)

# all SNPs to be removed from HD because potentially PAR
par_all_hd <- rbind(par_snps_hd, par_snps_hd_pot1, par_snps_hd_pot2) %>% 
           write_delim("data/Ram_PAR_SNPs_HD_ALL.txt", col_names = FALSE)

# check again
df %>% filter(!(snp.name %in% par_all$snp.name)) %>% 
  ggplot(aes(position, num_het)) + geom_point()
# which snps were unmapped previously?




############### Same exercise for LD Chip #########################

# plink name
#sheep_plink_name <- "merged_sheep_geno_ram"
sheep_plink_name <- "../sheep/data/SNP_chip/ramb_mapping/Plates_1to87_QC4_ram"
#sheep_plink_name <- "../sheep/data/SNP_chip/Plates_1-2_HD_QC2"

# read merged plink data
sheep_bed <- paste0(sheep_plink_name, ".bed")
sheep_bim <- paste0(sheep_plink_name, ".bim")
sheep_fam <- paste0(sheep_plink_name, ".fam")
full_sample_ld <- read.plink(sheep_bed, sheep_bim, sheep_fam)
inds_ld <- full_sample_ld$fam
snps_on_x <- full_sample_ld$map %>% filter(chromosome == 27) 

geno_x <- as_tibble(as(full_sample_ld$genotypes[, snps_on_x$snp.name], Class = "numeric"),
                    rownames = "id") %>% 
  # filter only males
  filter(inds_ld$sex == 1)

geno_x_mat <- geno_x %>% 
  select(-id) %>% 
  as.matrix() 

geno_x_mat[!(geno_x_mat == 1)] <- NA
all_het <- colSums(geno_x_mat, na.rm = TRUE)

df <- tibble(snp.name = names(all_het), num_het = all_het) %>% 
  left_join(snps_on_x, by = "snp.name") 
  #left_join(par_snps) %>% 
 # mutate(par_snp = replace_na(as.character(par_snp), 0))

# whats marked as par_snp already and where should the line be?  
df %>% ggplot(aes(position, num_het)) + geom_point() +
  xlim(c(0, 10000000)) +
  geom_vline(xintercept = 7700000) # probably at 7700000

# Here we do several filtering steps to not include PAR SNPs in the imputation
# even when they are wrongly mapped 

# save new PAR snps
df %>% filter(position < 7700000) %>% select(snp.name) -> par_snps_ld
par_snps_ld %>% write_delim("data/Ram_PAR_SNPs_LD.txt", col_names = FALSE)

# find SNPs that were previously unmapped and which are now wrongly mapped
unmapped_snps <- read_delim("data/old_unmmapped_X_ld.txt", " ") %>% 
  rename(snp.name = snp)

df %>% left_join(unmapped_snps, by = "snp.name") %>% 
  ggplot(aes(position, num_het, color = factor(pos_old))) + geom_point()

# they are potentially par snps and need to be removed 
par_snps_ld_pot1 <- unmapped_snps %>% select(snp.name)
par_snps_ld_pot1 %>%  write_delim("data/Ram_PAR_SNPs_LD_pot1.txt", col_names = FALSE)

# find SNPs with lots of heterozygotes (not necessary here)
# df %>% filter(!(snp.name %in% par_snps_ld$snp.name)) %>% 
#   filter(!(snp.name %in% unmapped_snps$snp.name)) %>% 
#   filter(num_het > 50)  %>% select(snp.name) -> par_snps_ld_pot2
# par_snps_ld_pot2 %>% write_delim("data/Ram_PAR_SNPs_LD_pot2.txt", col_names = FALSE)

# all SNPs to be removed from HD because potentially PAR
par_all_ld <- rbind(par_snps_ld, par_snps_ld_pot1) %>% 
  write_delim("data/Ram_PAR_SNPs_LD_ALL.txt", col_names = FALSE)

# check again
df %>% filter(!(snp.name %in% par_all$snp.name)) %>% 
  ggplot(aes(position, num_het)) + geom_point()
# which snps were unmapped previously?




# combine LD and HD PAR snps into one

par_snps_ldhd <- rbind(par_snps_hd, par_snps_ld) %>% 
                  distinct() %>%  write_delim("data/Ram_PAR_SNPs_HDLD_all.txt", col_names = FALSE)












inds <- full_sample$fam %>% rename(id = member) 

wrong_snp_names <- read_delim("../sheep/data/sheep_genome/susie_ld_mapping_output/All_Chromosome_Common_SNPs_QCed.txt", 
                              delim = "\t") 

mult_snps <- read_delim("../sheep_imputation/data/multimapped_X.txt", delim = " ") 
mult_snps_df <- tibble(snp.name = unique(mult_snps$snp), mult = 1)
# pseudoautosomal region SNPs
#par_snps <- readLines("/exports/csce/eddie/biology/groups/pemberton/martin/plink_genotypes/Oar3.1_PAR_SNPs_HD.txt")
par_snps <- read_delim("../sheep_imputation/data/Oar3.1_PAR_SNPs_HD.txt", " ")[[1]] %>% 
    as_tibble() %>% 
    mutate(par_snp = as.factor(1)) %>% 
    rename(snp.name = value)

# check hh file plink
hh_file <- read_delim("../sheep/data/SNP_chip/ramb_mapping/Plates_1-2_HD_QC3_ram.hh", "\t",
                      col_names = FALSE) %>% 
    rename(id = X2, snp.name = X3)


check_X <- hh_file %>% inner_join(inds) %>% 
    # filter(sex == 1) %>% 
    right_join(snps_on_x) %>% 
    group_by(snp.name) %>% 
    summarise(num_ind_het = n()-1, pos = first(position)) %>% 
    left_join(par_snps) %>% 
    mutate(par_snp = ifelse(is.na(par_snp), 0, 1)) %>% 
    mutate(par_snp = as.factor(par_snp)) 

library(wesanderson)
options(scipen = 999)
p <- check_X %>% 
    ggplot(aes(pos, num_ind_het, fill = par_snp, alpha = par_snp)) + 
    geom_point(size = 3, shape = 21) +
    theme_classic() +
    ylab("number of male individuals with het. genotype") +
    xlab("position on X") +
    #facet_wrap(~snp_chip, scales = "free_y", nrow = 3) +
    scale_alpha_manual(values = c(0.5, 1)) +
    scale_fill_manual(values = wes_palette("Darjeeling2")[1:2]) + 
    ggtitle("SNPs from hh files \n(at least one heterozygote genotype which shouldnt be there)")
p
ggsave("figs/X_snps_hh_hd_chip.jpg",plot = p, width = 6, height = 5)


p_ld <- check_X %>% 
    ggplot(aes(pos, num_ind_het, color = par_snp)) + 
    geom_point(size = 2, alpha = 0.5) +
    theme_classic() +
    ylab("number of LD individuals with het. genotype") +
    xlab("position on X") 


hets_per_snp_hd <- hh_file %>% filter(X2 %in% hd_inds) %>% group_by(X3) %>% tally()
hist(hets_per_snp_hd$n)

hh_snps_in_par <- hh_file$X3 %in% par_snps
sum(!hh_snps_in_par)

snps_on_x <- full_sample$map %>% filter(chromosome == 27) 













# CHECK LD chip
# plink name
#sheep_plink_name <- "merged_sheep_geno_ram"
sheep_plink_name <- "Plates_1to87_QC5_ram"
# read merged plink data
sheep_bed <- paste0(plink_geno_path, sheep_plink_name, ".bed")
sheep_bim <- paste0(plink_geno_path, sheep_plink_name, ".bim")
sheep_fam <- paste0(plink_geno_path, sheep_plink_name, ".fam")
full_sample <- read.plink(sheep_bed, sheep_bim, sheep_fam)
full_sample$fam

snps_on_x <- full_sample$map %>% filter(chromosome == 27) 

inds <- full_sample$fam %>% rename(id = member) 

# pseudoautosomal region SNPs
#par_snps <- readLines("/exports/csce/eddie/biology/groups/pemberton/martin/plink_genotypes/Oar3.1_PAR_SNPs_HD.txt")
par_snps <- read_delim("../sheep_imputation/data/Oar3.1_PAR_SNPs_HD.txt", " ")[[1]] %>% 
    as_tibble() %>% 
    mutate(par_snp = as.factor(1)) %>% 
    rename(snp.name = value)

# check hh file plink
hh_file <- read_delim("../sheep/data/SNP_chip/ramb_mapping/Plates_1-2_HD_QC3_ram.hh", "\t",
                      col_names = FALSE) %>% 
    rename(id = X2, snp.name = X3)


check_X <- hh_file %>% inner_join(inds) %>% 
    # filter(sex == 1) %>% 
    right_join(snps_on_x) %>% 
    group_by(snp.name) %>% 
    summarise(num_ind_het = n()-1, pos = first(position)) %>% 
    left_join(par_snps) %>% 
    mutate(par_snp = ifelse(is.na(par_snp), 0, 1)) %>% 
    mutate(par_snp = as.factor(par_snp)) 

library(wesanderson)
options(scipen = 999)
p <- check_X %>% 
    ggplot(aes(pos, num_ind_het, fill = par_snp, alpha = par_snp)) + 
    geom_point(size = 3, shape = 21) +
    theme_classic() +
    ylab("number of male individuals with het. genotype") +
    xlab("position on X") +
    #facet_wrap(~snp_chip, scales = "free_y", nrow = 3) +
    scale_alpha_manual(values = c(0.5, 1)) +
    scale_fill_manual(values = wes_palette("Darjeeling2")[1:2]) + 
    ggtitle("SNPs from hh files \n(at least one heterozygote genotype which shouldnt be there)")
p
ggsave("figs/X_snps_hh_ld_chip.jpg",plot = p, width = 6, height = 5)

p_ld <- check_X %>% 
    ggplot(aes(pos, num_ind_het, color = par_snp)) + 
    geom_point(size = 2, alpha = 0.5) +
    theme_classic() +
    ylab("number of LD individuals with het. genotype") +
    xlab("position on X") 


hets_per_snp_hd <- hh_file %>% filter(X2 %in% hd_inds) %>% group_by(X3) %>% tally()
hist(hets_per_snp_hd$n)

hh_snps_in_par <- hh_file$X3 %in% par_snps
sum(!hh_snps_in_par)

snps_on_x <- full_sample$map %>% filter(chromosome == 27) 












# plink name
sheep_plink_name <- "merged_sheep_geno_ram"
#sheep_plink_name <- "merged_sheep_geno"
# read merged plink data
sheep_bed <- paste0(plink_geno_path, sheep_plink_name, ".bed")
sheep_bim <- paste0(plink_geno_path, sheep_plink_name, ".bim")
sheep_fam <- paste0(plink_geno_path, sheep_plink_name, ".fam")
full_sample <- read.plink(sheep_bed, sheep_bim, sheep_fam)
full_sample$fam

snps_on_x <- full_sample$map %>% filter(chromosome == 27) 

geno_x <- full_sample$genotypes
geno_sub <- as_tibble(as(full_sample$genotypes[, snps_on_x$snp.name], Class = "numeric"),
                      rownames = "id") %>% 
  filter(id %in% hh_file$id)

geno_sub_hd <- geno_sub %>% filter(id %in% hd_inds$id) %>% 
  select(-id) %>% 
  as.matrix() 

geno_sub_hd[!(geno_sub_hd == 1)] <- NA
all_het <- colSums(geno_sub_hd, na.rm = TRUE)
plot(1:length(all_het), all_het)

ld_inds <- full_sample$fam[-c(1:188), ] %>% rename(id = member)
hd_inds <- full_sample$fam[c(1:188), ] %>% rename(id = member)

inds <- full_sample$fam %>% rename(id = member) 
inds$snp_chip <- c(rep("hd", 188), rep("ld", nrow(inds)-188))
# pseudoautosomal region SNPs
#par_snps <- readLines("/exports/csce/eddie/biology/groups/pemberton/martin/plink_genotypes/Oar3.1_PAR_SNPs_HD.txt")
par_snps <- read_delim("../sheep_imputation/data/Oar3.1_PAR_SNPs_HD.txt", " ")[[1]] %>% 
  as_tibble() %>% 
  mutate(par_snp = as.factor(1)) %>% 
  rename(snp.name = value)

# check hh file plink
hh_file <- read_delim("../sheep/data/SNP_chip/ramb_mapping/merged_sheep_geno_ram.hh", "\t",
                      col_names = FALSE) %>% 
  rename(id = X2, snp.name = X3)
# hh_file <- read_delim("../sheep/data/SNP_chip/old_merged_and_imputed/merged_sheep_geno.hh", "\t",
#                       col_names = FALSE) %>% 
#                       rename(id = X2, snp.name = X3)

snps_on_x %>% 
  inner_join(hd_inds)


check_X <- hh_file %>% left_join(inds) %>% 
  filter(sex == 1) %>% # only males in hh file anyway
  left_join(snps_on_x) %>% 
  group_by(snp_chip, snp.name) %>% 
  summarise(num_ind_het = n(), pos = first(position)) %>% 
  left_join(par_snps) %>% 
  mutate(par_snp = ifelse(is.na(par_snp), 0, 1)) %>% 
  mutate(par_snp = as.factor(par_snp))# %>% 
#filter(!is.na(snp_chip))

check_X %>% filter(snp_chip == "hd")
check_X %>% filter(snp_chip == "ld")
library(wesanderson)
options(scipen = 999)
p <- check_X %>% 
  ggplot(aes(pos, num_ind_het, fill = par_snp, alpha = par_snp)) + 
  geom_point(size = 3, shape = 21) +
  theme_classic() +
  ylab("number of male individuals with het. genotype") +
  xlab("position on X") +
  facet_wrap(~snp_chip, scales = "free_y", nrow = 3) +
  #scale_alpha_manual(values = c(0.5, 1)) +
  scale_fill_manual(values = wes_palette("Darjeeling2")[1:2]) + 
  ggtitle("SNPs from hh files \n(at least one heterozygote genotype which shouldnt be there)")
p
ggsave("figs/X_snps_hh_old_mapping.jpg",plot = p, width = 6, height = 5)

p_ld <- check_X %>% 
  ggplot(aes(pos, num_ind_het, color = par_snp)) + 
  geom_point(size = 2, alpha = 0.5) +
  theme_classic() +
  ylab("number of LD individuals with het. genotype") +
  xlab("position on X") 


hets_per_snp_hd <- hh_file %>% filter(X2 %in% hd_inds) %>% group_by(X3) %>% tally()
hist(hets_per_snp_hd$n)

hh_snps_in_par <- hh_file$X3 %in% par_snps
sum(!hh_snps_in_par)

snps_on_x <- full_sample$map %>% filter(chromosome == 27) 


