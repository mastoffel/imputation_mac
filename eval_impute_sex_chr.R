# evaluate alpha impute
library(data.table)
library(purrr)
library(tidyverse)
library(snpStats)
source("../bottleneck/R/martin.R")

# folder with imputed genotypes (only masked individuals)
# this was extracted from the ImputedGenotypeProbabilities file
run_name <- "cv_full_1_5_sex_chr"
imp_res_dir <- paste0("results/", run_name, "/genos/") # genos
chrs <- 27

# imputed genotypes
geno_imp <- fread(paste0(imp_res_dir, "/genos_chr_", chrs, ".txt"), header = FALSE)
inds_imputed <- as.character(geno_imp[[1]])
# how many SNPs?
chr_lengths <- ncol(geno_imp) - 1

#### imputed genotypes ####

######## # extract true genotypes ########
to_be_imputed <- read_lines("data/to_be_imputed_chr27.txt")

# true genotypes 
# plink name
sheep_plink_name <- "../sheep/data/SNP_chip/ramb_mapping/merged_sheep_geno_ram"
# read merged plink data
sheep_bed <- paste0(sheep_plink_name, ".bed")
sheep_bim <- paste0(sheep_plink_name, ".bim")
sheep_fam <- paste0(sheep_plink_name, ".fam")
full_sample <- read.plink(sheep_bed, sheep_bim, sheep_fam)

par_snps <- read_lines("data/Ram_PAR_SNPs_HDLD_all.txt")

# these snps were imputed
#imp_snps_names <- afread("data/imputed_snps_chr27.txt", nrows = 0)
# filter names of snps on one chromosome
all_chr_snps <- full_sample$map %>% filter(chromosome == 27) %>% 
                filter(!(snp.name %in% par_snps)) %>% 
                .$snp.name

geno_org <- as(full_sample$genotypes[, all_chr_snps], Class = "numeric")

# subset geno_org for individuals present in geno_imp
geno_org <- geno_org[rownames(geno_org) %chin% inds_imputed, ]
geno_org <- setDT(as.data.frame(geno_org), keep.rownames = TRUE)

# transform imputed genotypes to long format for faster processing
geno_imp_t <- data.table::transpose(geno_imp)
colnames(geno_imp_t) <-
    paste0("imputed_", as.character(geno_imp_t[1,]))
geno_imp_t <- geno_imp_t[-1,]
setDT(geno_imp_t)
#geno_imp_t[, snp := names(geno_org)[-1]]


# transform true genotypes to long format for faster processing
geno_org_t <- data.table::transpose(geno_org)
colnames(geno_org_t) <-
    paste0("org_", as.character(geno_org_t[1,]))
geno_org_t <- geno_org_t[-1,]

# bind them together
full_df <- cbind(geno_org_t, geno_imp_t)

# if no dosages
df_trans <- full_df

# translate 9 to NA
setDT(df_trans)
for (j in seq_along(df_trans)) {
    set(
        df_trans,
        i = which(df_trans[[j]] == 9),
        j = j,
        value = NA
    )
}
df_trans


# add chromosome information
#chr <- unlist(map(chrs, function(x) rep(x, chr_lengths[x])))
df_trans$chr <- chrs
#df_trans <- df_trans[snp %chin% to_be_imputed]
setDF(df_trans)

# Oar3.1_PAR_SNPs_HD
head(df_trans)
#par_snps <- read_delim("data/Oar3.1_PAR_SNPs_HD.txt", delim = " ", col_names = FALSE)


# calculate accuracies and proportion of imputed SNPs
calc_imp_results_per_ind <- function(ind_id, df_full) {
    # assign names
    ind_org <- paste0("org_", ind_id)
    ind_imp <- paste0("imputed_", ind_id)
    setDT(df_full)
    # positions where both are not NA
    both_not_na <-
        (!is.na(df_full[, get(ind_org)])) &
        (!is.na(df_full[, get(ind_imp)]))
    
    # proportion of SNPs imputed per chromosome
    prop_imputed <- 
        df_full[, .(prop_imputed = sum(!is.na(get(ind_imp))) / (length(get(ind_imp)))), 
                by = chr]
    prop_imputed[, prop_imputed := ifelse(prop_imputed == 0, NA, prop_imputed)]
    
    # get df where imputed and original are not NA
    df_full <- df_full[both_not_na,]
    
    # calculate the proportion of correct SNPs and how many SNPs were wrongly imputed
    df_out <- 
        df_full[, .(
            prop_correct = sum(get(ind_org) == get(ind_imp)) / length(get(ind_org)),
            snps_wrong = sum(get(ind_org) != get(ind_imp))
        ), by = chr]
    
    # In case any chromosome is missing from the imputation runs
    if (any(!(chrs %in% df_out$chr))) {
        df_out <-
            rbind(df_out,
                  data.frame(
                      chr = which(!(chrs %in% df_out$chr)),
                      prop_correct = NA,
                      snps_wrong = NA
                  ))
    }
    
    #df_out$prop_imputed <- prop_imputed
    df_final <- merge(df_out, prop_imputed, by = "chr")
    df_final[, ind_id := ind_id]
    df_final
}

# map accuracy function over individuals
all_inds <- geno_imp[[1]]

imp_acc <- map(all_inds, calc_imp_results_per_ind, df_trans)
# create table with all accuracies
full_acc <-rbindlist(imp_acc) %>% setDF %>% mutate(run = run_name) 

full_acc %>% write_delim("results/summaries/sex_chr_noPAR_ram.txt")

# 









# sexes
ind_sex <- read_delim("data/sex_chr_impute_table.txt", delim = " ", col_names = FALSE)
names(ind_sex) <- c("ind_id", "sex")
head(df_trans)


seq(from = 1, to = 16143, by = 500)

acc_acr_chr <- function(start) {
    imp_acc <- map(all_inds, calc_imp_results_per_ind, df_trans[start:(start+500), ])
    # create table with all accuracies
    full_acc <-rbindlist(imp_acc) %>% setDF %>% mutate(run = run_name) %>% 
                    mutate(start = {{ start }})
}

all_accs <- map_df(seq(from = 1, to = 16143, by = 100), acc_acr_chr)

all_accs %>% 
    left_join(ind_sex) %>% 
    ggplot(aes(start, prop_correct)) +
    geom_point()
    geom_boxplot() +
    geom_jitter(size = 3, alpha = 0.4)

par_snps <- read_delim("data/Oar3.1_PAR_SNPs_HD.txt", delim = " ", col_names = FALSE)
df_trans_nopar <- df_trans[!(snp %in% par_snps$X1)]

imp_acc <- map(all_inds, calc_imp_results_per_ind, df_trans_nopar)
    # create table with all accuracies
full_acc <-rbindlist(imp_acc) %>% setDF %>% mutate(run = run_name)

full_acc %>% 
    left_join(ind_sex) %>% 
    group_by(sex) %>% 
    summarise(mean(prop_correct), mean(prop_imputed))





