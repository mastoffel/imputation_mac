# evaluate alpha impute
library(data.table)
library(purrr)
library(tidyverse)
library(snpStats)
library(vroom)

# folder with imputed genotypes (only masked individuals)
# this was extracted from the ImputedGenotypeProbabilities file
run_name <- "cv_full_1_5_2020"
imp_res_dir <- paste0("results/", run_name, "/")
    
inds_imputed <- read_lines("data/cv_inds.txt")[1:10]
chrs <- 1:26
    
#### imputed genotypes ####
read_imputed <- function(chr) {
    geno_imp <- fread(paste0(imp_res_dir, "/genos_chr_", chr, ".txt"), header = FALSE) #[V1 %in% inds_imputed]
    geno_imp <- geno_imp[as.character(geno_imp$V1) %in% inds_imputed, ]
    if ((chr != chrs[1])) geno_imp <- geno_imp[,-1]
    setDF(geno_imp)
}

# how many SNPs per chromosomes?
geno_imp_list <- map(chrs, read_imputed)
chr_lengths <- map_chr(geno_imp_list, function(x)
    ncol(x))
chr_lengths <- as.numeric(chr_lengths)
chr_lengths[1] <- chr_lengths[1] - 1 # minus ID column

# bind all imputed chromosomes together
geno_imp <- map(chrs, read_imputed)
geno_imp <- setDT(bind_cols(geno_imp))

######## # extract true genotypes ########
#to_be_imputed <- read_lines("data/to_be_imputed.txt")

plink_geno <- "../sheep/data/SNP_chip/oar31_mapping/merged_sheep_geno_oar31"

# read merged plink data
sheep_bed <- paste0(plink_geno, ".bed")
sheep_bim <- paste0(plink_geno, ".bim")
sheep_fam <- paste0(plink_geno, ".fam")
full_sample <- read.plink(sheep_bed, sheep_bim, sheep_fam)
snps_autosomes <- full_sample$map %>% filter(chromosome %in% 1:26) %>% .$snp.name
geno_org <- as(full_sample$genotypes[rownames(full_sample$genotypes) %in% inds_imputed, snps_autosomes], "numeric") %>% 
            as.data.table(keep.rownames = "V1")
snp_map <- full_sample$map %>% 
                rename(snp = snp.name, chr = chromosome) %>% 
                dplyr::select(snp, chr) %>% 
                as_tibble() %>% 
                filter(chr %in% 1:26)
geno_org <- geno_org[as.character(geno_org$V1) %in% inds_imputed, ]

# get snps to be imputed
snp_sum <- col.summary(full_sample$genotypes) 
to_be_imputed <- snp_map$snp[snp_map$snp %in% rownames(snp_sum)[snp_sum$Calls < 200]]

# snps with low call rate
low_snps <- read_lines("../sheep_ID/data/oar_imp_low_call_snps095.txt")
snp_low_ind <- which(names(geno_org) %in% low_snps)

geno_org <- geno_org[, !snp_low_ind, with=FALSE]
geno_imp <- geno_imp[, !snp_low_ind, with=FALSE]
# order both dt by id
setorder(geno_org, "V1")
setorder(geno_imp, "V1")

# transform imputed genotypes to long format for faster processing
geno_imp_t <- data.table::transpose(geno_imp)
colnames(geno_imp_t) <-
    paste0("imputed_", as.character(geno_imp_t[1,]))
geno_imp_t <- geno_imp_t[-1,]
setDT(geno_imp_t)
geno_imp_t[, snp := names(geno_org)[-1]]

# transform true genotypes to long format for faster processing
geno_org_t <- data.table::transpose(geno_org)
colnames(geno_org_t) <-
    paste0("org_", as.character(geno_org_t[1,]))
geno_org_t <- geno_org_t[-1,]

# bind them together
full_df <- cbind(geno_org_t, geno_imp_t)

# transform dosages to genotypes using strict criteria +- 0.02 to call genotype, else set to missing
# err <- 0
# df_trans <- full_df  %>%
#     mutate_at(.vars = vars(contains("_")),
#               list( ~ case_when(
#                   . <= (0 + err) ~ 0,
#                   (. >= (1 - err)) & (. <= (1 + err)) ~ 1,
#                   (. >= (2 - err)) & (. <= (2 + err)) ~ 2,
#                   TRUE ~ 9
#               )))
# if no dosages
df_trans <- full_df

# translate 9 to NA
setDT(df_trans)
for (j in seq_along(df_trans)) {
    set(
        df_trans,
        i = which(df_trans[[j]] == 9),
        j = j,
        value = NA_real_
    )
}
df_trans

# add chromosome information
df_trans <- df_trans %>% left_join(snp_map)
df_trans <- df_trans[snp %chin% to_be_imputed]
setDF(df_trans)

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

# per chromosome and individual accuracy
write_delim(full_acc, paste0("results/summaries/", run_name, "per_chr_and_ind.txt"))

# per individual accuracy
acc_per_ind <- full_acc %>% group_by(ind_id) %>% 
                            summarise(prop_correct = mean(prop_correct), 
                                      genome_wide_snps_wrong = sum(snps_wrong),
                                      prop_imputed = mean(prop_imputed)) %>% 
                            as.data.frame()

write_delim(acc_per_ind, paste0("results/summaries/", run_name, "per_ind.txt"))

# accuracy per SNP
df_trans_mat <- as.matrix(df_trans)[,sort(names(df_trans))[-c(1, 32)]]
test <- df_trans_mat[, 1:15] == df_trans_mat[, 16:30]
snps_correct <- rowSums(test, na.rm = TRUE)
hist(snps_correct)

df_trans_sort <- df_trans

test <- full_acc %>% 
    group_by(chr) %>% 
    summarise(ave = mean(prop_imputed)) 

sum((test$ave * chr_lengths)) / sum(chr_lengths)  
sum((test$ave * chr_lengths)) / sum(chr_lengths)  
# plot
library(ggplot2)
library(viridis)
source("theme_simple.R")
library(gt)
median(acc_per_ind$prop_correct)

# first a table
full_acc %>% 
    group_by(chr) %>% 
    summarise(prop_correct = mean(prop_correct),
              prop_imputed = mean(prop_imputed)) %>% 
    add_column(n_snps = chr_lengths) %>% 
    select(chr,  prop_correct, prop_imputed, n_snps) %>% 
    mutate(prop_imputed = round(prop_imputed * 100, 3),
           prop_correct = round(prop_correct * 100, 3)) %>% 
    setNames(c("Chromosome",  "% SNPs correct", "% SNPs imputed", "N (Snps)")) %>% 
    mutate(`N (Individuals)` = 10) %>% 
    gt() %>% 
    gtsave("../sheep_ID/figs/tables/imputation.png")

full_acc %>% 
    mutate(chr = fct_inorder(as.factor(chr))) %>% 
    mutate(ind_id = as.factor(ind_id)) %>% 
    ggplot(aes(chr, prop_correct, fill = ind_id)) + 
    geom_point(shape = 21, size = 3) +
    scale_fill_viridis("individual", discrete = TRUE) +
    theme_minimal() +
    ylab("Correctly imputed SNPs") +
    xlab("Chromosome") -> p1

full_acc %>% 
    mutate(chr = fct_inorder(as.factor(chr))) %>% 
    mutate(ind_id = as.factor(ind_id)) %>% 
    mutate(prop_correct = prop_correct * 100) %>% 
    ggplot(aes(chr, prop_correct)) + 
    geom_boxplot() + 
    geom_point(size = 1, alpha = 0.5) +
    scale_fill_viridis("individual", discrete = TRUE) +
    theme_simple(grid_lines = FALSE, axis_lines = TRUE) +
    ylab("Correctly imputed SNPs (%)") +
    xlab("Chromosome")

ggsave("results/figs/acc_ramb_2020.jpg", width = 6, height = 3)

full_acc %>% 
    mutate(chr = fct_inorder(as.factor(chr))) %>% 
    mutate(ind_id = as.factor(ind_id)) %>% 
    ggplot(aes(chr, prop_correct, fill = ind_id)) + 
    geom_point(shape = 21, size = 3) +
    scale_fill_viridis("individual", discrete = TRUE) +
    theme_minimal() +
    ylab("Proportion of SNPs correctly imputed") +
    xlab("Chromosome") -> p1
