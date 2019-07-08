# evaluate alpha impute
library(data.table)
library(purrr)
library(tidyverse)
library(snpStats)
source("../bottleneck/R/martin.R")
library(vroom)

# folder with imputed genotypes (only masked individuals)
# this was extracted from the ImputedGenotypeProbabilities file
run_name <- "cv_full_1_5_sex_chr"
imp_res_dir <- paste0("results/", run_name, "/genos")
    
inds_imputed <- read_lines("data/inds_for_cv_full_first15.txt")
chrs <- 27
    
#### imputed genotypes ####
read_imputed <- function(chr) {
    geno_imp <- fread(paste0(imp_res_dir, "/genos_chr_", chr, ".txt"), header = FALSE) #[V1 %in% inds_imputed]
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
geno_imp <- setDT(bind_cols(map(chrs, read_imputed)))

######## # extract true genotypes ########
to_be_imputed <- read_lines("data/to_be_imputed.txt")

# grep true genotypes from merged geno txt file
inds_regex <- paste0(" ",paste(inds_imputed, collapse = " | "), " ")
# geno_org <- fread(cmd = paste0("grep -E ", "'", inds_regex, "'", " data/hdld_geno_merged.txt"), sep = " ", na.strings="9")
geno_org <- fread("data/hdld_geno_merged.txt", sep = " ", na.strings = "9")
geno_org <- geno_org[ID %in% inds_imputed]

# subset full genotype data for only the chromosomes of interest
all_chr_lengths <- as.numeric(read_lines("data/chr_lengths.txt"))
chr_pos_df <- data.frame(chr = rep(1:26, times = all_chr_lengths))
geno_org <- geno_org[, c(TRUE, as.character(chr_pos_df$chr) %chin% as.character(chrs)), with = FALSE]


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
        value = NA
    )
}
df_trans

# add chromosome information
chr <- unlist(map(chrs, function(x) rep(x, all_chr_lengths[x])))
df_trans$chr <- chr
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
write_delim(format(full_acc, digits = 4), paste0("results/summaries/", run_name, "per_chr_and_ind.txt"))

# per individual accuracy
acc_per_ind <- full_acc %>% group_by(ind_id) %>% 
                            summarise(prop_correct = mean(prop_correct), 
                                      genome_wide_snps_wrong = sum(snps_wrong),
                                      prop_imputed = mean(prop_imputed)) %>% 
                            as.data.frame()

write_delim(format(acc_per_ind, digits = 4), paste0("results/summaries/", run_name, "per_ind.txt"))

# accuracy per SNP
df_trans_mat <- as.matrix(df_trans)[,sort(names(df_trans))[-c(1, 32)]]
test <- df_trans_mat[, 1:15] == df_trans_mat[, 16:30]
snps_correct <- rowSums(test, na.rm = TRUE)
hist(snps_correct)

df_trans_sort <- df_trans



