# evaluate alpha impute
library(data.table)
library(purrr)
library(tidyverse)
library(snpStats)
source("../bottleneck/R/martin.R")

# folder with imputed genotypes (only masked individuals)
# this was extracted from the ImputedGenotypeProbabilities file
run_name <- "cv_full_1_5_sex_chr"
imp_res_dir <- paste0("results/", run_name, "/genos") # genos

inds_imputed <- read_lines("data/cv_inds.txt")
chrs <- 27

#### imputed genotypes ####
geno_imp <- fread(paste0(imp_res_dir, "/genos_chr_", chrs, ".txt"), header = FALSE)
# how many SNPs?
chr_lengths <- ncol(geno_imp) - 1

#par_snps <- read_lines("data/par_snps_in_data.txt")
######## # extract true genotypes ########
to_be_imputed <- read_lines("data/to_be_imputed_chr27.txt")

# grep true genotypes from merged geno txt file
inds_regex <- paste0(" ",paste(inds_imputed, collapse = " | "), " ")
# geno_org <- fread(cmd = paste0("grep -E ", "'", inds_regex, "'", " data/hdld_geno_merged.txt"), sep = " ", na.strings="9")
geno_org <- fread("data/hdld_geno_merged_sex_chr.txt", sep = " ", na.strings = "9")

# subset geno_org for individuals present in geno_imp
geno_org <- geno_org[ID %in% geno_imp$V1] %>% 
                select(ID, !!par_snps)

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
#     mutate_at(.vars = vars(contains("imputed")),
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
#chr <- unlist(map(chrs, function(x) rep(x, chr_lengths[x])))
df_trans$chr <- chrs
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





