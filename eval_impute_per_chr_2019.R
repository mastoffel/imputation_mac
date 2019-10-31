# evaluate alpha impute
library(data.table)
library(purrr)
library(tidyverse)
library(snpStats)
source("../bottleneck/R/martin.R")

# folder with imputed genotypes (only masked individuals)
# this was extracted from the ImputedGenotypeProbabilities file
run_name <- "cv_3_hap_2019"

#run_eval_pipeline <- function(run_name){
    

imp_res_dir <- paste0("results/", run_name)

# 2069 6106 6593 1881 1258 2322 1097 4741 1101 2167
# HMM 10-30 run: c(2069,6106,6593,1881,1258)
# 3 ind runs: c(2322, 1258,1097)
#inds_imputed <- c(2322, 1258,1097, 6034, 3002, 17, 3076, 2709, 1880, 7359) # 2322, 1258,1097, 
inds_imputed <- c(2069, 6106, 6593, 1881, 1258, 2322, 1097, 4741, 1101, 2167)
chrs <- c(6,14,20,26)

#### imputed genotypes ####
read_imputed <- function(chr) {
    geno_imp <- fread(paste0(imp_res_dir, "/genos_chr_", chr, ".txt"))[V1 %in% inds_imputed]
    if ((chr != chrs[1])) geno_imp <- geno_imp[, -1]
    setDF(geno_imp)
}

# how many SNPs per chromosomes?
geno_imp_list <- map(chrs, read_imputed)
chr_lengths <- map_chr(geno_imp_list, function(x) ncol(x))
chr_lengths <- as.numeric(chr_lengths)
chr_lengths[1] <- chr_lengths[1] - 1 # minus ID column

# bind all imputed chromosomes together
geno_imp <- setDT(bind_cols(map(chrs, read_imputed)))

######## # extract true genotypes ########

# to_be_imputed <- read_lines("data/to_be_imputed.txt")
#### original genotypes ####
read_org <- function(chr) {
    geno_imp <- fread(paste0("data/genos_org_2019/Genotypes_chr_", chr, ".txt"))
    if ((chr != chrs[1])) geno_imp <- geno_imp[, -1]
    setDF(geno_imp)
}

# get original data for imputed chromosomes
geno_org <- bind_cols(map(chrs, read_org))

# index of snps to be imputed
na_dist <- colSums(as.matrix(geno_org) == 9)
to_be_imputed <- which(na_dist > 1000)

# filter only imputed individuals out
geno_org <- dplyr::filter(geno_org, V1 %in% inds_imputed)
# match individual sequence
geno_org <- setDT(geno_org[match(geno_imp[[1]], geno_org[[1]]), ])
# check whether sequence oif individuals is the same
geno_org[[1]] == geno_imp[[1]]
# get 

# transform imputed genotypes to long format for faster processing
geno_imp_t <- data.table::transpose(geno_imp)
colnames(geno_imp_t) <- paste0("imputed_", as.character(geno_imp_t[1, ]))
geno_imp_t <- geno_imp_t[-1, ]
setDT(geno_imp_t)

#geno_imp_t[, snp := call_rate$snp]

# transform true genotypes to long format for faster processing
geno_org_t <- data.table::transpose(geno_org)
colnames(geno_org_t) <- paste0("org_", as.character(geno_org_t[1, ]))
geno_org_t <- geno_org_t[-1, ]

# bind them together
full_df <- cbind(geno_org_t, geno_imp_t)

# transform dosages to genotypes using strict criteria +- 0.02 to call genotype, else set to missing
library(stringr)
# if (str_detect(run_name, "HMM")) {
#     df_trans <- full_df
# } else {
#     err <- 0
#     df_trans <- full_df  %>%
#         mutate_at(.vars = vars(contains("_")),
#                   list(~case_when(
#                       . <= (0 + err) ~ 0,
#                       (. >= (1-err)) & (. <= (1+err)) ~ 1,
#                       (. >= (2-err)) & (. <= (2+err)) ~ 2,
#                       TRUE ~ 9
#                   )))
# }

df_trans <- full_df
# translate 9 to NA
setDT(df_trans)
for(j in seq_along(df_trans)){
    set(df_trans, i=which(df_trans[[j]]==9), j=j, value=NA)
}
df_trans

# sort columns (probably not necessary)
#setDF(setcolorder(df_trans, sort(colnames(df_trans))))
#df_trans

# add chromosome information
chr <- unlist(map(1:length(chrs), function(x) rep(chrs[x], each = chr_lengths[x])))
df_trans$chr <- chr
df_trans <- df_trans[to_be_imputed, ] #
setDF(df_trans)

# calculate accuracies and proportion of imputed SNPs
calc_imp_results_per_ind <- function(ind_id, df_full) {
    ind_org <- paste0("org_", ind_id)
    ind_imp <- paste0("imputed_", ind_id)
    setDT(df_full)
    both_not_na <- (!is.na(df_full[, get(ind_org)])) & (!is.na(df_full[, get(ind_imp)]))
    
    prop_imputed <- df_full[, .(prop_imputed = sum(!is.na(get(ind_imp)))/(length(get(ind_imp)))), by = chr]
    prop_imputed[, prop_imputed := ifelse(prop_imputed==0, NA, prop_imputed)]
    
    df_full <- df_full[both_not_na, ]
    #prop_imputed = sum(!is.na(get(ind_imp)))/(length(get(ind_imp))),
    df_out <- df_full[, .(prop_correct = sum(get(ind_org) == get(ind_imp)) / length(get(ind_org)),
                          snps_wrong = sum(get(ind_org) != get(ind_imp))), by = chr]
    
    if (any(!(chrs %in% df_out$chr))) {
        df_out <- rbind(df_out, data.frame(chr = which(!(chrs %in% df_out$chr)), prop_correct = NA, snps_wrong = NA))
    }
    #df_out$prop_imputed <- prop_imputed
    df_final <- merge(df_out, prop_imputed, by = "chr")
    df_final[, ind_id := ind_id]
    df_final
}

# map accuracy function over individuals
all_inds <- geno_imp[[1]]
imp_acc <- map(all_inds, calc_imp_results_per_ind, df_trans)

# create large df
full_acc <- rbindlist(imp_acc) %>% setDF %>% mutate(run = run_name)
write_delim(full_acc, paste0("results/summaries2/", run_name, ".txt"))


library(viridis)
full_acc %>% 
    mutate(chr = as.factor(chr),
           ind_id = as.factor(ind_id)) %>% 
    pivot_longer(starts_with("prop"), names_to = "evaluation", values_to = "proportion") %>% 
ggplot( aes(chr, proportion, color = ind_id)) + 
    geom_point(size = 2.5, alpha = 0.8) +
    facet_wrap(.~evaluation, scales = "free_y") + 
    scale_color_brewer(type = "div") + 
    theme_classic() + 
    ylab("Proportion") +
    xlab("Chromosome") -> p_acc


ggsave("results/figs/acc_ramb_map_2019.jpg", p_acc, width = 7, height = 3)
# 
#}

old_acc <- read_delim("results/summaries2/cv_3_hap_2018.txt", ' ')
new_acc <- full_acc %>% filter(chr !=20 )

plot(old_acc$prop_correct, new_acc$prop_correct)
sum(round(new_acc$prop_correct,3) - old_acc$prop_correct)

run_names <- c("cv_1_5_hap_3chr", 
               "cv_2_10_hap_3chr", 
               "cv_10_30_hap_3chr", 
               "HMM_1_5_hap_3chr", 
               "HMM_2_10_hap_3chr",
               "HMM_10_30_hap_3chr")

# write all accuracies to files in results folder
map(run_names, run_eval_pipeline)

