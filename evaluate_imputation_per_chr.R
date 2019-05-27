# evaluate alpha impute
library(data.table)
library(purrr)
library(tidyverse)
source("../bottleneck/R/martin.R")

# folder with imputed genotypes (only masked individuals)
# this was extracted from the ImputedGenotypeProbabilities file
imp_res_dir <- "results/imputation_by_chr"

#### imputed genotypes ####
read_imputed <- function(chr) {
    geno_imp <- fread(paste0(imp_res_dir, "/genos_chr_", chr, ".txt"))
    if ((chr != 1) & (chr != 2)) geno_imp <- geno_imp[, -1]
    setDF(geno_imp)
}

# with column binding / here chr_1 hasnt finished for 1 individual 
geno_imp_chr1 <- setDT(read_imputed(1))
geno_imp <- setDT(bind_cols(map(2:26, read_imputed)))
geno_imp <- merge(geno_imp_chr1, geno_imp, all = TRUE, by = "V1")
########

#### original genotypes ####
org_dir <- "results/geno_full_by_chr/"
read_org <- function(chr) {
    geno_imp <- fread(paste0(org_dir , "/genos_chr_", chr, ".txt"))
    if (chr != 1) geno_imp <- geno_imp[, -1]
    setDF(geno_imp)
}

geno_org_list <- map(1:26, read_org)

# how long are the chromosomes
chr_lengths <- map_chr(geno_org_list, function(x) ncol(x))
chr_lengths <- as.numeric(chr_lengths)
chr_lengths[1] <- chr_lengths[1] - 1

geno_org <- bind_cols(map(1:26, read_org))

# filter 
inds <- geno_org[[1]]
geno_imp <- geno_imp[V1 %chin% inds]
########


# transform imputed genotypes to long format for faster processing
geno_imp_t <- data.table::transpose(geno_imp)
colnames(geno_imp_t) <- paste0("imputed_", as.character(geno_imp_t[1, ]))
geno_imp_t <- geno_imp_t[-1, ]
setDT(geno_imp_t)

# same for true genotypes, also subset to the individuals which were imputed
#geno_org <- geno[ID %chin% geno_imp[[1]]]
geno_org_t <- data.table::transpose(geno_org)
colnames(geno_org_t) <- paste0("org_", as.character(geno_org_t[1, ]))
geno_org_t <- geno_org_t[-1, ]

# bring them together
full_df <- cbind(geno_org_t, geno_imp_t)


# transform dosages to genotypes using strict criteria +- 0.05 to call genotype, else set to missing
# df_trans <- full_df %>%
#     mutate_at(.vars = vars(contains("imputed")),
#               funs(case_when(
#                   . <= 0.05 ~ 0,
#                   (. >= 0.95) & (. <= 1.05)  ~ 1,
#                   . >= 1.95 ~ 2,
#                   TRUE ~ 9
#               ))) %>%
#     setDT()

df_trans <- full_df  %>% 
    mutate_at(.vars = vars(contains("imputed")),
              funs(case_when(
                  . <= 0.02 ~ 0,
                  (. >= 0.98) & (. <= 1.02)  ~ 1,
                  (. >= 1.98) & (. <= 2.02) ~ 2,
                  TRUE ~ 9
              ))) 

# translate 9 to NA
setDT(df_trans)
for(j in seq_along(df_trans)){
    set(df_trans, i=which(df_trans[[j]]==9), j=j, value=NA)
}
df_trans

# sort columns (probably not necessary)
setDF(setcolorder(df_trans, sort(colnames(df_trans))))
df_trans
# add chromosome information
chr <- unlist(map(1:26, function(x) rep(x, chr_lengths[x])))
df_trans$chr <- chr

# split up again into imputed and true (original) genotypes
#df_imputed <- df_trans_imp[1:(ncol(df_trans_imp)/2)]
#df_org <- df_trans_imp[((ncol(df_trans_imp)/2)+1):ncol(df_trans_imp)]



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
    if (any(!(1:26 %in% df_out$chr))) {
        df_out <- rbind(df_out, data.frame(chr = which(!(1:26 %in% df_out$chr)), prop_correct = NA, snps_wrong = NA))
    }
    #df_out$prop_imputed <- prop_imputed
    df_final <- merge(df_out, prop_imputed, by = "chr")
    df_final[, ind_id := ind_id]
    df_final
}

all_inds <- geno_imp[[1]]
imp_acc <- map(all_inds, calc_imp_results_per_ind, df_trans)


# create large df
full_acc <- rbindlist(imp_acc) %>% setDF


library(ggplot2)
library(ggrepel)

p1 <- ggplot(full_acc, aes(as.factor(chr), prop_correct, col = as.factor(chr))) +
    geom_boxplot() +
    geom_point(size = 3, alpha = 0.5) +
    theme_martin() +
    scale_color_viridis_d() +
    #scale_colour_brewer(type = "qual") +
    theme(legend.position = "none") +
    ylab("proportion of SNPs correct") +
    xlab("chromosome") +
    ggtitle("Proportion of correctly imputed SNPs per individual and 
            chromosome\n-based on SNPs which have been imputed (with high probabilities)")
    # geom_text_repel(aes(label = ind_id),
    #                 size = 2)
    
p2 <- ggplot(full_acc, aes(as.factor(chr), prop_imputed, col = as.factor(chr))) +
    geom_boxplot() +
    geom_point(size = 3, alpha = 0.5) +
    theme_martin() +
    scale_color_viridis_d() +
    theme(legend.position = "none") +
    ylab("proportion of SNPs imputed") +
    xlab("chromosome") +
    ggtitle("Proportion of imputed SNPs per individual and chromosome")

library(patchwork)

p_per_chr <- p1 / p2
ggsave(plot = p_per_chr, filename = "imputation_accuracy.jpg", width = 8, height = 7)

 full_acc %>% group_by(ind_id) %>% 
    summarise(prop_correct_mean_across_chr = mean(prop_correct, na.rm = TRUE),
              prop_correct_lowest = range(prop_correct, na.rm = TRUE)[1],
              prop_correct_highest = range(prop_correct, na.rm = TRUE)[2],
              prop_imputed_mean_across_chr = mean(prop_imputed, na.rm = TRUE),
              prop_imputed_lowest = range(prop_imputed, na.rm = TRUE)[1],
              prop_imputed_highest = range(prop_imputed, na.rm = TRUE)[2]) -> acc_per_ind

WriteXLS::WriteXLS(acc_per_ind, "accuracy_per_ind.xls")



# make df for plotting
snp_names <- paste0("SNP_", 1:length(prop_wrong_per_SNP))
wrong_per_SNP <- tibble(prop_wrong = prop_wrong_per_SNP, SNP = snp_names, num_incorrect = incor_imputed_SNPs,
                        snp_group = cut_number(1:length(prop_wrong_per_SNP), 100))

p1 <- wrong_per_SNP %>% 
    group_by(snp_group) %>% 
    summarise(mean_wrong = mean(prop_wrong, na.rm = TRUE)) %>% 
    #summarise(abs_incor = sum(num_incorrect, na.rm = TRUE)) %>% 
    ggplot(aes(snp_group, mean_wrong)) +
    geom_col() +
    theme_martin() +
    xlab("SNP bin (each containing around 450 SNPs") +
    ylab("Mean proportion of incorrectly \nimputed individuals per SNP") +
    theme(axis.text.x = element_blank(),
          panel.grid = element_blank()) +
    ggtitle("AlphaImpute imputation inaccuracy across chromosome 1")

ggsave(filename = "accuracy_across_chromosome1.jpg", width = 7, height = 4)



# only first bin

