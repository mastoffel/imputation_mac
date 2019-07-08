library(tidyverse)
library(ggplot2)
library(ggrepel)
library(forcats)
run_names <- c("cv_1_5_hap_3chr",
               "cv_2_10_hap_3chr",
               "cv_10_30_hap_3chr",
               "HMM_1_5_hap_3chr",
               "HMM_2_10_hap_3chr",
               "HMM_10_30_hap_3chr")

#run_names <- c("cv_full_1_5")

read_results <- function(run_names){
    run_results <- read_delim(paste0("results/summaries2/",run_names,".txt"), delim = " ")
}
all_accs <- map_df(run_names, read_results) %>% 
                mutate(chr = as.factor(chr),
                       run = fct_inorder(run),
                       ind_id = as.factor(ind_id))
# all_accs$prop_correct[all_accs$prop_correct < 0.925] <- 0.93 # actually 0.88

run_list <- c(
    "cv_1_5_hap_3chr" = "heuristic_short_cores",
    "cv_2_10_hap_3chr" = "heuristic_medium_cores",
    "cv_10_30_hap_3chr" = "heuristic_long_cores",
    "HMM_1_5_hap_3chr" = "HMM_heuristic_short_cores",
    "HMM_2_10_hap_3chr" = "HMM_heuristic_medium_cores",
    "HMM_10_30_hap_3chr" = "HMM_heuristic_long_cores")

p1 <- ggplot(all_accs, aes(chr, prop_correct)) +
    geom_boxplot(color = "lightgrey", width = 0.5) +
    geom_point(size = 3, alpha = 0.7, aes(col = ind_id)) +
    #theme_martin() +
    ylab("Correctly imputed SNPs") +
    xlab("chromosome") +
    theme_classic() +
    theme(panel.grid.major = element_line()) +
    scale_color_viridis_d() +
    facet_wrap(~run, labeller= as_labeller(run_list)) +
    ggtitle("AlphaImpute comparison with only heuristic method and heuristic+HMM method
(1) First panel of 6 plots shows the proportion of correctly imputed SNPs 
out of those that have been imputed.
(2) Second panel of 6 plots shows the proportion of SNPs which have been imputed.
(3) Only genotypes where 100% of cores had the same allele were included.
(4) The following core lengths were used:
short_cores = 1 to 5% of chromosome, medium = 2 to 10%, long = 10 to 30%")
p1
p2 <- ggplot(all_accs, aes(chr, prop_imputed)) +
    geom_boxplot(color = "lightgrey", width = 0.5) +
    geom_point(size = 3, alpha = 0.7, aes(col = ind_id)) +
    #theme_martin() +
    ylab("Proportion of imputed SNPs") +
    xlab("chromosome") +
    theme_classic() +
    theme(panel.grid.major = element_line()) +
    scale_color_viridis_d() +
    facet_wrap(~run, labeller= as_labeller(run_list))
p2
library(patchwork)

p_out <- p1/p2
p_out
ggsave( "results/summaries/alphaimpute_comparison_ImputeGenotypesHMM.jpg", p_out, width = 10, height = 9)


all_accs %>% 
    group_by(ind_id, run) %>% 
    summarise(prop_correct = mean(prop_correct),
              prop_imputed = mean(prop_imputed)) %>% 
    mutate(run = run_list) -> accs_out

WriteXLS::WriteXLS(accs_out, "results/summaries/alphaimpute_6methods_10individuals_3chromosomes_0perc_err.xls")






ggplot(all_accs, aes(prop_correct, prop_imputed, col=run)) +
    geom_point(size = 5, alpha = 0.5)

    #scale_colour_brewer(type = "qual") +
    #theme(legend.position = "none") +
    ylab("proportion of SNPs correct") +
    xlab("chromosome") +
    ggtitle("Proportion of correctly imputed SNPs per individual and 
            chromosome\n-based on SNPs which have been imputed (with high probabilities)")
p1

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
ggsave(plot = p_per_chr, filename = "imputation_accuracy_1_5hap.jpg", width = 8, height = 7)

full_acc %>% group_by(ind_id) %>% 
    summarise(prop_correct_mean_across_chr = mean(prop_correct, na.rm = TRUE),
              prop_correct_lowest = range(prop_correct, na.rm = TRUE)[1],
              prop_correct_highest = range(prop_correct, na.rm = TRUE)[2],
              prop_imputed_mean_across_chr = mean(prop_imputed, na.rm = TRUE),
              prop_imputed_lowest = range(prop_imputed, na.rm = TRUE)[1],
              prop_imputed_highest = range(prop_imputed, na.rm = TRUE)[2]) -> acc_per_ind

WriteXLS::WriteXLS(acc_per_ind, "imputation_accuracy_1_5hap.xls")






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

