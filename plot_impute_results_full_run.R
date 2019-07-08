library(tidyverse)
library(ggplot2)
library(ggrepel)
library(forcats)
library(patchwork)

run_names <- c("cv_full_1_5")

read_results <- function(run_names){
    run_results <- read_delim(paste0("results/summaries/",run_names,"per_chr_and_ind.txt"), delim = " ")
}
all_accs <- map_df(run_names, read_results) %>% 
    mutate(chr = as.factor(chr),
           run = fct_inorder(run),
           ind_id = as.factor(ind_id))


p1 <- ggplot(all_accs, aes(chr, prop_correct)) +
    geom_boxplot(color = "lightgrey", width = 0.5) +
    geom_point(size = 3, alpha = 0.7, aes(col = ind_id)) +
    #theme_martin() +
    ylab("Correctly imputed SNPs") +
    xlab("chromosome") +
    theme_classic() +
    theme(panel.grid.major = element_line()) +
    scale_color_viridis_d() 
p1
p2 <- ggplot(all_accs, aes(chr, prop_imputed)) +
    geom_boxplot(color = "lightgrey", width = 0.5) +
    geom_point(size = 3, alpha = 0.7, aes(col = ind_id)) +
    #theme_martin() +
    ylab("Proportion of imputed SNPs") +
    xlab("chromosome") +
    theme_classic() +
    theme(panel.grid.major = element_line()) +
    scale_color_viridis_d() 
p2

p_out <- p1/p2
ggsave( "results/summaries/full_run_genos.jpg", p_out, width = 10, height = 9)


all_accs %>% 
    group_by(ind_id) %>% 
    summarise(prop_correct = mean(prop_correct),
              prop_imputed = mean(prop_imputed)) -> accs_out

WriteXLS::WriteXLS(accs_out, "results/summaries/full_run_1_5_genos.xls")



