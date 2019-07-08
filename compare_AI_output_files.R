# compare full cv output when taking imputeGenotypes vs imputeGenotypeProbabilities 
# output files from AI

library(tidyverse)
library(data.table)
imp_genos <- read_delim(paste0("results/summaries/cv_full_1_5genos.txt"), delim = " ")
imp_dosages <- read_delim(paste0("results/summaries/cv_full_1_5.txt"), delim = " ")

names(imp_genos)

plot(imp_genos$prop_correct, imp_dosages$prop_correct)

imp_genos %>% 
    group_by(ind_id) %>% 
    summarise(prop_correct = mean(prop_correct),
              prop_imputed = mean(prop_imputed))
imp_dosages %>% 
    group_by(ind_id) %>% 
    summarise(prop_correct = mean(prop_correct),
              prop_imputed = mean(prop_imputed)) 
