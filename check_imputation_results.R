# check imputation results
library(tidyverse)
accuracies <- read_delim("run_1/Miscellaneous/AccuracyPerAnimal.txt", delim = " ", 
                        col_names = FALSE)

# imputation accuracies
accs <- accuracies %>% 
    mutate(X1 = str_trim(X1, side = "left"),
           X2 = as.numeric(str_trim(X2, side = "left"))) %>% 
    select(-X3)
accs

# pedigree
pedigree <- read_delim("../sheep/data/SNP_chip/20190208_Full_Pedigree.txt", delim = "\t") %>% 
                mutate(ID = as.character(ID))
str(pedigree)

check_df <- accs %>% 
    filter(X2 != 1) %>% 
    rename(ID = X1, accuracies = X2) %>% 
    left_join(pedigree, by = "ID")
