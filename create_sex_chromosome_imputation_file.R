# create file with ID and sex for sex chromosome imputation

library(snpStats)
library(dplyr)
# on mac
plink_geno_path <- "../sheep/data/SNP_chip/"
# plink name
sheep_plink_name <- "merged_sheep_geno"
# read merged plink data
sheep_bed <- paste0(plink_geno_path, sheep_plink_name, ".bed")
sheep_bim <- paste0(plink_geno_path, sheep_plink_name, ".bim")
sheep_fam <- paste0(plink_geno_path, sheep_plink_name, ".fam")
full_sample <- read.plink(sheep_bed, sheep_bim, sheep_fam)

# 7374 sheep are imputed, these are their ids
sheep_ids <- data.frame(ID = readLines("data/sheep_ids_for_imputation.txt"))

full_sample$fam[duplicated(full_sample$fam$member), ]

# 1 is male, 2 is female
sex_file_out <- full_sample$fam %>% 
    group_by(member) %>% 
    summarise(mean_sex = mean(sex)) %>% 
    mutate(ID = as.character(member)) %>% 
    filter(ID %in% sheep_ids$ID) %>% 
    select(ID, mean_sex)
    
write_delim(sex_file_out, "data/sex_chr_impute_table.txt", delim = " ", 
            col_names = FALSE)
