library(data.table)
library(readr)
dat <- fread("main_files_cv/Genotypes_snps_in_row.txt", nThread = 1, header = TRUE )
sessionInfo()


dat <- read_table("main_files_cv/Genotypes_snps_in_row.txt")


?read_table


install.packages("data.table")
dat <- fread("main_files_cv/Genotypes.txt", sep = " ")
