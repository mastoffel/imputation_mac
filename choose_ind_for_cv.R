library(data.table)
library(tidyverse)
# already cved 2322 1258 1097

ids <- read_delim("../sheep_ID/data/hd_ids.txt", " ")
ids

# choose 50
set.seed(144)
first50 <- sample(ids$id, size = 50, replace = FALSE)
write_lines(first50, "data/ids_first50.txt")


# old


# already cved full run :
cv_done <- as.numeric(read_lines("data/inds_for_cv_full_first10.txt"))

ids <- fread("data/genotypes_imputed_2018.txt", select = 1, nrows = 188, data.table = FALSE, header = FALSE)
ids_vec <- unlist(ids)

ids_vec_sub <- ids_vec[!(ids_vec %in% cv_done)]

# choose 10 individuals for full run cv
write_lines(as.character(sample(ids_vec_sub, 15)), "data/inds_for_cv_full_second15.txt")


# 7 individuals not done from second 15


# for sex chromosome
# 3799 4532 1746