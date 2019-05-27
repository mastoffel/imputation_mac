library(data.table)
library(tidyverse)
# already cved 2322 1258 1097

?fread
ids <- fread("data/hdld_geno_merged.txt", select = 1, nrows = 188, data.table = FALSE)
ids_vec <- unlist(ids)
cv_done <- c(2322, 1258, 1097)

ids_vec_sub <- ids_vec[!(ids_vec %in% cv_done)]

# choose 7 more
sample(ids_vec_sub, 7)

write_lines(as.character(sample(ids_vec_sub, 7)), "data/new_ind_for_cv.txt")
