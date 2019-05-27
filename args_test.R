# capture command line test
library(purrr)
args <- commandArgs(trailingOnly=TRUE)
chr_num <- args[[1]]
if (purrr::is_empty(chr_num))  stop("Provide a chromosome number as argument")

# folder with all the files
main_files <- paste0("/exports/eddie/scratch/v1mstoff/all_chr_cv/chr_", chr_num, "/AI_main_files/")
print(main_files)
#print(as.numeric(args[[1]]))
#print(str(args))