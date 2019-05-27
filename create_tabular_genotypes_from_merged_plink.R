# on mac
plink_geno_path <- "../sheep/data/SNP_chip/"

# plink name
sheep_plink_name <- "merged_sheep_geno"
# read merged plink data
sheep_bed <- paste0(plink_geno_path, sheep_plink_name, ".bed")
sheep_bim <- paste0(plink_geno_path, sheep_plink_name, ".bim")
sheep_fam <- paste0(plink_geno_path, sheep_plink_name, ".fam")
full_sample <- read.plink(sheep_bed, sheep_bim, sheep_fam)

# filter names of snps on one chromosome
all_chr_snps <- full_sample$map %>% filter(chromosome %in% 1:26) %>% .$snp.name

# filter those snps from full dataset and coerce from raw to numeric
sheep_geno <- as(full_sample$genotypes[, all_chr_snps], Class = "numeric")

# plink puts double ids when merging, extract unique ids here
sheep_ids <- unlist(lapply(str_split(rownames(full_sample$genotypes), "\\.", 2), function(x) x[[2]]))
rownames(sheep_geno) <- sheep_ids

# clear some space
rm(full_sample)

# make tibble and put rownames as ID column
sheep_geno <- as_tibble(sheep_geno, rownames = "ID")

###### merge individuals on both ld and hd chip #####
#dup_ids <-  which(duplicated(sheep_ids))
setDT(sheep_geno)

# function to merge SNP data from the same individual, when it is both 
# in the HD and LD chip
# if genotypoe is missing on one chip, take the existing genotype
# if genotypes differ between chips, set NA

merge_geno <- function(vec) {
    # vec <- as.numeric(vec)
    if (length(vec) == 1) return(as.numeric(vec))
    if (sum(is.na(vec)) == 2) return(as.numeric(NA))
    if (sum(is.na(vec)) == 1) return(vec[!is.na(vec)])
    if (sum(is.na(vec)) == 0) {
        if (vec[1] == vec[2]){
            return(vec[1]) 
        } else {
            print(vec)
            return(as.numeric(NA))
        }
    }
}

# still needs to be merged properly here
#sheep_geno_dup <- sheep_geno[ID %chin% sheep_geno[["ID"]][1:188]]
#sheep_geno_merged <- sheep_geno_dup[, lapply(.SD, merge_geno), by=ID]

sheep_geno_filt <- sheep_geno[!duplicated(sheep_geno[["ID"]])]




# write to file with col names for masking script
fwrite(sheep_geno_filt, paste0("data/hdld_geno_merged.txt"), 
       sep = " ", col.names = TRUE, na = "9")

# grep -E '^6106|^2069|^6593|^1881|^1258|^2322|^1097|^4741|^1101|^2167'

