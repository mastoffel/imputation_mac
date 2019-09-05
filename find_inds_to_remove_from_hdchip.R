

# plink name
sheep_plink_name <- "../sheep/data/SNP_chip/Plates_1-2_HD_QC1"
# read merged plink data
sheep_bed <- paste0(sheep_plink_name, ".bed")
sheep_bim <- paste0(sheep_plink_name, ".bim")
sheep_fam <- paste0(sheep_plink_name, ".fam")
full_sample <- read.plink(sheep_bed, sheep_bim, sheep_fam)

# plink name
sheep_plink_name <- "../sheep/data/SNP_chip/older_genotypes/20140214_SheepHD_QC1_Polym"
# read merged plink data
sheep_bed <- paste0(sheep_plink_name, ".bed")
sheep_bim <- paste0(sheep_plink_name, ".bim")
sheep_fam <- paste0(sheep_plink_name, ".fam")
full_sample2 <- read.plink(sheep_bed, sheep_bim, sheep_fam)

nrow(full_sample$fam)
nrow(full_sample2$fam)

fam_id <- which(!(full_sample$fam$member %in% full_sample2$fam$member))
ids <- full_sample$fam$member[which(!(full_sample$fam$member %in% full_sample2$fam$member))]
df <- data.frame(ids = fam_id, ids2 = ids)
write_delim(df, "../sheep/data/SNP_chip/ids_to_remove", col_names = FALSE, delim = "\t")
