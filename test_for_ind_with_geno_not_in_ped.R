
new_ped <- read_delim("../sheep/data/SNP_chip/Pedigree_SoaySheep_2019-07-11_SNPonly_SEJ_Update.txt", 
                      delim = " ",
                      col_types = c(rep("ccccccccc")))
sheep_geno_merged$ID[which(!(sheep_geno_merged$ID %in% new_ped$id))]

new_ped2 <- read_csv("../sheep/data/SNP_chip/Pedigree_SoaySheep_2019-07-11_toDB_SEJ_Update.csv")
sheep_geno_merged$ID[which(!(sheep_geno_merged$ID %in% new_ped2$Code))]


old_ped <- read_delim("../sheep/data/SNP_chip/20190208_Full_Pedigree.txt", delim = "\t")

sheep_ped <- read_delim("")

new_ped$Code[which(!(new_ped$Mother %in% old_ped$MOTHER))]

sheep_geno_merged$ID[!(sheep_geno_merged$ID %in% new_ped$Code)]
