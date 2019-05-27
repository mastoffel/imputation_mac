ind_org <- fread("results/geno_full_by_chr/ind_6106_org2")
ind_org <- sheep_geno[ID == 6106]
ind_imp <- fread("results/geno_full_by_chr/ind_6106_imp")[-2]

ind <- data.table::transpose(rbindlist(list(ind_org, ind_imp)))

df_trans <- ind  %>% 
    mutate_at(.vars = vars(contains("V")),
              funs(case_when(
                  . <= 0.02 ~ 0,
                  (. >= 0.98) & (. <= 1.02)  ~ 1,
                  (. >= 1.98) & (. <= 2.02) ~ 2,
                  TRUE ~ 9
              ))) %>% 
    setDT()

df_trans <- ind  %>% 
    mutate_at(.vars = vars(contains("V")),
              funs(case_when(
                  . <= 0.0 ~ 0,
                  (. >= 1) & (. <= 1)  ~ 1,
                  (. >= 2) & (. <= 2) ~ 2,
                  TRUE ~ 9
              ))) %>% 
    setDT()

df <- df_trans %>% 
    mutate_all(na_if, 9) 

df %>% 
    mutate(correct = V1==V2) %>% 
    summarise(prop_correct = sum(correct, na.rm = TRUE)/sum(!is.na(correct)))

# 6106 2069