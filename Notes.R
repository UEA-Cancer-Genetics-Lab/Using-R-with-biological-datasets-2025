#Notes

dat <- read_delim("16s_Community.tsv")
dat_longer <- dat |> 
  filter(primaryCat %in% c("CB","H")) |>
  pivot_longer(2:158)
ttest_res <- dat_longer |> group_by(name) |> t_test(value~primaryCat) |>
  adjust_pvalue(method = "BH") |> filter(!is.na(p)) |> arrange(p)

dat_16S <- read_delim("16s_Community.tsv") |> 
  select(sample_id, primaryCat, D_5__Arcanobacterium) |> 
  filter(primaryCat %in% c("CB","H")) |>
  rename(Arcanobacterium = D_5__Arcanobacterium)
t.test(Arcanobacterium ~ primaryCat, dat = dat_16S)

write_csv(dat_16S, "Arcanobacterium_16S.csv")


dat_longer <- dat |> 
#  filter(primaryCat %in% c("CB","L","I","H")) |>
#  mutate(primaryCat = factor(primaryCat, levels = c("CB", "L", "I", "H"))) |>
  pivot_longer(2:158)
  
aov_res <- dat_longer |> group_by(name) |> anova_test(value~primaryCat) |>
  data.frame |> arrange(p)

dat_16S2 <- read_delim("16s_Community.tsv") |> 
  select(sample_id, primaryCat, D_5__Veillonella) |> 
  rename(Veillonella = D_5__Veillonella)

aov_res <- aov(Veillonella ~ primaryCat, dat = dat_16S2)
summary(aov_res)

write_csv(dat_16S, "Arcanobacterium_16S_full.csv")

library(janitor)
cow_diet <- read_csv("cow_diet.csv") |>
  select(Sample, Diet, `3-HP`, Glutamate, Proline, Uracil, Xanthine) |>
  clean_names() |>
  mutate(diet = fct_collapse(factor(diet), A = 0, B = c(15,30), C = 45))

write_rds(cow_diet, "cow_diet.rds")

aov_cow <- aov(uracil ~ diet, data = cow_diet)

load(url("https://raw.githubusercontent.com/UEA-Cancer-Genetics-Lab/MMB-masterclass-biomarkers/main/16s_PCa_data.RData"))
# Remove those values that are less than 5% and convert to presence/absence
s16_community <- s16_community %>% mutate_if(is.numeric, ~1 * (. > 5))

# Only select genera with more than 2 hits
s16_community <- s16_community %>% select_if(function(col) is.character(col) || (is.numeric(col) && sum(col) >2))

# Merge taxa and survival
s16_merge <-  s16_clin_data %>% left_join(s16_community, by = join_by(sample_ID))

abbs_genera <- c("Ezakiella","Peptoniphilus","Porphyromonas","Anaerococcus","Fusobacterium")
s16_merge$abbs <- rowSums(s16_merge[,colnames(s16_merge) %in% abbs_genera]) > 0

write_rds(s16_merge, "data/s16_merge.rds")
