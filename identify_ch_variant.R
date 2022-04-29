library(dplyr)

df <- read.csv("2022-03-28_16-11-18_MCD289_genotypes.csv")

df2 <- df %>%
  filter(Covered.Case == 3) %>%
  group_by(Variant.ID) %>%
  summarise(sample_name = Sample.Name) %>%
  mutate(child_in_variant = "GHART291b" %in% sample_name) %>%
  filter(child_in_variant == TRUE) %>%
  summarise(mother_in_variant = "GHART1506" %in% sample_name,
            father_in_variant = "GHART1542a" %in% sample_name) %>%
  mutate(genotype = case_when(mother_in_variant == FALSE & father_in_variant == FALSE ~ "0/0",
                          mother_in_variant == TRUE & father_in_variant == FALSE ~ "1/0",
                          mother_in_variant == FALSE & father_in_variant == TRUE ~ "0/1",
                          mother_in_variant == TRUE & father_in_variant == TRUE ~ "1/1")) %>%
  filter(genotype %in% c("1/0","0/1"))