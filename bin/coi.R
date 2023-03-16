# LIBRARIES ----
## Data Management ----
require(argparse)
require(tidyverse)
require(dplyr)
require(data.table)
require(zoo)
library(haven)
library(sjlabelled)
library(measurements)

## Visualization and Plotting ----
require(ggpubr)
require(ggtext)
require(ggrepel)
require(scales)
require(cowplot)
require(corrplot)
require(colorspace)
require(viridis)
require(wesanderson)
require(jtools)
require(plotly)

## ARGUMENTS ----
parser <- ArgumentParser()

parser$add_argument("--dir",
                    help = "Working directory", 
                    default = "~/Dropbox (UMass Medical School)/Projects/DarwinsArk/git/darwins-ark-geno")

parser$add_argument("--dog",
                    help = "CSV file with Darwin's Ark `dogs` data", 
                    default = "/dat/dog/DarwinsArk_20221120_dogs.csv")

parser$add_argument("--gen",
                    help = "PLINK .fam file", 
                    default = "/dat/gen/DarwinsArk_gp-0.70_biallelic-snps_maf-0.001_geno-0.05_hwe-1e-20-midp-keep-fewhet_N-3465.fam")

parser$add_argument("--adm",
                    help = "CSV file with global ancestry results from ADMIXTURE", 
                    default = "/dat/adm/DarwinsArk_gp-0.70_biallelic-snps_maf-0.001_geno-0.05_hwe-1e-20-midp-keep-fewhet_N-3465_global-adm_supervised_K-115.csv")

parser$add_argument("--coi",
                    help = "CSV file with inbreeding", 
                    default = "/dat/coi/DarwinsArk_2022-11-20_inbreeding.csv")

parser$add_argument("--dapcoi",
                    help = "CSV file with inbreeding", 
                    default = "/dat/coi/DogAgingProject_2023-01-01_inbreeding.csv")

args <- parser$parse_args()
print(args)

### Set working directory ----
workDir = args$dir
setwd(workDir)

## LOAD DATA ----
dog = read_csv(file = paste(workDir,
                            args$dog,
                            sep = ""))

gen = read_tsv(file = paste(workDir,
                            args$gen,
                            sep = ""),
               col_names = c("fam",
                             "dog",
                             "pat",
                             "mat",
                             "sex",
                             "phe"))

adm = read_csv(file = paste(workDir,
                            args$adm,
                            sep = ""))

coi = read_csv(file = paste(workDir,
                            args$coi,
                            sep = ""))

dapcoi = read_csv(file = paste(workDir,
                            args$dapcoi,
                            sep = ""))

# STATS ----
coi.summary.all = coi %>%
  select(-dog) %>%
  describe() %>%
  as.data.table() %>%
  mutate(var = colnames(coi)[-1])

coi %>%
  merge(dog, by = "dog") %>%
  group_by(purebred=="no" | is.na(purebred)) %>%
  summarize(mean(coi),
            sd(coi),
            n())


load("~/Dropbox (UMass Medical School)/Projects/DogAgingProject/DAP_2021_HLES_dog_owner_v1.0.RData")

# PLOT ----
coi.all = bind_rows((coi %>%
                       merge((dog %>%
                                select(dog,
                                       purebred,
                                       breed)),
                             by = "dog") %>%
                       mutate(cohort = "Darwin's Ark")),
                    (dapcoi %>%
                       merge((HLES_dog_owner %>%
                                mutate(purebred= as_character(dd_breed_pure_or_mixed),
                                       breed= if_else(!is.na(dd_breed_pure),
                                                      as_character(dd_breed_pure),
                                                      as_character(dd_breed_pure_non_akc)))  %>% 
                                select(dog= dog_id,
                                       purebred,
                                       breed) %>%
                                mutate(cohort = "Dog Aging Project") %>%
                                mutate(breed = tolower(breed)) %>%
                                mutate(breed = if_else(breed== "american pitbull terrier",
                                                       "american pit bull terrier",
                                                       breed)))))) %>%
  mutate(pedigreed = if_else(purebred == "yes" | purebred == "Purebred",
                             TRUE,
                             FALSE,
                             FALSE)) %>%
  mutate(breed = if_else(is.na(breed),
                         "no reported breed",
                          breed))


coi.all %>%
  group_by(cohort,pedigreed) %>%
  summarise(mean = mean(coi),
            sd = sd(coi),
            n = n())

top.breeds = coi.all %>%
  group_by(breed) %>%
  summarize(n=n()) %>%
  filter(n>50)

library(EnvStats)
p = coi.all %>%
  merge(top.breeds, by = "breed") %>%
  arrange(n) %>%
  ggplot(aes(x = reorder(breed,n), 
             y = coi)) +
  facet_grid(cohort~pedigreed) +
  geom_point(position = position_jitter(), alpha = 0.2) +
  geom_boxplot(outlier.shape = NA, 
               fill = "white", 
               alpha = 0.5) +
  coord_flip() +
  stat_n_text() + 
  theme_pubr()

ggsave(filename = "~/Dropbox (UMass Medical School)/Thesis/figures/inbreeding.pdf",
       plot = p,
       device = "pdf",
       units = "in",
       width = 6*1.5,
       height= 4*1.5)

t.test((test %>% filter(class == "confirmed purebred") %>% pull(F.roh)), (test %>% filter(class == "mutt") %>% pull(F.roh)))$p.value
