# PREAMBLE ----

## LIBRARIES ----

### Data Management ----
require(argparse)
require(tidyverse)
require(dplyr)
require(data.table)
require(zoo)

### Visualization and Plotting ----
require(ggpubr)
require(ggtext)
require(ggrepel)
require(scales)
require(cowplot)
require(corrplot)
require(colorspace)
require(viridis)
require(PerformanceAnalytics)
require(wesanderson)
require(jtools)

### Dates and Times ----
require(anytime)
require(lubridate)

### Free Text and Strings ----
require(stringr)
require(stringi)
require(stringdist)
require(english)
require(tm)
require(topicmodels)
require(lattice)
require(tidytext)
require(proxy)
require(corpus)
require(gender)
require(rxnorm)
require(RxNormR)
require(qdap)
data(stop_words)

### Locations ----
require(zipcodeR)

### Correlation and Regression ----
require(pspearman)
require(jtools)
require(caret)
require(ade4)
require(rstatix)

### Missing Data ----
require(naniar)
require(finalfit)

### Dimensional Reduction ----
require(psych)
require(nFactors)
require(FactoMineR)
require(factoextra)
require(corrr)
require(car)
require(BBmisc)
require(mice)
require(phateR)


## ARGUMENTS ----
parser <- ArgumentParser()
parser$add_argument("--dir",
                    help = "Working directory",
                    default = "~/Dropbox (UMass Medical School)/Projects/DarwinsArk/git/darwins-ark-data")

parser$add_argument("--date",
                    help = "Data freeze date",
                    default = "20221120")

parser$add_argument("--geno",
                    help = "Genetic dataset IDs",
                    default = "~/Dropbox (UMass Medical School)/Projects/DarwinsArk/git/darwins-ark-data/gwa/gen/DarwinsArk_gp-0.70_biallelic-snps_maf-0.001_geno-0.05_hwe-1e-20-midp-keep-fewhet_N-3465.fam")

args <- parser$parse_args()
print(args)

### Set working directory ----
workDir = args$dir
setwd(workDir)

### Set data freeze ----
freezeDate = args$date
freezeDir = paste(workDir, "data", "database", "data_freeze", freezeDate, sep = "/")

### Load data tables ----
dogs = read_csv(paste(
  workDir,
  "/dat/",
  "DarwinsArk_",
  as.character(freezeDate),
  "_dogs.csv",
  sep = ""
))

questions = read_csv(paste(
  workDir,
  "/dat/",
  "DarwinsArk_",
  as.character(freezeDate),
  "_questions.csv",
  sep = ""
))

answers = read_csv(paste(
  workDir,
  "/dat/",
  "DarwinsArk_",
  as.character(freezeDate),
  "_answers.csv",
  sep = ""
))

scores = read_csv(paste(
  workDir,
  "/dat/",
  "DarwinsArk_",
  as.character(freezeDate),
  "_factor-scores.csv",
  sep = ""
))

# GENETIC CORRELATIONS ----

## Load genetic correlations ----
input.gen.cor = "~/Dropbox (UMass Medical School)/Thesis/repo/dat/gen/DarwinsArk_2022-11-20_genetic-correlations.tsv"

gen.cor = read_tsv(input.gen.cor)

## Check basic stats and distribution ----
gen.cor %>%
  summarize(mean = mean(rG),
            sd = sd(rG))

gen.cor %>%
  ggplot(aes(x = rG)) +
  geom_density()

## Set phenotypes ----
gen.phe = sort(unique(c(gen.cor$phe_A,
                        gen.cor$phe_B)))

# PHENOTYPE CORRELATIONS ----

## Get fa.[]-filled-by-mean phenotypes ----
fa.all = scores %>%
  mutate(ind = paste(
    "fa.",
    str_pad(
      string = as.character(factor),
      width = 2,
      side = "left",
      pad = "0"
    ),
    "-filled-by-mean",
    sep = ""
  )) %>%
  mutate(phe = round(score_fill_mean_norm,
                     digits = 2)) %>%
  select(ind,
         dog,
         phe)

## Get bq.[] phenotypes ----
bq.all = answers %>%
  filter(question %in%
           c(1:110,
             144:155,
             165:175,
             177:182,
             184,185,
             195:240)) %>%
  select(question,dog,answer,option,age) %>%
  mutate(answer = if_else(option %in% c("I don't know",
                                        "I'm not sure",
                                        "Maybe"),
                          "",
                          answer)) %>%
  mutate(ind = paste("bq.",
                       str_pad(question,
                               3,
                               "left",
                               "0"), 
                       sep = "")) %>%
  mutate(phe = as.numeric(answer)) %>%
  select(ind,
         dog,
         phe)

## Get mq.[] and mp.[] phenotypes ----

### Get mq.242 ----
mq.242 = answers %>%
  filter(question == 242) %>%
  select(dog,
         answer,
         age) %>%
  mutate(answer = gsub("-lb",
                       " lb",
                       gsub("-kg",
                            " kg",
                            answer))) %>%
  separate(
    answer,
    into = c('value',
             'unit'),
    sep = " ",
    convert = T
  ) %>%
  mutate(value = abs(value)) %>%
  mutate(phe = if_else(unit == "lb",
                          value,
                          value * 2.2)) %>%
  mutate(ind = "mq.242") %>%
  select(ind,
         dog,
         phe)

### Get continuous mq.[] traits ----
mq.quant = answers %>%
  filter(question %in% c(121:128,
                         243:250)) %>%
  select(question,dog,option,answer) %>%
  mutate(answer = if_else(question == 125 & answer %in% c(5,6),
                          NA_character_,
                          answer)) %>%
  filter(!question %in% (questions %>% filter(style %in% c("MultiChoices",
                                                           "MultiChoicesWithoutNone",
                                                           "MultiChoicesTooltip",
                                                           "Colors")) %>% pull(id))) %>%
  filter(!option %in% c("I don't know",
                        "I'm not sure",
                        "Maybe")) %>%
  select(dog,question,answer) %>%
  unique() %>%
  mutate(answer = as.numeric(answer)) %>%
  mutate(ind = paste("mq.",
                     question,
                     sep = "")) %>%
  select(ind,
         dog,
         phe = answer)

### Get discrete mp.[] traits ----
#### mp.249.eye-color ----
mp_eye_color = answers %>%
  filter(question == 249) %>%
  as.data.table()

mp_eye_color[answer == 4, phe := 3]
mp_eye_color[answer == 3, phe := 3]
mp_eye_color[answer == 0, phe := 2]
mp_eye_color[answer == 2, phe := 1]
mp_eye_color[answer == 1, phe := 0]

mp_eye_color = mp_eye_color %>% 
  mutate(ind = "mp.249.eye-color") %>%
  select(ind,
         dog,
         phe) %>%
  group_by(dog) %>%
  slice_head(n=1) %>%
  ungroup()

#### mp.126.249.eye-heter ----
mp_eye_heter = bind_rows((answers %>%
                            filter(question == 249) %>%
                            select(dog,option) %>%
                            group_by(dog) %>%
                            mutate(n = length(unique(option))) %>%
                            group_by(dog) %>%
                            summarize(phe = n > 1) %>%
                            ungroup() %>%
                            unique()),
                         (answers %>%
                            filter(question == 126) %>%
                            select(dog,option) %>%
                            filter(option != "I don't know") %>%
                            mutate(phe = if_else(option == "Yes",
                                                 1,
                                                 0)) %>%
                            select(dog,phe) %>%
                            unique())) %>%
  mutate(ind = "mp.126.249.eye-heter") %>%
  select(ind,
         dog,
         phe) %>%
  ungroup() %>%
  group_by(dog) %>%
  slice_head(n=1) %>%
  ungroup() %>%
  unique()

#### mp.250.pad-[] ----
mp_pad_color = answers %>%
  filter(question == 250) %>%
  select(dog,answer) %>%
  mutate(value = 1) %>%
  pivot_wider(names_from = "answer",
              values_from = "value",
              values_fn = min, 
              values_fill = 0) %>%
  rename(`mp.250.pad-liver` = `1`,
         `mp.250.pad-black` = `0`,
         `mp.250.pad-grey` = `2`,
         `mp.250.pad-pink` = `3`) %>%
  pivot_longer(cols = -dog,
               names_to = "ind",
               values_to = "phe")

#### mp.248.[] ----
mp_features = answers %>%
  filter(question == 248) %>%
  select(dog,answer) %>%
  mutate(value = 1) %>%
  pivot_wider(names_from = "answer",
              values_from = "value",
              values_fn = min, 
              values_fill = 0) %>%
  rename(`mp.248.jowled` = `0`,
         `mp.248.wrinkly` = `1`,
         `mp.248.freckled` = `2`,
         `mp.248.full_white` = `3`,
         `mp.248.vitiligo` = `4`,
         `mp.248.ridged` = `5`,
         `mp.248.hairless` = `6`,
         `mp.248.underbite` = `7`,
         `mp.248.tongue_color` = `8`) %>%
  select(-`9`) %>%
  pivot_longer(cols = -dog,
               names_to = "ind",
               values_to = "phe")

#### mp.245.[] ----
mp_head_shape = answers %>%
  filter(question == 245) %>%
  select(dog,answer) %>%
  mutate(value = 1) %>%
  pivot_wider(names_from = "answer",
              values_from = "value",
              values_fn = min, 
              values_fill = 0) %>%
  mutate(`mp.245.brachycephaly` = if_else(`2`==1 | `4`==1 | `5`==1,
                                          1,
                                          0),
         `mp.245.dolichocephaly` = if_else(`0`==1 | `1`==1,
                                           1,
                                           0)) %>%
  rename(`mp.245.head-shape-A` = `0`,
         `mp.245.head-shape-B` = `1`,
         `mp.245.head-shape-C` = `2`,
         `mp.245.head-shape-D` = `3`,
         `mp.245.head-shape-E` = `4`,
         `mp.245.head-shape-F` = `5`) %>%
  pivot_longer(cols = -dog,
               names_to = "ind",
               values_to = "phe")

#### mp.243.coat-color-[] ----
mp_coat_colors_243 = answers %>%
  filter(question == 243) %>%
  select(dog,answer) %>%
  mutate(value = 1) %>%
  pivot_wider(names_from = "answer",
              values_from = "value",
              values_fn = min, 
              values_fill = 0) %>%
  rename(`mp.243.coat-color-black` = `0`,
         `mp.243.coat-color-brown` = `1`,
         `mp.243.coat-color-white` = `2`,
         `mp.243.coat-color-red` = `3`,
         `mp.243.coat-color-yellow` = `4`,
         `mp.243.coat-color-dilute-black` = `5`,
         `mp.243.coat-color-dilute-brown` = `6`,
         `mp.243.coat-color-skin` = `7`,
         `mp.243.coat-color-tan` = `8`,
         `mp.243.coat-color-cream` = `9`) %>%
  pivot_longer(cols = -dog,
               names_to = "ind",
               values_to = "phe")

#### mp.122.coat-color-[] and mp.122.coat-pattern-[] ----
mp_coat_colors_122 = answers %>%
  filter(question == 122) %>%
  select(dog,answer) %>%
  mutate(value = 1) %>%
  pivot_wider(names_from = "answer",
              values_from = "value",
              values_fn = min, 
              values_fill = 0) %>%
  rename(`mp.122.coat-color-white` = `0`,
         `mp.122.coat-color-red` = `1`,
         `mp.122.coat-color-yellow` = `2`,
         `mp.122.coat-color-gray` = `3`,
         `mp.122.coat-color-chocolate-brown` = `4`,
         `mp.122.coat-color-pure-black` = `5`,
         `mp.122.coat-pattern-merle` = `6`,
         `mp.122.coat-pattern-brindle` = `7`) %>%
  select(-`8`) %>%
  pivot_longer(cols = -dog,
               names_to = "ind",
               values_to = "phe")

#### mp.244.coat-pattern-[] ----
mp_coat_patterns = answers %>%
  filter(question == 244) %>%
  select(dog,answer) %>%
  mutate(value = 1) %>%
  pivot_wider(names_from = "answer",
              values_from = "value",
              values_fn = min, 
              values_fill = 0) %>%
  rename(`mp.244.coat-pattern-striping` = `0`,
         `mp.244.coat-pattern-patching` = `1`,
         `mp.244.coat-pattern-splotching` = `2`,
         `mp.244.coat-pattern-spotting` = `3`,
         `mp.244.coat-pattern-roaning` = `4`,
         `mp.244.coat-pattern-ticking` = `5`,
         `mp.244.coat-pattern-none` = `6`) %>%
  pivot_longer(cols = -dog,
               names_to = "ind",
               values_to = "phe")

#### mq.244.coat-pattern-spots ----
mp.coat_patterns_spots = answers %>%
  filter(question == 244) %>%
  select(dog,answer) %>%
  filter(answer %in% c(3,4,5)) %>%
  select(dog,answer) %>%
  as.data.table()
mp.coat_patterns_spots[answer==3, `mq.244.coat-pattern-spots` := 3]
mp.coat_patterns_spots[answer==4, `mq.244.coat-pattern-spots` := 2]
mp.coat_patterns_spots[answer==5, `mq.244.coat-pattern-spots` := 1]

mp.coat_patterns_spots = mp.coat_patterns_spots %>%
  group_by(dog) %>%
  summarize(`mq.244.coat-pattern-spots` = max(`mq.244.coat-pattern-spots`)) %>%
  bind_rows((answers %>%
               filter(question == 244) %>%
               filter(!dog %in% mp.coat_patterns_spots$dog) %>%
               select(dog) %>%
               unique() %>%
               mutate(`mq.244.coat-pattern-spots` = 0))) %>%
  arrange(dog) %>%
  pivot_longer(cols = -dog,
               names_to = "ind",
               values_to = "phe")

#### mp.125.ears-[] ----
mp.earshape = answers %>%
  filter(question == 125) %>%
  select(dog,answer) %>%
  mutate(value = 1) %>%
  pivot_wider(names_from = "answer",
              values_from = "value",
              values_fn = min, 
              values_fill = 0) %>%
  rename(`mp.125.ears-pendant` = `0`,
         `mp.125.ears-dropped` = `1`,
         `mp.125.ears-rose-button` = `2`,
         `mp.125.ears-pricked` = `3`,
         `mp.125.ears-cropped` = `4`,
         `mp.125.ears-other` = `5`) %>%
  pivot_longer(cols = -dog,
               names_to = "ind",
               values_to = "phe")

### Set mq.all ----
mq.all = bind_rows(mq.242,
                   mq.quant,
                   mp_eye_color,
                   mp_eye_heter,
                   mp_features,
                   mp_head_shape,
                   mp.earshape,
                   mp_pad_color,
                   mp_coat_colors_122,
                   mp_coat_colors_243,
                   mp.coat_patterns_spots,
                   mp_coat_patterns)

## Prepare phenotype correlation matrix: ----
phe.all = bind_rows(fa.all,
                    bq.all,
                    mq.all)

phe.cor = phe.all %>%
  pivot_wider(id_cols = dog,
              names_from = ind,
              values_from = phe) %>%
  select(-dog) %>%
  cor_mat(method = "pearson")

phe.cor.mat = phe.cor %>%
  select(-rowname) %>%
  as.matrix()
rownames(phe.cor.mat) = phe.cor$rowname

# SET MATRICES ----

## Get common phenotypes ----
com.phe = intersect(phe.cor$rowname,
                    gen.phe)
  
## Make phenotype correlation matrix ----
phe.cor.mat = phe.cor %>%
  filter(rowname %in% com.phe) %>%
  select(all_of(com.phe)) %>%
  as.matrix()
rownames(phe.cor.mat) = colnames(phe.cor.mat)

## Make genetic correlation matrix ----
gen.cor.mat = gen.cor %>%
  filter((phe_A %in% com.phe) & (phe_B %in% com.phe)) %>%
  mutate(phe1 = factor(x = phe_A,
                       levels = com.phe),
         phe2 = factor(x = phe_B,
                       levels = com.phe)) %>%
  select(phe1,
         phe2,
         gen.cor= rG) %>%
  complete(phe1,phe2) %>%
  arrange(phe1,phe2) %>%
  pivot_wider(id_cols = "phe1",
              names_from = "phe2",
              values_from = "gen.cor") %>%
  select(-phe1) %>%
  as.matrix()

colnames(gen.cor.mat) = com.phe
rownames(gen.cor.mat) = com.phe

# PERFORM MANTEL TEST ----
## Make gen.cor.mat symmetric ----
gen.cor.mat[lower.tri(gen.cor.mat)] = t(gen.cor.mat)[lower.tri(gen.cor.mat)]

gen.cor.mat.nona = ifelse(is.na(gen.cor.mat), 0, gen.cor.mat)
phe.cor.mat.nona = ifelse(is.na(phe.cor.mat), 0, phe.cor.mat)

## Convert to correlation distances (d = 1 - |r|) ----
gen.cor.mat.inv = matrix(1, 
                         nrow = nrow(gen.cor.mat.nona), 
                         ncol = ncol(gen.cor.mat.nona)) - abs(gen.cor.mat.nona)

phe.cor.mat.inv = matrix(1, nrow = nrow(phe.cor.mat.nona), 
                         ncol = ncol(phe.cor.mat.nona)) - abs(phe.cor.mat.nona)

mantel.rtest(m1 = as.dist(phe.cor.mat.inv),
             m2 = as.dist(gen.cor.mat.inv),
             nrepet = 100000)

# CORRELATION OF CORRELATIONS ----
cxc.cor = phe.cor %>%
  cor_gather() %>%
  merge(gen.cor,
        by.x = c("var1","var2"),
        by.y = c("phe_A","phe_B"))

cor.test(cxc.cor$cor,cxc.cor$rG) %>% View()
cor.test(cxc.cor$cor,cxc.cor$rG)$p.value

p= cxc.cor %>%
  mutate(pair = paste(var1,
                      "x",
                      var2)) %>%
  mutate(mpleio = (var1 %like% "mq." & var2 %like% "fa.") |
           (var1 %like% "fa." & var2 %like% "mq.") |
           (var1 %like% "mq." & var2 %like% "bq.") |
           (var1 %like% "bq." & var2 %like% "mq.") |
           (var1 %like% "mp." & var2 %like% "fa.") |
           (var1 %like% "fa." & var2 %like% "mp.") |
           (var1 %like% "mp." & var2 %like% "bq.") |
           (var1 %like% "bq." & var2 %like% "mp.")) %>%
  ggplot(aes(x = rG,
             y = cor,
             label = pair)) +
  geom_point(aes(color = mpleio),
             alpha = 0.5,
             shape = 1) +
  scale_x_continuous(limits = c(-1,1)) +
  scale_y_continuous(limits = c(-1,1)) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  geom_text_repel() +
  theme_pubr()

ggsave(filename = "~/Dropbox (UMass Medical School)/Thesis/figures/GEN_COR_PAIRS.pdf",
       plot = p,
       device = "pdf",
       units = "in",
       width = 6*2,
       height = 4*2)

# FACTORS AND ITEMS LOADED ----
fa.loadings = read_tsv("~/Dropbox (UMass Medical School)/Projects/DarwinsArk/git/darwins-ark-data/dat/DarwinsArk_20221120_factor-analysis_discovery-no-na_qn-182_dn-7883_fn-25.loadings.tsv")

fa.load.gen = fa.loadings %>%
  mutate(idf = paste("fa.",
                    str_pad(factor,2,"left","0"),
                    "-filled-by-mean",
                    sep = ""),
         idb = paste("bq.",
                     str_pad(question,3,"left","0"),
                     sep = "")) %>%
  merge(gen.cor, 
        by.x = c("idf","idb"),
        by.y = c("phe_B","phe_A"))


fa.load.gen = fa.loadings %>%
  mutate(idf = paste("fa.",
                     str_pad(factor,2,"left","0"),
                     "-filled-by-mean",
                     sep = ""),
         idb = paste("bq.",
                     str_pad(question,3,"left","0"),
                     sep = "")) %>%
  merge(gen.cor, 
        by.x = c("idf","idb"),
        by.y = c("phe_B","phe_A"),
        all.y = T) %>%
  filter(idf %like% "fa." & idb %like% "bq.") %>%
  mutate(factored = !is.na(pattern)) %>%
  mutate(factor = gsub("-filled-by-mean","",idf)) %>%
  select(factor,factored,idf,idb,pattern,structure,rG)

fa.load.gen %>%
  ggplot(aes(x = factored,
             y = abs(rG))) +
  geom_violin(alpha = 0.5) +
  geom_point(alpha = 0.25) +
  facet_wrap(.~factor) +
  theme_pubr()

p =fa.load.gen %>%
  ggplot(aes(x = abs(rG),
             y = abs(pattern),
             label = paste(idf,idb))) +
  geom_point()

# UNIQUENESS VS GENETIC CORRELATIONS ----
# Survey item uniqueness ----
fa.uniqueness = read_tsv("~/Dropbox (UMass Medical School)/Projects/DarwinsArk/git/darwins-ark-data/dat/DarwinsArk_20221120_factor-analysis_discovery-no-na_qn-182_dn-7883_fn-25.uniqueness.tsv")

fa.uniqueness.gen.cor = bind_rows((fa.uniqueness %>%
                                     mutate(id = paste("bq.",
                                                       str_pad(string = question,
                                                               width = 3,
                                                               pad = "0"),
                                                       sep = "")) %>%
                                     merge((gen.cor %>%
                                              select(id= phe_A,
                                                     to = phe_B,
                                                     rG)), 
                                           by = "id")),
                                  (fa.uniqueness %>%
                                     mutate(id = paste("bq.",
                                                       str_pad(string = question,
                                                               width = 3,
                                                               pad = "0"),
                                                       sep = "")) %>%
                                     merge((gen.cor %>%
                                              select(id= phe_B,
                                                     to = phe_A,
                                                     rG)), 
                                           by = "id")))

p = fa.uniqueness.gen.cor %>%
  group_by(id) %>%
  summarize(uniqueness,
            n = n(),
            mean = mean(abs(rG)),
            sd = sd(abs(rG))) %>%
  ggplot(aes(x = uniqueness,
             y = mean,
             ymin = mean-sd,
             ymax = mean+sd,
             label = id)) +
  geom_point(shape = 21) + 
  geom_errorbar(alpha = 0.3) +
  theme_pubr()

ggsave(filename = "~/Dropbox (UMass Medical School)/Thesis/figures/GEN_COR_UNIQUENESS.pdf",
       plot = p,
       device = "pdf",
       units = "in",
       width = 6,
       height = 4)

fa.uniqueness.gen.cor %>%
  filter(!to %like% "mq." & !to %like% "mp.") %>%
  ggplot(aes(x = uniqueness,
             y = rG)) +
  geom_point(alpha=0.25) +
  theme_pubr()

matrix_data <- spread((gen.cor %>%
                         select(phe_A,
                                phe_B,
                                rG)), 
                      phe_A, 
                      rG)

complete_matrix_data <- complete(matrix_data, 
                                 phe_B, 
                                 fill = list(rG = NA))


gen.cor.mat= gen.cor %>%
  select(phe_A,
         phe_B,
         rG) %>%
  pivot_wider(names_from = phe_A,
              values_from = rG) %>%
  column_to_rownames('phe_B') %>%
  as.matrix()

gen.cor.mat= gen.cor.mat[order(rownames(gen.cor.mat)), order(colnames(gen.cor.mat))]

gen.cor.mat[upper.tri(gen.cor.mat)] <- t(gen.cor.mat)[upper.tri(gen.cor.mat)]


# CORR PLOT ----

## Set phenotypes ----
phe = sort(unique(c(cx$phe1,cx$phe2)))

## Complete pairs of phenotypes ----
gx = gx %>%
  mutate(phe1 = factor(x = phe1,
                       levels = phe),
         phe2 = factor(x = phe2,
                       levels = phe)) %>%
  complete(phe1,phe2) %>%
  arrange(phe1,phe2) %>%
  pivot_wider(id_cols = "phe1",
              names_from = "phe2",
              values_from = "corr") %>%
  select(-phe1) %>%
  as.matrix()

colnames(gx) = phe
rownames(gx) = phe

cx = cx %>%
  mutate(phe1 = factor(x = phe1,
                       levels = phe),
         phe2 = factor(x = phe2,
                       levels = phe)) %>%
  complete(phe1,phe2) %>%
  arrange(phe1,phe2) %>%
  pivot_wider(id_cols = "phe1",
              names_from = "phe2",
              values_from = "corr") %>%
  select(-phe1) %>%
  as.matrix()

colnames(cx) = phe
rownames(cx) = phe

## Combine matrices ----
dx = gx
dx[lower.tri(dx, diag = FALSE)] = cx[lower.tri(cx, diag = FALSE)]
diag(dx) = 1

## Plot ----
pdf(out)
corrplot.mixed(dx,
               upper = "circle",
               lower = "circle",
               is.corr = F)
dev.off()
