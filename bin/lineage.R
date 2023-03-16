# LIBRARIES ----
## Data Management ----
require(argparse)
require(tidyverse)
require(dplyr)
require(data.table)
require(zoo)

## Dimensional Reduction ----
require(DelayedMatrixStats)
require(phateR)
require(slingshot)

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

parser$add_argument("--pca",
                    help = "PCA results", 
                    default = "/dat/pca/DarwinsArk_gp-0.70_biallelic-snps_maf-0.001_geno-0.05_hwe-1e-20-midp-keep-fewhet_N-3465_indep-pairwise_kb-250_r2-0.2.eigenvec")

parser$add_argument("--pca",
                    help = "PCA results", 
                    default = "/dat/pca/DarwinsArk_gp-0.70_biallelic-snps_maf-0.001_geno-0.05_hwe-1e-20-midp-keep-fewhet_N-3465_GlobalAncestrySNPs_RefMerge.eigenvec")

#DogAgingProject_gp-0.70_biallelic-snps_N-6358_maf-0.001_geno-0.05_hwe-1e-20-midp-keep-fewhet_indep-pairwise_kb-250_r2-0.2.eigenvec

parser$add_argument("--adm",
                    help = "CSV file with global ancestry results from ADMIXTURE", 
                    default = "/dat/adm/DarwinsArk_gp-0.70_biallelic-snps_maf-0.001_geno-0.05_hwe-1e-20-midp-keep-fewhet_N-3465_global-adm_supervised_K-115.csv")

parser$add_argument("--coi",
                    help = "CSV file with inbreeding", 
                    default = "/dat/coi/DarwinsArk_2022-11-20_inbreeding.csv")

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
pca = read_tsv(file = paste(workDir,
                            args$pca,
                            sep = ""))
coi = read_csv(file = paste(workDir,
                            args$coi,
                            sep = ""))
this_set = "DarwinsArk_with-context"
#this_set = "DarwinsArk_cohort-only"
#this_set = "DogAgingProject_cohort-only"
# INFER LINEAGES ----
## label reference data ----
#pca = pca %>%
 # mutate(reference = `#FID` != IID)

## convert PCA to matrix ----
pca.matrix = as.matrix(pca %>% select(-c("#FID","IID")))

## unsupervised ----
phate.output = phate(pca.matrix, 
                     ndim=2, 
                     knn=60, 
                     gamma=0, 
                     decay=100)

phate.clusters = cluster_phate(phate.output, k = 30) %>%
  as.character()

phate.output.df = as.data.frame(phate.output$embedding)
phate.output.df$cluster = phate.clusters
phate.output.df$dog = pca$IID

phate.output.df = phate.output.df %>% as_tibble()

p = ggplot(phate.output.df,
       aes(x = PHATE1,
           y = PHATE2,
           color = as.factor(cluster))) +
  geom_point(alpha = 0.2) +
  geom_point(size = 3, 
             color = "red", 
             x = median(phate.output.df$PHATE1), 
             y = median(phate.output.df$PHATE2)) +
  theme_pubr()

ggplotly(p)

start.clusters = c("29","10","25","15","27")
end.clusters = c(14,8,18,6,2,17,1,3,7,4,9,12,13)

phate.mat = as.matrix(phate.output$embedding)

phate.sling = slingshot(data = phate.mat,
                        clusterLabels = phate.clusters,
                        start.clus = start.clusters)
                        #end.clus = end.clusters)

phate.pseudotimes = slingPseudotime(x = phate.sling)

curves <- slingCurves(SlingshotDataSet(phate.sling), as.df = TRUE)

p = phate.output.df %>%
  ggplot(aes(x = PHATE1,
             y = PHATE2)) +
  geom_point(aes(text = dog),
             color = "#505050",
             alpha = 0.25) +
  geom_path(data = curves %>% 
              arrange(Order),
            aes(group = Lineage)) +
  scale_x_continuous(limits = c(-0.05,0.05)) +
  scale_y_continuous(limits = c(-0.05,0.05)) +
  theme_pubr()

i_knn=60

ggsave(filename = paste(workDir,
                        "/fig/",
                        "lineage_phate",
                        ".set-", this_set,
                        ".knn-", i_knn,
                        ".pdf",
                        sep = ""),
       plot = p,
       width = 2.5*2,
       height = 2*2,
       dpi = 300,
       device = "pdf")

## perform PHATE at range of nearest neighbors ----
phate = list()
phate.curves = list()
phate.ptvals = list()
i=1
for (i_knn in c(60,80,100)) {
  print(i_knn)
  phate.output = phate(pca_matrix, 
                       ndim=2, 
                       knn=i_knn, 
                       gamma=0, 
                       decay=100)

  phate.clusters = cluster_phate(phate.output, 
                                 k = 50)

  phate[[i]] = as.data.frame(phate.output$embedding)
  phate[[i]]$knn = i_knn
  phate[[i]]$cluster = phate.clusters
  
  phate.clusters.starts = phate[[i]] %>%
    mutate(med_PHATE1 = median(PHATE1),
           med_PHATE2 = median(PHATE2)) %>%
    filter(PHATE1 < med_PHATE1+0.01 & PHATE1 > med_PHATE1-0.01 & PHATE2 < med_PHATE2+0.01 & PHATE2 > med_PHATE2-0.01) %>% 
    pull(cluster) %>%
    unique()
    
  phate.mat = as.matrix(phate.output$embedding)
  
  phate.sling = slingshot(data = phate.mat,
                          clusterLabels = as.character(phate.clusters),
                          start.clus = as.character(phate.clusters.starts))
  
  phate.pseudotimes = slingPseudotime(x = phate.sling)
  
  curves <- slingCurves(SlingshotDataSet(phate.sling), as.df = TRUE)
  
  phate.ptvals[[i]] = as.data.frame(phate.pseudotimes)
  phate.ptvals[[i]]$knn = i_knn
  
  phate.curves[[i]] = curves
  phate.curves[[i]]$knn = i_knn
  
  p = ggplot(phate[[i]],
             aes(x = PHATE1,
                 y = PHATE2)) +
    geom_point(color = "#505050",
               alpha = 0.25) +
    geom_path(data = phate.curves[[i]] %>% 
                arrange(Order),
              aes(group = Lineage)) +
    theme_pubr()
  
  ggsave(filename = paste(workDir,
                          "/fig/",
                          "lineage_phate",
                          ".set-", this_set,
                          ".knn-", i_knn,
                          ".pdf",
                          sep = ""),
         plot = p,
         width = 2*2,
         height = 2*2,
         dpi = 300,
         device = "pdf")
  
  i=i+1
}

phate.all = rbindlist(phate)
phate.curves.all = rbindlist(phate.curves)
phate.ptvals.all = phate.ptvals[[2]]
phate.ptvals.all$dog = pca$IID
phate.ptvals.all = phate.ptvals.all %>%
  replace(is.na(.), 0) %>%
  select(dog,starts_with("Lineage"))
colnames(phate.ptvals.all) = c("dog",paste("lineage.int.",str_pad(as.character(seq(1,7)),width=2,side="left",pad="0")))
  

phate.all$dog = rep(x = pca$IID, times = length(c(40,60,80,100)))
phate.all$fid = rep(x = pca$`#FID`, times = length(c(40,60,80,100)))
phate.all$ref = rep(x = pca$reference, times = length(c(40,60,80,100)))



## set up plotting table ----
phate.plot = bind_rows((merge((dog %>%
                                 select(dog,purebred,breed)),
                              (adm %>%
                                 pivot_wider(id_cols = "dog",
                                             names_from = "pop",
                                             values_from = "pct",
                                             values_fill = 0)),
                              by = "dog") %>%
                          mutate(dog = as.character(dog))),
                       (pca %>% 
                          filter(reference == T) %>%
                          mutate(pct = 1,
                                 pop = `#FID`,
                                 dog = IID) %>%
                          select(dog,pop,pct) %>%
                          pivot_wider(id_cols = "dog",
                                      names_from = "pop",
                                      values_from = "pct",
                                      values_fill = 0) %>%
                          mutate(dog = as.character(dog)))) %>%
  merge(phate.all, by = "dog", all.x = T)

write_tsv(x = phate.all,
          file = paste(workDir,
                       args$pca,
                       ".phate-embeddings",
                       ".tsv",
                       sep = ""))

## plot results ----
p = ggplot(phate.all,
           aes(x = PHATE1,
               y = PHATE2)) +
  facet_wrap(.~knn,
             scales = "free",
             ncol = 4,
             nrow = 1) +
  geom_point(color = "#505050",
             alpha = 0.25) +
  geom_path(data = phate.curves.all %>% arrange(Order),
            aes(group = Lineage)) +
  theme_pubr()

ggsave(filename = paste(workDir,
                        "/fig/",
                        "lineage_phate",
                        ".set-", this_set,
                        ".pdf",
                        sep = ""),
       plot = p,
       width = 6*2,
       height = 2*2,
       dpi = 300,
       device = "pdf")

# INFER PSEUDOTIME VALUES ----
this_knn = 80

## cluster containing wild canids ----
cluster.starts = phate.all %>%
  filter(knn==this_knn) %>%
  filter(fid %like% "wolf_" | fid %like% "coyote" | fid %like% "village") %>%
  pull(cluster) %>%
  unique()

## perform trajectory inference with slingshot ----
phate.mat = phate.all %>% 
  filter(knn==this_knn) %>% 
  select(PHATE1,PHATE2) %>% 
  as.matrix()

phate.lab = phate.all %>% 
  filter(knn==this_knn) %>% 
  pull(cluster)

# with cluster starts...
phate.sling = slingshot(data = phate.mat,
                        clusterLabels = phate.lab,
                        start.clus = cluster.starts)

# without cluster starts...
phate.sling = slingshot(data = phate.mat,
                        clusterLabels = phate.lab)

## extract pseudotimes with slingshot ----
phate.pseudotimes = slingPseudotime(x = phate.sling)

## re-visualize lineages ----
curves <- slingCurves(SlingshotDataSet(phate.sling), as.df = TRUE)

p = phate.all %>%
  filter(knn==this_knn) %>% 
  bind_cols(phate.pseudotimes) %>%
  pivot_longer(cols = colnames(phate.pseudotimes),
               names_to = "lineage",
               values_to = "pseudotime") %>%
  arrange(-is.na(pseudotime)) %>%
  ggplot(aes(x = PHATE1,
             y = PHATE2)) +
  geom_point(aes(color = pseudotime),
             alpha = 0.25) +
  geom_path(data = curves %>% arrange(Order),
            aes(group = Lineage)) +
  facet_wrap(.~lineage) +
  theme_pubr()

ggsave(filename = paste(workDir,
                        "/fig/",
                        "lineage_pseudotimes",
                        ".knn-", this_knn,
                        ".set-", this_set,
                        ".pdf",
                        sep = ""),
       plot = p,
       device = "pdf",
       width = 12,
       height = 8,
       units = "in",
       dpi = 300)

curves <- slingCurves(SlingshotDataSet(phate.sling), as.df = TRUE)

p = phate.all %>%
  filter(knn==this_knn) %>% 
  bind_cols(phate.pseudotimes) %>%
  pivot_longer(cols = colnames(phate.pseudotimes),
               names_to = "lineage",
               values_to = "pseudotime") %>%
  filter(lineage == colnames(phate.pseudotimes)[1]) %>%
  ggplot(aes(x = PHATE1,
             y = PHATE2)) +
  geom_point(color = "#505050") +
  geom_path(data = curves %>% arrange(Order),
              aes(group = Lineage)) +
  scale_x_continuous(limits = c(-0.06,0.06)) +
  scale_y_continuous(limits = c(-0.06,0.06)) +
  theme_pubr()

ggsave(filename = paste(workDir,
                        "/fig/",
                        "lineage_pseudotimes",
                        ".knn-", this_knn,
                        ".set-", this_set,
                        ".base",
                        ".pdf",
                        sep = ""),
       plot = p,
       device = "pdf",
       width = 12,
       height = 12,
       units = "in",
       dpi = 300)

## extract as phenotypes ----
i=1
for (lineage in colnames(phate.pseudotimes)){
  
  p = phate.all %>%
    filter(knn==this_knn) %>% 
    bind_cols(phate.pseudotimes) %>%
    filter(ref == F) %>%
    select(dog,PHATE1,PHATE2,!!lineage) %>%
    pivot_longer(cols = lineage,
                 names_to = "lineage",
                 values_to = "pseudotime") %>%
    mutate(pseudotime = if_else(is.na(pseudotime),
                                0,
                                pseudotime)) %>%
    ggplot(aes(x = PHATE1,
               y = PHATE2,
               alpha = pseudotime)) +
    geom_point() +
    scale_x_continuous(limits = c(-0.06,0.06)) +
    scale_y_continuous(limits = c(-0.06,0.06)) +
    scale_alpha_continuous(guide="none") +
    theme_pubr()
  
  ggsave(filename = paste(workDir,
                          "/fig/",
                          "lineage_pseudotimes.",
                          "lineage.",
                          str_pad(as.character(i),width=2,side="left",pad="0"),
                          ".pdf",
                          sep = ""),
         plot = p,
         device = "pdf",
         width = 12,
         height = 12,
         units = "in",
         dpi = 300)
  
  write_tsv(x = (phate.all %>%
    filter(knn==this_knn) %>% 
    bind_cols(phate.pseudotimes) %>%
    filter(ref == F) %>%
    select(dog,!!lineage) %>%
    mutate(fid=dog,
           iid=dog,
           phe=get(lineage)) %>%
    mutate(phe = if_else(!is.na(phe),
                          phe,
                          0)) %>%
    select(fid,iid,phe)),
    file = paste(workDir,
                 "/gwa/phe/",
                 "lineage.",
                 str_pad(as.character(i),width=2,side="left",pad="0"),
                 ".tsv",
                 sep=""),
    col_names = F)
  i=i+1
}

lineage.all = phate.all %>%
  filter(knn==this_knn) %>% 
  bind_cols(phate.pseudotimes) %>%
  filter(ref==F) %>%
  select(dog,colnames(phate.pseudotimes)) %>%
  replace(is.na(.), 0)

colnames(lineage.all) = c("dog",paste("lineage.",str_pad(as.character(seq(1,12)),width=2,side="left",pad="0")))
write_tsv(x = lineage.all,
          file = paste(workDir,
                       "/gwa/phe/",
                       "lineage.all",
                       ".tsv",
                       sep=""))


write_tsv(x = phate.ptvals.all,
          file = paste(workDir,
                       "/gwa/phe/",
                       "lineage.int",
                       ".tsv",
                       sep=""))

lineage.all = read_tsv(paste(workDir,
                             "/gwa/phe/",
                             "lineage.all",
                             ".tsv",
                             sep=""))
colnames(lineage.all) = c("dog",paste("lineage.ext.",str_pad(as.character(seq(1,12)),width=2,side="left",pad="0")))
write_tsv(x = lineage.all,
          file = paste(workDir,
                       "/gwa/phe/",
                       "lineage.ext",
                       ".tsv",
                       sep=""))


# RELOAD ----

phate.all = read_tsv(paste(workDir,
                       args$pca,
                       ".phate-embeddings",
                       ".tsv",
                       sep = ""))

this_knn = 80

cluster.starts = phate.all %>%
  filter(knn==this_knn) %>%
  filter(fid %like% "wolf_" | fid %like% "coyote" | fid %like% "village") %>%
  pull(cluster) %>%
  unique()

phate.mat = phate.all %>% 
  filter(knn==this_knn) %>% 
  select(PHATE1,PHATE2) %>% 
  as.matrix()

phate.lab = phate.all %>% 
  filter(knn==this_knn) %>% 
  pull(cluster)

phate.sling = slingshot(data = phate.mat,
                        clusterLabels = phate.lab,
                        start.clus = cluster.starts)

phate.pseudotimes = slingPseudotime(x = phate.sling)

curves <- slingCurves(SlingshotDataSet(phate.sling), as.df = TRUE)

phate.coi = phate.all %>%
  filter(knn==this_knn) %>% 
  bind_cols(phate.pseudotimes) %>%
  merge(coi, by = "dog")

for (i in paste("Lineage",1:12,sep="")){
  print(i)
  print(cor.test(get(i,phate.coi), phate.coi$coi))
}

phate.all %>%
  filter(knn==this_knn) %>% 
  bind_cols(phate.pseudotimes) %>%
  merge(coi, by = "dog") %>%
  pivot_longer(cols = colnames(phate.pseudotimes),
               names_to = "lineage",
               values_to = "pseudotime") %>%
  mutate(lineage = as.integer(gsub("Lineage","",lineage))) %>% filter(lineage==1) %>% arrange(-pseudotime) %>% select(dog,purebred,breed,greyhound,coi,pseudotime,lineage) %>% filter(!is.na(pseudotime)) %>% View()

p = phate.all %>%
  filter(knn==this_knn) %>% 
  bind_cols(phate.pseudotimes) %>%
  merge(coi, by = "dog") %>%
  pivot_longer(cols = colnames(phate.pseudotimes),
               names_to = "lineage",
               values_to = "pseudotime") %>%
  mutate(lineage = as.integer(gsub("Lineage","",lineage))) %>%
  ggplot(aes(x = pseudotime,
             y = coi)) +
  geom_point(alpha = 0.5) +
  stat_cor(p.accuracy = 1e-200) +
  facet_wrap(.~lineage, scales = "free") +
  theme_pubr()
  
ggsave(filename = paste(workDir,
                        "/fig/",
                        "lineage_pseudotimes",
                        ".knn-", this_knn,
                        ".set-", this_set,
                        ".vs-coi",
                        ".pdf",
                        sep = ""),
       plot = p,
       device = "pdf",
       width = 6*2,
       height = 4*2,
       units = "in",
       dpi = 300)


p = phate.all %>%
  filter(knn==this_knn) %>% 
  bind_cols(phate.pseudotimes) %>%
  merge(coi, by = "dog") %>%
  pivot_longer(cols = colnames(phate.pseudotimes),
               names_to = "lineage",
               values_to = "pseudotime") %>%
  mutate(lineage = as.integer(gsub("Lineage","",lineage))) %>%
  filter(lineage == 6) %>%
  arrange(coi) %>%
  ggplot(aes(x = PHATE1,
             y = PHATE2)) +
  geom_point(aes(color = pseudotime)) +
  geom_path(data = curves %>% arrange(Order),
            aes(group = Lineage)) +
  scale_x_continuous(limits = c(-0.06,0.06)) +
  scale_y_continuous(limits = c(-0.06,0.06)) +
  theme_pubr()

ggsave(filename = paste(workDir,
                        "/fig/",
                        "lineage_pseudotimes",
                        ".knn-", this_knn,
                        ".set-", this_set,
                        ".base",
                        ".pdf",
                        sep = ""),
       plot = p,
       device = "pdf",
       width = 12,
       height = 12,
       units = "in",
       dpi = 300)