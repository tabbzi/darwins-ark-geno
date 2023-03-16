library(stringr)
# Set the directory path
dir_path <- "~/Dropbox (UMass Medical School)/Thesis/repo/dat/gsa/"

# Get a list of files in the directory
file_list <- list.files(path = dir_path, full.names = TRUE)

# Test list
file_list = file_list[1:10]

# Read in all files in the list using map_dfr()
read_table_set = function(x){ read_table(file = x,
                                     comment = "#") }

data_list <- map_dfr(file_list, 
                     read_table_set, 
                     .id = "file")

# Get file names
data_list = data_list %>% 
  merge((data.frame(filename=file_list) %>% 
           mutate(file=row_number())), by= "file")

# Split file name for info
data_list = data_list %>%
  mutate(filename = basename(filename)) %>%
  separate(filename,
           into = c("cohort",
                    "gp",
                    "var",
                    "maf",
                    "geno",
                    "hwe",
                    "samp",
                    "phe",
                    "dcov",
                    "qcov",
                    "extra"),
           sep = "_",
           remove = T, extra = "merge") %>%
  mutate(set = str_extract(extra, "(genes\\..+)")) %>%
  mutate(Padj = p.adjust(P, method = "BH")) %>%
  mutate(trait = if_else(phe %like% "bq.",
                         "survey item",
                         if_else(phe %like% "fa.",
                                 "factor",
                                 if_else(phe %like% "mq." | phe %like% "mp.",
                                         "physical trait",
                                         "other")))) %>%
  mutate(subset = if_else(!is.na(FULL_NAME),
                          FULL_NAME,
                          VARIABLE))

unique(data_list$set)

data = as.data.frame(data_list) %>%
  filter(!phe %like% "lineage") %>%
  filter(!set %in% c("genes.mousebraincell.top500.set.gsa.out",
                     "genes.brainspan.genes.top500.set.gsa.out",
                     "genes.canine-pseudo-phenotypes.set.gsa.out",
                     "genes.canine-pseudo-phenotypes.gsa.out",
                     "genes.gtex.gsa.out")) %>%
  group_by(set) %>%
  mutate(Padj = p.adjust(P, method = "BH")) %>%
  filter(P<0.1)

# Plot!
p = ggplot(data = data,
           aes(x = phe,
               y = subset)) +
  geom_point(data = (data %>% filter(Padj>=0.05)),
             aes(size = -log10(P),
                 alpha = -log10(P)),
             color = "#525252",
             shape = 16) +
  geom_point(data = (data %>% filter(Padj<0.05)),
             aes(size=-log10(P)),
             color="#cb181d",
             shape=16) +
  geom_point(data = (data %>% filter(P<=0.005)),
             aes(size=-log10(P)),
             color="black",
             shape=21) +
  facet_grid(set~trait,
             scales="free",
             space="free") +
  scale_y_discrete(limits=rev) +
  scale_alpha_continuous(range=c(0.1,0.75)) +
  scale_size(range=c(1,3)) +
  scale_fill_manual(values=c("black","red","black")) +
  scale_color_manual(values=c("black","red","black")) +
  theme_minimal() +
  theme(legend.position = "bottom",
        panel.grid.minor=element_blank(),
        plot.title = element_text(size=7, 
                                  face="bold"),
        plot.subtitle = element_text(size=6),
        axis.text.y = element_text(size=6),
        axis.text.x = element_text(angle=90,
                                   size=6,
                                   hjust=1,
                                   vjust=0.5),
        axis.title.y=element_blank(),
        axis.title.x=element_blank(),
        legend.title = element_blank(),
        legend.text=element_text(hjust=0.5,
                                 size=6),
        legend.key.size = unit(0.1, 'in'),
        strip.text.y=element_text(size=5,
                                  vjust=0,
                                  hjust=0.5),
        strip.text.x=element_text(size=5,
                                  vjust=0,
                                  hjust=0.5))

ggsave(plot=p,
       filename="~/Dropbox (UMass Medical School)/Thesis/figures/MAGMA.all_circles.all-sets.pdf",
       width=6*2,
       height=5*2) 


data %>%
  filter(P<0.005) %>%
  ggplot(aes(x = phe,
             y = BETA,
             fill = Padj<0.05)) +
  geom_bar(stat = "identity") +
  geom_errorbar(aes(ymin= BETA-BETA_STD,
                    ymax= BETA+BETA_STD)) +
  coord_flip() +
  facet_wrap(.~subset)
