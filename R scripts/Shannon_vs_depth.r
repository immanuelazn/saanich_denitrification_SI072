#Absolute abdunace across depth for each gene

#Note that the number of genes classified by TreeSAPP was used as an estimate 
#of absolute abundance 

library(dplyr)
library(ggplot2)
library(tidyverse)
# Process taxonomic information
taxa_ranks <- c("Root", "Domain", "Phylum", "Class",
                "Order", "Family", "Genus", "Species")

# Load the geochemical data 
geochem_df <- read.table("Saanich_TimeSeries_Chemical_DATA.csv",
                         header=TRUE, sep=',') %>% 
  filter(Cruise == 72) %>% 
  select(Cruise, Depth,
         CTD_O2, NO3,
         Mean_NH4, Mean_NO2, Mean_H2S, Mean_CH4, Mean_O2, Mean_N2O) %>% 
  mutate(Mean_CH4 = Mean_CH4*1E-3) %>% 
  mutate(across(.fns=as.numeric))

# functions ---------------------------------------------------------------

format_data <- function(df) {
  
  df <- df %>%
    separate(col = Taxonomy,
             into = taxa_ranks,
             sep = "; ", fill = "right", remove = T)
  
  # Distinguishing metagenomes from metatranscriptomes
  df <- df %>% 
    mutate(Sample = gsub("_pe", "_MetaG", Sample)) %>% 
    separate(col = Sample,
             into = c("Cruise", "Depth", "SeqType"),
             extra = "drop")
  
  # Processing depth data 
  df <- df %>%
    mutate(Depth.m = as.numeric(gsub('m', '', Depth))) %>% 
    mutate(Cruise = as.numeric(gsub('SI0', '', Cruise)))
  
  # Merge the geochemical data for cruise 72 with the TreeSAPP 
  # classifications for NirK by Sample
  
  df <- left_join(df,
                  geochem_df,
                  by=c("Cruise",
                  "Depth.m" = "Depth"))
  
  
  df <- df %>%
    mutate(across(Domain:Species, gsub, pattern = ".__", replacement = ""))
}


# Shannon diversity given a treesapp gene df
shannon_diversity <- function(df_gene, n, g) {
  test_diversity <- df_gene %>% filter(Abundance != 0) %>%
    group_by(Depth, SeqType) %>%
    summarise(total_N = sum(Abundance),
              H = sum((Abundance / total_N)*-log(Abundance / total_N))) %>% 
    add_column(Name = n) %>% add_column(G = g) %>%
    select(Name, Depth, SeqType, H, G) %>% filter(H > 0)
  return(test_diversity) 
}
# load data ---------------------------------------------------------------

#Nosz
nosz <- read.table(file="SI072_NosZ_layered_classifications.tsv",
                          header=TRUE,
                          sep="\t")
nosz <- format_data(nosz)
#NorB
norb <- read.table(file="SI072_NorB_layered_classifications.tsv",
                  header=TRUE,
                  sep="\t")
norb <- format_data(norb)

#NarI
nari <- read.table(file="SI072_NarI_layered_classifications.tsv",
                  header=TRUE,
                  sep="\t")
nari <- format_data(nari)


#NirK
nirk <- read.table(file="SI072_NirK_layered_classifications.tsv",
                  header=TRUE,
                  sep="\t")
nirk <- format_data(nirk)

#NapA
napa <- read.table(file="SI072_NapA_layered_classifications.tsv",
                  header=TRUE,
                  sep="\t")
napa <- format_data(napa)

# # Calculate diversity ----------------------------------------------------
# http://vm-trypanocyc.toulouse.inra.fr/META/NEW-IMAGE?type=REACTION&object=RXN-12130
shannon_nosz <- shannon_diversity(nosz, "nosZ", 34.26947) 
# https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=NITRIC-OXIDE-REDUCTASE-RXN
shannon_norb <- shannon_diversity(norb, "norB", 8.969482)
shannon_nari <- shannon_diversity(nari, "narI", -69.334)
# https://www.sciencedirect.com/science/article/pii/S0005272806001046
shannon_nirk <- shannon_diversity(nirk, "nirK", -45.3)
# http://vm-trypanocyc.toulouse.inra.fr/META/NEW-IMAGE?type=REACTION&object=NITRATE-REDUCTASE-CYTOCHROME-RXN&detail-level=2&orgids=LEISH
shannon_napa <- shannon_diversity(napa, "napA", -74.12053)
total_data <- rbind(shannon_nosz, shannon_norb, shannon_nari, shannon_nirk, shannon_napa)

# # Model diversity ----------------------------------------------------
# Remove outliers and subset data
total_data_test <- total_data 
metaG_data <- filter(total_data, SeqType == "MetaG")
metaT_data <- filter(total_data, SeqType == "MetaT") %>% filter(H > 0.2)

summary(lm(data=total_data_test, H ~ G*SeqType + I(G^2)))

# diagnostic plots 
plot(lm(data=total_data_test, H ~ G*SeqType + I(G^2)))


# Plot diversity across MAGs and Metatranscriptomic data -----------------
ggplot(data = metaG_data, aes(x=G, y=H)) +
  geom_point() + theme_classic() + xlab("Standard Delta G") +
  ggtitle("Standard Delta G vs Shannon Diversity in Metagenomic Data") +
  stat_smooth(method = 'lm', formula = y ~ poly(x, 2), se = TRUE)

ggplot(data = metaT_data, aes(x=G, y=H)) +
  geom_point() + theme_classic() + xlab("Standard Delta G") +
  ggtitle("Standard Delta G vs Shannon Diversity in Metatranscriptomic Data") +
  stat_smooth(method = 'lm', formula = y ~ poly(x, 2), se = TRUE)


# Shannon vs Depth --------------------------------------------------------

total_data_test$Depth = str_sub(total_data_test$Depth, 1,
                                   nchar(total_data_test$Depth)-1)

total_data_test <- total_data_test %>% 
  mutate(Depth = as.numeric(Depth))

#Plotting Shannon Index vs Depth

ggplot(data = metaG_data, aes(x=Depth, y=H)) +
  geom_line(linetype = "dashed", aes(group = 1))+
  geom_point()+
  theme_classic() + xlab("Depth") + ylab("Shannon Diversity Index")+
  ggtitle("Depth vs Shannon Diversity in Metagenomic Data") +
  facet_wrap(~Name)+

ggplot(data = metaT_data, aes(x=Depth, y=H, color = Name)) +
  geom_line(linetype = "dashed", aes(group = 1))+
  geom_point()+
  theme_classic() + xlab("Depth") + ylab("Shannon Diversity Index")+
  ggtitle("Depth vs Shannon Diversity in Metatranscriptomic Data") +
  facet_wrap(~Name)

#Overlap of MetaT and MetaG data
ggplot(data = total_data_test, aes(x=Depth, y=H, colour = factor(SeqType)))+
  geom_line(linetype = "dashed")+
  geom_point()+
  labs(colour = "Omic Type")+
  xlab("Depth") + ylab("Shannon Diversity Index")+
  ggtitle("Depth vs Shannon Diversity") +
  facet_wrap(~Name)

ggsave("depth_vs_shannon.png",
       width = 2950,
       height = 1700,
       units = "px")

