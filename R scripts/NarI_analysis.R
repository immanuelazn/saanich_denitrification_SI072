# Load required packages 
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggpubr)
library(grid)

# Load metagenomic and metatranscriptomic data 
NarI_ts_dat <- read.table(file="SI072_NarI_layered_classifications.tsv",
                          header=TRUE,
                          sep="\t")

# Filter for only NarI
# (not really necessary bc it's already all NorB)
NarI_ts_df <- NarI_ts_dat %>%
  filter(Marker == "NarI")

# Process taxonomic information
taxa_ranks <- c("Root", "Domain", "Phylum", "Class",
                "Order", "Family", "Genus", "Species")

NarI_ts_df <- NarI_ts_df %>%
  separate(col = Taxonomy,
           into = taxa_ranks,
           sep = "; ", fill = "right", remove = T)

# Distinguishing metagenomes from metatranscriptomes
NarI_ts_df <- NarI_ts_df %>% 
  mutate(Sample = gsub("_pe", "_MetaG", Sample)) %>% 
  separate(col = Sample,
           into = c("Cruise", "Depth", "SeqType"),
           extra = "drop")


# Processing depth data 
NarI_ts_df <- NarI_ts_df %>%
  mutate(Depth.m = as.numeric(gsub('m', '', Depth))) %>% 
  mutate(Cruise = as.numeric(gsub('SI0', '', Cruise)))

# Load the geochemical data 
geochem_df <- read.table("Saanich_TimeSeries_Chemical_DATA.csv",
                         header=TRUE, sep=',') %>% 
  filter(Cruise == 72) %>% 
  select(Cruise, Depth,
         CTD_O2, NO3,
         Mean_NH4, Mean_NO2, Mean_H2S, Mean_CH4, Mean_O2, Mean_N2O) %>% 
  mutate(Mean_CH4 = Mean_CH4*1E-3) %>% 
  mutate(across(.fns=as.numeric))

# Merge the geochemical data for cruise 72 with the TreeSAPP 
# classifications for NorB by Sample

NarI_ts_geo_df <- left_join(NarI_ts_df,
                            geochem_df,
                            by=c("Cruise",
                                 "Depth.m" = "Depth"))

# ----------------- Data visualization --------------------

# Stacked barplot to show the relative change in proportions at 
# different depths

NarI_ts_geo_df %>% 
  group_by(Depth, Depth.m, SeqType) %>% 
  mutate(Proportion = Abundance/sum(Abundance)) %>% 
  group_by(Depth, Depth.m, Phylum, SeqType) %>% 
  summarise(sum = sum(Proportion)) %>% 
  ungroup() %>% 
  mutate(Depth = reorder(Depth, Depth.m)) %>% 
  ggplot(aes(x=Depth, y=sum, fill=Phylum)) +
  geom_bar(stat = "identity") +
  facet_wrap(~SeqType) +
  labs(x="Depth",
       y="Relative abundance") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(legend.text=element_text(size=9)) +
  ggtitle("Abundance of NarI-containing Phyla (metagenome)
          by depth in the Saanich Inlet")

ggsave("NarI_taxonomic_phyla_barplot.png")








