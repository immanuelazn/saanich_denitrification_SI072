# NosZ R analysis 

# -------------- Manipulate data files -------------

# Set working directory

# Load required packages 
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggpubr)

# Load metagenomic and metatranscriptomic data 
NosZ_ts_dat <- read.table(file="SI072_NosZ_layered_classifications.tsv",
                          header=TRUE,
                          sep="\t")

# Process taxonomic information
taxa_ranks <- c("Root", "Domain", "Phylum", "Class",
                "Order", "Family", "Genus", "Species")

NosZ_ts_df <- NosZ_ts_dat %>%
  separate(col = Taxonomy,
           into = taxa_ranks,
           sep = "; ", fill = "right", remove = T)

# Distinguishing metagenomes from metatranscriptomes
NosZ_ts_df <- NosZ_ts_df %>% 
  mutate(Sample = gsub("_pe", "_MetaG", Sample)) %>% 
  separate(col = Sample,
           into = c("Cruise", "Depth", "SeqType"),
           extra = "drop")


# Processing depth data 
NosZ_ts_df <- NosZ_ts_df %>%
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
# classifications for NirK by Sample

NosZ_ts_geo_df <- left_join(NosZ_ts_df,
                       geochem_df,
                       by=c("Cruise",
                            "Depth.m" = "Depth"))


NosZ_ts_geo_df <- NosZ_ts_geo_df %>%
  mutate(across(Domain:Species, gsub, pattern = ".__", replacement = ""))


# Investigating proportions --------------------------------------------

proportions_df <- NosZ_ts_geo_df %>% 
  group_by(Depth, Phylum,SeqType) %>% 
  summarise(ab_sum = sum(Abundance))

proportions_df <- proportions_df %>% 
  group_by(Depth,SeqType) %>%
  mutate(rel_ab = ab_sum/sum(ab_sum))

max_proportions_df <- proportions_df %>% 
  group_by(Phylum,SeqType) %>% 
  summarise(max(rel_ab))

#N2O
NosZ_ts_geo_df %>% 
  select(Depth, Mean_N2O) %>%
  group_by(Depth) %>% 
  summarize(max(Mean_N2O))


# Order 
proportions_df_o <- NosZ_ts_geo_df %>% 
  group_by(Depth, Family,SeqType) %>% 
  summarise(ab_sum = sum(Abundance))

proportions_df_o <- proportions_df_o %>% 
  group_by(Depth,SeqType) %>%
  mutate(rel_ab = ab_sum/sum(ab_sum))

# ----------------- Data visualization --------------------

# Line plot to show the distribution of abundances of taxonomic orders 
# at different depths
NosZ_ts_geo_df %>%
  filter(SeqType == "MetaG") %>% 
  group_by(Depth.m, Phylum) %>%
  summarise(Sum = sum(Abundance)) %>%
  ungroup() %>%
  ggplot(aes(x=Depth.m, y=Sum, colour=Phylum)) +
  geom_line() +
  coord_flip() +
  scale_x_reverse() +
  labs(x="Depth (m)",
       y="Relative abundance (TPM)") +
  theme(legend.text=element_text(size=9)) +
  ggtitle("NosZ-expressing Phyla (metagenome) 
          by depth in the Saanich Inlet")

#ggsave("NosZ_taxonomic_order_lineplot.png")

ggsave("NosZ_taxonomic_order_lineplot_Phylum_mg.png",
       width = 2950,
       height = 1700,
       units = "px")

NosZ_ts_geo_df %>%
  filter(SeqType == "MetaT") %>% 
  group_by(Depth.m, Phylum) %>%
  summarise(Sum = sum(Abundance)) %>%
  ungroup() %>%
  ggplot(aes(x=Depth.m, y=Sum, colour=Phylum)) +
  geom_line() +
  coord_flip() +
  scale_x_reverse() +
  labs(x="Depth (m)",
       y="Relative abundance (TPM)") +
  theme(legend.text=element_text(size=9)) +
  ggtitle("NosZ-expressing Phyla (metatranscriptome) 
          by depth in the Saanich Inlet")

ggsave("NosZ_taxonomic_order_lineplot_Phylum_MT.png",
       width = 2950,
       height = 1700,
       units = "px")

NosZ_ts_geo_df %>%
  filter(SeqType == "MetaT") %>% 
  group_by(Depth.m, Genus) %>%
  summarise(Sum = sum(Abundance)) %>%
  ungroup() %>%
  ggplot(aes(x=Depth.m, y=Sum, colour=Genus)) +
  geom_line() +
  coord_flip() +
  scale_x_reverse() +
  labs(x="Depth (m)",
       y="Relative abundance (TPM)") +
  theme(legend.text=element_text(size=9)) +
  ggtitle("NosZ-expressing Orders (metatranscriptome) 
          by depth in the Saanich Inlet")

ggsave("NosZ_taxonomic_order_lineplot_MT.png",
       width = 2950,
       height = 1700,
       units = "px")

NosZ_ts_geo_df %>%
  filter(SeqType == "MetaG") %>% 
  group_by(Depth.m, Order) %>%
  summarise(Sum = sum(Abundance)) %>%
  ungroup() %>%
  ggplot(aes(x=Depth.m, y=Sum, colour=Order)) +
  geom_line() +
  coord_flip() +
  scale_x_reverse() +
  labs(x="Depth (m)",
       y="Relative abundance (TPM)") +
  theme(legend.text=element_text(size=9)) +
  ggtitle("NosZ-expressing Orders (metagenome) 
          by depth in the Saanich Inlet")


# Stacked barplot to show the relative change in proportions at 
# different depths
NosZ_ts_geo_df %>% 
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
  ggtitle("Abundance of NosZ-containing Phylum (metagenome)
          by depth in the Saanich Inlet")
  
ggsave("NosZ_taxonomic_order_barplot_phylum.png",
       width = 2950,
       height = 1700,
       units = "px")


# Bubble-plot where the bubbles at each depth-order intersection are 
# scaled by the TPM value (Abundance). 
NosZ_ts_geo_df %>% 
  group_by(Depth, Depth.m, SeqType, Order) %>% 
  summarise(Sum = sum(Abundance)) %>% 
  ungroup() %>%
  mutate(Depth = reorder(Depth, desc(Depth.m))) %>% 
  ggplot(aes(x=Depth, y=Order)) +
  geom_point(aes(size=Sum)) +
  facet_wrap(~SeqType) +
  coord_flip() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("NosZ_functional_bubbleplot.png")


# Line plot to visualize chemical profiles for each depth
NosZ_ts_geo_df %>% 
  pivot_longer(cols=c(starts_with("Mean"), NO3, CTD_O2),
               values_to = "Value.uM",
               names_to = "Molecule") %>%
  ggplot(aes(x=Depth.m, y=Value.uM)) +
  geom_point() +
  geom_line() +
  facet_wrap(~Molecule, scales = "free_x") +
  coord_flip() +
  scale_x_reverse()

ggsave("NosZ_geochemical_lineplot.png")





