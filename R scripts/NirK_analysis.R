# NirK R analysis 

# -------------- Manipulate data files -------------

# Set working directory
setwd("/Users/KristiMacBookAir2015/Desktop/saanich/capstone")

# Load required packages 
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggpubr)

# Load metagenomic and metatranscriptomic data 
NirK_ts_dat <- read.table(file="SI072_NirK_layered_classifications.tsv",
                          header=TRUE,
                          sep="\t")

# Filter for only Nirk
# (not really necessary bc it's already all NirK)
NirK_ts_df <- NirK_ts_dat %>%
  filter(Marker == "NirK")

# Process taxonomic information
taxa_ranks <- c("Root", "Domain", "Phylum", "Class",
                "Order", "Family", "Genus", "Species")

NirK_ts_df <- NirK_ts_df %>%
  separate(col = Taxonomy,
           into = taxa_ranks,
           sep = "; ", fill = "right", remove = T)

# Distinguishing metagenomes from metatranscriptomes
NirK_ts_df <- NirK_ts_df %>% 
  mutate(Sample = gsub("_pe", "_MetaG", Sample)) %>% 
  separate(col = Sample,
           into = c("Cruise", "Depth", "SeqType"),
           extra = "drop")


# Processing depth data 
NirK_ts_df <- NirK_ts_df %>%
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

NirK_ts_geo_df <- left_join(NirK_ts_df,
                       geochem_df,
                       by=c("Cruise",
                            "Depth.m" = "Depth"))

# -------------- Data visualization: ORDER level ---------------

# Line plot to show the distribution of abundances of taxonomic orders 
# at different depths
NirK_ts_geo_df %>%
  mutate(across(Domain:Species, gsub, pattern = ".__", replacement = "")) %>% 
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
  ggtitle("NirK-expressing Orders (metatranscriptome) 
          by depth in the Saanich Inlet")

ggsave("NirK_taxonomic_order_lineplot.png")


# Stacked barplot to show the relative change in proportions at 
# different depths
NirK_ts_geo_df %>% 
  mutate(across(Domain:Species, gsub, pattern = ".__", replacement = "")) %>% 
  group_by(Depth, Depth.m, SeqType) %>% 
  mutate(Proportion = Abundance/sum(Abundance)) %>% 
  group_by(Depth, Depth.m, Order, SeqType) %>% 
  summarise(sum = sum(Proportion)) %>% 
  ungroup() %>% 
  mutate(Depth = reorder(Depth, Depth.m)) %>% 
  ggplot(aes(x=Depth, y=sum, fill=Order)) +
  geom_bar(stat = "identity") +
  facet_wrap(~SeqType) +
  labs(x="Depth",
       y="Relative abundance") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(legend.text=element_text(size=9)) +
  ggtitle("Abundance of NirK-containing Orders (metagenome)
          by depth in the Saanich Inlet")
  
ggsave("NirK_taxonomic_order_barplot.png")


# Bubble-plot where the bubbles at each depth-order intersection are 
# scaled by the TPM value (Abundance). 
NirK_ts_geo_df %>% 
  mutate(across(Domain:Species, gsub, pattern = ".__", replacement = "")) %>% 
  group_by(Depth, Depth.m, SeqType, Order) %>% 
  summarise(Sum = sum(Abundance)) %>% 
  ungroup() %>%
  mutate(Depth = reorder(Depth, desc(Depth.m))) %>% 
  ggplot(aes(x=Depth, y=Order)) +
  geom_point(aes(size=Sum)) +
  facet_wrap(~SeqType) +
  coord_flip() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("NirK_functional_bubbleplot.png")


# Line plot to visualize chemical profiles for each depth
NirK_ts_geo_df %>% 
  pivot_longer(cols=c(starts_with("Mean"), NO3, CTD_O2),
               values_to = "Value.uM",
               names_to = "Molecule") %>%
  ggplot(aes(x=Depth.m, y=Value.uM)) +
  geom_point() +
  geom_line() +
  facet_wrap(~Molecule, scales = "free_x") +
  coord_flip() +
  scale_x_reverse()

ggsave("NirK_geochemical_lineplot.png")


# Combine the NO2 concentration data with the TPM values for Class.
# Colour by mean NO2 concentration because NirK catalyzes NO2 -> NO.
NirK_ts_geo_df %>%
  mutate(across(Domain:Species, gsub, pattern = ".__", replacement = "")) %>% 
  filter(SeqType=="MetaT") %>%
  filter(!is.na(Class)) %>% 
  mutate(Depth = reorder(Depth, Depth.m)) %>% 
  ggplot(aes(x=Depth, y=Abundance, colour=Mean_NO2)) +
  geom_boxplot() +
  geom_point() +
  facet_wrap(~Order) +
  ggtitle("NirK by depth and mean NO2 (too many Orders)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("NirK_taxonomic_order_boxplot.png")

# Too many Classes with 0 abundance. Filter to only include  
# Orders Nitrosomonadales, Nitrososphaerales, and SAR324.

NirK_select_Order <- NirK_ts_geo_df %>% 
  filter(Order %in% c("Nitrososphaerales",
                      "Nitrosomonadales",
                      "SAR324"))

NirK_select_Order %>%
  mutate(across(Domain:Species, gsub, pattern = ".__", replacement = "")) %>% 
  filter(SeqType=="MetaT") %>%
  filter(!is.na(Class)) %>% 
  mutate(Depth = reorder(Depth, Depth.m)) %>% 
  ggplot(aes(x=Depth, y=Abundance, colour=Mean_NO2)) +
  geom_boxplot() +
  geom_point() +
  facet_wrap(~Order) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("Orders expressing (metatranscriptome) highest 
          levels of NirK, by mean NO2 and depth")

ggsave("NirK_taxonomic_order_boxplot_select_Orders.png")


# Literature search says that NirK-expressing microbes are involved
# in ammonia oxidation.
# Create the same plots, but this time colour by mean NH3/NH4 concentration. 

NirK_select_Order %>%
  mutate(across(Domain:Species, gsub, pattern = ".__", replacement = "")) %>% 
  filter(SeqType=="MetaT") %>%
  filter(!is.na(Class)) %>% 
  mutate(Depth = reorder(Depth, Depth.m)) %>% 
  ggplot(aes(x=Depth, y=Abundance, colour=Mean_NH4)) +
  geom_boxplot() +
  geom_point() +
  facet_wrap(~Order) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("Orders expressing (metatranscriptome) highest 
          levels of NirK, by mean NH4 and depth")


# -------------- Data visualization: PHYLUM level ---------------

# Line plot to show the distribution of abundances of taxonomic orders 
# at different depths
NirK_ts_geo_df %>%
  mutate(across(Domain:Species, gsub, pattern = ".__", replacement = "")) %>% 
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
  ggtitle("NirK-expressing Phyla (metatranscriptome) 
          by depth in the Saanich Inlet")

ggsave("NirK_taxonomic_order_lineplot.png")


# Stacked barplot to show the relative change in proportions at 
# different depths
NirK_ts_geo_df %>% 
  mutate(across(Domain:Species, gsub, pattern = ".__", replacement = "")) %>% 
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
  ggtitle("Abundance of Phyla containing and expressing NirK
          by depth in the Saanich Inlet")

ggsave("NirK_taxonomic_order_barplot.png")


# Bubble-plot where the bubbles at each depth-order intersection are 
# scaled by the TPM value (Abundance). 
NirK_ts_geo_df %>% 
  mutate(across(Domain:Species, gsub, pattern = ".__", replacement = "")) %>% 
  group_by(Depth, Depth.m, SeqType, Phylum) %>% 
  summarise(Sum = sum(Abundance)) %>% 
  ungroup() %>%
  mutate(Depth = reorder(Depth, desc(Depth.m))) %>% 
  ggplot(aes(x=Depth, y=Phylum)) +
  geom_point(aes(size=Sum)) +
  facet_wrap(~SeqType) +
  coord_flip() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("NirK_functional_bubbleplot.png")


# Line plot to visualize chemical profiles for each depth
NirK_ts_geo_df %>% 
  pivot_longer(cols=c(starts_with("Mean"), NO3, CTD_O2),
               values_to = "Value.uM",
               names_to = "Molecule") %>%
  ggplot(aes(x=Depth.m, y=Value.uM)) +
  geom_point() +
  geom_line() +
  facet_wrap(~Molecule, scales = "free_x") +
  coord_flip() +
  scale_x_reverse()

ggsave("NirK_geochemical_lineplot.png")


# Combine the NO2 concentration data with the TPM values for Phylum.
# Colour by mean NO2 concentration because NirK catalyzes NO2 -> NO.
NirK_ts_geo_df %>%
  mutate(across(Domain:Species, gsub, pattern = ".__", replacement = "")) %>% 
  filter(SeqType=="MetaT") %>%
  filter(!is.na(Phylum)) %>% 
  mutate(Depth = reorder(Depth, Depth.m)) %>% 
  ggplot(aes(x=Depth, y=Abundance, colour=Mean_NO2)) +
  geom_boxplot() +
  geom_point() +
  facet_wrap(~Phylum) +
  ggtitle("NirK-expressing Phyla by depth and mean NO2") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("NirK_taxonomic_order_boxplot.png")

# Too many Phyla with near 0 abundance. Filter to only include  
# Phyla ______.

NirK_select_Phylum <- NirK_ts_geo_df %>% 
  filter(Phylum %in% c("p__Proteobacteria",
                      "p__Thaumarchaeota",
                      "p__SAR324"))

NirK_select_Phylum %>%
  mutate(across(Domain:Species, gsub, pattern = ".__", replacement = "")) %>% 
  filter(SeqType=="MetaT") %>%
  filter(!is.na(Phylum)) %>% 
  mutate(Depth = reorder(Depth, Depth.m)) %>% 
  ggplot(aes(x=Depth, y=Abundance, colour=Mean_NO2)) +
  geom_boxplot() +
  geom_point() +
  facet_wrap(~Phylum) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("NirK-expressing Phyla by depth and mean NO2")

ggsave("NirK_taxonomic_order_boxplot_select_Orders.png")


# Literature search says that NirK-expressing microbes are involved
# in ammonia oxidation.
# Create the same plots, but this time colour by mean NH3/NH4 concentration. 

NirK_select_Phylum %>%
  mutate(across(Domain:Species, gsub, pattern = ".__", replacement = "")) %>% 
  filter(SeqType=="MetaT") %>%
  filter(!is.na(Phylum)) %>% 
  mutate(Depth = reorder(Depth, Depth.m)) %>% 
  ggplot(aes(x=Depth, y=Abundance, colour=Mean_NH4)) +
  geom_boxplot() +
  geom_point() +
  facet_wrap(~Phylum) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("NirK-expressing Phyla by depth and mean NH4")


