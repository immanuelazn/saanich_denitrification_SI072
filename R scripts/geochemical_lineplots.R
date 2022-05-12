# Geochemical profile line plot

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
         Mean_NH4, Mean_NO2, Mean_H2S, Mean_CH4, Mean_N2O) %>% 
  mutate(Mean_CH4 = Mean_CH4*1E-3) %>% 
  mutate(across(.fns=as.numeric))

# Merge the geochemical data for cruise 72 with the TreeSAPP 
# classifications for NirK by Sample

NirK_ts_geo_df <- left_join(NirK_ts_df,
                            geochem_df,
                            by=c("Cruise",
                                 "Depth.m" = "Depth"))


# Faceted line plot to visualize chemical profiles for each depth

geochem_facet_labels <- as_labeller(c("CTD_O2" = "O2", "Mean_H2S" = "H2S", 
                                      "Mean_N2O" = "N2O", "Mean_NH4" = "NH4",
                                      "Mean_NO2" = "NO2", "NO3" = "NO3"))

NirK_ts_geo_df %>% 
  pivot_longer(cols=c(CTD_O2, Mean_H2S, Mean_N2O, Mean_NH4, Mean_NO2, NO3),
               values_to = "Value.uM",
               names_to = "Molecule") %>%
  ggplot(aes(x=Depth.m, y=Value.uM)) +
  geom_point() +
  geom_line() +
  facet_wrap(~Molecule, scales = "free_x", 
             labeller = geochem_facet_labels) +
  coord_flip() +
  scale_x_reverse() +
  xlab("Depth (m)") +
  ylab("Concentration (µM)")

ggsave("NirK_geochemical_lineplot.png")

