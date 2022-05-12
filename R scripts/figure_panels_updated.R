# NorB R analysis 

# -------------- Manipulate data files -------------

# Set working directory
setwd("/Users/david/Desktop/Denitrification")

# Load required packages 
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggpubr)
library(grid)

# Load metagenomic and metatranscriptomic data for each gene 
NapA_ts_df <- read.table(file="SI072_NapA_layered_classifications.tsv",
                          header=TRUE,
                          sep="\t")

NarI_ts_df <- read.table(file="SI072_NarI_layered_classifications.tsv",
                          header=TRUE,
                          sep="\t")

NirK_ts_df <- read.table(file="SI072_NirK_layered_classifications.tsv",
                          header=TRUE,
                          sep="\t")

NorB_ts_df <- read.table(file="SI072_NorB_layered_classifications.tsv",
                          header=TRUE,
                          sep="\t")

NosZ_ts_df <- read.table(file="SI072_NosZ_layered_classifications.tsv",
                          header=TRUE,
                          sep="\t")

# Process taxonomic information
taxa_ranks <- c("Root", "Domain", "Phylum", "Class",
                "Order", "Family", "Genus", "Species")

NapA_ts_df <- NapA_ts_df %>%
  separate(col = Taxonomy,
           into = taxa_ranks,
           sep = "; ", fill = "right", remove = T)

NarI_ts_df <- NarI_ts_df %>%
  separate(col = Taxonomy,
           into = taxa_ranks,
           sep = "; ", fill = "right", remove = T)

NirK_ts_df <- NirK_ts_df %>%
  separate(col = Taxonomy,
           into = taxa_ranks,
           sep = "; ", fill = "right", remove = T)

NorB_ts_df <- NorB_ts_df %>%
  separate(col = Taxonomy,
           into = taxa_ranks,
           sep = "; ", fill = "right", remove = T)

NosZ_ts_df <- NosZ_ts_df %>%
  separate(col = Taxonomy,
           into = taxa_ranks,
           sep = "; ", fill = "right", remove = T)

# Distinguishing metagenomes from metatranscriptomes
NapA_ts_df <- NapA_ts_df %>% 
  mutate(Sample = gsub("_pe", "_MetaG", Sample)) %>% 
  separate(col = Sample,
           into = c("Cruise", "Depth", "SeqType"),
           extra = "drop")

NarI_ts_df <- NarI_ts_df %>% 
  mutate(Sample = gsub("_pe", "_MetaG", Sample)) %>% 
  separate(col = Sample,
           into = c("Cruise", "Depth", "SeqType"),
           extra = "drop")

NirK_ts_df <- NirK_ts_df %>% 
  mutate(Sample = gsub("_pe", "_MetaG", Sample)) %>% 
  separate(col = Sample,
           into = c("Cruise", "Depth", "SeqType"),
           extra = "drop")

NorB_ts_df <- NorB_ts_df %>% 
  mutate(Sample = gsub("_pe", "_MetaG", Sample)) %>% 
  separate(col = Sample,
           into = c("Cruise", "Depth", "SeqType"),
           extra = "drop")

NosZ_ts_df <- NosZ_ts_df %>% 
  mutate(Sample = gsub("_pe", "_MetaG", Sample)) %>% 
  separate(col = Sample,
           into = c("Cruise", "Depth", "SeqType"),
           extra = "drop")

# Processing depth data 
NapA_ts_df <- NapA_ts_df %>%
  mutate(Depth.m = as.numeric(gsub('m', '', Depth))) %>% 
  mutate(Cruise = as.numeric(gsub('SI0', '', Cruise)))

NarI_ts_df <- NarI_ts_df %>%
  mutate(Depth.m = as.numeric(gsub('m', '', Depth))) %>% 
  mutate(Cruise = as.numeric(gsub('SI0', '', Cruise)))

NirK_ts_df <- NirK_ts_df %>%
  mutate(Depth.m = as.numeric(gsub('m', '', Depth))) %>% 
  mutate(Cruise = as.numeric(gsub('SI0', '', Cruise)))

NorB_ts_df <- NorB_ts_df %>%
  mutate(Depth.m = as.numeric(gsub('m', '', Depth))) %>% 
  mutate(Cruise = as.numeric(gsub('SI0', '', Cruise)))

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

# Merge the geochemical data for cruise 72 with the TreeSAPP classifications by Sample
NapA_ts_geo_df <- left_join(NapA_ts_df,
                            geochem_df,
                            by=c("Cruise",
                                 "Depth.m" = "Depth"))

NarI_ts_geo_df <- left_join(NarI_ts_df,
                            geochem_df,
                            by=c("Cruise",
                                 "Depth.m" = "Depth"))

NirK_ts_geo_df <- left_join(NirK_ts_df,
                            geochem_df,
                            by=c("Cruise",
                                 "Depth.m" = "Depth"))

NorB_ts_geo_df <- left_join(NorB_ts_df,
                            geochem_df,
                            by=c("Cruise",
                                 "Depth.m" = "Depth"))

NosZ_ts_geo_df <- left_join(NosZ_ts_df,
                            geochem_df,
                            by=c("Cruise",
                                 "Depth.m" = "Depth"))

#rename taxa
NapA_ts_geo_df <- NapA_ts_geo_df %>% mutate(across(Domain:Species, gsub, pattern = ".__", replacement = ""))

NarI_ts_geo_df <- NarI_ts_geo_df %>% mutate(across(Domain:Species, gsub, pattern = ".__", replacement = ""))

NirK_ts_geo_df <- NirK_ts_geo_df %>% mutate(across(Domain:Species, gsub, pattern = ".__", replacement = ""))

NorB_ts_geo_df <- NorB_ts_geo_df %>% mutate(across(Domain:Species, gsub, pattern = ".__", replacement = ""))

NosZ_ts_geo_df <- NosZ_ts_geo_df %>% mutate(across(Domain:Species, gsub, pattern = ".__", replacement = ""))

# Combine the dataframes using rbind
combined <- rbind(NapA_ts_geo_df, NarI_ts_geo_df, NirK_ts_geo_df, NorB_ts_geo_df, NosZ_ts_geo_df)
combinedMetaG <- combined %>% filter(SeqType == "MetaG")
combinedMetaT <- combined %>% filter(SeqType == "MetaT")


#Apply threshold to aggregate non-significant phyla under "Other" for barplot and line plot

combinedMetaG2 <- combinedMetaG %>%
  mutate(Proportion2 = Abundance/sum(Abundance))

combinedMetaG2$Phylum[combinedMetaG2$Proportion2 < .00005 & is.na(combinedMetaG2$Phylum) == FALSE] <- "z Other"

combinedMetaT2 <- combinedMetaT %>%
  mutate(Proportion2 = Abundance/sum(Abundance))

combinedMetaT2$Phylum[combinedMetaT2$Proportion2 < .0000037 & is.na(combinedMetaT2$Phylum) == FALSE] <- "z Other"


# ----------------- Data visualization --------------------

### Make the relative abundance barplot panel 

MetaG_barplot <- combinedMetaG2 %>%
  group_by(Depth, Depth.m, Marker) %>% 
  mutate(Proportion = Abundance/sum(Abundance)) %>% 
  group_by(Depth, Depth.m, Phylum, Marker) %>% 
  summarise(sum = sum(Proportion)) %>%
  ungroup() %>%
  mutate(Depth = reorder(Depth, Depth.m)) %>% 
  ggplot(aes(x=Depth, y=sum, fill=Phylum)) +
  geom_bar(stat = "identity") +
  facet_wrap(~Marker, ncol=1) +
  labs(x = "Depth", y=NULL) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(legend.text=element_text(size=9))

MetaT_barplot <- combinedMetaT2 %>%
  group_by(Depth, Depth.m, Marker) %>% 
  mutate(Proportion = Abundance/sum(Abundance)) %>% 
  group_by(Depth, Depth.m, Phylum, Marker) %>% 
  summarise(sum = sum(Proportion)) %>% 
  ungroup() %>% 
  mutate(Depth = reorder(Depth, Depth.m)) %>% 
  ggplot(aes(x=Depth, y=sum, fill=Phylum)) +
  geom_bar(stat = "identity") +
  facet_wrap(~Marker, ncol=1) +
  labs(x = "Depth", y=NULL) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(legend.text=element_text(size=9))

MetaG_MetaT_barplot <- ggarrange(MetaG_barplot, MetaT_barplot, labels = c("A", "B"),  ncol = 2, common.legend = TRUE, legend = "right")

MetaG_MetaT_barplot <- annotate_figure(MetaG_MetaT_barplot, left = textGrob("Relative Abundance", rot = 90))

### Make the relative abundance lineplot panel

MetaG_lineplot <- combinedMetaG2 %>%
  group_by(Depth.m, Phylum, Marker) %>%
  summarise(Sum = sum(Abundance)) %>%
  ungroup() %>%
  ggplot(aes(x=Depth.m, y=Sum, colour=Phylum)) +
  geom_line() +
  coord_flip() +
  scale_x_reverse() +
  facet_wrap(~Marker, scales ="free", ncol=1) +
  labs(x=NULL,
       y="Abundance (TPM)") +
  theme_gray(base_size = 14) +
  theme(legend.text=element_text(size=11))

MetaT_lineplot <- combinedMetaT2 %>%
  group_by(Depth.m, Phylum, Marker) %>%
  summarise(Sum = sum(Abundance)) %>%
  ungroup() %>%
  ggplot(aes(x=Depth.m, y=Sum, colour=Phylum)) +
  geom_line() +
  coord_flip() +
  scale_x_reverse() +
  facet_wrap(~Marker, scales = "free", ncol=1) +
  labs(x=NULL,
       y="Abundance (TPM)") +
  theme_gray(base_size = 14) +
  theme(legend.text=element_text(size=11))

MetaG_MetaT_lineplot <- ggarrange(MetaG_lineplot, MetaT_lineplot, labels = c("A", "B"),  ncol = 2, common.legend = TRUE, legend = "right")

MetaG_MetaT_lineplot <- annotate_figure(MetaG_MetaT_lineplot, left = textGrob("Depth (m)", rot = 90, gp = gpar(fontsize = 15)))

