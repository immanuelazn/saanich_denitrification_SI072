### MICB 425 NapA Figures from TreeSAPP Outputs
## April 23, 2022
# Jerry He


# Load packages -----------------------------------------------------------
library(tidyverse)


# Code --------------------------------------------------------------------
ts_dat <- read_tsv(file = "./../local/SI072_layered_classifications.tsv")

ts_df <- ts_dat %>% 
  filter(Marker == "NapA")

taxa_ranks <- c("Root", "Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")

ts_df <- ts_df %>% 
  separate(col = Taxonomy,
           into = taxa_ranks,
           sep = "; ", fill = "right", remove = TRUE)

ts_df <- ts_df %>% 
  mutate(Sample = gsub("_pe", "_MetaG", Sample)) %>% 
  separate(col = Sample,
           into = c("Cruise", "Depth", "SeqType"), 
           extra = "drop")

ts_df <- ts_df %>% 
  mutate(Depth.m = as.numeric(gsub("m", "", Depth))) %>% 
  mutate(Cruise = as.numeric(gsub("SI0", "", Cruise)))

geochem_df <- read_csv("Saanich_TimeSeries_Chemical_DATA.csv") %>% 
  filter(Cruise == 72) %>% 
  select(Cruise, Depth, CTD_O2, NO3, Mean_NH4, Mean_NO2, Mean_H2S, Mean_CH4) %>% 
  mutate(Mean_CH4 = Mean_CH4*1E-3) %>% 
  mutate(across(.fns = as.numeric))

NapA_ts_geo_df <- left_join(ts_df,
                            geochem_df,
                            by = c("Cruise", 
                                   "Depth.m" = "Depth"))
NapA_ts_geo_df <- NapA_ts_geo_df %>% 
  mutate(across(Domain:Species, gsub, pattern = ".__", replacement = ""))



# Plots at Order level ----------------------------------------------------

NapA_ts_geo_df %>% 
  filter(SeqType == "MetaT") %>% 
  group_by(Depth.m, Order) %>% 
  summarize(Sum = sum(Abundance)) %>% 
  ungroup() %>% 
  ggplot(aes(x = Depth.m, y = Sum, colour = Order)) +
  geom_line() +
  theme_classic() +
  coord_flip() +
  scale_x_reverse() +
  labs(x = "Depth (m)",
       y = "Relative abundance (TPM)") +
  ggtitle("NapA-expressing Orders (metatranscriptome) by depth in the Saanich Inlet") +
  theme_bw() 

NapA_ts_geo_df %>% 
  filter(SeqType == "MetaG") %>% 
  group_by(Depth.m, Order) %>% 
  summarize(Sum = sum(Abundance)) %>% 
  ungroup() %>% 
  ggplot(aes(x = Depth.m, y = Sum, colour = Order)) +
  geom_line() +
  theme_classic() +
  coord_flip() +
  scale_x_reverse() +
  labs(x = "Depth (m)",
       y = "Relative abundance (TPM)") +
  ggtitle("NapA-containing Orders (metagenome) by depth in the Saanich Inlet ") +
  theme_bw()   
  

NapA_ts_geo_df %>% 
  group_by(Depth, Depth.m, SeqType) %>% 
  mutate(Proportion = Abundance/sum(Abundance)) %>% 
  group_by(Depth, Depth.m, Order, SeqType) %>% 
  summarize(sum = sum(Proportion)) %>% 
  ungroup() %>% 
  mutate(Depth = reorder(Depth, Depth.m)) %>% 
  ggplot(aes(x = Depth, y = sum, fill = Order)) +
  geom_bar(stat = "identity") +
  facet_wrap(~SeqType) +
  labs(x = "Depth", 
       y = "Relative abundance") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("Abundance of NapA-containing Orders (metagenome) by depth in the Saanich Inlet")
  theme_bw()
  

# Plots at Phylum level --------------------------------------------------
NapA_ts_geo_df %>% 
  filter(SeqType == "MetaT") %>% 
  group_by(Depth.m, Phylum) %>% 
  summarize(Sum = sum(Abundance)) %>% 
  ungroup() %>% 
  ggplot(aes(x = Depth.m, y = Sum, colour = Phylum)) +
  geom_line() +
  theme_classic() +
  coord_flip() +
  scale_x_reverse() +
  ggtitle("NapA-expressing Phyla (metatranscriptome) by depth in the Saanich Inlet") +
  labs(x = "Depth (m)",
       y = "Relative abundance (TPM)") +
  theme_bw() 

NapA_ts_geo_df %>% 
  filter(SeqType == "MetaG") %>% 
  group_by(Depth.m, Phylum) %>% 
  summarize(Sum = sum(Abundance)) %>% 
  ungroup() %>% 
  ggplot(aes(x = Depth.m, y = Sum, colour = Phylum)) +
  geom_line() +
  theme_classic() +
  coord_flip() +
  scale_x_reverse() +
  ggtitle("NapA-containing Phyla (metagenome) by depth in the Saanich Inlet") +
  labs(x = "Depth (m)",
       y = "Relative abundance (TPM)") +
  theme_bw() 

NapA_ts_geo_df %>% 
  group_by(Depth, Depth.m, SeqType) %>% 
  mutate(Proportion = Abundance/sum(Abundance)) %>% 
  group_by(Depth, Depth.m, Phylum, SeqType) %>% 
  summarize(sum = sum(Proportion)) %>% 
  ungroup() %>% 
  mutate(Depth = reorder(Depth, Depth.m)) %>% 
  ggplot(aes(x = Depth, y = sum, fill = Phylum)) +
  geom_bar(stat = "identity") +
  facet_wrap(~SeqType) +
  labs(x = "Depth", 
       y = "Relative abundance") +
  ggtitle("Abundance of NapA-containing Phyla (metagenome) by depth in the Saanich Inlet") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 


# Plots at Class level ----------------------------------------------------
NapA_ts_geo_df %>% 
  filter(SeqType == "MetaT") %>% 
  group_by(Depth.m, Class) %>% 
  summarize(Sum = sum(Abundance)) %>% 
  ungroup() %>% 
  ggplot(aes(x = Depth.m, y = Sum, colour = Class)) +
  geom_line() +
  theme_classic() +
  coord_flip() +
  scale_x_reverse() +
  labs(x = "Depth (m)",
       y = "Relative abundance (TPM)") +
  ggtitle("NapA-expressing Classes (metatranscriptome) by depth in the Saanich Inlet") +
  theme_bw() 

NapA_ts_geo_df %>% 
  filter(SeqType == "MetaG") %>% 
  group_by(Depth.m, Class) %>% 
  summarize(Sum = sum(Abundance)) %>% 
  ungroup() %>% 
  ggplot(aes(x = Depth.m, y = Sum, colour = Class)) +
  geom_line() +
  theme_classic() +
  coord_flip() +
  scale_x_reverse() +
  labs(x = "Depth (m)",
       y = "Relative abundance (TPM)") +
  ggtitle("NapA-containing Classes (metagenome) by depth in the Saanich Inlet ") +
  theme_bw()   


NapA_ts_geo_df %>% 
  group_by(Depth, Depth.m, SeqType) %>% 
  mutate(Proportion = Abundance/sum(Abundance)) %>% 
  group_by(Depth, Depth.m, Class, SeqType) %>% 
  summarize(sum = sum(Proportion)) %>% 
  ungroup() %>% 
  mutate(Depth = reorder(Depth, Depth.m)) %>% 
  ggplot(aes(x = Depth, y = sum, fill = Class)) +
  geom_bar(stat = "identity") +
  facet_wrap(~SeqType) +
  labs(x = "Depth", 
       y = "Relative abundance") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("Abundance of NapA-containing Classes (metagenome) by depth in the Saanich Inlet")
theme_bw()



