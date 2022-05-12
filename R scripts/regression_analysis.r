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
ggplotRegression <- function (fit) {
  
  require(ggplot2)
  
  ggplot(fit$model, aes_string(x = names(fit$model)[2], y = names(fit$model)[1])) + 
    geom_point() +
    stat_smooth(method = "lm", col = "red") +
    labs(title = paste("Adj R2 = ",signif(summary(fit)$adj.r.squared, 5),
                       "Intercept =",signif(fit$coef[[1]],5 ),
                       " Slope =",signif(fit$coef[[2]], 5),
                       " P =",signif(summary(fit)$coef[2,4], 5)))
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
shannon_nosz <- shannon_diversity(nosz, "nosz", 34.26947) 
# https://biocyc.org/META/NEW-IMAGE?type=REACTION&object=NITRIC-OXIDE-REDUCTASE-RXN
shannon_norb <- shannon_diversity(norb, "norb", 8.969482)
shannon_nari <- shannon_diversity(nari, "nari", -69.334)
# https://www.sciencedirect.com/science/article/pii/S0005272806001046
shannon_nirk <- shannon_diversity(nirk, "nirk", -45.3)
# http://vm-trypanocyc.toulouse.inra.fr/META/NEW-IMAGE?type=REACTION&object=NITRATE-REDUCTASE-CYTOCHROME-RXN&detail-level=2&orgids=LEISH
shannon_napa <- shannon_diversity(napa, "napa", -74.12053)
total_data <- rbind(shannon_nosz, shannon_norb, shannon_nari, shannon_nirk, shannon_napa)

# # Model diversity ----------------------------------------------------
# Remove outliers and subset data
total_data_test <- total_data %>% filter(H>0.2)


summary(lm(data=total_data_test, H ~ G*SeqType + I(G^2) ))

# diagnostic plots 
plot(lm(data=total_data_test, H ~ G*SeqType + I(G^2)))
#Save model and use for stat_smooth (doesnt work)
model <- lm(data=total_data_test, H ~ G*SeqType + I(G^2))

metaG_data <- filter(total_data, SeqType == "MetaG") %>% filter(H>0.2)
metaT_data <- filter(total_data, SeqType == "MetaT") %>% filter(H > 0.2)

#Do this because I do not know how to plot the dummy multivariable lm correctly
metaG_curve <- total_data_test %>% filter(H>0.2)
metaG_curve$H <- metaG_curve$H + .485
metaT_curve <- total_data_test %>% filter(H>0.2)
metaT_curve$H <- metaT_curve$H - .515

# Plot diversity across MAGs and Metatranscriptomic data -----------------
ggplot(data = total_data_test, aes(x=G, y=H)) + 
  stat_smooth(se = FALSE, data= metaG_curve, method = lm, formula= y ~ poly(x, 2), colour = "red") +
  geom_point(alpha=0.5, size=2, aes(color=SeqType)) + xlab("Standard Delta G (kcal/mol)") + ylab("H'")
  stat_smooth(se = FALSE, data= metaT_curve, method = lm, formula= y ~ poly(x, 2), colour = "turquoise2")
  




