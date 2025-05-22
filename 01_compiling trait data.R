### Combining and plotting trait data
###
### Author: Kevin Wilcox
### March 27, 2025

### set up workspace
rm(list=ls())
setwd("C:\\Users\\K_WILCOX\\Dropbox\\Kruger Compound Extremes\\Data\\2023 Data\\")

library(tidyverse)

### read in individual cleaned data
leaf_traits <- read.csv("Leaves\\leaf_calculations_2-12.csv") %>%
  group_by(block, plot, plant_tag, species) %>%
  summarize_at(vars(total_area:LDMC), .funs=c(mean), na.rm=T) %>%
  ungroup()

bunch_traits <- read.csv("Pop demo\\cache\\popdemo traits_for merging_March2025.csv") %>%
  dplyr::select(-year, -month)

aboveground_traits <-  bunch_traits %>%
  full_join(leaf_traits, by=c("block", "plot", "plant_tag", "species"))

leaf_traits_nondom <- read.csv("Traits\\NExS2023_traits_nondoms_cleaned_20250410.csv") %>%
  rename(plant_tag=tag_num, species=species_id) %>%
  mutate(species=replace(species, species=="ari_con", "ari con")) %>%
  mutate(species=replace(species, species=="the_tri", "the tri")) %>%
  mutate(species=replace(species, species=="bot_rad", "bot rad")) %>%
  mutate(species=replace(species, species=="pan_col", "pan col")) %>%
  mutate(species=replace(species, species=="pan_max", "pan max")) %>%
  group_by(block, plot, plant_tag, species) %>%
  summarize_at(vars(tiller_number:SLA), .funs=c(mean), na.rm=T) %>%
  ungroup() %>%
  mutate(LDMC=LDMC*1000) %>% # fix unit inconsistancy
  dplyr::select(block:species, veg_height, total_area, thickness, SLA, LDMC) %>%
  rename(thick=thickness)

traits_dom_and_nondom <- trait_df %>%
  dplyr::select(block:species, veg_height, total_area, thick, SLA, LDMC) %>%
  bind_rows(leaf_traits_nondom)
  

### cleaning and combining srl and root diameter data
tags_to_remove <- c(26,21,22,74,183) # these are tags that we are missing either scans or samples for

srl_df <- read.csv("..\\root_data\\NExS_root_traits_SRL_2023_cleaned_2.csv") %>%
  mutate(plant_tag=replace(plant_tag, plant_tag %in%c("22B","22A"),"22")) %>%
  rename(SRL=srl) %>%
  dplyr::select(block:species, SRL) %>%
  mutate(plant_tag=as.numeric(plant_tag)) %>%
  filter(!plant_tag %in% tags_to_remove)
  
rootDiam_df <- read.csv("..\\root_data\\NExS_root_traits_extra_2023_cleaned.csv") %>%
  dplyr::select(block:species, Average.Diameter.mm) %>%
  group_by(block, plot, plant_tag, species) %>%
  summarize(Average.Diameter.mm = mean(Average.Diameter.mm,na.rm=T)) %>%
  ungroup() %>%
  rename(root_diam_mm = Average.Diameter.mm) %>%
  mutate(species=replace(species, species=="aricon", "ari con")) %>%
  mutate(species=replace(species, species=="thetri", "the tri")) %>%
  mutate(species=replace(species, species=="botrad", "bot rad")) %>%
  mutate(species=replace(species, species=="pancol", "pan col")) %>%
  mutate(species=replace(species, species=="panmax", "pan max")) %>%
  filter(!plant_tag %in% tags_to_remove)

roots_df <- srl_df %>%
  full_join(rootDiam_df, by=c("block", "plot", "plant_tag", "species"))

write.csv(roots_df, file="Roots\\NExS_srl and diamter traits_working_2023.csv", row.names=F)
### Bring it all together
trait_df <- bunch_traits %>%
  full_join(leaf_traits, by=c("block", "plot", "plant_tag", "species")) %>%
  full_join(rootDiam_df, by=c("block", "plot", "plant_tag", "species")) %>%
  full_join(srl_df, by=c("block", "plot", "plant_tag", "species")) %>%
  dplyr::select(block:species,
                veg_height,
                basal_arial_ratio,
                tillerdensity,
                total_area,
                thick,
                SLA,
                LDMC,
                root_diam_mm,
                SRL)
corrplot(trait_df[,5:13])

biplot(trait_df)
corPlot(trait_df[,5:18])
install.packages("corPlot")
install.packages("GGally")
install.packages("corrplot")
install.packages("ggfortify")
library(corrplot)
library(GGally)
library(ggfortify)

?ggcorr
ggcorr(subset(trait_df, tillerdensity<11&thick<1&LDMC<1000&root_diam_mm<1.2)[,5:13], label=T, label_alpha=T)
ggpairs(subset(trait_df, tillerdensity<11&thick<1&LDMC<1000&root_diam_mm<1.2), columns = 5:13,
        ggplot2::aes(col=species))

trait_data <- trait_df %>% 
  na.omit(.) %>%
  dplyr::select(veg_height:SRL)

trait_env <- trait_df %>% 
  na.omit(traits_data) %>%
  dplyr::select(block:species)

pca_trait <- prcomp(trait_data, scale.=T)

autoplot(pca_trait, data=trait_env, colour="species",
         loadings = TRUE, loadings.colour = 'blue',
         loadings.label = TRUE, loadings.label.size = 3,
         size=3)

### non-doms only
trait_nondom_data <- leaf_traits_nondom %>% 
  na.omit(.) %>%
  dplyr::select(veg_height:LDMC)

trait_nondom_env <- leaf_traits_nondom %>% 
  na.omit(.) %>%
  dplyr::select(block:species)

pca_nondom_trait <- prcomp(trait_nondom_data, scale.=T)

autoplot(pca_nondom_trait, data=trait_nondom_env, colour="species",
         loadings = TRUE, loadings.colour = 'blue',
         loadings.label = TRUE, loadings.label.size = 3,
         size=3, frame=T)+
  scale_color_manual(values=c25[1:21]) +
  scale_fill_manual(values=c25[1:21])

### dominant and non doms together
trait_allsp_data <- traits_dom_and_nondom %>% 
  na.omit(.) %>%
  dplyr::select(veg_height:LDMC)

trait_allsp_env <- traits_dom_and_nondom %>% 
  na.omit(.) %>%
  dplyr::select(block:species)

pca_allsp_trait <- prcomp(trait_allsp_data, scale.=T)

autoplot(pca_allsp_trait, data=trait_allsp_env, colour="species",
         loadings = TRUE, loadings.colour = 'blue',
         loadings.label = TRUE, loadings.label.size = 3,
         size=3, frame=T)+
         scale_color_manual(values=c25[1:21])

ggpairs(subset(traits_dom_and_nondom, LDMC<1000 & !species %in% c("cen_cil","com_car","ehr_rig","eri_spp","her_gla","ind_vic","oci_ame","pav_bur","pel_leu","sch_pap","tep_pur")), columns = 5:9)
ggpairs(subset(traits_dom_and_nondom, LDMC<1000 & !species %in% c("cen_cil","com_car","ehr_rig","eri_spp","her_gla","ind_vic","oci_ame","pav_bur","pel_leu","sch_pap","tep_pur")), columns = 5:9,
  ggplot2::aes(col=species))


ggpairs(subset(leaf_traits_nondom, LDMC<1000 & !species %in% c("cen_cil","com_car","ehr_rig","eri_spp","her_gla","ind_vic","oci_ame","pav_bur","pel_leu","sch_pap","tep_pur")), columns = 5:9,
        ggplot2::aes(col=species))

with(leaf_traits_nondom, table(species))


c25 <- c(
  "dodgerblue2", "#E31A1C", # red
  "green4",
  "#6A3D9A", # purple
  "#FF7F00", # orange
  "black", "gold1",
  "skyblue2", "#FB9A99", # lt pink
  "palegreen2",
  "#CAB2D6", # lt purple
  "#FDBF6F", # lt orange
  "gray70", "khaki2",
  "maroon", "orchid1", "deeppink1", "blue1", "steelblue4",
  "darkturquoise", "green1", "yellow4", "yellow3",
  "darkorange4", "brown"
)
