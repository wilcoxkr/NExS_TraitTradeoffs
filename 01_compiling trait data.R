### Combining and plotting trait data
###
### Author: Kevin Wilcox (k_wilcox@uncg.edu)
### created: March 27, 2025; Last updated: July 22, 2025

### set up workspace
rm(list=ls())
setwd("C:\\Users\\K_WILCOX\\Dropbox\\Kruger Compound Extremes\\Data\\2023 Data\\")

library(tidyverse)
library(readr)

#Standard error function
sefxn = function(x, na.rm=na.rm) {
  SE = sd(x, na.rm=TRUE)/sqrt(length(x[!is.na(x)]))
  return(SE)
}

###
### Clean leaf data (Modified from Shelby's script)
###
{
  leaf_path <- "C:\\Users\\K_WILCOX\\Dropbox\\Kruger Compound Extremes\\Data\\2023 Data\\Leaves\\"
  ### Read in data
  leaftraits_noarea_raw <- read.csv(paste0(leaf_path, "NExS_main-trait-data_2023_sw_JA.csv"))
  leaftraits_area_raw <- read.csv(paste0(leaf_path, "NExS_Leaf_area_2023_Report_cleaned_sw_2-12.csv"))
  
  dom_sp_to_keep <- c("ari con",
                      "the tri",
                      "bot rad",
                      "pan col",
                      "pan max")
  
  ## Leaf area data
  ##
  leaftraits_area_df <- leaftraits_area_raw %>% 
    filter(QC_flag != 1) %>%
    dplyr::select(block, plot, plant_tag, species, leaf, total_area) %>% 
    mutate(leaf = as.numeric(leaf)) %>% 
    mutate(block = as.factor(block)) %>% 
    drop_na(plant_tag) %>% #drop non-doms
    mutate(species = replace(species, species=="pan col ", "pan col")) %>% # Pesky space at the end was removing pan col data from tag 165 -- found during merge issue below
    filter(species %in% dom_sp_to_keep) %>%
    filter(!plant_tag %in% 16:20)
  
  # Check leaf area data frame for missing or weird data -- I fixed these in the NExS_Leaf_area_2023_Report_cleaned_sw_2-12.csv
  plot(with(leaftraits_area_df, table(plant_tag)))
  with(leaftraits_area_df, table(block, plot))
  filter(leaftraits_area_df, plot==2) 
    # NOTE there are 30 leaves from tagged individuals in plot 2. This is because we tagged 
    # Pan col as a non dom in this plot (and its of course in our 'dom' list, just not for block A)
    # I am keeping these data in here, but if we do any type of dom vs non-dom analysis, that is
    # agnostic to species identity, we should put the pan col here into the 'non dom' category
    
    # NOTE 2: we actually have pan col tagged as a non-dom (different tag numbers, which I think were updated)
    # so I will remove pan col from this data frame (a few lines above) to avoid double counting
  leaftraits_area_df %>% # looks good
    ggplot(aes(x=species, y=total_area)) + geom_boxplot() + geom_jitter()
    
  ## All other leaf traits
  ##
  
  # Check main trait data for missing or weird data - looks good
  with(leaftraits_noarea_raw, table(plant_tag))
  
  leaftraits_noarea_df <- leaftraits_noarea_raw %>% 
    select(!(tiller_num:diameter_t3)) %>% 
    select(!(year:month)) %>% 
    pivot_longer(
      cols = !c(block:species), 
      names_to = c(".value", "leaf"), 
      names_sep = "_l",
      values_drop_na = TRUE ## NOTE: There are some missing values in here, but these were all field or processing errors, so we should just drop them and move on
    ) %>% 
    mutate(leaf = as.numeric(leaf)) %>%
    mutate(block = as.factor(block))
  
  with(leaftraits_noarea_df, table(plot))
  
  leaftraits_full_df <- leaftraits_area_df %>% 
    full_join(leaftraits_noarea_df, by=c("block","plot","species","plant_tag","leaf")) %>% #change join to only keep tagged plants... KW note -- the left join masked some data errors found when full joining -- we should always full join then drop NAs after examination
    filter(
      !(block=="A" & plot==5 & plant_tag==31 & species=="ari con" & leaf==2) & # REMOVE: A5 tag 31 aricon leaf 2 is missing area and dry mass for leaf 2, missing from spreadsheet
      !(block=="A" & plot==5 & plant_tag==35 & species=="ari con" & leaf==2) & # REMOVE A5 tag 35 aricon leaf 2: two leaf#2 values for area so we need to drop both
      !(block=="A" & plot==6 & plant_tag==37 & species=="ari con" & leaf==2) & # REMOVE A6 tag 37 leaf 2  – missing area data point in raw data
      !(block=="A" & plot==7 & plant_tag==44 & species=="ari con" & leaf==3) & #	REMOVE A7 tag 44 leaf 3 – two leaves labeled 3 so have to remove both
      !(block=="A" & plot==7 & plant_tag==45 & species=="ari con" & leaf==1) & #	REMOVE A7 tag 45 leaf 1 – no area in raw spreadsheet or dry mass
      !(block=="A" & plot==7 & plant_tag==45 & species=="ari con" & leaf==3) & #	REMOVE A7 tag 45 leaf 3 – no area in raw spreadsheet.
      !(block=="F" & plot==41 & plant_tag==215 & species=="pan max" & leaf==1) & #	REMOVE: F41 tag 215 leaf 1&3 – there were two leaves labeled #3 in area spreadsheet, so we have to remove both of these since we don’t know which is which.
      !(block=="F" & plot==41 & plant_tag==215 & species=="pan max" & leaf==3) & #	REMOVE: F41 tag 215 leaf 1&3 – there were two leaves labeled #3 in area spreadsheet, so we have to remove both of these since we don’t know which is which.
      !(block=="D" & plot==28 & plant_tag==148 & species=="pan col") &           #	REMOVE: D28 tag 148 leaf 1,2,&3 – Just missing the whole plant for mass ... bummer.
      !(block=="C" & plot==24 & plant_tag==128 & species=="bot rad" & leaf==3) & #	REMOVE: C24 plant tag 128 bot rad leaf 3, missing thickness and mass from spreadsheet
      !(block=="B" & plot==10 & plant_tag==56 & species=="the tri" & leaf==1) & #	REMOVE: missing area data from spreadsheet
      !(block=="B" & plot==10 & plant_tag==57 & species=="the tri" & leaf==1) & #	REMOVE:  missing area data from spreadsheet
      !(block=="B" & plot==10 & plant_tag==57 & species=="the tri" & leaf==3) & #	REMOVE:  missing area data from spreadsheet
      !(block=="B" & plot==10 & plant_tag==58 & species=="the tri" & leaf==2) & #	REMOVE:  missing area data from spreadsheet
      !(block=="B" & plot==10 & plant_tag==59 & species=="the tri" & leaf==1) & #	REMOVE:  missing area data from spreadsheet
      !(block=="B" & plot==10 & plant_tag==59 & species=="the tri" & leaf==3) & #	REMOVE:  missing area data from spreadsheet
      !(block=="B" & plot==10 & plant_tag==60 & species=="the tri" & leaf==2) &  #	REMOVE:  missing mass data from spreadsheet
      !(block=="E" & plot==39 & plant_tag==204 & species=="pan col" & leaf==2) &  #	REMOVE:  missing mass data from spreadsheet
      !(block=="D" & plot==31 & plant_tag==161 & species=="pan col" & leaf==3) &  #	REMOVE:  missing mass data from spreadsheet
      !(block=="D" & plot==27 & plant_tag==145 & species=="pan col" & leaf==3) & #	REMOVE:  missing thickness data from spreadsheet, removing for completeness
      !(block=="A" & plot==4 & plant_tag==30 & species=="ari con" & leaf==2) & # REMOVE: Dry mass is greater than wet_mass -- likely an instrument error since mass values were very small (.006 g dry mass)
      !(block=="A" & plot==8 & plant_tag==47 & species=="ari con" & leaf==3) # REMOVE: Dry mass is greater than wet_mass -- likely an instrument error since mass values were very small
      ) %>%
    mutate(SLA = total_area/dry_mass) %>% #SLA
    mutate(LDMC = dry_mass/wet_mass) #Check LDMC calculations
  
  # check to see if we have any missing data
  leaftraits_full_df %>% filter(is.na(total_area))
  leaftraits_full_df %>% filter(is.na(thick))
  leaftraits_full_df %>% filter(is.na(dry_mass))
  leaftraits_full_df %>% filter(is.na(wet_mass)) ## Woohoo, no more missing data!
  
  
  ## Look for outliers and check data
  ##
  
  # LDMC
  leaftraits_full_df %>%
    ggplot(aes(x=species, y=LDMC, label=plant_tag)) + geom_jitter(position=position_jitter(seed=1)) +
    geom_text(position=position_jitter(seed=1))
  
  leaftraits_full_df %>% filter(plant_tag==30) # TO REMOVE: Dry mass is greater than wet_mass -- likely an instrument error since mass values were very small (.006 g dry mass)
  leaftraits_full_df %>% filter(plant_tag==47) # TO REMOVE: Dry mass is greater than wet_mass -- likely an instrument error since mass values were very small
  leaftraits_full_df %>% filter(plant_tag==28) # KEEP: Dry mass is just close to wet_mass -- 
  leaftraits_full_df %>% filter(plant_tag==121) # FIXED: Dry mass is greater than wet_mass -- This was a data entry error, the data sheet scan showed 0.040 instead of 0.004 -- I fixed this in the NExS_main-trait-data_2023_sw_JA.csv spreadsheet
  leaftraits_full_df %>% filter(plant_tag==248) # KEEP: Dry mass is just close to wet mass
  
  # SLA
  leaftraits_full_df %>%
    ggplot(aes(x=species, y=SLA, label=plant_tag)) + geom_jitter(position=position_jitter(seed=1)) +
    geom_text(position=position_jitter(seed=1))
  
  leaftraits_full_df %>% filter(plant_tag==94) # KEEP: this seems like a real value -- total area and dry mass pretty in line with other values
  leaftraits_full_df %>% filter(plant_tag==250) # KEEP: this also seems like a real value -- total area and dry mass pretty in line with other values
  
  ## Remove data errors and replot (I actually brought this up to the main cleaning step of leaftraits_full_df)
  ## 
  # leaftraits_full_df <- leaftraits_full_df %>%
  #   filter(
  #     !(block=="A" & plot==4 & plant_tag==30 & species=="ari con" & leaf==2) & # REMOVE: Dry mass is greater than wet_mass -- likely an instrument error since mass values were very small (.006 g dry mass)
  #     !(block=="A" & plot==8 & plant_tag==47 & species=="ari con" & leaf==3) # REMOVE: Dry mass is greater than wet_mass -- likely an instrument error since mass values were very small
  #   )
    

  write.csv(leaftraits_full_df, "Traits\\NExS 2023_leaf area thickness drymass wetmass SLA LDMC_dom species_CLEANED_kw_20250722.csv", row.names=T)
  
  leaf_calculations %>% 
    #filter(SLA<600) %>% 
    ggplot(aes(x=SLA, group=Species, fill=Species)) +
    geom_density(adjust=1.5, alpha=.4) +
    theme(text = element_text(size = 30)) +
    theme(legend.text = element_text(size = 30))
  
  leaf_calculations %>% 
    #filter(LDMC < 800) %>% 
    ggplot(aes(x=LDMC, group=Species, fill=Species)) +
    geom_density(adjust=1.5, alpha=.4) +
    theme(text = element_text(size = 30)) +
    theme(legend.text = element_text(size = 30)) 

}


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
  dplyr::select(block:species, tiller_number, veg_height, tiller_diameter, total_area, thick, SLA, LDMC) %>%
  bind_rows(leaf_traits_nondom)
  

### cleaning and combining srl and root diameter data
tags_to_remove <- c(26,21,22,74,183) # these are tags that we are missing either scans or samples for

srl_df <- read.csv("..\\root_data\\NExS_root_traits_SRL_2023_cleaned_2.csv") %>%
  mutate(plant_tag=replace(plant_tag, plant_tag %in% c("22B","22A"),"22")) %>%
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
