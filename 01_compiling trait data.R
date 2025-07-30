### Combining and plotting trait data
###
### Author: Kevin Wilcox (k_wilcox@uncg.edu)
### created: March 27, 2025; Last updated: July 28, 2025

###
### Set up workspace
###
{
rm(list=ls())
setwd("C:\\Users\\K_WILCOX\\Dropbox\\Kruger Compound Extremes\\Data\\2023 Data\\")

library(tidyverse)
library(readr)
library(GGally)


#Standard error function
sefxn = function(x, na.rm=na.rm) {
  SE = sd(x, na.rm=TRUE)/sqrt(length(x[!is.na(x)]))
  return(SE)
}

##
## Plot Key
plot_key <- read.csv("..//NExS_plot_key.csv") %>%
  rename(block=Block, 
         plot=Plot, 
         drought=Drought, 
         grazing=Grazing, 
         surface_sensor=Surface.Sensor,
         sensor_install_date=Install.Date)


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
  # NOTE: I found an error in E36 plant tag 186, where leaf thickness was 16 in one of the leaves. This should have been 0.16 -- I fixed this in the NExS_main-trait-data_2023_sw_JA
  with(leaftraits_noarea_raw, table(plant_tag))
  
  leaftraits_noarea_df <- leaftraits_noarea_raw %>% 
    rename(diameter_l1=diameter_t1,diameter_l2=diameter_t2,diameter_l3=diameter_t3) %>%
    select(!(tiller_num:height_t3)) %>% 
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
      !(block=="A" & plot==8 & plant_tag==47 & species=="ari con" & leaf==3) & # REMOVE: Dry mass is greater than wet_mass -- likely an instrument error since mass values were very small
      !(block=="B" & plot==9 & plant_tag==53 & species=="the tri") # REMOVE: missing all leaf data, except tiller diameter
      ) %>%
    mutate(SLA = total_area/dry_mass) %>% #SLA
    mutate(LDMC = dry_mass/wet_mass) %>% #Check LDMC calculations
    mutate(species=replace(species, species=="ari con", "ari_con")) %>%
    mutate(species=replace(species, species=="the tri", "the_tri")) %>%
    mutate(species=replace(species, species=="bot rad", "bot_rad")) %>%
    mutate(species=replace(species, species=="pan col", "pan_col")) %>%
    mutate(species=replace(species, species=="pan max", "pan_max"))
    
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
  
  leaftraits_full_df %>% filter(plant_tag==30) # TO REMOVE: Dry mass is greater than wet_mass -- likely an instrument error since mass values were very small (.006 g dry mass) - removed in data step above
  leaftraits_full_df %>% filter(plant_tag==47) # TO REMOVE: Dry mass is greater than wet_mass -- likely an instrument error since mass values were very small - removed in data step above
  leaftraits_full_df %>% filter(plant_tag==28) # KEEP: Dry mass is just close to wet_mass -- 
  leaftraits_full_df %>% filter(plant_tag==121) # FIXED: Dry mass is greater than wet_mass -- This was a data entry error, the data sheet scan showed 0.040 instead of 0.004 -- I fixed this in the NExS_main-trait-data_2023_sw_JA.csv spreadsheet
  leaftraits_full_df %>% filter(plant_tag==248) # KEEP: Dry mass is just close to wet mass
  
  # SLA
  leaftraits_full_df %>%
    ggplot(aes(x=species, y=SLA, label=plant_tag)) + geom_jitter(position=position_jitter(seed=1)) +
    geom_text(position=position_jitter(seed=1))
  
  leaftraits_full_df %>% filter(plant_tag==94) # KEEP: this seems like a real value -- total area and dry mass pretty in line with other values
  leaftraits_full_df %>% filter(plant_tag==250) # KEEP: this also seems like a real value -- total area and dry mass pretty in line with other values

  # Tiller diameter
  leaftraits_full_df %>%
    ggplot(aes(x=species, y=diameter, label=plant_tag)) + geom_jitter(position=position_jitter(seed=1)) +
    geom_text(position=position_jitter(seed=1)) # they all look realistic
  
  # Tiller diameter
  leaftraits_full_df %>%
    ggplot(aes(x=species, y=total_area, label=plant_tag)) + geom_jitter(position=position_jitter(seed=1)) +
    geom_text(position=position_jitter(seed=1)) # they all look realistic
  
  ##
  ## Write cleaned data to file and clean up other data objects
  
    write.csv(leaftraits_full_df, "Traits\\NExS 2023_leaf area thickness drymass wetmass SLA LDMC tillerDiam_dom species_CLEANED_kw_20250729.csv", row.names=T)
    rm(leaftraits_area_df, 
       leaftraits_area_raw, 
       leaftraits_noarea_df,
       leaftraits_noarea_raw)
}

###
### Read in other aboveground traits, check, and combine with non-dom data
###
{
### read in individual cleaned data
leaftraits_df <- leaftraits_full_df %>%
  group_by(block, plot, plant_tag, species) %>%
  summarize_at(vars(c(total_area, thick, dry_mass, wet_mass, SLA, LDMC, diameter)), .funs=c(mean), na.rm=T) %>%
  ungroup() %>%
  rename(leaf_thick=thick,
         tiller_diam=diameter)
  

with(leaftraits_df, table(block)) # looks good, just a couple of plots missing data from individuals
with(leaftraits_df, table(block, plot)) # looks good, just a couple of plots missing data from individuals
rm(leaftraits_full_df)

bunchtraits_df <- read.csv("..\\01_PopDemo\\04_data_clean\\NExS_PopDemo_March2023_clean.csv") %>%
  full_join(dplyr::select(plot_key, block, plot, drought, grazing, Fire),
            by=c("block", "plot")) %>%
  mutate(basal1 = replace(basal1, plant_tag==97, 1.0)) %>% # I'm pretty sure this is a decimal error, but we should find the raw data to confirm
  mutate(basal2 = replace(basal2, plant_tag==97, 1.0)) %>% # I'm pretty sure this is a decimal error, but we should find the raw data to confirm
  mutate(basal2 = replace(basal2, plant_tag==68, 5.0)) %>% # I'm pretty sure this is a decimal error, but we should find the raw data to confirm
  mutate(basal_area = pi*(basal1/2)*(basal2/2)) %>% 
#  mutate(aerialarea = pi*(aerial1/2)*(aerial2/2)) %>% # Not using this for trait tradeoffs
#  mutate(conevolume = (1/3)*pi* veg_height*(((basal1+basal2)/2)^2)+((basal1+basal2)/2)* (((aerial1+aerial2)/2)^2)+(((aerial1+aerial2)/2))) %>% # not using this for trait tradeoffs
  dplyr::select(-comments, -year, -month, -drought, -grazing, -Fire)

#bunchtraits_df %>% filter(basal2<1)

# Quick look at values of some key traits
bunchtraits_df %>% ggplot(aes(x=block, y=veg_height)) + geom_jitter(width=0.2) + geom_boxplot(alpha=0.5)
bunchtraits_df %>% ggplot(aes(x=block, y=tiller_num)) + geom_jitter(width=0.2) + geom_boxplot(alpha=0.5)
bunchtraits_df %>% ggplot(aes(x=block, y=basal_area)) + geom_jitter(width=0.2) + geom_boxplot(alpha=0.5)

leafandbunchtraits_df <-  bunchtraits_df %>%
  full_join(leaftraits_df, by=c("block", "plot", "plant_tag", "species")) %>%
  mutate(tiller_pack = tiller_diam*tiller_num/basal_area) %>%
  mutate(dom_tag = 1) %>%
  rename(leaf_area=total_area) %>%
  dplyr::select(block:species, dom_tag, SLA, leaf_area, leaf_thick, veg_height, LDMC, tiller_pack, tiller_diam)
### NOTE: Tag #80 and 129 are missing tiller diameter data on raw data sheet. We will keep this in there for now, but they will have to be removed for any multivariate analyses

## quick look at some of these values
leafandbunchtraits_df %>%
  ggplot(aes(species, tiller_pack, col=species)) + geom_jitter(width=0.1)

leafandbunchtraits_df %>%
  ggplot(aes(species, tiller_diam, col=species)) + geom_jitter(width=0.1) 


leaftraits_nondom <- read.csv("Traits\\NExS2023_traits_nondoms_cleaned_20250424.csv") %>%
  rename(plant_tag=tag_num, species=species_id) %>%
  group_by(block, plot, plant_tag, species) %>%
  summarize_at(vars(tiller_number:SLA), .funs=c(mean), na.rm=T) %>%
  ungroup() %>%
  mutate(dom_tag=0) %>%
  mutate(tiller_pack=NA) %>%
#  mutate(LDMC=LDMC*1000) %>% # fix unit inconsistancy
  rename(leaf_thick=thickness, tiller_diam=tiller_diameter, leaf_area=total_area) %>%
  dplyr::select(block:species, dom_tag, SLA, leaf_area, leaf_thick, veg_height, LDMC, tiller_pack, tiller_diam) %>%
  mutate(species = replace(species, species=="com_car", "com_pil")) %>% # KW NOTE: I think this species name was mistakenly put in as com_car (I think for commicarpus), but it should be com_pil for commicarpus pilosus
  mutate(species = replace(species, species=="evo_nut", "evo_als")) # KW NOTE: I think this species name was mistakenly labeled evo_nut (for the Konza species) -- it should be evo_als (evolvulus alsinoides)

leaftraits_dom_and_nondom <- leafandbunchtraits_df %>%
  bind_rows(leaftraits_nondom)
}

###
### Root traits (srl and root diameter)
###
{
tags_to_remove <- c(26,21,22,74,183) # these are tags that we are missing either scans or samples for

srl_df <- read.csv("..\\root_data\\NExS_root_traits_SRL_2023_cleaned_2.csv") %>%
  mutate(plant_tag=replace(plant_tag, plant_tag %in% c("22B","22A"),"22")) %>%
  rename(SRL=srl) %>%
  dplyr::select(block:species, SRL) %>%
  mutate(plant_tag=as.numeric(plant_tag)) %>%
  filter(!plant_tag %in% tags_to_remove) %>%
  mutate(species=replace(species, species=="ari con", "ari_con")) %>%
  mutate(species=replace(species, species=="the tri", "the_tri")) %>%
  mutate(species=replace(species, species=="bot rad", "bot_rad")) %>%
  mutate(species=replace(species, species=="pan col", "pan_col")) %>%
  mutate(species=replace(species, species=="pan max", "pan_max"))

  
rootDiam_df <- read.csv("..\\root_data\\NExS_root_traits_extra_2023_cleaned.csv") %>%
  dplyr::select(block:species, Average.Diameter.mm) %>%
  group_by(block, plot, plant_tag, species) %>%
  summarize(Average.Diameter.mm = mean(Average.Diameter.mm,na.rm=T)) %>%
  ungroup() %>%
  rename(root_diam = Average.Diameter.mm) %>%
  mutate(species=replace(species, species%in%c("aricon","ari con"), "ari_con")) %>%
  mutate(species=replace(species, species%in%c("thetri","the tri"), "the_tri")) %>%
  mutate(species=replace(species, species%in%c("botrad","bot rad"), "bot_rad")) %>%
  mutate(species=replace(species, species%in%c("pancol","pan col"), "pan_col")) %>%
  mutate(species=replace(species, species%in%c("panmax","pan max"), "pan_max")) %>%
  filter(!plant_tag %in% tags_to_remove)

roottraits_df <- srl_df %>%
  full_join(rootDiam_df, by=c("block", "plot", "plant_tag", "species")) %>%
#  mutate(dom_tag=1) %>%
  dplyr::select(block:species, SRL, root_diam)

roottraits_df %>% ggplot(aes(species, SRL)) + geom_jitter(width=0.1)
roottraits_df %>% ggplot(aes(species, root_diam)) + geom_jitter(width=0.1) ## three high values that might be worth looking into

#write.csv(roots_df, file="Roots\\NExS_srl and diamter traits_working_2023.csv", row.names=F)

}

###
### Ecophysiological traits
###
{
  ecophys_c3nondoms <- read.csv("..//02_Ecophysiology//04_data_clean//NExS_ACi_C3nondoms_outliers_marked.csv") %>%
    filter(fgroup=="Forb") %>%
    filter(species %in% c("COMPIL","EHRRIG","EVOALS","HELSTU","TEPPUR")) %>%
    mutate(species=replace(species, species=="COMPIL","com_pil")) %>%
    mutate(species=replace(species, species=="EHRRIG","ehr_rig")) %>%
    mutate(species=replace(species, species=="EVOALS","evo_als")) %>%
    mutate(species=replace(species, species=="HELSTU","hel_stu")) %>%
    mutate(species=replace(species, species=="TEPPUR","tep_pur")) %>%
    mutate(block="outside",
           plot=999,
           plant_tag=1:nrow(.)+10000) %>%
    dplyr::select(block, plot, species, plant_tag, Jmax_at_25, Vcmax_at_25) %>%
    rename(Jmax=Jmax_at_25, Vcmax=Vcmax_at_25)

  ecophys_c4nondoms <- read.csv("..//02_Ecophysiology//04_data_clean//NExS_ACi_C4nondoms_outliers_marked.csv") %>%
    mutate(species=replace(species, species=="ERASUP","era_sup")) %>%
    mutate(species=replace(species, species=="HETCON","het_con")) %>%
    mutate(species=replace(species, species=="UROMOS","uro_mos")) %>%
    mutate(block="outside",
           plot=999,
           plant_tag=1:nrow(.)+10100) %>%
    dplyr::select(block, plot, species, plant_tag, Light_Jmax, Rubisco_Vcmax) %>%
    rename(Jmax=Light_Jmax, Vcmax=Rubisco_Vcmax)

  ecophys_doms <- read.csv("..//02_Ecophysiology//04_data_clean//NExS_ACi_C4doms_pretreat_outliers_marked.csv") %>%
    mutate(species=replace(species, species=="ARICON","ari_con")) %>%
    mutate(species=replace(species, species=="BOTRAD","bot_rad")) %>%
    mutate(species=replace(species, species=="PANCOL","pan_col")) %>%
    mutate(species=replace(species, species=="PANMAX","pan_max")) %>%
    mutate(species=replace(species, species=="THETRI","the_tri")) %>%
    filter(Outlier==FALSE) %>% # A4 tag 30, D29 tag 155, and F48 tag 250 are all identified as outliers, so I remove these here
    dplyr::select(block, plot, species, plant_tag, Light_Jmax, Rubisco_Vcmax) %>%
    rename(Jmax=Light_Jmax, Vcmax=Rubisco_Vcmax)
  
  ecophystraits_df <- ecophys_c3nondoms %>%
    bind_rows(
      ecophys_c4nondoms,
      ecophys_doms
    ) %>%
    filter(!plant_tag %in% c(37,10112)) # both of these have very large Jmax values, likely outliers
  
  ecophystraits_df %>% 
    ggplot(aes(species, Jmax, label=plant_tag)) + 
    geom_text() +
    geom_boxplot(alpha=0.5)
  
  ecophystraits_df %>% 
    ggplot(aes(species, Vcmax, label=plant_tag)) + 
    geom_text() +
    geom_boxplot(alpha=0.5)
  
 }

###
### Bring all traits together
###
{
  # library(corrplot)
  # library(ggfortify)
  
  ### species key with functional groups, life history, and photosynthetic pathway
  species_key <- data.frame(
    species = unique(alltraits_long$species),
    fxn_grp = c("G","G","G","G","G","W","G","F","G",
                "G","G","F","F","F","F","F","F","F",
                "F","G","F"),
    life_hist = c("A/B","P","P","P","P","P","P","P","P",
                  "P","P","P","P","P","A/B","P","P","P",
                  "A/B","P","P"),
    ps_path = c("C4","C4","C4","C4","C4","C3","C4","C3","C4", # largely redundant with fxn group
                "C4","C4","C3","C3","C3","C3","C3","C3","C3",
                "C3","C4","C3")
  )
  
alltraits_df <- leaftraits_dom_and_nondom %>%
  full_join(roottraits_df, by=c("block", "plot", "plant_tag", "species")) %>%
  full_join(ecophystraits_df, by=c("block", "plot", "plant_tag", "species")) %>%
  full_join(species_key, by="species")

with(alltraits_df, table(species))
corrplot(trait_df[,6:16])

# Take a quick look at the biplots
ggpairs(alltraits_df, columns = 6:16,
        ggplot2::aes(col=species))

# put trait df into long format and calculate means
alltraits_long <- alltraits_df %>%
  pivot_longer(cols=SLA:Vcmax, names_to="trait_name", values_to="trait_value")

alltraits_long %>%
  ggplot(aes(x=species, y=trait_value, col=species)) +
  geom_jitter(width=0.1) + facet_wrap(~trait_name, scales="free_y")
}

###
### Transformations for trait values
###
{
  # SLA
  alltraits_long %>%
    filter(trait_name=="SLA") %>%
    ggplot(aes(x=sqrt(trait_value))) + geom_histogram()

    qqPlot(subset(alltraits_long, trait_name=="SLA")$trait_value)
    qqPlot(log(subset(alltraits_long, trait_name=="SLA")$trait_value))
    qqPlot(log10(subset(alltraits_long, trait_name=="SLA")$trait_value))
    qqPlot(subset(alltraits_long, trait_name=="SLA")$trait_value^(1/2))
    qqPlot(subset(alltraits_long, trait_name=="SLA")$trait_value^(1/3))
    qqPlot(subset(alltraits_long, trait_name=="SLA")$trait_value^(1/4))
    
    shapiro.test(subset(alltraits_long, trait_name=="SLA")$trait_value)
    shapiro.test(log(subset(alltraits_long, trait_name=="SLA")$trait_value))
    shapiro.test(log10(subset(alltraits_long, trait_name=="SLA")$trait_value))
    shapiro.test(subset(alltraits_long, trait_name=="SLA")$trait_value^(1/2))
    shapiro.test(subset(alltraits_long, trait_name=="SLA")$trait_value^(1/3))
    shapiro.test(subset(alltraits_long, trait_name=="SLA")$trait_value^(1/4))
    shapiro.test(log10(subset(alltraits_long, trait_name=="SLA")$trait_value))
  ## Square root transform is best
    # leaf_area
    alltraits_long %>% filter(trait_name=="leaf_area") %>% ggplot(aes(x=trait_value)) + geom_histogram()
    alltraits_long %>% filter(trait_name=="leaf_area") %>% ggplot(aes(x=log(trait_value))) + geom_histogram()
    alltraits_long %>% filter(trait_name=="leaf_area") %>% ggplot(aes(x=log10(trait_value))) + geom_histogram()
    alltraits_long %>% filter(trait_name=="leaf_area") %>% ggplot(aes(x=trait_value^(1/2))) + geom_histogram()
    alltraits_long %>% filter(trait_name=="leaf_area") %>% ggplot(aes(x=trait_value^(1/3))) + geom_histogram()
    alltraits_long %>% filter(trait_name=="leaf_area") %>% ggplot(aes(x=trait_value^(1/4))) + geom_histogram()
    
    qqPlot(subset(alltraits_long, trait_name=="leaf_area")$trait_value)
    qqPlot(log(subset(alltraits_long, trait_name=="leaf_area")$trait_value))
    qqPlot(log10(subset(alltraits_long, trait_name=="leaf_area")$trait_value))
    qqPlot(subset(alltraits_long, trait_name=="leaf_area")$trait_value^(1/2))
    qqPlot(subset(alltraits_long, trait_name=="leaf_area")$trait_value^(1/3))
    qqPlot(subset(alltraits_long, trait_name=="leaf_area")$trait_value^(1/4))
    
    shapiro.test(subset(alltraits_long, trait_name=="leaf_area")$trait_value)
    shapiro.test(log(subset(alltraits_long, trait_name=="leaf_area")$trait_value))
    shapiro.test(log10(subset(alltraits_long, trait_name=="leaf_area")$trait_value))
    shapiro.test(subset(alltraits_long, trait_name=="leaf_area")$trait_value^(1/2))
    shapiro.test(subset(alltraits_long, trait_name=="leaf_area")$trait_value^(1/3))
    shapiro.test(subset(alltraits_long, trait_name=="leaf_area")$trait_value^(1/4))
    # Try Box-Cox transformation
    leaf_area_vec <- alltraits_long %>% filter(trait_name=="leaf_area") %>% pull(trait_value)
    leaf_area_bc <- boxcox(leaf_area_vec ~ 1, lambda = seq(-2,2,.1))
    leaf_area_lambda <- leaf_area_bc$x[which.max(leaf_area_bc$y)] ## optimal lambda is -0.7070707
    
    if (leaf_area_lambda == 0) {
      leaf_area_boxcox <- log(leaf_area)
    } else {
      leaf_area_boxcox <- (leaf_area_vec^leaf_area_lambda - 1) / leaf_area_lambda
    }
    hist(leaf_area_boxcox)
    qqPlot(leaf_area_boxcox)
    shapiro.test(leaf_area_boxcox)
    ## quadroot (^1/4) is best
    
    # leaf_thick
    alltraits_long %>% filter(trait_name=="leaf_thick") %>% ggplot(aes(x=trait_value)) + geom_histogram()
    alltraits_long %>% filter(trait_name=="leaf_thick") %>% ggplot(aes(x=log(trait_value))) + geom_histogram()
    alltraits_long %>% filter(trait_name=="leaf_thick") %>% ggplot(aes(x=log10(trait_value))) + geom_histogram()
    alltraits_long %>% filter(trait_name=="leaf_thick") %>% ggplot(aes(x=trait_value^(1/2))) + geom_histogram()
    alltraits_long %>% filter(trait_name=="leaf_thick") %>% ggplot(aes(x=trait_value^(1/3))) + geom_histogram()
    alltraits_long %>% filter(trait_name=="leaf_thick") %>% ggplot(aes(x=trait_value^(1/4))) + geom_histogram()

    qqPlot(subset(alltraits_long, trait_name=="leaf_thick")$trait_value)
    qqPlot(log(subset(alltraits_long, trait_name=="leaf_thick")$trait_value))
    qqPlot(log10(subset(alltraits_long, trait_name=="leaf_thick")$trait_value))
    qqPlot(subset(alltraits_long, trait_name=="leaf_thick")$trait_value^(1/2))
    qqPlot(subset(alltraits_long, trait_name=="leaf_thick")$trait_value^(1/3))
    qqPlot(subset(alltraits_long, trait_name=="leaf_thick")$trait_value^(1/4))
    
    shapiro.test(subset(alltraits_long, trait_name=="leaf_thick")$trait_value)
    shapiro.test(log(subset(alltraits_long, trait_name=="leaf_thick")$trait_value))
    shapiro.test(log10(subset(alltraits_long, trait_name=="leaf_thick")$trait_value))
    shapiro.test(subset(alltraits_long, trait_name=="leaf_thick")$trait_value^(1/2))
    shapiro.test(subset(alltraits_long, trait_name=="leaf_thick")$trait_value^(1/3))
    shapiro.test(subset(alltraits_long, trait_name=="leaf_thick")$trait_value^(1/4))
    shapiro.test(log10(subset(alltraits_long, trait_name=="leaf_thick")$trait_value))
    
    # Try Box-Cox transformation
    leaf_thick_vec <- alltraits_long %>% filter(trait_name=="leaf_thick") %>% pull(trait_value)
    leaf_thick_bc <- boxcox(leaf_thick_vec ~ 1, lambda = seq(-2,2,.1))
    leaf_thick_lambda <- leaf_thick_bc$x[which.max(leaf_thick_bc$y)] ## optimal lambda is -0.7070707
    
    if (leaf_thick_lambda == 0) {
      leaf_thick_boxcox <- log(leaf_thick)
    } else {
      leaf_thick_boxcox <- (leaf_thick_vec^leaf_thick_lambda - 1) / leaf_thick_lambda
    }
    hist(leaf_thick_boxcox)
    qqPlot(leaf_thick_boxcox)
    shapiro.test(leaf_thick_boxcox)
    
    ## Box-Cox transform is best
    
    # veg_height
    alltraits_long %>% filter(trait_name=="veg_height") %>% ggplot(aes(x=trait_value)) + geom_histogram()
    alltraits_long %>% filter(trait_name=="veg_height") %>% ggplot(aes(x=log(trait_value))) + geom_histogram()
    alltraits_long %>% filter(trait_name=="veg_height") %>% ggplot(aes(x=log10(trait_value))) + geom_histogram()
    alltraits_long %>% filter(trait_name=="veg_height") %>% ggplot(aes(x=trait_value^(1/2))) + geom_histogram()
    alltraits_long %>% filter(trait_name=="veg_height") %>% ggplot(aes(x=trait_value^(1/3))) + geom_histogram()
    alltraits_long %>% filter(trait_name=="veg_height") %>% ggplot(aes(x=trait_value^(1/4))) + geom_histogram()
    
    qqPlot(subset(alltraits_long, trait_name=="veg_height")$trait_value)
    qqPlot(log(subset(alltraits_long, trait_name=="veg_height")$trait_value))
    qqPlot(log10(subset(alltraits_long, trait_name=="veg_height")$trait_value))
    qqPlot(subset(alltraits_long, trait_name=="veg_height")$trait_value^(1/2))
    qqPlot(subset(alltraits_long, trait_name=="veg_height")$trait_value^(1/3))
    qqPlot(subset(alltraits_long, trait_name=="veg_height")$trait_value^(1/4))
    
    shapiro.test(subset(alltraits_long, trait_name=="veg_height")$trait_value)
    shapiro.test(log(subset(alltraits_long, trait_name=="veg_height")$trait_value))
    shapiro.test(log10(subset(alltraits_long, trait_name=="veg_height")$trait_value))
    shapiro.test(subset(alltraits_long, trait_name=="veg_height")$trait_value^(1/2))
    shapiro.test(subset(alltraits_long, trait_name=="veg_height")$trait_value^(1/3))
    shapiro.test(subset(alltraits_long, trait_name=="veg_height")$trait_value^(1/4))
    shapiro.test(log10(subset(alltraits_long, trait_name=="veg_height")$trait_value))
    ## Square root transformation is best
    
    # LDMC
    alltraits_long %>% filter(trait_name=="LDMC") %>% ggplot(aes(x=trait_value)) + geom_histogram()
    alltraits_long %>% filter(trait_name=="LDMC") %>% ggplot(aes(x=log(trait_value))) + geom_histogram()
    alltraits_long %>% filter(trait_name=="LDMC") %>% ggplot(aes(x=log10(trait_value))) + geom_histogram()
    alltraits_long %>% filter(trait_name=="LDMC") %>% ggplot(aes(x=trait_value^(1/2))) + geom_histogram()
    alltraits_long %>% filter(trait_name=="LDMC") %>% ggplot(aes(x=trait_value^(1/3))) + geom_histogram()
    alltraits_long %>% filter(trait_name=="LDMC") %>% ggplot(aes(x=trait_value^(1/4))) + geom_histogram()
    alltraits_long %>% filter(trait_name=="LDMC") %>% ggplot(aes(x=asin(sqrt(trait_value)))) + geom_histogram()
    
    qqPlot(subset(alltraits_long, trait_name=="LDMC")$trait_value)
    qqPlot(log(subset(alltraits_long, trait_name=="LDMC")$trait_value))
    qqPlot(log10(subset(alltraits_long, trait_name=="LDMC")$trait_value))
    qqPlot(subset(alltraits_long, trait_name=="LDMC")$trait_value^(1/2))
    qqPlot(subset(alltraits_long, trait_name=="LDMC")$trait_value^(1/3))
    qqPlot(subset(alltraits_long, trait_name=="LDMC")$trait_value^(1/4))
    qqPlot(asin(sqrt(subset(alltraits_long, trait_name=="LDMC")$trait_value))[-37])
    
    shapiro.test(subset(alltraits_long, trait_name=="LDMC")$trait_value)
    shapiro.test(log(subset(alltraits_long, trait_name=="LDMC")$trait_value))
    shapiro.test(log10(subset(alltraits_long, trait_name=="LDMC")$trait_value))
    shapiro.test(subset(alltraits_long, trait_name=="LDMC")$trait_value^(1/2))
    shapiro.test(subset(alltraits_long, trait_name=="LDMC")$trait_value^(1/3))
    shapiro.test(subset(alltraits_long, trait_name=="LDMC")$trait_value^(1/4))
    shapiro.test(asin(sqrt(subset(alltraits_long, trait_name=="LDMC")$trait_value))[-37])
    ## untransformed is best -- and it passes shapiro wilks test when row 37 is removed
    
    # tiller_pack
    alltraits_long %>% filter(trait_name=="tiller_pack") %>% ggplot(aes(x=trait_value)) + geom_histogram()
    alltraits_long %>% filter(trait_name=="tiller_pack") %>% ggplot(aes(x=log(trait_value))) + geom_histogram()
    alltraits_long %>% filter(trait_name=="tiller_pack") %>% ggplot(aes(x=log10(trait_value))) + geom_histogram()
    alltraits_long %>% filter(trait_name=="tiller_pack") %>% ggplot(aes(x=trait_value^(1/2))) + geom_histogram()
    alltraits_long %>% filter(trait_name=="tiller_pack") %>% ggplot(aes(x=trait_value^(1/3))) + geom_histogram()
    alltraits_long %>% filter(trait_name=="tiller_pack") %>% ggplot(aes(x=trait_value^(1/4))) + geom_histogram()
    
    qqPlot(subset(alltraits_long, trait_name=="tiller_pack")$trait_value)
    qqPlot(log(subset(alltraits_long, trait_name=="tiller_pack")$trait_value))
    qqPlot(log10(subset(alltraits_long, trait_name=="tiller_pack")$trait_value))
    qqPlot(subset(alltraits_long, trait_name=="tiller_pack")$trait_value^(1/2))
    qqPlot(subset(alltraits_long, trait_name=="tiller_pack")$trait_value^(1/3))
    qqPlot(subset(alltraits_long, trait_name=="tiller_pack")$trait_value^(1/4))
    
    shapiro.test(subset(alltraits_long, trait_name=="tiller_pack")$trait_value)
    shapiro.test(log(subset(alltraits_long, trait_name=="tiller_pack")$trait_value))
    shapiro.test(log10(subset(alltraits_long, trait_name=="tiller_pack")$trait_value))
    shapiro.test(subset(alltraits_long, trait_name=="tiller_pack")$trait_value^(1/2))
    shapiro.test(subset(alltraits_long, trait_name=="tiller_pack")$trait_value^(1/3))
    shapiro.test(subset(alltraits_long, trait_name=="tiller_pack")$trait_value^(1/4))
    # natural log transformation is best
    
    # tiller_diam
    alltraits_long %>% filter(trait_name=="tiller_diam"&trait_value<0.5) %>% ggplot(aes(x=trait_value)) + geom_histogram()
    alltraits_long %>% filter(trait_name=="tiller_diam"&trait_value<0.5) %>% ggplot(aes(x=log(trait_value))) + geom_histogram()
    alltraits_long %>% filter(trait_name=="tiller_diam"&trait_value<0.5) %>% ggplot(aes(x=log10(trait_value))) + geom_histogram()
    alltraits_long %>% filter(trait_name=="tiller_diam"&trait_value<0.5) %>% ggplot(aes(x=trait_value^(1/2))) + geom_histogram()
    alltraits_long %>% filter(trait_name=="tiller_diam"&trait_value<0.5) %>% ggplot(aes(x=trait_value^(1/3))) + geom_histogram()
    alltraits_long %>% filter(trait_name=="tiller_diam"&trait_value<0.5) %>% ggplot(aes(x=trait_value^(1/4))) + geom_histogram()
    alltraits_long %>% filter(trait_name=="tiller_diam"&trait_value<0.5) %>% ggplot(aes(x=asin(trait_value^(1/2)))) + geom_histogram()
    
    qqPlot(subset(alltraits_long, trait_name=="tiller_diam")$trait_value[-c(312,307)])
    qqPlot(log(subset(alltraits_long, trait_name=="tiller_diam")$trait_value[-c(312,307)]))
    qqPlot(log10(subset(alltraits_long, trait_name=="tiller_diam")$trait_value[-c(312,307)]))
    qqPlot(subset(alltraits_long, trait_name=="tiller_diam")$trait_value[-c(312,307)]^(1/2))
    qqPlot(subset(alltraits_long, trait_name=="tiller_diam")$trait_value[-c(312,307)]^(1/3))
    qqPlot(subset(alltraits_long, trait_name=="tiller_diam")$trait_value[-c(312,307)]^(1/4))
    qqPlot(asin(subset(alltraits_long, trait_name=="tiller_diam")$trait_value[-c(312,307)]^(1/2)))
    
    shapiro.test(subset(alltraits_long, trait_name=="tiller_diam")$trait_value[-c(312,307)])
    shapiro.test(log(subset(alltraits_long, trait_name=="tiller_diam")$trait_value[-c(312,307)]))
    shapiro.test(log10(subset(alltraits_long, trait_name=="tiller_diam")$trait_value[-c(312,307)]))
    shapiro.test(subset(alltraits_long, trait_name=="tiller_diam")$trait_value[-c(312,307)]^(1/2))
    shapiro.test(subset(alltraits_long, trait_name=="tiller_diam")$trait_value[-c(312,307)]^(1/3))
    shapiro.test(subset(alltraits_long, trait_name=="tiller_diam")$trait_value[-c(312,307)]^(1/4))
    shapiro.test(asin(subset(alltraits_long, trait_name=="tiller_diam")$trait_value[-c(312,307)]^(1/2)))
    # Try Box-Cox transformation
    tiller_diam_vec <- alltraits_long %>% filter(trait_name=="tiller_diam") %>% pull(trait_value)
    tiller_diam_bc <- boxcox(tiller_diam_vec ~ 1, lambda = seq(-2,2,.1))
    tiller_diam_lambda <- tiller_diam_bc$x[which.max(tiller_diam_bc$y)] ## optimal lambda is -0.7070707
    
    if (tiller_diam_lambda == 0) {
      tiller_diam_boxcox <- log(tiller_diam)
    } else {
      tiller_diam_boxcox <- (tiller_diam_vec^tiller_diam_lambda - 1) / tiller_diam_lambda
    }
    hist(tiller_diam_boxcox)
    qqPlot(tiller_diam_boxcox)
    shapiro.test(tiller_diam_boxcox)
    ## I can't really get this to be normal but log transformed is the best-- too many 0.1 values, so maybe we remove it?
    
    # SRL
    alltraits_long %>% filter(trait_name=="SRL") %>% ggplot(aes(x=trait_value)) + geom_histogram()
    alltraits_long %>% filter(trait_name=="SRL") %>% ggplot(aes(x=log(trait_value))) + geom_histogram()
    alltraits_long %>% filter(trait_name=="SRL") %>% ggplot(aes(x=log10(trait_value))) + geom_histogram()
    alltraits_long %>% filter(trait_name=="SRL") %>% ggplot(aes(x=trait_value^(1/2))) + geom_histogram()
    alltraits_long %>% filter(trait_name=="SRL") %>% ggplot(aes(x=trait_value^(1/3))) + geom_histogram()
    alltraits_long %>% filter(trait_name=="SRL") %>% ggplot(aes(x=trait_value^(1/4))) + geom_histogram()
    
    qqPlot(subset(alltraits_long, trait_name=="SRL")$trait_value)
    qqPlot(log(subset(alltraits_long, trait_name=="SRL")$trait_value))
    qqPlot(log10(subset(alltraits_long, trait_name=="SRL")$trait_value))
    qqPlot(subset(alltraits_long, trait_name=="SRL")$trait_value^(1/2))
    qqPlot(subset(alltraits_long, trait_name=="SRL")$trait_value^(1/3))
    qqPlot(subset(alltraits_long, trait_name=="SRL")$trait_value^(1/4))
    
    shapiro.test(subset(alltraits_long, trait_name=="SRL")$trait_value)
    shapiro.test(log(subset(alltraits_long, trait_name=="SRL")$trait_value))
    shapiro.test(log10(subset(alltraits_long, trait_name=="SRL")$trait_value))
    shapiro.test(subset(alltraits_long, trait_name=="SRL")$trait_value^(1/2))
    shapiro.test(subset(alltraits_long, trait_name=="SRL")$trait_value^(1/3))
    shapiro.test(subset(alltraits_long, trait_name=="SRL")$trait_value^(1/4))
    shapiro.test(log10(subset(alltraits_long, trait_name=="SRL")$trait_value))
    ## Square root transformation is best
    
    # root_diam
    alltraits_long %>% filter(trait_name=="root_diam") %>% ggplot(aes(x=trait_value)) + geom_histogram()
    alltraits_long %>% filter(trait_name=="root_diam") %>% ggplot(aes(x=log(trait_value))) + geom_histogram()
    alltraits_long %>% filter(trait_name=="root_diam") %>% ggplot(aes(x=log10(trait_value))) + geom_histogram()
    alltraits_long %>% filter(trait_name=="root_diam") %>% ggplot(aes(x=trait_value^(1/2))) + geom_histogram()
    alltraits_long %>% filter(trait_name=="root_diam") %>% ggplot(aes(x=trait_value^(1/3))) + geom_histogram()
    alltraits_long %>% filter(trait_name=="root_diam") %>% ggplot(aes(x=trait_value^(1/4))) + geom_histogram()
    
    qqPlot(subset(alltraits_long, trait_name=="root_diam")$trait_value)
    qqPlot(log(subset(alltraits_long, trait_name=="root_diam")$trait_value[-c(124,111,82)]))
    qqPlot(log10(subset(alltraits_long, trait_name=="root_diam")$trait_value))
    qqPlot(subset(alltraits_long, trait_name=="root_diam")$trait_value^(1/2))
    qqPlot(subset(alltraits_long, trait_name=="root_diam")$trait_value^(1/3))
    qqPlot(subset(alltraits_long, trait_name=="root_diam")$trait_value^(1/4))
    
    shapiro.test(subset(alltraits_long, trait_name=="root_diam")$trait_value)
    shapiro.test(log(subset(alltraits_long, trait_name=="root_diam")$trait_value[-c(124,111,82)]))
    shapiro.test(log10(subset(alltraits_long, trait_name=="root_diam")$trait_value))
    shapiro.test(subset(alltraits_long, trait_name=="root_diam")$trait_value^(1/2))
    shapiro.test(subset(alltraits_long, trait_name=="root_diam")$trait_value^(1/3))
    shapiro.test(subset(alltraits_long, trait_name=="root_diam")$trait_value^(1/4))
    shapiro.test(log10(subset(alltraits_long, trait_name=="root_diam")$trait_value))
    ## Log transformed is best -- potential outliers rows 124, 111, 82    
    # Jmax
    alltraits_long %>% filter(trait_name=="Jmax") %>% ggplot(aes(x=trait_value)) + geom_histogram()
    alltraits_long %>% filter(trait_name=="Jmax") %>% ggplot(aes(x=log(trait_value))) + geom_histogram()
    alltraits_long %>% filter(trait_name=="Jmax") %>% ggplot(aes(x=log10(trait_value))) + geom_histogram()
    alltraits_long %>% filter(trait_name=="Jmax") %>% ggplot(aes(x=trait_value^(1/2))) + geom_histogram()
    alltraits_long %>% filter(trait_name=="Jmax") %>% ggplot(aes(x=trait_value^(1/3))) + geom_histogram()
    alltraits_long %>% filter(trait_name=="Jmax") %>% ggplot(aes(x=trait_value^(1/4))) + geom_histogram()
    
    qqPlot(subset(alltraits_long, trait_name=="Jmax")$trait_value)
    qqPlot(log(subset(alltraits_long, trait_name=="Jmax")$trait_value))
    qqPlot(log10(subset(alltraits_long, trait_name=="Jmax")$trait_value))
    qqPlot(subset(alltraits_long, trait_name=="Jmax")$trait_value^(1/2))
    qqPlot(subset(alltraits_long, trait_name=="Jmax")$trait_value^(1/3))
    qqPlot(subset(alltraits_long, trait_name=="Jmax")$trait_value^(1/4))
    
    shapiro.test(subset(alltraits_long, trait_name=="Jmax")$trait_value)
    shapiro.test(log(subset(alltraits_long, trait_name=="Jmax")$trait_value))
    shapiro.test(log10(subset(alltraits_long, trait_name=="Jmax")$trait_value))
    shapiro.test(subset(alltraits_long, trait_name=="Jmax")$trait_value^(1/2))
    shapiro.test(subset(alltraits_long, trait_name=="Jmax")$trait_value^(1/3))
    shapiro.test(subset(alltraits_long, trait_name=="Jmax")$trait_value^(1/4))
    # Log transformation is best
    
    # Vcmax
    alltraits_long %>% filter(trait_name=="Vcmax") %>% ggplot(aes(x=trait_value)) + geom_histogram()
    alltraits_long %>% filter(trait_name=="Vcmax") %>% ggplot(aes(x=log(trait_value))) + geom_histogram()
    alltraits_long %>% filter(trait_name=="Vcmax") %>% ggplot(aes(x=log10(trait_value))) + geom_histogram()
    alltraits_long %>% filter(trait_name=="Vcmax") %>% ggplot(aes(x=trait_value^(1/2))) + geom_histogram()
    alltraits_long %>% filter(trait_name=="Vcmax") %>% ggplot(aes(x=trait_value^(1/3))) + geom_histogram()
    alltraits_long %>% filter(trait_name=="Vcmax") %>% ggplot(aes(x=trait_value^(1/4))) + geom_histogram()
    
    qqPlot(subset(alltraits_long, trait_name=="Vcmax")$trait_value)
    qqPlot(log(subset(alltraits_long, trait_name=="Vcmax")$trait_value))
    qqPlot(log10(subset(alltraits_long, trait_name=="Vcmax")$trait_value))
    qqPlot(subset(alltraits_long, trait_name=="Vcmax")$trait_value^(1/2))
    qqPlot(subset(alltraits_long, trait_name=="Vcmax")$trait_value^(1/3))
    qqPlot(subset(alltraits_long, trait_name=="Vcmax")$trait_value^(1/4))
    
    shapiro.test(subset(alltraits_long, trait_name=="Vcmax")$trait_value)
    shapiro.test(log(subset(alltraits_long, trait_name=="Vcmax")$trait_value))
    shapiro.test(log10(subset(alltraits_long, trait_name=="Vcmax")$trait_value))
    shapiro.test(subset(alltraits_long, trait_name=="Vcmax")$trait_value^(1/2))
    shapiro.test(subset(alltraits_long, trait_name=="Vcmax")$trait_value^(1/3))
    shapiro.test(subset(alltraits_long, trait_name=="Vcmax")$trait_value^(1/4))
    ## cube or quad-root is the best... I'll stick with cube-root because it's more common
    
    
    ##
    ## Perform transformations
    alltraits_long <- alltraits_long %>%
      mutate(trait_value_tf = case_when(
               trait_name=="SLA" ~           sqrt(trait_value),
               trait_name=="leaf_area" ~     trait_value^(1/4),
               trait_name=="leaf_thick" ~    (trait_value^(-.70707)-1)/(-.70707),
               trait_name=="veg_height" ~    sqrt(trait_value),             
               trait_name=="LDMC" ~          trait_value,             
               trait_name=="tiller_pack" ~   log(trait_value),
               trait_name=="tiller_diam" ~   log(trait_value),
               trait_name=="SRL" ~           sqrt(trait_value),             
               trait_name=="root_diam" ~     log(trait_value),
               trait_name=="Jmax" ~          log(trait_value),
               trait_name=="Vcmax" ~         trait_value^(1/3),
               TRUE ~ trait_value
            ),
            tf_type = case_when(
              trait_name=="SLA" ~           "sqrt",
              trait_name=="leaf_area"  ~    "1/4root",
              trait_name=="leaf_thick" ~    "boxcox",
              trait_name=="veg_height" ~    "sqrt",             
              trait_name=="LDMC" ~          "none",             
              trait_name=="tiller_pack" ~   "ln",
              trait_name=="tiller_diam" ~   "ln",
              trait_name=="SRL" ~           "sqrt",             
              trait_name=="root_diam" ~     "ln",
              trait_name=="Jmax" ~          "ln",
              trait_name=="Vcmax" ~         "cuberoot",
              TRUE ~ NA          
          )
        )
       

}



###
### Plotting traits
###
{
alltraits_means <- alltraits_long %>%
  group_by(species, fxn_grp, life_hist, ps_path, trait_name) %>%
  summarize_at(.vars=vars(trait_value, trait_value_tf), .funs=c(mean=mean, se=sefxn), na.rm=T) %>%
  ungroup()

# plot up means and se's
# alltraits_means %>%
#   ggplot(aes(x=species, y=mean, ymin=mean-se, ymax=mean+se, fill=species)) +
#   geom_errorbar(width=0.1) + geom_point(pch=21, col="black") + facet_wrap(~trait_name, scales="free_y") +
#   theme(axis.text.x=element_text(angle=45, hjust=1))
  
unique(alltraits_means$trait_name)

alltraits_means_wide <- alltraits_means %>%
  pivot_wider(names_from=trait_name, values_from=c(trait_value_mean, 
                                                   trait_value_se, 
                                                   trait_value_tf_mean,
                                                   trait_value_tf_se))
names_temp <- c(colnames(alltraits_means_wide[1:4]),
                gsub("trait_value_","",colnames(alltraits_means_wide[5:ncol(alltraits_means_wide)])))

colnames(alltraits_means_wide) <- names_temp
rm(names_temp)
# alltraits_means_wide %>%
#   ggplot(aes(x=mean_Jmax, y=mean_LDMC, 
#              xmin=mean_Jmax-se_Jmax, xmax=mean_Jmax+se_Jmax, 
#              ymin=mean_LDMC-se_LDMC, ymax=mean_LDMC+se_LDMC)) +
#   geom_errorbar(width=0) +
#   geom_errorbarh(height=0) +
#   geom_point(size=2)

##
## Loop that goes through trait names and plots all biplots separately
species_n <- as.data.frame(with(alltraits_df, table(species)))
high_n_sp <- as.character(species_n %>% filter(Freq>10) %>% pull(species))
trait_name_vec <- unique(alltraits_long$trait_name)

for(TRAIT_X in 1:(length(trait_name_vec)-1)){
  for(TRAIT_Y in (TRAIT_X+1):length(trait_name_vec)){
    df_temp <- data.frame(
    x_mean = alltraits_means_wide %>% pull(paste0("tf_mean_",trait_name_vec[TRAIT_X])),
    x_se = alltraits_means_wide %>% pull(paste0("tf_se_",trait_name_vec[TRAIT_X])),
    y_mean = alltraits_means_wide %>% pull(paste0("tf_mean_",trait_name_vec[TRAIT_Y])),
    y_se = alltraits_means_wide %>% pull(paste0("tf_se_",trait_name_vec[TRAIT_Y])),
    species = alltraits_means_wide %>% pull(species),
    fxn_grp = alltraits_means_wide %>% pull(fxn_grp),
    life_hist = alltraits_means_wide %>% pull(life_hist),
    ps_path = alltraits_means_wide %>% pull(ps_path)
    ) %>%
    na.omit() %>%
    filter(fxn_grp=="G" & species %in% high_n_sp)

    print(
    ggplot(df_temp, aes(x=x_mean, y=y_mean,
               xmin=x_mean-x_se, xmax=x_mean+x_se, 
               ymin=y_mean-y_se, ymax=y_mean+y_se,
               fill=fxn_grp)) +
    geom_errorbar(width=0, linetype=1) +
    geom_errorbarh(height=0, linetype=1) +
    geom_point(pch=21, size=3) +
    geom_smooth(method="lm",se=F, aes(col=fxn_grp)) +
    scale_color_manual(values=c("coral","darkgreen")) +
    xlab(paste0(trait_name_vec[TRAIT_X],"_tf")) + ylab(paste0(trait_name_vec[TRAIT_Y],"_tf"))
    )
  }
}

alltraits_long %>% dplyr::select(trait_name, tf_type) %>% unique(.)
with(alltraits_long, unique(c(trait_name, tf_type)))
trait_names <- colnames(alltraits_means_wide[27:36])

label_map <- setNames(
  nm=trait_names,
  object = gsub("trait_value_transformed_mean_", "", trait_names)
)


ggpairs(alltraits_means_wide %>% filter(fxn_grp=="G" & species %in% high_n_sp),
        columns = 27:36,
        lower = list(continuous = wrap("smooth", method = "lm", se = T)),
        labeller = as_labeller(label_map))

ggpairs(alltraits_means_wide, columns = 2:11,
        ggplot2::aes(col=species))

ggplot(alltraits_means_wide, aes(trait_value_mean_Jmax, trait_value_mean_Vcmax, label=species)) +
  geom_point() + geom_text()
ggplot(alltraits_means_wide, aes(trait_value_transformed_mean_Jmax, trait_value_transformed_mean_Vcmax, label=species)) +
  geom_point() + geom_text()


#  [1] "Jmax"        "LDMC"        "SLA"         "SRL"         "Vcmax"       "leaf_thick" 
#  [7] "root_diam"   "tiller_diam" "tiller_pack" "veg_height"

}



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
