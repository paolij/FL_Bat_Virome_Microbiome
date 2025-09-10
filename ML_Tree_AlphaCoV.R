## R script for visualizing maximum likelihood tree 

# Load libraries
library(ggtree)
library(ggimage)
library(phytools)
library(ggplot2)
library(ape)
library(ggbreak)
library(dplyr)
library(readr)


# Set working directory
setwd("/Users/julia/Dropbox/MAVIAN_Lab/Manuscripts/Julia_FL_Bat/Figures_plus_data_scripts/AlphaCoV/Concatenated_tree/")

# Load the tree
ml_tree <- read.tree("Alphacov_concatenate_fixed.nwk")

# Load metadata and prepare image paths
metadata <- read_csv("alphacov_md.csv") %>%
  mutate(
    animal = tolower(animal),
    image = ifelse(is.na(animal) | animal == "", NA, paste0(animal, ".png"))
  )

# Get tip coordinates and join with metadata
tip_positions <- ggtree::fortify(ml_tree) %>% 
  filter(isTip) %>%
  select(label, x, y)

img_data <- left_join(tip_positions, metadata, by = c("label" = "full"))

#Finding Nodes
which(ml_tree$tip.label == "Ma_Frozen_3|AlphaCoV-Ma3-FL-2021|Myotis_austroriparius|USA")


# Visualize the rooted tree with bootstrap values as text labels
q <- ggtree(ml_tree, ladderize = TRUE) +
  geom_nodepoint(aes(subset = !is.na(as.numeric(label)) & as.numeric(label) > 90), 
                 size = 13, shape = 18, color = "black") +
  geom_treescale(x=0, y=25, color='black', fontsize = 15, linesize = 3) +
  geom_tiplab(aes(subset = node %in% c(42), label = label), size = 10, align = TRUE, offset = 0.15, fontface = 'bold', color = "black") +  
  geom_tiplab(aes(subset = !node %in% c(42), label = label), size = 10, align = TRUE, offset = 0.15) + # Regular for others
  coord_cartesian(xlim = c(0, max(nodeHeights(ml_tree)) * 2.0))
q

# Get max x position of tips
tip_x_max <- max(img_data$x, na.rm = TRUE)

# Compute aligned x-position (same as geom_tiplab offset)
aligned_x <- tip_x_max + 0.06

# Add animal icons aligned to labels
q <- q + geom_image(
  data = img_data,
  aes(x = aligned_x, y = y, image = image),  # fixed x for alignment
  size = 0.03,
  inherit.aes = FALSE
)


# Highlight a clade by MRCA
TADBR_ACoV2 <- getMRCA(ml_tree, c("Ma_Frozen_3|AlphaCoV-Ma3-FL-2021|Myotis_austroriparius|USA", 
                                  "OP715780.1|Tadarida_brasiliensis_bat_alphacoronavirus_2|Tadarida_brasiliensis|Argentina"))

Luchacovirus <- getMRCA(ml_tree, c("NC_032730.1|Lucheng_rat_coronavirus|Rattus_norvegicus|China", 
                                   "NC_034972.1|Coronavirus_AcCoV-JC34|Apodemus_chevrieri|China"))

Sunacovirus <- getMRCA(ml_tree, c("NC_035191.1|Wencheng_shrew_coronavirus|Suncus_murinus|China", 
                                  "NC_048211.1|Wencheng_shrew_coronavirus|Suncus_murinus|China"))

Rhinacovirus <- getMRCA(ml_tree, c("NC_009988.1|Rhinolophus_bat_coronavirus_HKU2|Rhinolophus_sinicus|China", 
                                   "NC_028824.1|BtRf-AlphaCoV/YN2012|Rhinolophus_ferrumequinum|China"))

Tegacovirus <- getMRCA(ml_tree, c("KC175340.1|Canine_coronavirus|Canis_lupus_familiaris|USA", 
                                  "DQ848678.1|Feline_coronavirus|Felis_catus|United_Kingdom"))


Minacovirus <- getMRCA(ml_tree, c("NC_030292.1|Ferret_coronavirus|Mustela_putorius|Netherlands",
                                  "NC_023760.1|Mink_coronavirus_strain_WD1127|Neogale_vison|USA"))


Setracovirus <- getMRCA(ml_tree, c("NC_032107.1|NL63-related_bat_coronavirus|Triaenops_afer|Kenya",
                                   "NC_048216.1|NL63-related_bat_coronavirus|Triaenops_afer|Kenya"))


Duvinacovirus <- getMRCA(ml_tree, c("NC_028752.1|Camel_alphacoronavirus|Camelus_dromedarius|Saudi_Arabia", 
                                    "NC_002645.1|Human_coronavirus_229E|Homo_sapiens|NA"))


Minunacovirus <- getMRCA(ml_tree, c("NC_010437.1|Bat_coronavirus_1A|Miniopterus_sp|Hong_Kong", 
                                    "NC_010438.1|Miniopterus_bat_coronavirus_HKU8|Miniopterus_sp|Hong_Kong"))

Nyctacovirus <- getMRCA(ml_tree, c("NC_046964.1|Alphacoronavirus_Bat-CoV/P.kuhlii/Italy/3398-19/2015|Pipistrellus_kuhlii|Italy", 
                                   "MK720944.1|Tylonycteris_bat_coronavirus_HKU33|Tylonycteris_robustula|China"))


Decacovirus <- getMRCA(ml_tree, c("NC_028814.1|BtRf-AlphaCoV/HuB2013|Rhinolophus_ferrumequinum|China", 
                                  "NC_018871.1|Rousettus_bat_coronavirus_HKU10|Rousettus_leschenaulti|China"))


Pedacovirus <- getMRCA(ml_tree, c("NC_003436.1|Porcine_epidemic_diarrhea_virus|Sus_scrofa|NA", 
                                 "NC_009657.1|Scotophilus_bat_coronavirus_512|Scotophilus_kuhlii|China"))

Unclassified <- getMRCA(ml_tree, c("MZ081397.1|Alphacoronavirus_bat/Yunnan/MlYN20/2020|Myotis_laniger|China", 
                                   "OL410607.1|Eptesicus_bat_coronavirus|Eptesicus_fuscus|USA"))


# Get node numbers for single tips
Soracovirus    <- which(ml_tree$tip.label == "KY370053.1|Common_shrew_coronavirus_Tibet-2014|Sorex_araneus|China")
Myotacovirus   <- which(ml_tree$tip.label == "NC_028811.1|BtMr-AlphaCoV/SAX2011|Myotis_pilosus|China")
Amalacovirus   <- which(ml_tree$tip.label == "MT663548.1|Bat_alphacoronavirus_AMA_L_F|Desmodus_rotundus|Peru")
Colacovirus    <- which(ml_tree$tip.label == "NC_022103.1|Bat_coronavirus_CDPHE15/USA/2006|Myotis_lucifugus|USA")
Unclassified2  <- which(ml_tree$tip.label == "MZ081383.1|Alphacoronavirus_bat/Yunnan/CpYN11/2019|Chaerephon_plicatus|China")


# Add clade-style labels with visible bars
q <- q +
  geom_cladelabel(node = Soracovirus, label = "Soracovirus",
                  align = TRUE, offset = 1.8, extend = 0.5,
                  barsize = 3, fontsize = 11) +
  geom_cladelabel(node = Myotacovirus, label = "Myotacovirus",
                  align = TRUE, offset = 1.8, extend = 0.5,
                  barsize = 3, fontsize = 11) +
  geom_cladelabel(node = Amalacovirus, label = "Amalacovirus",
                  align = TRUE, offset = 1.8, extend = 0.5,
                  barsize = 3, fontsize = 11) +
  geom_cladelabel(node = Colacovirus, label = "Colacovirus",
                  align = TRUE, offset = 1.8, extend = 0.5,
                  barsize = 3, fontsize = 11) +
  geom_cladelabel(node = Unclassified2, label = "Unclassified Bat AlphaCoV",
                  align = TRUE, offset = 1.8, extend = 0.5,
                  barsize = 3, fontsize = 11)


# Highlight and label
q <- q +
  geom_hilight(node = TADBR_ACoV2, fill = "#4DBBD5", alpha = 0.4) +
  geom_cladelabel(node = TADBR_ACoV2, label = "Unclassified Bat AlphaCoV", align = TRUE, offset = 1.8, barsize = 3, fontsize = 11) +
  geom_cladelabel(node = Luchacovirus, label = "Luchacovirus", align = TRUE, offset = 1.8, barsize = 3, fontsize = 11) +
  geom_cladelabel(node = Setracovirus, label = "Setracovirus", align = TRUE, offset = 1.8, barsize = 3, fontsize = 11) +
  geom_cladelabel(node = Minacovirus, label = "Minacovirus", align = TRUE, offset = 1.8, barsize = 3, fontsize = 11) +
  geom_cladelabel(node = Rhinacovirus, label = "Rhinacovirus", align = TRUE, offset = 1.8, barsize = 3, fontsize = 11) +
  geom_cladelabel(node = Tegacovirus, label = "Tegacovirus", align = TRUE, offset = 1.8, barsize = 3, fontsize = 11) +
  geom_cladelabel(node = Duvinacovirus, label = "Duvinacovirus", align = TRUE, offset = 1.8, barsize = 3, fontsize = 11) +
  geom_cladelabel(node = Decacovirus, label = "Decacovirus", align = TRUE, offset = 1.8, barsize = 3, fontsize = 11) +
  geom_cladelabel(node = Minunacovirus, label = "Minunacovirus", align = TRUE, offset = 1.8, barsize = 3, fontsize = 11) +
  geom_cladelabel(node = Nyctacovirus, label = "Nyctacovirus", align = TRUE, offset = 1.8, barsize = 3, fontsize = 11) +
  geom_cladelabel(node = Pedacovirus, label = "Pedacovirus", align = TRUE, offset = 1.8, barsize = 3, fontsize = 11) +
  geom_cladelabel(node = Unclassified, label = "Unclassified Bat AlphaCoV", align = TRUE, offset = 1.8, barsize = 3, fontsize = 11) +
  geom_cladelabel(node = Sunacovirus, label = "Sunacovirus", align = TRUE, offset = 1.8, barsize = 3, fontsize = 11) 



# Save final plot
ggsave("AlphaCoV_concatenate.tiff", plot = q,
       width = 55, height = 35, units = "in", dpi = 600, limitsize = FALSE)
