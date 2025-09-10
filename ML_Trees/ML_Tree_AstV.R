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
setwd("/Users/julia/Dropbox/MAVIAN_Lab/Manuscripts/Julia_FL_Bat/Figures_plus_data_scripts/AstV/Concatenated/")

# Load the tree
ml_tree <- read.tree("astrovirus_concatenate_fixed.nwk")

# Load metadata and prepare image paths
metadata <- read_csv("mamastrovirus_md.csv") %>%
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
which(ml_tree$tip.label == "Ma_Frozen_4|AstV-Ma4-FL-2021|Myotis_austroriparius|USA")

# Visualize the rooted tree with bootstrap values as text labels
q <- ggtree(ml_tree, ladderize = TRUE) +
  geom_nodepoint(aes(subset = !is.na(as.numeric(label)) & as.numeric(label) > 90), 
                 size = 13, shape = 18, color = "black") +
  geom_treescale(x=0, y=25, color='black', fontsize = 15, linesize = 3) +
  geom_tiplab(aes(subset = node %in% c(28), label = label), size = 10, align = TRUE, offset = 0.15, fontface = 'bold', color = "black") +  
  geom_tiplab(aes(subset = !node %in% c(28), label = label), size = 10, align = TRUE, offset = 0.15) + # Regular for others
  coord_cartesian(xlim = c(0, max(nodeHeights(ml_tree)) * 1.9))
q

# MRCA nodes
Bat_AstV <- getMRCA(ml_tree, c("MZ218054.1|Bat_astrovirus|Myotis_daubentonii|Denmark", 
                                "Ma_Frozen_4|AstV-Ma4-FL-2021|Myotis_austroriparius|USA"))

GII <- getMRCA(ml_tree, c("NC_013060.1|Mamastrovirus_9|Homo_sapiens|USA", 
                          "NC_038368.1|Mamastrovirus_17|Hipposideros_pomona|China"))

GI <- getMRCA(ml_tree, c("NC_030922.1|Mamastrovirus_1|Homo_sapiens|India", 
                         "NC_033792.1|Qinghai_Himalayan_marmot_astrovirus_1|Marmota_himalayana|China"))

# Highlight and label
q <- q +
  geom_hilight(node = Bat_AstV, fill = "#4DBBD5", alpha = 0.6) +
  #geom_hilight(node = subclade1, fill = NA, color = "#E64B35", alpha = 1, linewidth = 1) +
  geom_hilight(node = GII, fill = "white", alpha = 0) +
  geom_cladelabel(node = GII, label = "GII", align = TRUE, offset = 2.1, barsize = 3, fontsize = 11) +
  geom_hilight(node = GI, fill = "white", alpha = 0) +
  geom_cladelabel(node = GI, label = "GI", align = TRUE, offset = 2.1, barsize = 3, fontsize = 11)

# Get max x position of tips
tip_x_max <- max(img_data$x, na.rm = TRUE)

# Compute aligned x-position (same as geom_tiplab offset)
aligned_x <- tip_x_max + 0.01

# Add animal icons aligned to labels
q <- q + geom_image(
  data = img_data,
  aes(x = aligned_x, y = y, image = image),  # fixed x for alignment
  size = 0.03,
  inherit.aes = FALSE
)

# Save final plot
ggsave("AstV_BatCoV_FL_with_images_concatenated.pdf", plot = q,
       width = 55, height = 35, units = "in", dpi = 600, limitsize = FALSE)

