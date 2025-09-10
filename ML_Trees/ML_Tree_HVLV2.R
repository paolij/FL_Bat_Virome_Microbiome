setwd("/Users/julia/Dropbox/MAVIAN_Lab/Manuscripts/Julia_FL_Bat/Figures_plus_data_scripts/HVLV2/ML/")
library("ggtree")
library("phytools")
library("ggplot2")
library("ape")
library("ggbreak")

# Load the tree
ml_tree <- read.tree("HVLV2_concatenate_fixed.nwk")

#show node numberss
ggtree(ml_tree,ladderize = TRUE) +
  geom_text(aes(label = node))

#Finding Nodes
which(ml_tree$tip.label == "HVLV2-Tb15-FL-2024|USA|2024")

# Visualize the rooted tree with bootstrap values as text labels
q <- ggtree(ml_tree, ladderize = TRUE) +
  geom_nodepoint(aes(subset = !is.na(as.numeric(label)) & as.numeric(label) > 90), 
                 size = 13, shape = 18, color = "black") +
  geom_treescale(x=0.07, y=35, color='black', fontsize = 15, linesize = 3) +
  geom_tiplab(aes(subset = node %in% c(4), label = label), size = 10, align = TRUE, offset = 0.0005, fontface = 'bold', color = "black") +  
  geom_tiplab(aes(subset = !node %in% c(4), label = label), size = 10, align = TRUE, offset = 0.0005) + # Regular for others
  coord_cartesian(xlim = c(0, max(nodeHeights(ml_tree)) * 1.9))
q

# Load metadata
metadata <- read.csv('HVLV2_metadata.csv', header = TRUE)
unique(metadata$Accession)
unique(metadata$Country)
unique(metadata$Host)
unique(metadata$Date)

# Define colors for bathost
mycolor <- c("Tadarida brasiliensis" = "#AD002A", "Culex tarsalis" = "#42B540", "Culex erythrothorax" = "#925E9F", "Culex sp." = "#FDAF91", "Culicidae sp." = "#FFDC91" )


# Plot with metadata
p <- q %<+% metadata + 
  geom_tippoint(aes(color = Host), size = 10.0) + 
  scale_color_manual(values = mycolor) + 
  theme(legend.position = 'left', 
        legend.text = element_text(size = 40),  # Adjust text size
        legend.title = element_text(size = 40)) + # Adjust title size
  guides(color = guide_legend(override.aes = list(size = 10))) +
  coord_cartesian(xlim = c(0, max(nodeHeights(ml_tree)) * 1.9))


# Save the plot
ggsave("HVLV2_ML_concatenate.pdf", width = 55, height = 15, units = "in", dpi = 600, limitsize = FALSE)
