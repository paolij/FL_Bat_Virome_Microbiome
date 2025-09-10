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
  #geom_tippoint(size = 15.0, color = "black") +
  coord_cartesian(xlim = c(0, max(nodeHeights(ml_tree)) * 1.9))
q


# MRCA nodes
#subclade1 <- getMRCA(ml_tree, c("OM817543.1|USA|2020", "KX883772.1|China|2013"))


# # Highlight and label
# q <- q +
#   geom_hilight(node = subclade1, fill = "#4DBBD5", alpha = 0.4) 

# Load metadata
metadata <- read.csv('metadata.csv', header = TRUE)
unique(metadata$Accession_short)
unique(metadata$Accession)
unique(metadata$Country)
unique(metadata$host)
unique(metadata$date)

# Define colors for bathost
mycolor <- c("Tadarida brasiliensis" = "#AD002A", "Culex tarsalis" = "#42B540", "Culex erythrothorax" = "#925E9F", "Culex sp." = "#FDAF91", "Culicidae sp." = "#FFDC91" )


# Plot with metadata
p <- q %<+% metadata + 
  geom_tippoint(aes(color = host), size = 10.0) + 
  scale_color_manual(values = mycolor) + 
  theme(legend.position = 'left', 
        legend.text = element_text(size = 40),  # Adjust text size
        legend.title = element_text(size = 40)) + # Adjust title size
  guides(color = guide_legend(override.aes = list(size = 10))) +
  coord_cartesian(xlim = c(0, max(nodeHeights(ml_tree)) * 1.9))


# Save the plot
ggsave("Hubei_no_outgroup_concatenate.pdf", width = 55, height = 15, units = "in", dpi = 600, limitsize = FALSE)