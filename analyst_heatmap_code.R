### code to create heatmap

library(gplots)
#Create a heatmap of the gene-expression of each gene across all samples
# Load metadata
metadata <- read.csv('/project/bf528/project_1/doc/proj_metadata.csv')
metadata <- subset(metadata, select = c(geo_accession, cit.coloncancermolecularsubtype))
colors <- ifelse(metadata$cit.coloncancermolecularsubtype == 'C3', 'red', 'blue')

# Save heatmap as a file
png('/projectnb/bf528/users/glass_bottom_boat/project_1/analyst/heatmap_2.png', width=1920, height=1080, res=100)
# Vector of colors for annotation: “red” if the sample belongs to the C3 subtype and “blue” otherwise
heatmap.2(as.matrix(var_filter.2), xlab='Patient Sample', ylab='Gene', main='Gene Expression Across Samples',ColSideColors = colors)

# Legend for Molecular Subtype
legend(x=100,y=150,legend=c("C3","C4"),fill = c("red","other"))
dev.off()
#xpd=TRUE
