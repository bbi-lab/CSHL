# Chapter 1. Basic scRNA-seq workflow in Monocle3, some ggplot, and why you must 
# explore your data! 
# Mary B. O'Neill & Anh Vo, BBI BAT-Lab, Single Cell Genomics
# CSHL Computational Genomics Course 2024

################################################################################
###### Setup ###################################################################
################################################################################
# Set your working directory to where the data is downloaded
setwd("/Users/oneillmb/Documents/CourseMaterials/CSHL") #CHANGE ME!!!

library(monocle3)
library(ggplot2)
library(tidyverse)
library(viridis)
library(data.table)
library(scuttle)
library(cowplot)
library(edgeR)
library(scales)

################################################################################
###### Let's start with a small, publicly available dataset ####################
###### from 10X Genomics website; 1K PBMCs #####################################
################################################################################
# Read the data
pbmc <- load_cellranger_data("data/inputs/10X_1kPBMCs")
?load_cellranger_data
?load_mm_data

## An alternative method of loading the data
#pbmc2 <- load_mm_data(mat_path = "10X_1kPBMCs/outs/filtered_feature_bc_matrix/matrix.mtx.gz", 
#                    feature_anno_path = "10X_1kPBMCs/outs/filtered_feature_bc_matrix/features.tsv.gz", 
#                    cell_anno_path = "10X_1kPBMCs/outs/filtered_feature_bc_matrix/barcodes.tsv.gz")

# learn about the cell_data_set class
class(pbmc) #cell_data_set
?cell_data_set
pbmc #dim: 33538 1222 

# multiple ways to access the column data
head(pData(pbmc))
head(colData(pbmc))

# multiple ways to access the row data
head(fData(pbmc))
head(rowData(pbmc))

# multiple ways to access the assays
counts(pbmc)[1:10,1:10]
assay(pbmc, 'counts')[1:10,1:10]

# Pre-process the data
?preprocess_cds
pbmc <- preprocess_cds(pbmc, num_dim = 50)
pbmc #reducedDimNames(1): PCA

# Plot PC variance
plot_pc_variance_explained(pbmc)

# Reduce dimensionality
?reduce_dimension
pbmc <- reduce_dimension(pbmc)
pbmc #reducedDimNames(2): PCA UMAP

#save in object for easier plotting later
pbmc$UMAP1 <- reducedDim(pbmc, "UMAP")[,1]
pbmc$UMAP2 <- reducedDim(pbmc, "UMAP")[,2]

# Plot within Monocle3
?plot_cells
plot_cells(pbmc)
plot_cells(pbmc, genes=c("CD14", "LYZ", # CD14+ Mono
                        "IL7R", "CCR7", # Naive CD4+ T
                        "S100A4", #	Memory CD4+
                        "MS4A1", # B
                        "CD8A", # CD8+ T
                        "FCGR3A", "MS4A7", # FCGR3A+ Mono
                        "GNLY", "NKG7", #	NK
                        "FCER1A", "CST3", #	DC
                        "PPBP" # Platelet
                        ))


# Run an unsupervised clustering
?cluster_cells
pbmc <- cluster_cells(pbmc)
pbmc$clusters <- clusters(pbmc) #save this in colData for easier access
plot_cells(pbmc)

# Marker gene detection - data driven
?top_markers
marker_test_res <- top_markers(pbmc, group_cells_by="cluster")

?saveRDS
saveRDS(marker_test_res, file="data/outputs/PBMC_marker_test_results.RDS") #save for later use

str(marker_test_res) #get in the habit of looking at data structures
head(marker_test_res) #show the first 6 lines

# filter the list to retain 3 top markers for each group
top_specific_markers <- marker_test_res %>%
  filter(fraction_expressing >= 0.10) %>%
  group_by(cell_group) %>%
  top_n(3, pseudo_R2)

# reference the gene_id
top_specific_marker_ids <- unique(top_specific_markers %>% pull(gene_id))

# make a dot plot of data-driven top markers and some cannonical marker genes
?plot_genes_by_group
plot_genes_by_group(pbmc,
                    c(top_specific_marker_ids, #data-driven markers
                      "CD14", "LYZ", # CD14+ Mono
                      "IL7R", "CCR7", # Naive CD4+ T
                      "S100A4", #	Memory CD4+
                      "MS4A1", # B
                      "CD8A", # CD8+ T
                      "FCGR3A", "MS4A7", # FCGR3A+ Mono
                      "GNLY", "NKG7", #	NK
                      "FCER1A", "CST3", #	DC
                      "PPBP"),
                    group_cells_by="cluster",
                    ordering_type="maximal_on_diag",
                    max.size=3)

################################################################################
###### Quality Control #########################################################
################################################################################

# Quality Control (QC) is a trade off. Being too aggressive can results in the 
# loss of rare cell populations, being too permissive can complicate cell 
# annotation and differential expression analyses. EVERY DATASET IS UNIQUE!

# The heart of scRNA-seq data is the feature-by-cell matrix
pbmc
counts(pbmc)[1:10,1:10] #sparse matrix, "." = 0; single-cell data is zero inflated! 

?detect_genes
pbmc <- detect_genes(pbmc)
head(fData(pbmc))
head(pData(pbmc))

expressed <- data.frame(rowData(pbmc)) %>% arrange(desc(num_cells_expressed))
counts(pbmc)[head(expressed$id, n=10), 1:10]

# Calculate mito content
head(fData(pbmc))

# How might we identify mitochondrial genes?
fData(pbmc)$MT <- grepl("^MT-", rowData(pbmc)$gene_short_name)
table(fData(pbmc)$MT) #sanity check

pData(pbmc)$MT_reads <- Matrix::colSums(exprs(pbmc)[fData(pbmc)$MT,])
pData(pbmc)$MT_perc <- pData(pbmc)$MT_reads/Matrix::colSums(exprs(pbmc))*100
summary(pbmc$MT_perc)

pData(pbmc)$n.umi <- Matrix::colSums(exprs(pbmc))

head(pData(pbmc))

# Lets examine some common QC metrics & learn about ggplot, a very flexible
# R package for plotting

# Your data must be a dataframe if you want to plot it using ggplot2
# We achieve this by forcing the column data to a dataframe
# The aes function specifies how variables in your dataframe map to features on your plot
# Geoms specify how your data should be represented 

# Most basic
ggplot(data.frame(colData(pbmc)), aes(x=n.umi, y=num_genes_expressed, color=MT_perc)) +
  geom_point()

# Let's spruce this up a bit
ggplot(data.frame(colData(pbmc)), aes(x=n.umi, y=num_genes_expressed, color=MT_perc)) +
  geom_point() +
  theme_light() + #there are several preset themes, you can also customize 
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) + #fancy log scale
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) #fancy log scale

# And color by cluster  
ggplot(data.frame(colData(pbmc)), aes(x=n.umi, y=num_genes_expressed, color=clusters)) +
  geom_point() +
  theme_light() +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)))


# We can use other types of plots (and layer!): violin, boxplot
ggplot(data.frame(colData(pbmc)), aes(x=clusters, y=n.umi, fill=clusters)) +
  geom_violin() + 
  geom_boxplot(width=0.2, fill="white", alpha=0.3) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  theme_light()

summary(pbmc$n.umi)

# And the % mitochondrial
ggplot(data.frame(colData(pbmc)), aes(y=MT_perc)) +
  geom_density(fill="salmon") + #yet another type of plot
  coord_flip() +
  theme_light()

summary(pbmc$MT_perc)

# Let's ensure there are no cells with high mito & high umi
ggplot(data.frame(colData(pbmc)), aes(x=n.umi, y=MT_perc)) +
  geom_point(alpha=0.3) + 
  theme_light() + 
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)))

# This is usually the case, but sometimes a certain cell type or cell state
# is associated with high mitochondrial content in a biologically meaningful
# context and we do not want to throw those cells out.

# What thresholds should be used? One approach is to use adaptive thresholds.
# If you can reasonably assume that most cells are of acceptable quality,
# identifying and removing outliers may be a solid approach.

?isOutlier

# The number of UMIs are not normally distributed so we use the log. Play with
# the number of MADs and see how this changes.
lib <- isOutlier(pbmc$n.umi, 
          log=TRUE, 
          nmads=3, #change me
          type=c("both"))
attr(lib, "thresholds") #1736 UMIs at 3 MADs

genes <- isOutlier(pbmc$num_genes_expressed, 
                 log=TRUE,
                 nmads=3, #change me
                 type=c("both"))
attr(genes, "thresholds") #680 genes at 3 MADs

high.mito <- isOutlier(pbmc$MT_perc, 
                      nmads=3, #change me
                      type=c("higher"))
attr(high.mito, "thresholds") #24% at 5 MADs

# Let's plot those thresholds!
ggplot(data.frame(colData(pbmc)), aes(x=n.umi, y=num_genes_expressed, color=MT_perc)) + #also try coloring by clusters
  geom_point(alpha=0.3) +
  theme_light() +
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  geom_vline(xintercept=attr(lib, "thresholds")[1], linetype="dotted", color="red") +
  geom_vline(xintercept=attr(lib, "thresholds")[2], linetype="dotted", color="red") +
  geom_hline(yintercept=attr(genes, "thresholds")[1], linetype="dotted", color="blue") +
  geom_hline(yintercept=attr(genes, "thresholds")[2], linetype="dotted", color="blue")
  
# flag the cells failing QC 
pbmc$passQC <- ifelse(pbmc$n.umi > attr(lib, "thresholds")[2] |
                             pbmc$n.umi < attr(lib, "thresholds")[1] |
                             pbmc$num_genes_expressed > attr(genes, "thresholds")[2] |
                             pbmc$num_genes_expressed < attr(genes, "thresholds")[1] |
                           pbmc$MT_perc > 25, #setting manually
                           FALSE, TRUE)
?table
table(pbmc$passQC, pbmc$clusters)

# For this particular dataset, the low quality cells all clustered together.
# Filtering based on hard or adaptive thresholds, or based on cluster membership
# both are valid options.

# Save!
?saveRDS
saveRDS(pbmc, file="data/outputs/processed_pbmc_cds.RDS") #save for use on Monday


################################################################################
###### Why are there are no rules or gold standards? ###########################
################################################################################
# Read in the metadata (column data) from 50K randomly subsampled barcodes from
# a benchmarking experiment our group did. These were performed with different
# technologies, on fixed samples, and much more shallowly sequenced.
df <- readRDS("data/inputs/benchmarking_metadata_50Kbarcodes_df.RDS")
str(df)

# Plot the distribution of UMIs and the threshold we determined above
# Does the red dotted line seem appropriate to you? HOPEFULLY NOT!
# Would you use the same threshold for all of these sample types? I SURE WOULD NOT!
ggplot(df, aes(x=paste(sample, assay), y=n.umi, fill=assay)) +
  facet_wrap(~samptype, scales="free_x", drop=T) +
  geom_violin() + 
  geom_boxplot(width=0.2, fill="white", alpha=0.3) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_fill_manual(values=c("purple", "blue")) +
  theme_light() +
  theme(axis.text.x = element_blank()) +
  geom_hline(yintercept=attr(lib, "thresholds")[1], linetype="dotted", color="red", size=1) #thresholds from the 10X PBMC data
  
# Plot the distribution of genes detected and the threshold we determined above
# Again, does this feel appropriate?
ggplot(df, aes(x=paste(sample, assay), y=num_genes_expressed, fill=assay)) +
  facet_wrap(~samptype, scales="free_x", drop=T) +
  geom_violin() + 
  geom_boxplot(width=0.2, fill="white", alpha=0.3) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_fill_manual(values=c("purple", "blue")) +
  theme_light() +
  theme(axis.text.x = element_blank()) +
  geom_hline(yintercept=attr(genes, "thresholds")[1], linetype="dotted", color="red", size=1) #thresholds from the 10X PBMC data

# Plot the distribution of MT and the threshold we determined above
# Again, does this feel appropriate?

# Because the dataset has both human and mouse derived samples, create a single
# column that we can plot on
df$MT_per <- ifelse(df$species == "hs", df$human_MT_reads/df$human_reads*100,
                    df$mouse_MT_reads/df$mouse_reads*100)

ggplot(df, aes(x=paste(sample, assay), y=MT_per, fill=assay)) +
  facet_wrap(~samptype, scales="free_x", drop=T) +
  geom_violin() + 
  geom_boxplot(width=0.2, fill="white", alpha=0.3, outlier.shape=NA) +
  scale_fill_manual(values=c("purple", "blue")) +
  theme_light() +
  theme(axis.text.x = element_blank()) +
  geom_hline(yintercept=attr(high.mito, "thresholds")[2], linetype="dotted", color="red", size=1) #thresholds from the 10X PBMC data

# Notice anything? Discuss.

################################################################################
###### Homework ################################################################
################################################################################      
# Play with thresholds
# Pick another dataset, download it, and play around with the above code. Or
# use your own data!
# Future tutorials will go more in depth, including multiple samples from 
# different runs, etc.


