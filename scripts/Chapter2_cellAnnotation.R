# Chapter 2. The bottleneck of cell annotation.
# Mary B. O'Neill & Anh Vo, BBI BAT-Lab, Single Cell Genomics
# CSHL Computational Genomics Course 2024

################################################################################
###### Setup ###################################################################
################################################################################
# Set your working directory to where the data is downloaded
setwd("/Users/oneillmb/Documents/CourseMaterials/CSHL") #CHANGE ME!!!

library(monocle3)
library(Matrix)
library(ggplot2)
library(tidyverse)
library(clusterProfiler)
library("org.Hs.eg.db", character.only = TRUE)
library(celldex)
library("SingleR")
library(pheatmap)
library(garnett)


################################################################################
###### Let's start with the PBMC dataset you processed on Friday ###############
################################################################################
# load the cds object you generated on Friday
pbmc <- readRDS("data/outputs/processed_pbmc_cds.RDS")

# Remind yourself what it looked like
plot_cells(pbmc) #built-in monocle UMAP plotting

# We can already get a good sense for cell types based on common marker genes
# if you have a little domain knowledge
plot_genes_by_group(pbmc,
                    c("CD14", "LYZ", # CD14+ Mono
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

# We can alternatively plot these on the UMAP
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

################################################################################
###### Marker gene enrichment analysis #########################################
###### Code modified from Wei Yang, SASC workshop 2024-09-23 ###################
################################################################################
# Load the organism whose GO terms we will be looking at
organism = "org.Hs.eg.db" #human 

# Read in the results of the marker gene test you did on Friday
marker_test_res <- readRDS("data/outputs/PBMC_marker_test_results.RDS")

# Filter based on significance
df = marker_test_res %>% filter(marker_test_q_value < 0.05)

# Convert gene names to ENTREZID (GO terms associated to ENTREZID)
# Prepare for the marker gene sets
gene.df <- select(org.Hs.eg.db,
                  keys = df$gene_id,  
                  columns = c("ENTREZID", "ENSEMBL"),  
                  keytype = "ENSEMBL")

# Drop duplicate by selecting the first entrezid for each gene name
gene.df <- gene.df[!duplicated(gene.df$ENSEMBL), ]
df_entrez = left_join(df, gene.df, by=c('gene_id'='ENSEMBL')) %>% drop_na()

# Prepare background gene sets (all genes expressed in pbmc)
all.gene.df <- select(org.Hs.eg.db,
                      keys = rowData(pbmc)$id,
                      columns = c("ENTREZID", "ENSEMBL"),
                      keytype = "ENSEMBL")

all.gene.df <- all.gene.df[!duplicated(all.gene.df$ENSEMBL), ] %>% drop_na()

# KEGG pathway
?compareCluster
compKEGG = compareCluster(geneClusters=ENTREZID~cell_group, data = df_entrez, fun = "enrichKEGG",
                          organism='hsa',universe = all.gene.df$ENTREZID)

dotplot(compKEGG, title = "KEGG analysis")
setReadable(compKEGG, OrgDb = org.Hs.eg.db, keyType="ENTREZID")

# Gene ontology (BP, MF, CC)
compGO =compareCluster(geneClusters=ENTREZID~cell_group, data = df_entrez, fun = "enrichGO",OrgDb=organism,
                       ont='BP',universe = all.gene.df$ENTREZID)
dev.off()
dotplot(compGO, title = "GO BP analysis")

# Are these pathways consistent with our cell type predictions?
# Next, let's see how our inclinations stack against an automated cell type prediction

################################################################################
###### SingleR - built-in reference ############################################
################################################################################
# Cedex contains a number of pre-built reference datasets that contain immune cell 
# circulating in the blood. SingleR takes a long time when the search space is large
# so we are going to use the smallest reference as a demonstration and load the 
# results from one others to demonstrate how different references affect predictions.

?fetchReference
hpca <- fetchReference("hpca", "2024-02-26") # Human primary cell atlas (HPCA)
hpca #dim: 19363 713, SummarizedExperiment class
?SummarizedExperiment
assay(hpca, "logcounts")[1:5, 1:5] #these are bulk, not single-cell!

blueprint <- fetchReference("blueprint_encode", "2024-02-26") # Blueprint/ENCODE
blueprint #dim: 19859 259 

# Prepare our data
countmat <- normalized_counts(pbmc) #SingleR expects the normalized count matrix
?normalized_counts

# Do gene nomenclatures match? NO! 
identical(rownames(countmat), rowData(hpca)$id)

# Modify our gene names
rownames(countmat) <- rowData(pbmc)$gene_short_name
countmat[1:5, 1:5]

# Run SingleR
# Since this is very slow, I have commented out the code and you will simply
# read in the results
#hpca.pred <- SingleR(test = countmat, ref = assay(hpca, "logcounts"), assay.type.test=1,
#                          labels = hpca$label.main)
#saveRDS(hpca.pred, file="data/inputs/PBMC_SingleR_predictions_based_on_HPCA_reference.RDS")

# Load the pre-run results
hpca.pred <- readRDS("data/inputs/PBMC_SingleR_predictions_based_on_HPCA_reference.RDS")
table(hpca.pred$pruned.labels)

# Append to our cds object
pbmc$hpca_predictions <- hpca.pred$pruned.labels

# One way to visualize this...
# Using a large pseudo-count for a smoother color transition
# between 0 and 1 cell in each 'tab'.
dev.off()

tab <- table(paste("Cluster", pbmc$clusters), 
             paste("SingleR", hpca.pred$pruned.labels))

hpcaPlot <- pheatmap(log10(tab+10), main="HPCA Reference",
                 color=viridis::viridis(100), silent=TRUE)

hpcaPlot

# How about another reference?
# This reference dataset is much smaller and thus faster, so we will actually run
# this one.
blueprint.pred <- SingleR(test = countmat, ref = assay(blueprint, "logcounts"), assay.type.test=1,
                     labels = blueprint$label.main)

table(blueprint.pred$pruned.labels)

# Append to our cds object
pbmc$blueprint_predictions <- blueprint.pred$pruned.labels

# And plot
dev.off() #clear the plot area
tab <- table(paste("Cluster", pbmc$clusters), 
             paste("SingleR", blueprint.pred$pruned.labels))

blueprintPlot <- pheatmap(log10(tab+10), main="Blueprint Reference",
                 color=viridis::viridis(100), silent=TRUE)
blueprintPlot

# On the UMAP
ggplot(data.frame(colData(pbmc)), aes(x=UMAP1, y=UMAP2, color=blueprint_predictions)) +
  geom_point()

# It looks like our initial clustering actually wasn't fine-grained enough. All
# the T-cells were in one cluster and we clearly have the power to resolve them!
ggplot(data.frame(colData(pbmc)), aes(x=UMAP1, y=UMAP2, color=clusters)) +
  geom_point()

# We can tweak the resolution of clustering by altering the 'resolution' or 'k'
# parameters. 
?cluster_cells
pbmc <- cluster_cells(pbmc, k = 10)
plot_cells(pbmc) # We see separation between presumed CD4 and CD8 T-cells!
pbmc$clusters_k10 <- clusters(pbmc) #save this in colData for easier access

# But what about classical vs non-classical monocytes? Based on the cananonical
# marker genes expression levels it seemed they resolved as well.
pbmc <- cluster_cells(pbmc, k = 5)
plot_cells(pbmc)
pbmc$clusters_k5 <- clusters(pbmc) #save this in colData for easier access

# We can also use a different clustering method
pbmc <- cluster_cells(pbmc, cluster_method = "louvain")
plot_cells(pbmc)
pbmc$louvain_clusters <- clusters(pbmc) #save this in colData for easier access

# Where do you stop? Truth is, this is often an interative process requiring lots
# of exploration and domain-specific knowledge. Don't get attached to your first 
# clustering analsis!

################################################################################
###### Let's try another program: Garnett ######################################
################################################################################
# Load a pre-trained Garnett classifier for PBMCs
classifier <- readRDS("data/inputs/hsPBMC_20191017.RDS")

pbmc <- classify_cells(pbmc, classifier,
                       db = "none",
                       cluster_extend = TRUE)
# It returns a cds object with two new columns!

# And plot!
ggplot(data.frame(colData(pbmc)), aes(x=UMAP1, y=UMAP2, color=cluster_ext_type)) +
  geom_point()

# Or this visualization
dev.off()
tab <- table(paste("Cluster", pbmc$clusters_k5), 
             paste("Garnett", pbmc$cluster_ext_type))
garnPlot <- pheatmap(log10(tab+10), main="Pre-trained Garnett Reference",
                          color=viridis::viridis(100), silent=TRUE)
garnPlot

# The resolution between CD4 and CD8 T cells does not seem to work as well in this
# particular case. Hopefully this demonstrates the value of trying multiple 
# methods and references. 
# This afternoon we will work with a dataset for which the Garnett classifier 
# worked exceedingly well on.

################################################################################
###### How about another approach? #############################################
###### Projecting query data onto a reference ##################################
################################################################################
# What if there aren't good pre-trained reference datasets? 
# If you can find appropriate reference datasets, you can train your own
# classifiers for both Garnett and SingleR. 
# An alternative is to project your dataset onto another well annotated one,
# and transfer the labels. We will demonstrate this approach with a 
# whole mouse embryo dataset that we have generated at our platform.

# Load the query data set.
cds_qry <- readRDS("data/inputs/Scale_Embryo-mm-mix_68886genes_4410barcodes_cds.RDS")

cds_qry #dim: 68886 4410

head(colData(cds_qry))
table(cds_qry$donor) # mix of B6 and F1 embryos

# Load the well-annotated reference dataset. For this exercise, we will
# load a sub-sample of the 12.4M cell atlas data generated in: 
# "A single-cell time-lapse of mouse prenatal development from gastrula to birth"
# https://doi.org/10.1038/s41586-024-07069-w

# Load the downsampled dataset (65000 nuclei)
cds_ref <- readRDS("data/inputs/cx_mouse_combined_cds_sub5000_celltype.RDS")
names(rowData(cds_ref)) <- c("id", "gene_short_name", "chr", "chr_start", "chr_end", "chr_strand", "gene_type")

cds_ref #dim: 49585 65000 

# Let's peak at the annotations
head(colData(cds_ref))
table(cds_ref$experimental_id) # 5 runs
table(cds_ref$embryo_id) # many embryos
table(cds_ref$embryo_sex) # both sexes represented
table(cds_ref$major_trajectory) # we will use these for cell type assignment

### Remove genes that are not in both data sets
# Genes in reference.
genes_ref <- row.names(cds_ref)
head(genes_ref)

# Genes in query.
genes_qry <- row.names(cds_qry)
tail(genes_qry)
# notice something?

# Shared genes.
genes_shared <- intersect(genes_ref, genes_qry)
length(genes_shared) # no overlap because the reference naming conventions are different!

### Unify gene names
head(rowData(cds_ref))
head(rowData(cds_qry))
tail(rowData(cds_qry))

# Somewhat commonly you need to strip everything after a decimal point... here 
# is one way to approach that. It is not necessary for these datasets, however.
#rownames(cds_ref) <- sapply(strsplit(as.character(rownames(cds_ref)), "\\."), `[`, 1) 
#rowData(cds_ref)$gene_id <- sapply(strsplit(as.character(rowData(cds_ref)$gene_id), "\\."), `[`, 1) 

# Remove the assembly name and underscores in front of the ensembl gene ids in the query gene names
rownames(cds_qry) <- substring(rownames(cds_qry), 8) 
rowData(cds_qry)$id <- substring(rowData(cds_qry)$id, 8) 
rowData(cds_qry)$gene_short_name <- substring(rowData(cds_qry)$gene_short_name, 8) 

# Shared genes... try 2
genes_ref <- row.names(cds_ref)
head(genes_ref)

# Genes in query.
genes_qry <- row.names(cds_qry)
tail(genes_qry)

#Overlap
genes_shared <- intersect(genes_ref, genes_qry)
length(genes_shared) # 29621, that's more like it!

# Remove non-shared genes.
cds_ref <- cds_ref[genes_shared,]
cds_qry <- cds_qry[genes_shared,]

### Estimate size factors
cds_ref <- estimate_size_factors(cds_ref)
cds_qry <- estimate_size_factors(cds_qry)

### Process the reference data set
cds_ref <- preprocess_cds(cds_ref) #takes a bit
cds_ref <- reduce_dimension(cds_ref, build_nn_index=TRUE) #takes a bit

# Save the PCA and UMAP transform models for use with projection.
save_transform_models(cds_ref, 'data/outputs/cds_ref_test_models')

### Project the query data set into the reference space
# Load the reference transform models into the query cds.
cds_qry <- load_transform_models(cds_qry, 'data/outputs/cds_ref_test_models')

# Apply the reference transform models to the query cds.
cds_qry <- preprocess_transform(cds_qry)
cds_qry <- reduce_dimension_transform(cds_qry)

# Plot the combined data sets
# Note that we weighted our subsampling of the reference dataset to include 
# more dense representation of rare cell types.
plot_cells(cds_ref, color_cells_by = "major_trajectory", label_cell_groups = TRUE)
plot_cells(cds_qry, color_cells_by = "donor")

# Label the data sets.
colData(cds_ref)[['data_set']] <- 'reference'
colData(cds_qry)[['data_set']] <- 'query'

# Combine the reference and query data sets.
cds_combined <- combine_cds(list(cds_ref, cds_qry),  keep_all_genes=TRUE, cell_names_unique=TRUE, keep_reduced_dims=TRUE)
plot_cells(cds_combined, color_cells_by='data_set')

?transfer_cell_labels
cds_qry <- transfer_cell_labels(cds_qry, reduction_method='UMAP', ref_coldata=colData(cds_ref), ref_column_name='major_trajectory', query_column_name='cell_type_xfr', transform_models_dir='cds_ref_test_models')
table(is.na(cds_qry$cell_type_xfr))

?fix_missing_cell_labels
cds_qry <- fix_missing_cell_labels(cds_qry, reduction_method='UMAP', from_column_name='cell_type_xfr', to_column_name='cell_type_fix')
table(is.na(cds_qry$cell_type_fix))

plot_cells(cds_qry, color_cells_by = "cell_type_fix")

################################################################################
###### Takeaways ###############################################################
################################################################################      
# Clustering and annotation are a major bottleneck in single-cell analyses.
# As with many things single-cell related, every dataset is unique and what may
# work best for one dataset will not necessarily be the best for yours.
# It is best to think of cell type annotation as an iterative process.
# You will likely refine your annotations many times.

################################################################################
###### Homework ################################################################
################################################################################      
# Pick another dataset, download it, and play around with one or more of the
# annotation strategies outlined above. Or use your own data!

