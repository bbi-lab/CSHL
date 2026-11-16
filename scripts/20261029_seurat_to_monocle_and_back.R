##################### Seurat <-> Monocle3 conversion ###########################

# General process is :
#  1. Extract count matrix
#  2. Extract metadata
#  3. Construct object to the package's specifications
#  4. Copy over reductions (PCA, UMAP) and clusters

################################################################################

  library(Seurat)        
  library(SeuratObject)
  library(monocle3)     

################################################################################

# -------------------- Seurat -> cell_data_set --------------------------------

data <- readRDS("data/inputs/Seurat_data_example.rds")

# Cell_data_set is the format used in monocle3

## 1. Count matrix ----

# Monocle3 uses it's own normalization , so we're taking the raw "counts" layer.

counts <- LayerData(data, assay = "RNA", layer = "counts")

## 2. Cell metadata ----

cell_metadata <- data@meta.data

## 3. Gene metadata ----

# Monocle3 uses gene metadata, stored in rowData(), and expects it to have a 
# column called "gene_short_name". plot_cells(genes = ...) and other functions l
# ook genes up by this column, so it needs to hold the gene symbols.

gene_metadata <- data.frame(gene_short_name = rownames(counts),
                            row.names       = rownames(counts))

## 4. Build the cds ----

cds <- new_cell_data_set(
  expression_data = counts,
  cell_metadata   = cell_metadata,
  gene_metadata   = gene_metadata)

## 5. Transferring reductions and clustering ----

# This step is optional, you can run dimensionality reduction and clustering
# in the cell_data_set, but the results will probably be a bit different
# from Seurat, so if you want to keep the same UMAP and clusters you need
# to manually transfer them.

# Reductions live in reducedDims(cds). Monocle names them in capitals
# ("PCA", "UMAP"), while Seurat uses lowercase ("pca", "umap").
# Embeddings() returns a cells x dimensions matrix with the same cell order
# as the counts, so the embeddings line up with the cds.

# Check the name of your reductions with:
Reductions(data)

reducedDim(cds, "PCA")  <- Embeddings(data, "pca")
reducedDim(cds, "UMAP") <- Embeddings(data, "umap")

# Monocle stores clusters in cds@clusters[["UMAP"]], vs Seurat identities (Idents)
# Monocle also stores a lower resolution clustering, partitions. As we don't have 
# information in Seurat, we're settiong partition to "1".

cluster_ids   <- Idents(data)
partition_ids <- factor(setNames(rep("1", ncol(cds)), colnames(cds)))

cds@clusters[["UMAP"]] <- list(cluster_result = list(),
                               partitions     = partition_ids,
                               clusters       = cluster_ids)


# Let's check that everything worked with a plot
plot_cells(cds,
           group_label_size = 5)


## Standard monocle pipeline ----

# Alternatively, once you have built the cell_data_set object, you can run the 
# monocle3 pipeline to calculate the reductions and clusters.

cds <- preprocess_cds(cds)
set.seed(1000)
cds <- reduce_dimension(cds,
                        reduction_method = "UMAP",
                        preprocess_method = "PCA",
                        umap.min_dist = min_dist,
                        umap.n_neighbors = n_neigh) 

cds$UMAP1 <- reducedDim(cds, "UMAP")[,1]
cds$UMAP2 <- reducedDim(cds, "UMAP")[,2]

set.seed(1000)
cds <- cluster_cells(cds, 
                     k = k,
                     cluster_method = "leiden")
cds$leiden_clusters <- clusters(cds)


### Trajectory analysis ----

# The most popular use of monocle3 is performing a pseudotime trajectory analysis.
# more information in: https://cole-trapnell-lab.github.io/monocle3/docs/trajectories/

cds <- learn_graph(cds)
cds <- order_cells(cds)
plot_cells(cds, color_cells_by = "pseudotime")




# ----------------------- cell_data_set -> Seurat -----------------------------

# Using the processed dataset we created in the first tutorial
cds <-  load_monocle_objects("data/outputs/processed_pbmc_cds")


## 1. Count matrix ----

# counts() gives the raw counts stored in the cds.

counts <- counts(cds)

## 2. Metadata ----

# colData() is the cds equivalent of meta.data. 

cell_metadata <- as.data.frame(colData(cds))

## 3. Build the Seurat object ----

seurat <- CreateSeuratObject(counts    = counts,
                             assay     = "RNA",
                             meta.data = cell_metadata)

# A new Seurat object only has counts, so we normalize to fill the "data" layer.
seurat <- NormalizeData(seurat)

## 4. Transferring reductions and clustering ----

# Seurat stores each reduction as a DimReduc object. Each needs a "key": a
# prefix ending in "_" that names the columns (UMAP_1, UMAP_2, ...). The
# column names of the embedding matrix must match that key.

umap <- reducedDim(cds, "UMAP")
colnames(umap) <- paste0("UMAP_", seq_len(ncol(umap)))
seurat[["umap"]] <- CreateDimReducObject(embeddings = umap,
                                         key        = "UMAP_",
                                         assay      = "RNA")

pca <- reducedDim(cds, "PCA")
colnames(pca) <- paste0("PC_", seq_len(ncol(pca)))
seurat[["pca"]] <- CreateDimReducObject(embeddings = pca,
                                        key        = "PC_",
                                        assay      = "RNA")

# Monocle results go into meta.data as regular columns.
seurat$monocle_clusters   <- clusters(cds)
seurat$monocle_partitions <- partitions(cds)
seurat$pseudotime         <- pseudotime(cds)

# Optional: use the Monocle clusters as the active identities
Idents(seurat) <- "monocle_clusters"

# Checking everything worked with a plot

DimPlot(seurat, label = TRUE)


## Standard Seurat pipeline ----

# Alternatively, once you have built the Seurat object, you can run the 
# Seurat pipeline to calculate the reductions and clusters.

seurat <- NormalizeData(seurat)
seurat <- FindVariableFeatures(seurat)
seurat <- ScaleData(seurat)
seurat <- RunPCA(seurat, npcs = 50)
seurat <- FindNeighbors(seurat, dims = 1:50)
seurat <- FindClusters(seurat, resolution = 0.1)
seurat <- RunUMAP(seurat, dims = 1:50, reduction = "pca")


