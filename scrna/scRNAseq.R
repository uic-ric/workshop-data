# =============================================================================
# Single Cell RNA-seq
# Research Informatics Core
# Fall 2026
#
# Generated from scRNAseq.Rmd by purl_handouts.R -- do not edit by hand.
# =============================================================================

# =============================================================================
# DAY 1
# =============================================================================

# -----------------------------------------------------------------------------
# PREPARATION: Installing R packages
# -----------------------------------------------------------------------------

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install('multtest', update=F)
BiocManager::install("ComplexHeatmap", update=F)
BiocManager::install("circlize", update=F)
BiocManager::install("edgeR", update=F)
install.packages('Seurat')
install.packages('Matrix')
install.packages('tidyr')
install.packages('dplyr')
install.packages('ggplot2')
install.packages('devtools')
devtools::install_github('immunogenomics/presto')
library(Seurat)
library(multtest)
library(circlize)
library(ComplexHeatmap)
library(dplyr)
library(tidyr)
library(edgeR)
library(ggplot2)
library(Matrix)
library(presto)

# -----------------------------------------------------------------------------
# Read single-cell data into R
# -----------------------------------------------------------------------------

library(Seurat)
library(Matrix)

# --- Sparse Matrix (10X data) -----------------------------------------------

    sparse_10X <- Read10X("~/Downloads/10X_test/filtered_feature_bc_matrix/",
                          gene.column=1)
    sparse_10X <- Read10X("~/../Downloads/10X_test/filtered_feature_bc_matrix/",
                          gene.column=1)
seurat_10X <- CreateSeuratObject(counts=sparse_10X, project="10X_sample_data")
seurat_10X

# --- Full matrix data -------------------------------------------------------

counts_sparse <- Matrix(as.matrix(counts_full),sparse=T)
format(object.size(counts_full), units="auto")
format(object.size(counts_sparse), units="auto")
seurat_indrop <- CreateSeuratObject(counts=counts_sparse, project="indrop_sample_data")
seurat_indrop

# --- Read in h5 file (TAKE-HOME) --------------------------------------------

h5_10X <- Read10X_h5("...path to file.../10X_h5test.h5", use.names=F)
seurat_from_h5 <- CreateSeuratObject(h5_10X, project="data_from_h5")

# -----------------------------------------------------------------------------
# Quality control and filtering
# -----------------------------------------------------------------------------

# --- Download the data and read into R (TAKE HOME: Skip for now, try later at home.) ---

library(Seurat)
V035_F031 <- Read10X("~/Downloads/GSE138867/sub_V035_F031", gene.column=1)
V039_F093 <- Read10X("~/Downloads/GSE138867/sub_V039_F093", gene.column=1)
V040_F198 <- Read10X("~/Downloads/GSE138867/sub_V040_F198", gene.column=1)
V041_F043 <- Read10X("~/Downloads/GSE138867/sub_V041_F043", gene.column=1)
V043_F047 <- Read10X("~/Downloads/GSE138867/sub_V043_F047", gene.column=1)
V044_F037 <- Read10X("~/Downloads/GSE138867/sub_V044_F037", gene.column=1)
V045_F120 <- Read10X("~/Downloads/GSE138867/sub_V045_F120", gene.column=1)
V054_F209 <- Read10X("~/Downloads/GSE138867/sub_V054_F209", gene.column=1)
V035_F031 <- Read10X("~/../Downloads/GSE138867/sub_V035_F031", gene.column=1)
V039_F093 <- Read10X("~/../Downloads/GSE138867/sub_V039_F093", gene.column=1)
V040_F198 <- Read10X("~/../Downloads/GSE138867/sub_V040_F198", gene.column=1)
V041_F043 <- Read10X("~/../Downloads/GSE138867/sub_V041_F043", gene.column=1)
V043_F047 <- Read10X("~/../Downloads/GSE138867/sub_V043_F047", gene.column=1)
V044_F037 <- Read10X("~/../Downloads/GSE138867/sub_V044_F037", gene.column=1)
V045_F120 <- Read10X("~/../Downloads/GSE138867/sub_V045_F120", gene.column=1)
V054_F209 <- Read10X("~/../Downloads/GSE138867/sub_V054_F209", gene.column=1)
scV035_F031 <- CreateSeuratObject(counts=V035_F031, project="V035_F031")
scV039_F093 <- CreateSeuratObject(counts=V039_F093, project="V039_F093")
scV040_F198 <- CreateSeuratObject(counts=V040_F198, project="V040_F198")
scV041_F043 <- CreateSeuratObject(counts=V041_F043, project="V041_F043")
scV043_F047 <- CreateSeuratObject(counts=V043_F047, project="V043_F047")
scV044_F037 <- CreateSeuratObject(counts=V044_F037, project="V044_F037")
scV045_F120 <- CreateSeuratObject(counts=V045_F120, project="V045_F120")
scV054_F209 <- CreateSeuratObject(counts=V054_F209, project="V054_F209")
sc_data <- merge(scV035_F031, 
  y=c(scV039_F093,scV040_F198,scV041_F043,scV043_F047,
      scV044_F037,scV045_F120,scV054_F209),
  add.cell.ids=c("s1","s2","s3","s4","s5","s6","s7","s8"))

# --- Cell QC steps (Start at this point in the workshop) --------------------

sc_data <- readRDS(url("https://wd.cri.uic.edu/scrna/GSE138867.rds"))
sc_data <- JoinLayers(sc_data)
head(sc_data@meta.data)
table(sc_data$orig.ident)
# add in sample group metadata
metadata <- read.delim("https://wd.cri.uic.edu/scrna/GSE138867_metadata.txt",row.names=1)
sc_data@meta.data$group <- metadata[sc_data@meta.data$orig.ident,"Group"]
head(sc_data@meta.data)
genes_to_ids <- read.table("https://wd.cri.uic.edu/scrna/features.tsv",
  row.names=1, col.names=c("","gene.symbol"))
head(genes_to_ids)
sc_data[["RNA"]] <- AddMetaData(sc_data[["RNA"]],
  genes_to_ids[Features(sc_data),,drop=F])

# --- Basic quality control checks -------------------------------------------

genes <- read.table("https://wd.cri.uic.edu/scrna/hg19_genes.bed")
dim(genes)
# column 4 have the gene IDs, and column 1 the chromosome
head(genes)
mt.genes <- as.character(genes[genes[,1]=="chrM", 4])
length(mt.genes)
sc_data[["percent.MT"]] <- PercentageFeatureSet(sc_data, features=mt.genes)
VlnPlot(sc_data, features = c("nFeature_RNA","nCount_RNA","percent.MT"), pt.size=0.2)
FeatureScatter(sc_data, feature1="nCount_RNA", feature2="percent.MT")

# --- Filtering cells for downstream analysis --------------------------------

sc_subset <- subset(sc_data, 
  nFeature_RNA>500 & nCount_RNA>2000 & percent.MT < 10)
sc_subset

# --- Doublet analysis (TAKE HOME: Skip for now, try later at home. Jump to the next sub-section.) ---

library(scds)
library(SingleCellExperiment)
samples <- unique(sc_subset@meta.data$orig.ident)
doublets_df <- data.frame()
for(s in samples){
  # subset to this sample
  this_sc <- subset(sc_subset, orig.ident == s)
  # get counts
  this_counts <- GetAssayData(this_sc, layer="counts")
  # set up SingleCellExperiment object
  this_sce <- SingleCellExperiment(list(counts=this_counts))
  # run hybrid method from scds
  cxds <- cxds_bcds_hybrid(this_sce)
  # make data frame of results
  this_df <- data.frame("Barcode"=names(cxds$hybrid_score),
    "Score" = cxds$hybrid_score,
    "Sample" = s)
  doublets_df <- rbind(doublets_df, this_df)
}

# --- Doublet removal --------------------------------------------------------

doublets_df <- read.delim("https://wd.cri.uic.edu/scrna/doublet_analysis.txt")
head(doublets_df)
library(dplyr)
doublet_cells <- doublets_df %>% 
  group_by(Sample) %>%
  slice_max(order_by=Score, prop=0.05)
# count how many cells we selected per sample
table(doublet_cells$Sample)
# add a doublet column to the metadata and filter
sc_subset@meta.data$IsDoublet <- rownames(sc_subset@meta.data) %in% doublet_cells$Barcode
head(sc_subset@meta.data)
sc_subset <- subset(sc_subset, IsDoublet == F)

# --- Count the final fraction of cells remaining ----------------------------

# counts before filtering
orig.counts <- table(sc_data$orig.ident)
# counts after filtering
subset.counts <- table(sc_subset$orig.ident)
# combine and compare counts
cell_stats <- cbind("Starting Cells" = orig.counts,
  "Retained Cells" = subset.counts,
  "Fraction" = subset.counts/orig.counts)
cell_stats

# -----------------------------------------------------------------------------
# Gene feature selection and scaling
# -----------------------------------------------------------------------------

sc_subset <- readRDS(url("https://wd.cri.uic.edu/scrna/postQC.rds"))
sc_subset <- NormalizeData(sc_subset)
sc_subset <- FindVariableFeatures(sc_subset, nfeatures=4000)
VariableFeaturePlot(sc_subset)
sc_subset <- ScaleData(sc_subset)
sc_subset <- RunPCA(sc_subset, npcs = 50, verbose=F)
ElbowPlot(sc_subset, ndims=50)
DimHeatmap(sc_subset, dims=1:50, cells=300, balanced=T)

# --- JackStraw (TAKE HOME) --------------------------------------------------

sc_subset <- JackStraw(sc_subset, num.replicate = 100, dims = 50)
sc_subset <- ScoreJackStraw(sc_subset, dims = 1:50)
JackStrawPlot(sc_subset, dims = 1:50)

# -----------------------------------------------------------------------------
# Clustering
# -----------------------------------------------------------------------------

# --- Run clustering on the top PCs to filter out noise ----------------------

pca.dims <- 1:21
sc_subset <- FindNeighbors(sc_subset, dims=pca.dims)
sc_subset <- FindClusters(sc_subset, resolution=0.5)
# see where clusters are stored
head(sc_subset@meta.data)
sc_subset <- RunUMAP(sc_subset, dims=pca.dims)
# see where UMAP coordinates are stored
head(sc_subset@reductions$umap@cell.embeddings)
DimPlot(sc_subset, reduction='umap', label=T)
saveRDS(sc_subset,"clustered.rds")

# =============================================================================
# DAY 2
# =============================================================================

library(Seurat)
sc_subset <- readRDS("clustered.rds")
library(Seurat)
sc_subset <- readRDS(url("https://wd.cri.uic.edu/scrna/clustered.rds"))

# -----------------------------------------------------------------------------
# Check for cell cycle and batch biases
# -----------------------------------------------------------------------------

# --- Cell cycle analysis ----------------------------------------------------

# load gene symbols for cell cycle
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
head(s.genes)
head(g2m.genes)
genes_to_ids <- data.frame(Features(sc_subset),
  sc_subset[["RNA"]]@meta.data["gene.symbol"])
## or ##
genes_to_ids <- read.table("https://wd.cri.uic.edu/scrna/features.tsv")
s.genes_ids <- merge(s.genes, genes_to_ids, by.x=1, by.y=2)
g2m.genes_ids <- merge(g2m.genes, genes_to_ids, by.x=1, by.y=2)
head(s.genes_ids)
# cell cycle scoring
sc_subset <- CellCycleScoring(sc_subset,
  s.features=s.genes_ids[,2], g2m.features=g2m.genes_ids[,2])
head(sc_subset@meta.data)
# visualize cell cycle in UMAP
DimPlot(sc_subset, group.by="Phase")
# visualize scores versus clusters in boxplot
library(ggplot2)
ggplot(sc_subset@meta.data, aes(x=RNA_snn_res.0.5, y=S.Score)) +
  geom_boxplot()
ggplot(sc_subset@meta.data, aes(x=RNA_snn_res.0.5, y=G2M.Score)) +
  geom_boxplot()
sc_subset2 <- ScaleData(sc_subset, vars.to.regress = c("S.Score", "G2M.Score"))

# --- Batch effect analysis --------------------------------------------------

# UMAP plot by cluster and sample
DimPlot(sc_subset, group.by="seurat_clusters",split.by="orig.ident")
# Number of cells per sample per cluster
table(sc_subset@meta.data[,c("orig.ident","seurat_clusters")])
# split object based on the orig.ident, or 
# whatever batch factor is appropriate
sc_split <- SplitObject(sc_subset, split.by="orig.ident")
# run normalization and variant feature selection separately for each sample
sc_split <- lapply(sc_split, function(x){
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x)})
# find integration features
integration_genes <- SelectIntegrationFeatures(sc_split)
# run integration
integration_anchors <- FindIntegrationAnchors(sc_split, anchor.features=integration_genes)
sc_integrated <- IntegrateData(anchorset = integration_anchors)
# this creates a NEW assay, "integrated"
# use this instead of the "RNA" assay
DefaultAssay(sc_integrated) <- "integrated"

# -----------------------------------------------------------------------------
# Putative cell type identification
# -----------------------------------------------------------------------------

sc_subset <- readRDS(url("https://wd.cri.uic.edu/scrna/clustered.rds"))

# --- Marker expression ------------------------------------------------------

markers <- read.delim("https://wd.cri.uic.edu/scrna/blood_markers.txt")
markers
genenames <- data.frame(Features(sc_subset),
  sc_subset[["RNA"]]@meta.data["gene.symbol"])
## or
genenames <- read.table("https://wd.cri.uic.edu/scrna/features.tsv")
gene_to_id <- genenames[,1]
names(gene_to_id) <- make.unique(genenames[,2])
markers$ID <- gene_to_id[markers$Marker]
markers
# get unique list of genes
DotPlot(sc_subset, features=markers$ID, group.by="RNA_snn_res.0.5" ) + RotatedAxis()
# format the marker names as "gene: cell type"
#  - this will make a named vector
#  - the values are the labels to plot
#  - the names are the gene IDs
markers_description = paste0(markers$Marker,": ", markers$Cell.Type)
names(markers_description) = markers$ID
markers_description
# Add custom x labels to our plot
library(ggplot2)
DotPlot(sc_subset, features=markers$ID, group.by="RNA_snn_res.0.5" ) + 
  RotatedAxis() +
  scale_x_discrete(labels=function(x) markers_description)

# --- Label transfer from gold standard --------------------------------------

# read in Seurat object with data
sc_ref <- readRDS(url("https://wd.cri.uic.edu/scrna/sc_reference.rds"))
# run normalization and find variable features
sc_ref <- NormalizeData(sc_ref)
sc_ref <- FindVariableFeatures(sc_ref, nfeatures=4000)
# confirm that the gene identifiers match, and subset to genes in common
head(rownames(sc_ref))
head(rownames(sc_subset))
common.genes <- intersect(rownames(sc_ref),rownames(sc_subset))
length(common.genes)
# cell types we can transfer
table(sc_ref@meta.data$rough_annot)
table(sc_ref@meta.data$fine_annot)
# run the label transfer (transfer anchors) with Seurat
anchors <- FindTransferAnchors(reference=sc_ref, query=sc_subset, dims=1:30)
predictions <- TransferData(anchorset=anchors, refdata=sc_ref@meta.data$rough_annot)
# compute the majority predicted cell type per cluster
predictions$Cluster <- sc_subset@meta.data$RNA_snn_res.0.5
prediction_summary <- table(predictions[,c("Cluster","predicted.id")])
prediction_summary
prediction_type <- apply(prediction_summary,1,function(x){
  colnames(prediction_summary)[which.max(x)]})
prediction_type
# assign these cell types in the Seurat object
sc_subset@meta.data$CellType <- prediction_type[sc_subset@meta.data$RNA_snn_res.0.5]
DimPlot(sc_subset, group.by = "CellType")
# compare to the clusters
DimPlot(sc_subset, reduction='umap', label=T)

# -----------------------------------------------------------------------------
# Differential statistics
# -----------------------------------------------------------------------------

# --- Differential expression between clusters -------------------------------

degs <- FindAllMarkers(sc_subset, test.use="wilcox", min.pct=0.2)
head(degs)
genenames <- sc_subset[["RNA"]]@meta.data["gene.symbol"]
rownames(genenames) = Features(sc_subset)
## or
genenames <- read.table("https://wd.cri.uic.edu/scrna/features.tsv",row.names=1)
degs$symbol <- genenames[degs$gene, 1]
head(degs)
table(degs[degs$p_val_adj<0.05 & degs$avg_log2FC>1, "cluster"])
library(dplyr)
top5 <- degs %>% group_by(cluster) %>% slice_max(n=5, order_by=avg_log2FC)
top5
top5_cluster5 <- top5[top5$cluster==5,]
VlnPlot(sc_subset, features=top5_cluster5$gene, pt.size=F)
DoHeatmap(sc_subset, features=top5$gene)

# --- Differential expression between samples, single cell level -------------

# example analyzing for cluster 1
sc_cluster1 = subset(sc_subset, subset = seurat_clusters == 1)
Idents(sc_cluster1) = "group"
table(Idents(sc_cluster1))
cluster1_stats <- FindMarkers(sc_cluster1, test.use="wilcox",
  ident.1="Smoker", ident.2="Non-Smoker", logfc.threshold=0)
cluster1_stats$p_val_adj <- p.adjust(cluster1_stats$p_val, method="fdr")
cluster1_stats$symbol <- genenames[rownames(cluster1_stats), 1]
head(cluster1_stats)
# count number of significant changes
table(cluster1_stats$p_val_adj < 0.05)
# heatmap visualization
degs1 <- rownames(subset(cluster1_stats, p_val_adj<0.05))
norm <- GetAssayData(sc_cluster1, layer="data")
norm <- as.matrix(norm[degs1,])
norm <- t(scale(t(norm)))
library(ComplexHeatmap)
library(circlize)
col_annot <- HeatmapAnnotation(df=sc_cluster1@meta.data["group"])
col_fun <- colorRamp2(c(-2,0,2),c("blue","white","red"))
Heatmap(norm, name="Z-score", top_annotation=col_annot, col=col_fun,
  show_column_names=F, row_labels=genenames[rownames(norm),1])

# --- Differential cell type abundance ---------------------------------------

cell_counts <- table(sc_subset@meta.data[,c("seurat_clusters","orig.ident")])
# obtain metadata and match column order
metadata <- read.delim("https://wd.cri.uic.edu/scrna/GSE138867_metadata.txt",row.names=1)
metadata <- metadata[colnames(cell_counts),]
# analysis in edgeR
library(edgeR)
counts_dge <- DGEList(cell_counts, group=metadata)
counts_dge <- estimateDisp(counts_dge)
# check the BCV: it is almost always much higher for cell type abundance comparisons
sqrt(counts_dge$common.disp)
counts_results <- exactTest(counts_dge)$table
counts_results$QValue <- p.adjust(counts_results$PValue)
counts_results
percent <- cpm(counts_dge)
# convert to percent (divide CPM by a million)
percent <- percent/1000000
# look at significant changes, add group info
percent_signif <- data.frame(t(percent[counts_results$QValue<0.1,]),
  group=metadata,check.names=F)
percent_signif
# make boxplot
library(ggplot2)
library(tidyr)
percent_signif_long <- pivot_longer(percent_signif,
  cols=-group, names_to="cluster", values_to="abundance")
percent_signif_long
ggplot(percent_signif_long, aes(x=cluster, y=abundance, fill=group)) +
  geom_boxplot()

# -----------------------------------------------------------------------------
# Incorporating Feature Barcoding (FBC) Data
# -----------------------------------------------------------------------------

# --- Loading 10X data into Seurat (INSTRUCTOR DEMONSTRATION) ----------------

data_w_fbc <- Read10X("Sample_A_matrix/", gene.column=1)
names(data_w_fbc)
sc_w_fbc <- CreateSeuratObject(counts = data_w_fbc[['Gene Expression']])
sc_w_fbc[["Protein"]] <- CreateAssayObject(counts = data_w_fbc[["Antibody Capture"]])

# --- Basic visualization of FBC data ----------------------------------------

sc_w_fbc <- readRDS(url("https://wd.cri.uic.edu/scrna/sc_fbc.rds"))
DimPlot(sc_w_fbc, reduction='umap', label=T)
row.names(sc_w_fbc[["Protein"]])
VlnPlot(sc_w_fbc, 
	features = c("nFeature_RNA","nCount_RNA", "nCount_Protein"), 
	pt.size=0.2)
sc_w_fbc <- NormalizeData(sc_w_fbc, assay="Protein", normalization.method="CLR")
DefaultAssay(sc_w_fbc)
DefaultAssay(sc_w_fbc) <- "Protein"
DefaultAssay(sc_w_fbc)
row.names(sc_w_fbc[["Protein"]])
sel_ab <- c("CD3", "CD4", "CD8a", "IgG1")
FeaturePlot(sc_w_fbc, features=sel_ab, label=T)
DotPlot(sc_w_fbc, features=sel_ab)
VlnPlot(sc_w_fbc, features=sel_ab)
    sc_w_fbc <- ScaleData(sc_w_fbc, assay="Protein") 
    DoHeatmap(sc_w_fbc, features=sel_ab)
    DoHeatmap(sc_w_fbc, features=row.names(sc_w_fbc[["Protein"]]))

# Plot both RNA (gene expression) and FBC features together

sel_features <- c("protein_CD4", "rna_ENSG00000010610", 
                  "protein_CD8a", "rna_ENSG00000153563")
FeaturePlot(sc_w_fbc, features = sel_features, label=T)
DotPlot(sc_w_fbc, features = sel_features) + RotatedAxis()
VlnPlot(sc_w_fbc, features = sel_features, ncol=2)

# --- Computing basic statistics of features (_TAKE HOME EXERCISE_) ----------

# Generate cell counts for each cluster with antibody counts above a threshold. (Can use violin plots to determine cutoff)

fbc_presence <- t(sc_w_fbc[["Protein"]]@data > 1)
cluster_ids <- data.frame(Cluster=sc_w_fbc$seurat_clusters)

fbc_presence <- merge(cluster_ids, fbc_presence, by.x=0, by.y=0)
cd4_counts <- table(fbc_presence[,c("Cluster", "CD4")])

cd4_counts

# Compile a table of CD4+ cell counts and a fraction of cells in each cluster.

cd4_counts_df <- as.data.frame(cd4_counts)
cd4_counts <- merge(subset(cd4_counts_df, CD4 == "TRUE", select=c(Cluster, Freq)),
                    subset(cd4_counts_df, CD4 == "FALSE", select=c(Cluster, Freq)),
                    by.x=1, by.y=1)
colnames(cd4_counts)[2:3] <- c("Present", "Absent")
cd4_counts$Fraction <- cd4_counts$Present / ( cd4_counts$Present + cd4_counts$Absent )
cd4_counts

# -----------------------------------------------------------------------------
# Incorporating V(D)J data
# -----------------------------------------------------------------------------

# --- Get basic TCR chain statistics for each cluster. -----------------------

sc_vdj_contigs <- read.csv("https://wd.cri.uic.edu/scrna/tcr_contigs.csv")
sc_vdj_contigs$barcode <- sub('-[0-9]+$', '', sc_vdj_contigs$barcode)
# start with all of the cells
all_cells <- sc_w_fbc@meta.data["seurat_clusters"]
all_cells$TRA <- rownames(all_cells) %in% sc_vdj_contigs$barcode[sc_vdj_contigs$chain=="TRA"]
all_cells$TRB <- rownames(all_cells) %in% sc_vdj_contigs$barcode[sc_vdj_contigs$chain=="TRB"]
all_cells$TRAB <- all_cells$TRA & all_cells$TRB
head(all_cells)
table(all_cells[,c("seurat_clusters","TRAB")])

# --- Visualize TCR chain presence/absence on UMAP plot. ---------------------

umap_coords <- as.data.frame(Embeddings(sc_w_fbc, reduction="umap"))
# Take a peek at the results
head(umap_coords)
plot_data <- merge(all_cells, umap_coords, by=0)
# Take a peek at the results
head(plot_data)
library(ggplot2)
ggplot(plot_data, aes(x=umap_1, y=umap_2, color=TRAB)) +
        geom_point(size=0.5) + theme_classic() +
        guides(color=guide_legend(override.aes = list(size=3)))
ggplot(plot_data, aes(x=umap_1, y=umap_2, color=seurat_clusters)) +
        geom_point(size=0.5) + theme_classic() +
        guides(color=guide_legend(override.aes = list(size=3)))

# =============================================================================
# EXTRA (TAKE HOME) EXERCISES
# =============================================================================

BiocManager::install("scds", update=F)
BiocManager::install("SingleCellExperiment", update=F)
BiocManager::install("monocle", update=F)
install.packages('fossil')
install.packages('cowplot')
install.packages('vegan')
if ( ! requireNamespace("devtools"))
   install.packages("devtools")
devtools::install_local("PATH/TO/DOWNLOAD/CytoTRACE_0.3.3.tar.gz")
BiocManager::install("sva", update=F)
library(cowplot)
library(CytoTRACE)
library(fossil)
library(monocle)
library(scds)
library(SingleCellExperiment)
library(vegan)

# -----------------------------------------------------------------------------
# Cluster comparisons
# -----------------------------------------------------------------------------

# --- First, let's rerun clustering with two other resolutions ---------------

sc_subset <- FindClusters(sc_subset, resolution=0.1)
sc_subset <- FindClusters(sc_subset, resolution=1)
table(sc_subset@meta.data$RNA_snn_res.0.1)
table(sc_subset@meta.data$RNA_snn_res.0.5)
table(sc_subset@meta.data$RNA_snn_res.1)

# --- Compare the similarity of the clustering results at a high level with the adjusted Rand index ---

library(fossil)
clust0.1 <- sc_subset@meta.data$RNA_snn_res.0.1
clust0.5 <- sc_subset@meta.data$RNA_snn_res.0.5
clust1 <- sc_subset@meta.data$RNA_snn_res.1
adj.rand.index(clust0.1, clust0.5)
adj.rand.index(clust0.1, clust1)
adj.rand.index(clust0.5, clust1)
# building a data frame for our 3 clustering results
cluster.df <- data.frame(clust0.1, clust0.5, clust1)
# define a function to make a similarity matrix
cluster_similarity <- function( clust_df ){
   # number of columns we're comparing
   columns <- ncol(clust_df)
   # set up a similarity matrix
   rand.sim <- matrix( nrow=columns, ncol=columns )
   rownames(rand.sim) = colnames(clust_df)
   colnames(rand.sim) = colnames(clust_df)
   # set all values to be 1 first
   # this will keep the diagonal entries 1 after we do the other calculations
   rand.sim[,] <- 1
   # loop over all pairs of columns
   for( i in 1:(columns-1) ){
      for( j in (i+1):columns ){
         # compute the similarity and store it at positions i,j and j,i
         sim <- adj.rand.index( clust_df[,i], clust_df[,j] )
         rand.sim[i,j] <- sim
         rand.sim[j,i] <- sim
      }
   }
   return(rand.sim)
}
# run the function
cluster.sim <- cluster_similarity(cluster.df)
cluster.sim

# --- Detailed comparison between clustering at resolutions 0.25 and 1 using the overlap index ---

DimPlot(sc_subset, reduction='umap', group.by="RNA_snn_res.0.1", label=T)
DimPlot(sc_subset, reduction='umap', group.by="RNA_snn_res.1", label=T)
source("https://wd.cri.uic.edu/scrna/overlap_index.R")
# function for overlap index
overlap_index <- function( cluster1, cluster2 ){
   # make factor vectors from clusters
   c1 = as.factor(cluster1)
   c2 = as.factor(cluster2)
   # define a set of "names" for our cells, which we'll compare between clusters
   # the names are arbitraty, so we'll just number them
   names = 1:length(c1)
   # make distance matrix to store our calculations in
   my_dist = matrix(nrow=length(levels(c1)),ncol=length(levels(c2)))
   rownames(my_dist) = levels(c1)
   colnames(my_dist) = levels(c2)
   for(i in 1:length(levels(c1))){
      for(j in 1:length(levels(c2))){
         # get the list of cells in each cluster
         c1.sub = names[c1==levels(c1)[i]]
         c2.sub = names[c2==levels(c2)[j]]
         # get the intersection and min cluster size
         int = length(intersect(c1.sub,c2.sub))
         minsize = min(length(c1.sub),length(c2.sub))
         # overlap index
         overlap = int/minsize
         # store in matrix
         my_dist[i,j] = overlap
      }
   }
   return(my_dist)
}
res_0.1vs1 <- overlap_index(clust0.1, clust1)
res_0.1vs1
library(ComplexHeatmap)
Heatmap(res_0.1vs1, 
        col = c("white","yellow","red"), 
        name = "Overlap Index",
        column_title = "Resolution 1", 
        row_title = "Resolution 0.1")

# -----------------------------------------------------------------------------
# Additional cell typing exercises
# -----------------------------------------------------------------------------

# --- Feature plots for gene markers -----------------------------------------

library(cowplot)
# run FeaturePlot with combine=F so we get all separate plots
# 'myplots' will be a list of ggplot objects
myplots <- FeaturePlot(sc_subset, markers$ID, combine=F)
# for each plot in the list, add a new title with the description
for(i in 1:length(myplots)){
  myplots[[i]] = myplots[[i]] + ggtitle(markers_description[i])
}
# then make a new combined plot with plot_grid
plot_grid(plotlist=myplots, ncol=3)

# --- Marker-based automated cell typing with scType -------------------------

# source the scType functions from github
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")
# prepare gene sets using the marker list above
celltypes <- unique(markers$Cell.Type)
celltypes_genes <- sapply(celltypes, function(x){
  subset(markers, Cell.Type==x)$ID})
# this is a named list of gene markers for each cell type
celltypes_genes
# also include CD3 genes as negative markers for NK cells
negative_markers <- list("NK cells"=markers$ID[grep("^CD3",markers$Marker)])
# obtain normalized expression for these genes
norm <- GetAssayData(sc_subset, layer="data")
norm <- as.matrix(norm[markers$ID,])
# run scType
# scaled=F means that we are supplying normalized, but not z-scored data
#   and scType will run z-scoring then
# the "gs2" parameter is for negative markers (optional in general)
sctype_result <- sctype_score( norm, scaled=F,
  gs=celltypes_genes, gs2=negative_markers )
# transpose and peek at results
sctype_result <- t(sctype_result)
head(sctype_result)
# assign a cell type per cell based on the max score
celltypes <- apply(sctype_result, 1, function(x){
  colnames(sctype_result)[which.max(x)]})
# count the number of each cell type per cluster
cluster_types <- data.frame(
  Cluster = sc_subset@meta.data$RNA_snn_res.0.5,
  Type = celltypes)
cluster_summary <- table(cluster_types)
cluster_summary
# obtain the majority cell type per cluster
cluster_type <- apply(cluster_summary,1,function(x){
  colnames(cluster_summary)[which.max(x)]})
cluster_type
# assign these cell types in the Seurat object
sc_subset@meta.data$CellType <- cluster_type[sc_subset@meta.data$RNA_snn_res.0.5]
DimPlot(sc_subset, group.by = "CellType")
# compare to the clusters
DimPlot(sc_subset, reduction='umap', label=T)

# -----------------------------------------------------------------------------
# Additional cluster-to-cluster differential tests and visualizations
# -----------------------------------------------------------------------------

# --- Visualizations for DEGs ------------------------------------------------

library(cowplot)
myplots <- VlnPlot(sc_subset, features=top5_cluster5$gene, combine=F)
for(i in 1:length(myplots)){
  myplots[[i]] = myplots[[i]] + ggtitle(top5_cluster5$symbol[i])
}
plot_grid(plotlist=myplots)
DotPlot(sc_subset, features=top5_cluster5$gene ) +
  scale_x_discrete(labels=function(x) genenames[x,1])
DoHeatmap(sc_subset, features=top5$gene) +
  scale_y_discrete(labels=function(x) genenames[x,1])

# --- Other ways to run differential expression ------------------------------

cluster1vs2 <- FindMarkers(sc_subset, ident.1=1, ident.2=2, group.by="seurat_clusters")
cluster1vs2and3 <- FindMarkers(sc_subset, ident.1=1, ident.2=c(2,3),
  group.by="seurat_clusters")
# check we've assigned cell types from the label transfer analysis
sc_subset@meta.data$CellType <- prediction_type[sc_subset@meta.data$RNA_snn_res.0.5]
Idents(sc_subset) <- "CellType"
stats_between_types <- FindAllMarkers(sc_subset, test.use="wilcox")

# -----------------------------------------------------------------------------
# Differential stats as pseudo-bulk
# -----------------------------------------------------------------------------

library(tidyr)
sc_cluster1 <- subset(sc_subset, subset = seurat_clusters == 1)
samples <- unique(sc_cluster1$orig.ident)
counts <- GetAssayData(sc_cluster1,layer="counts")
cluster1_counts <- sapply(samples, function(x){
  # boolean vector indicating which cells are in sample x
  sample.cells <- sc_cluster1$orig.ident == x
  # sum counts over all of these cells
  sample.counts <- Matrix::rowSums(counts[,sample.cells])
  return(sample.counts)
})
# this is the pseudo-bulk counts table for this cluster
head(cluster1_counts)
# now run differential expression statistics in edgeR
library(edgeR)
# obtain metadata and match column order
metadata <- read.delim("https://wd.cri.uic.edu/scrna/GSE138867_metadata.txt",row.names=1)
metadata <- metadata[colnames(cluster1_counts),]
# subset genes to those expressed
cluster1_counts_use = cluster1_counts[rowSums(cluster1_counts)>20,]
# prepare edgeR object
cluster1_dge <- DGEList(cluster1_counts_use, group=metadata)
# normalization and dispersion
cluster1_dge <- calcNormFactors(cluster1_dge)
cluster1_dge <- estimateDisp(cluster1_dge)
# check the BCV level
sqrt(cluster1_dge$common.dispersion)
# run differential stats
cluster1_result <- exactTest(cluster1_dge)$table
cluster1_result$QValue <- p.adjust(cluster1_result$PValue)
head(cluster1_result)
# count how many are significant
table(cluster1_result$QValue < 0.05)
# comparison of DEGs here and from single-cell comparison
colnames(cluster1_result) <- paste0("Pseudobulk_",colnames(cluster1_result))
deg_comparison <- merge(cluster1_result, cluster1_stats, by=0)
subset(deg_comparison, Pseudobulk_QValue < 0.05)

# -----------------------------------------------------------------------------
# Pseudotime
# -----------------------------------------------------------------------------

# --- Load the libraries and data (if necessary) -----------------------------

# load both libraries, if not done already
library(Seurat)
library(monocle)
sc_subset <- readRDS(url("https://wd.cri.uic.edu/scrna/clustered.rds"))

# --- Get Monocle object from Seurat -----------------------------------------

# source the R script
source("https://wd.cri.uic.edu/scrna/importCDS2.R")
# use the importCDS2 function (our function) instead of importCDS (monocle function)
monocle_data <- importCDS2(sc_subset)
monocle_data

# --- Start processing in Monocle --------------------------------------------

monocle_data <- estimateSizeFactors(monocle_data)
# get variable genes from Seurat object
top.genes <- VariableFeatures(sc_subset)
monocle_data <- setOrderingFilter(monocle_data, top.genes)
monocle_data <- reduceDimension(monocle_data, max_components=2, method='DDRTree')
monocle_data <- orderCells(monocle_data)

# --- Try different visualizations -------------------------------------------

plot_cell_trajectory(monocle_data)
plot_cell_trajectory(monocle_data, color_by = "Pseudotime")
plot_cell_trajectory(monocle_data, color_by = "RNA_snn_res.0.5")

# --- Further exploration of the data. ---------------------------------------

head(pData(monocle_data))
# Get the top marker gene for cluster 7
topgene = degs$gene[degs$cluster==7][1]
# plot using gene expression as size of dot
plot_cell_trajectory(monocle_data, markers = topgene, use_color_gradient=T)
head(t(monocle_data@reducedDimS))

# -----------------------------------------------------------------------------
# Statistics with Pseudotime
# -----------------------------------------------------------------------------

# --- Look for genes that are highly correlated with pseudotime --------------

to_test <- top.genes[1:20]
monocle_subset <- monocle_data[to_test,]
pseudotime_test <- differentialGeneTest(monocle_subset, 
   fullModelFormulaStr="~sm.ns(Pseudotime)")
head(pseudotime_test)
pseudotime <- pData(monocle_data)$Pseudotime
# for time purposes, we'll do the first 20 genes again
expression <- GetAssayData(sc_subset,layer="data")[to_test,]
# run test with apply statement, making our own short function to return the 
# correlation coefficient and p-value from cor.test
# we're going to transpose the result because otherwise the genes are going
# along the columns, and keep it as a data frame
spearman.corrs <- data.frame(t(apply(expression,1,function(x){
   r=cor.test(x,pseudotime,method="spearman"); return(c(r$estimate, r$p.value))})))
colnames(spearman.corrs) = c("rho","p.value")
# add FDR correction
q.value <- p.adjust(spearman.corrs$p.value,method="fdr")
spearman.corrs <- cbind(spearman.corrs,q.value)
# sort by p-value
spearman.corrs <- spearman.corrs[order(spearman.corrs$p.value),]
# see what the final table looks like
head(spearman.corrs)

# --- Compare our cell clusters with pseudotime ------------------------------

clusters <- sc_subset@meta.data$RNA_snn_res.0.5
pseudotime <- pData(monocle_data)$Pseudotime
cluster_vs_time <- data.frame(clusters, pseudotime)
ggplot(cluster_vs_time, aes(x=clusters, y=pseudotime)) +
  geom_boxplot()
kruskal.test(pseudotime ~ clusters)

# -----------------------------------------------------------------------------
# CytoTRACE
# -----------------------------------------------------------------------------

# --- Installing CytoTRACE ---------------------------------------------------

if (! requireNamespace("BiocManager"))
  install.packages("BiocManager")

BiocManager::install("sva", update=F)
if ( ! requireNamespace("devtools"))
  install.packages("devtools")

devtools::install_local("PATH/TO/DOWNLOAD/CytoTRACE_0.3.3.tar.gz")

# --- Running CytoTRACE ------------------------------------------------------

library(Seurat)
library(CytoTRACE,verbose=F)
sc_obj <- readRDS("clustered.rds")
sc_obj.1_2 <- subset(sc_obj, subset = RNA_snn_res.0.5 == 1 | RNA_snn_res.0.5 == 2)
sc_counts <- GetAssayData(sc_obj,layer="counts")
cyto_results <- CytoTRACE(as.matrix(sc_counts), enableFast=F)

# --- Analyzing the CytoTRACE results ----------------------------------------

sc_obj@meta.data$CytoTRACE <- cyto_results$CytoTRACE
FeaturePlot(sc_obj, reduction='umap', features='CytoTRACE')
# Get the count of cells expressing the gene.   
# The code sum(x>0) is a shortcut to get a count of items greater than zero (0)
counts_sum = apply(sc_counts, 1, function(x) sum(x>0) )
# Get the genes in which the count is greater than 25% of the number of cells.
genes_keep = counts_sum > ( 0.25 * ncol(sc_counts) )
# Subset the gene counts data to just those genes.
genes_data = GetAssayData(sc_obj,layer="data")[genes_keep,]
# Take a peek at the number of genes (rows) in the filtered data.
dim(genes_data)
# Use an apply function to run the Spearman test (correlation) for each row.
# We transpose the results as the apply function will combine the output of each 
# iteration column wise and we want to have it by row and save as a data.frame
cytotrace_corr <- data.frame(t(apply(genes_data, 1, function (x) {
  res <- cor.test(x, sc_obj@meta.data$CytoTRACE, method="spearman")
  return(c(res$estimate, res$p.value))
})))

# Set the column names
colnames(cytotrace_corr) <- c("estimate", "PValue")

# Add FDR corrected PValue
cytotrace_corr$QValue <- p.adjust(cytotrace_corr$PValue, method="BH")

# Sort by absolute value of the correlation estimate
cytotrace_corr <- cytotrace_corr[order(abs(cytotrace_corr$estimate), decreasing = T), ]

# Peek at the first few results
head(cytotrace_corr)
library(dplyr)
sc_obj@meta.data %>% 
  group_by(RNA_snn_res.0.5, group) %>%
  summarize(mean=mean(CytoTRACE), sd=sd(CytoTRACE))
library(ggplot2)
ggplot(sc_obj@meta.data, aes(x=orig.ident, y= CytoTRACE, fill=group)) + 
  theme_classic() + 
  geom_boxplot() + 
  facet_wrap(~ RNA_snn_res.0.5) +
  labs(x="", y="CytoTRACE time", fill="Sample") +
  theme(axis.text.x=element_text(angle=90))
    kruskal.test(CytoTRACE ~ RNA_snn_res.0.5, data = sc_obj@meta.data)
    cluster_stats <- sc_obj@meta.data %>% group_by(RNA_snn_res.0.5) %>%
      do(w=wilcox.test(CytoTRACE ~ group, data=.)) %>%
      summarize(ClusterID=RNA_snn_res.0.5, PValue=w$p.value)
    head(cluster_stats)

# -----------------------------------------------------------------------------
# Incorporate TCR data into Seurat object
# -----------------------------------------------------------------------------

# --- Convert and load the TCR data into the Seurat object -------------------

sc_vdj_contigs <- read.csv("https://wd.cri.uic.edu/scrna/tcr_contigs.csv")
sc_vdj_contigs$barcode <- sub('-[0-9]+$', '', sc_vdj_contigs$barcode)
# First get a list of the detected chains
table(sc_vdj_contigs$chain)

# Then compute the chain counts per cell
vdj_chain_counts <- as.data.frame(table(sc_vdj_contigs[, c("barcode", "chain")]))
vdj_chain_counts$presence <- vdj_chain_counts$Freq > 0

# Take a peek at the results
head(vdj_chain_counts)
cluster_ids <- data.frame(cluster=Idents(sc_w_fbc))

# Take a peek at the results
head(cluster_ids)
library(tidyr)
vdj_cell_chains <- pivot_wider(vdj_chain_counts[, c("barcode", "chain", "Freq")],
  names_from="chain", values_from="Freq")

dim(vdj_cell_chains)

head(vdj_cell_chains)
vdj_cell_chains <- merge(cluster_ids, vdj_cell_chains, 
                         by.x=0, by.y=1, all.x=T)

dim(vdj_cell_chains)
row.names(vdj_cell_chains) <- vdj_cell_chains[,1]
vdj_cell_chains <- as.matrix(vdj_cell_chains[, c(-1, -2)])

head(vdj_cell_chains)
vdj_cell_chains[is.na(vdj_cell_chains)] <- 0
sc_w_fbc[["TCR"]] <- CreateAssayObject(counts=t(vdj_cell_chains))
saveRDS(sc_w_fbc, file="seurat_obj_w_tcr.rds")

# --- Create TCR feature plot in Seurat. -------------------------------------

DefaultAssay(sc_w_fbc) <- "TCR"
FeaturePlot(sc_w_fbc, features = c("TRA", "TRB"), label=T)
DotPlot(sc_w_fbc, features=c("TRA", "TRB"))
DotPlot(sc_w_fbc, features=c("tcr_TRA", "tcr_TRB", "protein_CD4", "rna_ENSG00000010610")) +
  RotatedAxis()

# -----------------------------------------------------------------------------
# Diversity analysis of V(D)J data
# -----------------------------------------------------------------------------

sc_vdj_contigs <- read.delim("https://wd.cri.uic.edu/scrna/filtered_contig_annotations.txt")
sc_vdj_contigs$barcode <- sub('-[0-9]+$', '', sc_vdj_contigs$barcode)
cluster_ids <- read.delim("https://wd.cri.uic.edu/scrna/cell_info.txt", row.names=1)
head(cluster_ids)
vdj_data <- subset(sc_vdj_contigs, select=c(barcode, cdr3))
head(vdj_data)
vdj_data <- subset(sc_vdj_contigs, chain=="TRA", select=c(barcode, cdr3))
head(vdj_data)
vdj_data <- data.frame(barcode=sc_vdj_contigs$barcode, 
		vdj_annotation=paste(sc_vdj_contigs$v_gene, sc_vdj_contigs$d_gene, 
			sc_vdj_contigs$j_gene, sc_vdj_contigs$c_gene, sep="."))
head(vdj_data)
tra_data <- subset(sc_vdj_contigs, chain=="TRA")
vdj_data <- data.frame(barcode=tra_data$barcode, 
		vdj_annotation=paste(tra_data$v_gene, tra_data$d_gene, 
			tra_data$j_gene, tra_data$c_gene, sep="."))
head(vdj_data)
vdj_data <- merge(cluster_ids, vdj_data, by.x=0, by.y=1)
head(vdj_data)
library(dplyr)

vdj_counts <- vdj_data %>% group_by(sample, vdj_annotation) %>% summarize(counts=n())

head(vdj_counts)
library(tidyr)
vdj_counts_mat <- pivot_wider(vdj_counts, names_from="vdj_annotation", values_from="counts", values_fill=0)

head(vdj_counts_mat)[,1:4]

# The first column has the sample IDs.
# So, set as the row names and drop the first column when converting to matrix.
# Convert to a data frame first so we can add rownames
vdj_counts_mat <- as.data.frame(vdj_counts_mat)
row.names(vdj_counts_mat) <- vdj_counts_mat[,1]
vdj_counts_mat <- as.matrix(vdj_counts_mat[,-1, drop=F])

head(vdj_counts_mat)[,1:3]
library(vegan)
vdj_diversity <- data.frame(sample=row.names(vdj_counts_mat), 
			S=diversity(vdj_counts_mat, index="shannon"),
			N=specnumber(vdj_counts_mat))

head(vdj_diversity)

# --- Comparing and visualization of diversity indices. ----------------------

mapping <- read.delim("https://wd.cri.uic.edu/scrna/mapping.txt")
head(mapping)
vdj_diversity <- merge(vdj_diversity, mapping, by.x=1, by.y=1)

# Take a peek at the results
head(vdj_diversity)
# Comparing difference in computed Shannon entropy (S) with group
kruskal.test(S ~ Group, vdj_diversity)

# Comparing difference in richness (N) with group
kruskal.test(N ~ Group, vdj_diversity)
boxplot(S ~ Group, vdj_diversity)

boxplot(N ~ Group, vdj_diversity)
library(ggplot2)

ggplot(vdj_diversity, aes(x=Group, y=S)) + geom_boxplot()

ggplot(vdj_diversity, aes(x=Group, y=N)) + geom_boxplot()
vdj_counts_200 <- rrarefy(vdj_counts_mat, 200)

# Filter out any samples that were not subsampled to a depth of 200 (have less than 200 counts)
vdj_counts_200 <- vdj_counts_200[ rowSums(vdj_counts_200) == 200, ]

# Compute the diversity indices for each sample
# Will compute Shannon entropy (S) and richness, a.k.a. number of unique entities (N)
vdj_diversity <- data.frame(sample=row.names(vdj_counts_200), 
			S=diversity(vdj_counts_200, index="shannon"),
			N=specnumber(vdj_counts_200))

head(vdj_diversity)

