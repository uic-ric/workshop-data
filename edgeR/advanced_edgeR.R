# =============================================================================
# Advanced Statistics in EdgeR
# Research Informatics Core
# Fall 2026
#
# Generated from advanced_edgeR.Rmd by purl_handouts.R -- do not edit by hand.
# =============================================================================

# =============================================================================
# DAY 1
# =============================================================================

# -----------------------------------------------------------------------------
# Prepare counts and metadata tables from RNA-seq results files
# -----------------------------------------------------------------------------

# --- Make a list of all of the files ending in "counts.txt". ----------------

fc.files <- list.files("~/../Downloads/featurecounts",
   pattern="*.counts.txt", full.names=T)
fc.files
fc.files <- list.files("~/Downloads/featurecounts",
   pattern="*.counts.txt", full.names=T)
fc.files

# --- Read each file into R --------------------------------------------------

fc.data <- lapply(fc.files,
   function(x) { read.delim(x, row.names=1, comment.char="#")[, 6, drop=F] })

# --- Create a table of the gene ID and gene counts for each sample ----------

# Confirm every file has the same genes in the same order before cbind
stopifnot(all(sapply(fc.data, function(x) identical(rownames(x), rownames(fc.data[[1]])))))

counts.table <- Reduce(cbind, fc.data)
head(counts.table, n=3)
colnames(counts.table) <- gsub("\\.bam$", "", colnames(counts.table))
head(counts.table, n=3)

# --- Write this to a new file -----------------------------------------------

counts.table <- cbind(Gene=rownames(counts.table), counts.table)
write.table(counts.table, "combined_counts.txt", sep="\t", row.names=F, quote=F)

# -----------------------------------------------------------------------------
# Pairwise differential analysis in edgeR
# -----------------------------------------------------------------------------

# --- Load the edgeR library -------------------------------------------------

library(edgeR)

# --- Load the read counts and metadata tables -------------------------------

data <- read.delim("https://wd.cri.uic.edu/edgeR/pairwise_counts.txt", row.names=1)
dim(data)
head(data, n=3)
metadata <- read.delim("https://wd.cri.uic.edu/edgeR/pairwise_metadata.txt", row.names=1)
metadata

# --- Verify the metadata row order matches the data column order ------------

data <- data[, rownames(metadata)]
dim(data)
head(data, n=3)

# --- Subset counts to genes with >20 total counts ---------------------------

data_subset <- data[rowSums(data) > 20,]
dim(data_subset)

# --- Create edgeR object ----------------------------------------------------

genes <- DGEList(counts=data_subset, group=metadata[, 1])
genes

# --- Calculate TMM factors --------------------------------------------------

genes <- calcNormFactors(genes)
genes

# --- Estimate dispersion ----------------------------------------------------

genes <- estimateDisp(genes)
genes

# --- Calculate Biological Coefficient of Variation (BCV) --------------------

sqrt(genes$common.dispersion)
plotBCV(genes)

# --- Run differential analysis ----------------------------------------------

stats <- exactTest(genes)
stats

# --- Run False Discovery Rate (FDR) correction ------------------------------

stats$table$QValue <- p.adjust(stats$table$PValue, method="BH")
head(stats$table, n=3)
sum(stats$table$QValue < 0.05)

# --- Calculate normalized expression ----------------------------------------

norm <- cpm(genes)
head(norm, n=3)

# -----------------------------------------------------------------------------
# PCA plot
# -----------------------------------------------------------------------------

norm <- read.delim("https://wd.cri.uic.edu/edgeR/pairwise_norm.txt", row.names=1)
metadata <- read.delim("https://wd.cri.uic.edu/edgeR/pairwise_metadata.txt", row.names=1)

# --- Log-scale the data first -----------------------------------------------

norm.log <- log2(norm + 0.1)

# --- Run Principal Component Analysis (PCA) ---------------------------------

pca <- prcomp(t(norm.log))

# --- Examine the summary of the PCA -----------------------------------------

summary(pca)

# --- Create labels for the X and Y axes -------------------------------------

importance <- summary(pca)$importance[2,]
xlabel <- sprintf("PC1 (%.2f%%)", 100 * importance[1])
ylabel <- sprintf("PC2 (%.2f%%)", 100 * importance[2])

# --- Format the PCA data for use with ggplot2 -------------------------------

pca.data <- cbind(pca$x, metadata, Sample=rownames(metadata))
head(pca.data, n=3)

# --- Create a PCA plot using ggplot2 ----------------------------------------

library(ggplot2)
library(ggrepel)
ggplot(pca.data, aes(x=PC1, y=PC2, label=Sample, color=Group)) +
  geom_point() + geom_text_repel(show.legend=F, color="black", seed=42) +
  labs(x=xlabel, y=ylabel)

# --- Create a Screeplot using ggplot2 (optional) ----------------------------

importance.df <- data.frame(PC=names(importance), Variance=importance)
ggplot(importance.df, aes(x=PC, y=Variance)) + geom_col() +
  labs(x="Principal Component", y="Variance Explained")

# -----------------------------------------------------------------------------
# Heatmap
# -----------------------------------------------------------------------------

stats <- readRDS(url("https://wd.cri.uic.edu/edgeR/stats.rds"))
norm.log <- read.delim("https://wd.cri.uic.edu/edgeR/pairwise_lognorm.txt", row.names=1)

# --- Load the ComplexHeatmap and circlize libraries -------------------------

library(ComplexHeatmap)
library(circlize)

# --- Get the log-scaled CPMs of differentially expressed genes --------------

degs.norm <- norm.log[stats$table$QValue < 0.05,]
dim(degs.norm)

# --- Z-score across the genes -----------------------------------------------

degs.norm.z <- t(scale(t(degs.norm)))

# --- Setup the colors we will use -------------------------------------------

col_fun <- colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))

# --- Create the heatmap -----------------------------------------------------

ht <- draw(Heatmap(degs.norm.z, show_row_names=F, name="Z-scored log2 CPM",
  col=col_fun, row_split=2, column_split=2))

# --- Get the order of genes in the heatmap plot -----------------------------

genes_order <- row_order(ht)

# Get gene names for the first cluster

cluster1_genes <- rownames(degs.norm.z)[ genes_order[[1]] ]
head(cluster1_genes, n=3)

# Use an apply statement to get the genes from all clusters

gene_clusters <- sapply(genes_order, function(x) rownames(degs.norm.z)[x])
head(gene_clusters[[1]], n=3)
head(gene_clusters[[2]], n=3)

# --- Create a tall heatmap with gene names displayed (optional) -------------

pdf("pairwise_long_heatmap.pdf", height=100)
Heatmap(degs.norm.z, show_row_names=T, name="Z-scored log2 CPM",
  col=col_fun, row_split=2, column_split=2)
dev.off()

# -----------------------------------------------------------------------------
# One-way ANOVA
# -----------------------------------------------------------------------------

# --- Load the edgeR library -------------------------------------------------

library(edgeR)

# --- Read in counts table ---------------------------------------------------

data <- read.delim("https://wd.cri.uic.edu/edgeR/anova_counts.txt", row.names=1)
dim(data)
head(data, n=3)

# --- Read in metadata table -------------------------------------------------

metadata <- read.delim("https://wd.cri.uic.edu/edgeR/anova_metadata.txt", row.names=1)
head(metadata, n=3)

# --- Make sure the metadata row order matches the data column order ---------

data <- data[, rownames(metadata)]
dim(data)

# --- Subset counts to genes with >20 total counts ---------------------------

data_subset <- data[rowSums(data) > 20,]
dim(data_subset)

# --- Define our factors -----------------------------------------------------

treat <- factor(metadata[, 1])
treat

# --- Define our model matrix ------------------------------------------------

model <- model.matrix(~treat)
model

# --- Create edgeR object ----------------------------------------------------

genes <- DGEList(counts=data_subset)

# --- Calculate TMM factors --------------------------------------------------

genes <- calcNormFactors(genes)

# --- Estimate dispersion ----------------------------------------------------

genes <- estimateDisp(genes, model)

# --- Calculate Biological Coefficient of Variation (BCV) --------------------

sqrt(genes$common.dispersion)
plotBCV(genes)

# --- Fit the model ----------------------------------------------------------

fit <- glmQLFit(genes, model)

# --- Test if either coefficient 2 or 3 has an effect ------------------------

qlf <- glmQLFTest(fit, coef=2:3)

# --- Run False Discovery Rate (FDR) correction ------------------------------

qlf$table$QValue <- p.adjust(qlf$table$PValue, method="BH")
head(qlf$table, n=3)
sum(qlf$table$QValue < 0.05)

# --- Normalize expression ---------------------------------------------------

norm <- cpm(genes, log=T)

# --- Z-score across the genes -----------------------------------------------

norm.z <- t(scale(t(norm)))

# --- Subset the z-scored table based on differentially expressed genes ------

degs.norm.z <- norm.z[qlf$table$QValue < 0.05,]

# --- Create a heatmap -------------------------------------------------------

Heatmap(degs.norm.z, show_row_names=F, name="One-way ANOVA, Z-scored log2 CPM",
        col=col_fun)

# -----------------------------------------------------------------------------
# Contrasts
# -----------------------------------------------------------------------------

# --- Create the model -------------------------------------------------------

model.alt <- model.matrix(~0 + treat)
model.alt

# --- Create a new edgeR object ----------------------------------------------

genes.alt <- DGEList(counts=data_subset)

# --- calculate TMM factors --------------------------------------------------

genes.alt <- calcNormFactors(genes.alt)

# --- Estimate dispersion ----------------------------------------------------

genes.alt <- estimateDisp(genes.alt, model.alt)

# --- Calculate BCV from the common dispersion -------------------------------

sqrt(genes.alt$common.dispersion)

# --- Fit the model ----------------------------------------------------------

fit.alt <- glmQLFit(genes.alt, model.alt)

# --- Try some contrasts -----------------------------------------------------

comp1 <- makeContrasts(treatModel1 - treatModel2, levels=model.alt)
comp2 <- makeContrasts(treatControl - (treatModel1 + treatModel2)/2, levels=model.alt)

# --- Run qlf tests for these comparisons ------------------------------------

qlf1 <- glmQLFTest(fit.alt, contrast=comp1)
qlf2 <- glmQLFTest(fit.alt, contrast=comp2)

# --- Do B-H p-value adjustment for each -------------------------------------

qlf1$table$QValue <- p.adjust(qlf1$table$PValue, method="BH")
qlf2$table$QValue <- p.adjust(qlf2$table$PValue, method="BH")

# --- Subset the z-scored data for the DEGs ----------------------------------

degs.norm.z1 <- norm.z[qlf1$table$QValue < 0.05,]
degs.norm.z2 <- norm.z[qlf2$table$QValue < 0.05,]
length(intersect(rownames(degs.norm.z1), rownames(degs.norm.z2)))

# --- Create the heatmaps ----------------------------------------------------

Heatmap(degs.norm.z1, show_row_names=F, name="Contrast 1, Z-scored log2 CPM",
        col=col_fun)
Heatmap(degs.norm.z2, show_row_names=F, name="Contrast 2, Z-scored log2 CPM",
        col=col_fun)

# =============================================================================
# DAY 2
# =============================================================================

# -----------------------------------------------------------------------------
# Two-way ANOVA
# -----------------------------------------------------------------------------

# --- Load the edgeR library -------------------------------------------------

library(edgeR)

# --- Read in counts table ---------------------------------------------------

data <- read.delim("https://wd.cri.uic.edu/edgeR/twofactor_counts.txt", row.names=1)
dim(data)
head(data, n=3)

# --- Read in metadata table -------------------------------------------------

metadata <- read.delim("https://wd.cri.uic.edu/edgeR/twofactor_metadata.txt", row.names=1)
head(metadata, n=3)

# --- Make sure the metadata row order matches the data column order ---------

data <- data[, rownames(metadata)]
dim(data)

# --- Subset counts to genes with >20 total counts ---------------------------

data_subset <- data[rowSums(data) > 20,]
dim(data_subset)

# --- Define our factors and model matrix ------------------------------------

disease <- factor(metadata[, 1])
tissue <- factor(metadata[, 2])
model <- model.matrix(~ disease + tissue + disease:tissue)

# --- Create edgeR object ----------------------------------------------------

genes <- DGEList(counts=data_subset)

# --- Calculate TMM factors --------------------------------------------------

genes <- calcNormFactors(genes)

# --- Estimate dispersion ----------------------------------------------------

genes <- estimateDisp(genes, model)

# --- Calculate BCV from the common dispersion -------------------------------

sqrt(genes$common.dispersion)
plotBCV(genes)

# --- Fit the model ----------------------------------------------------------

fit <- glmQLFit(genes, model)

# --- Test each factor in turn -----------------------------------------------

qlf.disease <- glmQLFTest(fit, coef=2:3)
qlf.tissue <- glmQLFTest(fit, coef=4)
qlf.inter <- glmQLFTest(fit, coef=5:6)

# --- Do B-H p-value adjustment ----------------------------------------------

qlf.disease$table$QValue <- p.adjust(qlf.disease$table$PValue, method="BH")
qlf.tissue$table$QValue <- p.adjust(qlf.tissue$table$PValue, method="BH")
qlf.inter$table$QValue <- p.adjust(qlf.inter$table$PValue, method="BH")
sum(qlf.disease$table$QValue < 0.05)
sum(qlf.tissue$table$QValue < 0.05)
sum(qlf.inter$table$QValue < 0.05)

# -----------------------------------------------------------------------------
# Filtering genes from two-way ANOVA
# -----------------------------------------------------------------------------

# --- Load libraries ---------------------------------------------------------

library(ggplot2)
library(ComplexHeatmap)
library(circlize)

# --- Generate log-scaled normalized expression ------------------------------

norm <- cpm(genes, log=T)

# --- Run PCA ----------------------------------------------------------------

pca <- prcomp(t(norm))

# --- Setup dataframe for ggplot ---------------------------------------------

pca.data <- cbind(pca$x, metadata, Sample=rownames(metadata))

# --- Create labels for the X and Y axes -------------------------------------

importance <- summary(pca)$importance[2,]
xlabel <- sprintf("PC1 (%.2f%%)", 100 * importance[1])
ylabel <- sprintf("PC2 (%.2f%%)", 100 * importance[2])

# --- Create combined PCA plot -----------------------------------------------

library(ggrepel)
ggplot(pca.data, aes(x=PC1, y=PC2, label=Sample, color=Disease, shape=Tissue)) + 
  geom_point() + geom_text_repel(show.legend=F, color="black", seed=42, max.overlaps=Inf) +
  labs(x=xlabel, y=ylabel)

# --- Create a Screeplot (optional) ------------------------------------------

importance.pc <- factor(names(importance), levels=names(importance))
importance.df <- data.frame(PC=importance.pc, Variance=importance)
ggplot(importance.df, aes(x=PC, y=Variance)) + geom_col() +
labs(x="Principal Component", y="Variance Explained")

# --- Z-score the normalized expression --------------------------------------

norm.z <- t(scale(t(norm)))

# --- Make a Boolean selection vector for effects due to disease or interaction term ---

disease.sel <- qlf.disease$table$QValue < 0.05
inter.sel <- qlf.inter$table$QValue < 0.05

# --- Find genes with tissue-specific effects --------------------------------

different_effect <- norm.z[disease.sel & inter.sel,]
nrow(different_effect)

# --- Find genes with disease effect, but no interaction ---------------------

shared_effect <- norm.z[disease.sel & !inter.sel,]
nrow(shared_effect)

# --- Prepare color scales and heatmap annotations ---------------------------

col_fun <- colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))
col_labels <- HeatmapAnnotation(df=metadata)

# --- Make heatmaps of each of these -----------------------------------------

Heatmap(different_effect, column_title="Different effect", name="Z-scored log2 CPM",
  col=col_fun, top_annotation=col_labels, show_row_names=F, show_column_names=F,
  column_split=tissue)

Heatmap(shared_effect, column_title="Shared effect", name="Z-scored log2 CPM",
  col=col_fun, top_annotation=col_labels, show_row_names=F, show_column_names=F,
  column_split=tissue)

# --- Make separate norm tables for each tissue ------------------------------

norm.drg <- norm[, tissue == "DRG"]
norm.sciatic <- norm[, tissue == "Sciatic_nerve"]

# --- Z-score each table -----------------------------------------------------

norm.drg.z <- t(scale(t(norm.drg)))
norm.sciatic.z <- t(scale(t(norm.sciatic)))

# --- Combine the two tables back together -----------------------------------

norm.z <- cbind(norm.drg.z, norm.sciatic.z)

# --- Reorder the columns to match the original sample order -----------------

norm.z <- norm.z[, colnames(norm)]

# --- Filter for effects -----------------------------------------------------

different_effect <- norm.z[disease.sel & inter.sel,]
shared_effect <- norm.z[disease.sel & !inter.sel,]

# --- Create heatmaps --------------------------------------------------------

Heatmap(different_effect, column_title="Different effect", name="Z-scored log2 CPM",
  col=col_fun, top_annotation=col_labels, show_row_names=F, show_column_names=F,
  column_split=tissue)

Heatmap(shared_effect, column_title="Shared effect", name="Z-scored log2 CPM",
  col=col_fun, top_annotation=col_labels, show_row_names=F, show_column_names=F,
  column_split=tissue)

# -----------------------------------------------------------------------------
# Repeated measures differential
# -----------------------------------------------------------------------------

# --- Read in counts data ----------------------------------------------------

data <- read.delim("https://wd.cri.uic.edu/edgeR/16S_counts.txt", row.names=1)
dim(data)
head(data, n=3)

# --- Read in metadata table -------------------------------------------------

metadata <- read.delim("https://wd.cri.uic.edu/edgeR/16S_metadata.txt", row.names=1)
head(metadata, n=3)
table(metadata)
nrow(metadata)
ncol(data)
colnames(data)[!colnames(data) %in% rownames(metadata)]

# --- Subset the data to match the metadata ----------------------------------

data <- data[, rownames(metadata)]
dim(data)

# --- Subset counts to taxa with >500 total counts ---------------------------

data_subset <- data[rowSums(data) > 500,]
dim(data_subset)

# --- Define factors and model matrix ----------------------------------------

treat <- factor(metadata[, 1])
subject <- factor(metadata[, 2])
model <- model.matrix(~ 0 + treat + subject)

# --- Create edgeR object and fit the model ----------------------------------

taxa <- DGEList(counts=data_subset)
taxa <- calcNormFactors(taxa)
taxa <- estimateDisp(taxa, model)
sqrt(taxa$common.dispersion)
fit <- glmQLFit(taxa, model)

# --- Run contrasts ----------------------------------------------------------

colnames(model)
contrast <- makeContrasts(treatDay5.A - treatPre, levels=model)
qlf <- glmQLFTest(fit, contrast=contrast)
qlf$table$QValue <- p.adjust(qlf$table$PValue, method="BH")

# --- Computed normalized values adjusted for repeated measures effects ------

taxa.norm <- cpm(taxa, log=T)
taxa.norm <- removeBatchEffect(taxa.norm, subject)

# --- Check the number of significant taxa. ----------------------------------

taxa.signif <- rownames(qlf$table)[qlf$table$QValue < 0.05]
length(taxa.signif)

# --- Get the normalized values for significant taxa. ------------------------

norm.signif <- taxa.norm[taxa.signif,]

# --- Z-scale normalized values ----------------------------------------------

norm.signif <- t(scale(t(norm.signif)))

# --- Taxa names are too long for display, let's fix that. -------------------

rownames(norm.signif) <- gsub("(;Other)+$", "", rownames(norm.signif))
rownames(norm.signif) <- gsub(".*;", "", rownames(norm.signif))

# --- Create a color bar for the top of the heatmap based on treatment group names. ---

group.labels <- HeatmapAnnotation(Group=treat)

# --- Create the heatmap -----------------------------------------------------

Heatmap(norm.signif, col=col_fun, name="Z-scored log2 CPM", column_split=treat,
        top_annotation=group.labels, column_title=NULL, show_column_dend=F, show_column_names=F)

# --- Create a heatmap for group 5A vs pre -----------------------------------

contrast.sel <- treat == "Day5.A" | treat == "Pre"
norm.signif.subset <- norm.signif[, contrast.sel]
subset.labels <- HeatmapAnnotation(Group=treat[contrast.sel])
Heatmap(norm.signif.subset, col=col_fun, name="Z-scored log2 CPM",
        column_split=treat[contrast.sel], top_annotation=subset.labels,
        column_title=NULL, show_column_dend=F, show_column_names=F)

# -----------------------------------------------------------------------------
# Batch effect correction
# -----------------------------------------------------------------------------

library(edgeR)

# --- Read in counts data ----------------------------------------------------

data <- read.delim("https://wd.cri.uic.edu/edgeR/batch_counts.txt", row.names=1)
dim(data)
head(data, n=3)

# --- Read in metadata table -------------------------------------------------

metadata <- read.delim("https://wd.cri.uic.edu/edgeR/batch_metadata.txt", row.names=1)
head(metadata, n=3)

# --- Make sure the metadata row order matches the data column order ---------

data <- data[, rownames(metadata)]
dim(data)

# --- Subset counts to genes with >20 total counts ---------------------------

data_subset <- data[rowSums(data) > 20,]
dim(data_subset)

# --- Define our factors -----------------------------------------------------

geno <- factor(metadata[, 1])
batch <- factor(metadata[, 2])

# --- Create a model ---------------------------------------------------------

model <- model.matrix(~ geno)
head(model, n=3)

# --- Create the edgeR object, calculate TMM factors, and estimate dispersion ---

genes <- DGEList(counts=data_subset)
genes <- calcNormFactors(genes)
genes <- estimateDisp(genes, model)

# --- Calculate BCV from the common dispersion -------------------------------

sqrt(genes$common.dispersion)

# --- Fit the model ----------------------------------------------------------

fit <- glmQLFit(genes, model)

# --- Run ANOVA: test all coefficients other than the intercept --------------

qlf <- glmQLFTest(fit, coef=2:5)
qlf$table$QValue <- p.adjust(qlf$table$PValue, method="BH")
sum(qlf$table$QValue < 0.05)

# --- Generate log-scaled normalized expression ------------------------------

norm <- cpm(genes, log=T)

# --- Run PCA ----------------------------------------------------------------

pca <- prcomp(t(norm))

# --- Setup the data for ggplot ----------------------------------------------

pca.data <- cbind(pca$x, metadata, Sample=rownames(metadata))
pca.data$Batch <- factor(pca.data$Batch)

# --- Create the plot --------------------------------------------------------

library(ggrepel)
ggplot(pca.data, aes(x=PC1, y=PC2, shape=Batch, color=Genotype, label=Sample)) +
  geom_point() + geom_text_repel(show.legend=F, color="black", seed=42)

# --- Create a new model with the batch effect correction --------------------

batch.model <- model.matrix(~ geno + batch)
head(batch.model, n=3)

# --- Create edgeR object, calculate TMM factors, and estimate dispersion ----

batch.genes <- DGEList(counts=data_subset)
batch.genes <- calcNormFactors(batch.genes)
batch.genes <- estimateDisp(batch.genes, batch.model)

# --- Calculate BCV from the common dispersion -------------------------------

sqrt(batch.genes$common.dispersion)

# --- Fit the new model ------------------------------------------------------

batch.fit <- glmQLFit(batch.genes, batch.model)

# --- Test again for the genotype effect -------------------------------------

batch.qlf <- glmQLFTest(batch.fit, coef=2:5)
batch.qlf$table$QValue <- p.adjust(batch.qlf$table$PValue, method="BH")
sum(batch.qlf$table$QValue < 0.05)

# --- Correct for batch effect -----------------------------------------------

batch.norm <- removeBatchEffect(norm, batch)

# --- Run PCA ----------------------------------------------------------------

batch.pca <- prcomp(t(batch.norm))

# --- Setup the data for ggplot ----------------------------------------------

batch.pca.data <- cbind(batch.pca$x, metadata, Sample=rownames(metadata))
batch.pca.data$Batch <- factor(batch.pca.data$Batch)

# --- Create the plot --------------------------------------------------------

ggplot(batch.pca.data, aes(x=PC1, y=PC2, shape=Batch, color=Genotype, label=Sample)) +
  geom_point() + geom_text_repel(show.legend=F, color="black", seed=42)

# -----------------------------------------------------------------------------
# Analysis with a continuous variable
# -----------------------------------------------------------------------------

# --- Read in counts table ---------------------------------------------------

data <- read.delim("https://wd.cri.uic.edu/edgeR/cont_counts.txt", row.names=1)
dim(data)

# --- Read in metadata table -------------------------------------------------

metadata <- read.delim("https://wd.cri.uic.edu/edgeR/cont_metadata.txt", row.names=1)
head(metadata, n=3)

# --- Make sure the metadata row order matches the data column order ---------

data <- data[, rownames(metadata)]
dim(data)

# --- Subset counts to genes with >20 total counts ---------------------------

data_subset <- data[rowSums(data) > 20,]
dim(data_subset)

# --- Create a histogram of the Phenotype Score ------------------------------

breaks <- seq(0, max(metadata$Phenotype_Score) + 20, 20)
ggplot(metadata, aes(x=Phenotype_Score)) +
       geom_histogram(breaks=breaks, fill="white", col="black")

# --- Check a boxplot showing any association between score and disease ------

ggplot(metadata, aes(x=Disease, y=Phenotype_Score)) +
       geom_boxplot()

# --- Run the Kruskal-Wallis rank sum test -----------------------------------

kruskal.test(Phenotype_Score ~ Disease, data=metadata)

# --- Define the variable/factor and model matrix ----------------------------

score <- as.numeric(metadata[, 1])
disease <- factor(metadata[, 2])
model <- model.matrix(~ score + disease + score:disease)
head(model, n=3)

# --- Create the edgeR object, calculate TMM factors, and estimate dispersion ---

genes <- DGEList(counts=data_subset)
genes <- calcNormFactors(genes)
genes <- estimateDisp(genes, model)

# --- Calculate BCV from the common dispersion -------------------------------

sqrt(genes$common.dispersion)

# --- Fit the model ----------------------------------------------------------

fit <- glmQLFit(genes, model)
head(fit$coefficients, n=3)

# --- Now we test each variable/factor in turn -------------------------------

qlf.score <- glmQLFTest(fit, coef=2)
qlf.disease <- glmQLFTest(fit, coef=3)
qlf.inter <- glmQLFTest(fit, coef=4)

# --- Do B-H p-value adjustment for each -------------------------------------

qlf.score$table$QValue <- p.adjust(qlf.score$table$PValue, method="BH")
qlf.disease$table$QValue <- p.adjust(qlf.disease$table$PValue, method="BH")
qlf.inter$table$QValue <- p.adjust(qlf.inter$table$PValue, method="BH")
sum(qlf.score$table$QValue < 0.05)
sum(qlf.disease$table$QValue < 0.05)
sum(qlf.inter$table$QValue < 0.05)

# =============================================================================
# EXTRA (TAKE HOME) EXERCISES
# =============================================================================

# -----------------------------------------------------------------------------
# Run all pairwise comparisons
# -----------------------------------------------------------------------------

# read in, reconicle, and filter data as usual
data <- read.delim("https://wd.cri.uic.edu/edgeR/anova_counts.txt", row.names=1)
metadata <- read.delim("https://wd.cri.uic.edu/edgeR/anova_metadata.txt", row.names=1)
data <- data[, rownames(metadata)]
data_subset <- data[rowSums(data) > 20,]
dim(data_subset)
treat <- factor(metadata[, 1])
# start a new data frame to store the results in
pw_stats <- data.frame(Gene = rownames(data_subset))
# get a list of the levels, and the number
levs <- levels(treat)
nlevs <- length(levs)
# loop over all pairs of levels
for(i in 1:(nlevs-1)){
	for(j in (i + 1):nlevs){
		# names for these groups
		A <- levs[i]
		B <- levs[j]
		# get the logical vectors for these groups
		groupA <- A == treat
		groupB <- B == treat
		# subset the data for groupA or groupB data
		pw_metadata <- treat[as.logical(groupA + groupB)]
		pw_counts <- data_subset[, as.logical(groupA + groupB)]
		# run our edgeR commands
		pw_genes <- DGEList(counts=pw_counts, group=pw_metadata)
		pw_genes <- calcNormFactors(pw_genes)
		pw_genes <- estimateDisp(pw_genes)
		pw_test <- exactTest(pw_genes)
		pw_test$table$QValue <- p.adjust(pw_test$table$PValue, method="BH")
		# add the name of this comparison to the column names
		# edgeR does the comparisons as groupB / groupA, so
		# we list groupB first
		pw_name <- paste(B, A, sep="/")
		colnames(pw_test$table) = paste(pw_name, ":", colnames(pw_test$table))
		# add results to the pairwise list
		pw_stats <- cbind(pw_stats, pw_test$table)
	}
}
head(pw_stats, n=3)
pw_tests <- function(counts, groups){
  # start a new data frame to store the results in
  pw_stats <- data.frame(Gene = rownames(counts))
  # get a list of the levels, and the number
  levs <- levels(groups)
  nlevs <- length(levs)
  # loop over all pairs of levels
  for (i in 1:(nlevs - 1)) {
    for (j in (i + 1):nlevs) {
      # names for these groups
      A <- levs[i]
      B <- levs[j]
      # get the logical vectors for these groups
      groupA <- A == groups
      groupB <- B == groups
      # subset the data for groupA or groupB data
      pw_metadata <- treat[as.logical(groupA + groupB)]
      pw_counts <- data_subset[, as.logical(groupA + groupB)]
      # run our edgeR commands
      pw_genes <- DGEList(counts = pw_counts, group = pw_metadata)
      pw_genes <- calcNormFactors(pw_genes)
      pw_genes <- estimateDisp(pw_genes)
      pw_test <- exactTest(pw_genes)
      pw_test$table$QValue <- p.adjust(pw_test$table$PValue, method = "BH")
      # add the name of this comparison to the column names
      # edgeR does the comparisons as groupB / groupA, so
      # we list groupB first
      pw_name <- paste(B, A, sep = "/")
      colnames(pw_test$table) = paste(pw_name, ":", colnames(pw_test$table))
      # add results to the pairwise list
      pw_stats <- cbind(pw_stats, pw_test$table)
    }
  }
  return(pw_stats)
}
all_pairwise <- pw_tests(data_subset, treat)
head(all_pairwise, n=3)

# -----------------------------------------------------------------------------
# Analysis with no replicates
# -----------------------------------------------------------------------------

# --- Read in counts table ---------------------------------------------------

data <- read.delim("https://wd.cri.uic.edu/edgeR/norep_counts.txt", row.names=1)
dim(data)

# --- Read in metadata table -------------------------------------------------

metadata <- read.delim("https://wd.cri.uic.edu/edgeR/norep_metadata.txt", row.names=1)
metadata

# --- Make sure the metadata row order matches the data column order ---------

data <- data[, rownames(metadata)]
dim(data)

# --- Subset counts to genes with >20 total counts ---------------------------

data_subset <- data[rowSums(data) > 20,]
dim(data_subset)

# --- Create the edgeR object and calculate TMM factors ----------------------

genes <- DGEList(counts=data_subset, group=metadata[, 1])
genes <- calcNormFactors(genes)

# --- Estimate dispersion ----------------------------------------------------

bcv <- 0.2

# --- Run differential analysis, specifying the dispersion to use ------------

stats <- exactTest(genes, dispersion=bcv^2)
stats$table$QValue <- p.adjust(stats$table$PValue, method="BH")

# --- Sort by q-value to prioritize top genes --------------------------------

head(stats$table[order(stats$table$QValue),], n=3)

# -----------------------------------------------------------------------------
# Gene ID conversion with biomaRt
# -----------------------------------------------------------------------------

if (!require("BiocManager", quietly=T))
    install.packages("BiocManager")
BiocManager::install("biomaRt")
library(biomaRt)

# --- Read in the data table -------------------------------------------------

data <- read.delim("https://wd.cri.uic.edu/edgeR/mouse_ensembl_rnaseq.txt", row.names=1)
head(data, n=3)

# --- Get rid of version numbers in transcript names -------------------------

head(rownames(data), n=3)
newnames <- gsub("\\.[0-9]*", "", rownames(data))
head(newnames, n=3)

# --- Get the database of mouse gene annotations -----------------------------

mart <- useDataset("mmusculus_gene_ensembl", useMart("ensembl"))
gene_list <- getBM(filters="ensembl_gene_id", 
   attributes=c("ensembl_gene_id", "external_gene_name"),
   values=newnames, mart=mart)
gene_list = read.delim("gene_list.txt")
head(gene_list, n=3)

# --- Note that not ALL IDs may get matched ----------------------------------

length(newnames)
nrow(gene_list)

# --- Also note that gene symbols are NOT unique relative to the gene IDs ----

length(unique(gene_list[, 1]))
length(unique(gene_list[, 2]))

# --- Replace the row names with the IDs after stripping off the version numbers ---

rownames(data) = newnames

# --- Merge based on ensembl_gene_id and rownames ----------------------------

data.named <- merge(data, gene_list, by.x="row.names", by.y="ensembl_gene_id", all=T)
head(data.named, n=3)

# -----------------------------------------------------------------------------
# Two-way ANOVA: interaction but no independent effect
# -----------------------------------------------------------------------------

any_inter <- norm.z[inter.sel,]
inter_only <- norm.z[inter.sel & !disease.sel,]
Heatmap(inter_only, col=col_fun, name="Z-scored log2 CPM",
        top_annotation=col_labels, show_row_names=F, show_column_names=F,
        column_split = tissue)
samples.drg <- tissue=="DRG"
samples.sn <- tissue=="Sciatic_nerve"
data.drg <- data_subset[, samples.drg]
data.sn <- data_subset[, samples.sn]
disease.drg <- disease[samples.drg]
disease.sn <- disease[samples.sn]
model.drg <- model.matrix(~disease.drg)
model.sn <- model.matrix(~disease.sn)
genes.drg <- DGEList(data.drg)
genes.sn <- DGEList(data.sn)
genes.drg <- calcNormFactors(genes.drg)
genes.sn <- calcNormFactors(genes.sn)
genes.drg <- estimateDisp(genes.drg, model.drg)
genes.sn <- estimateDisp(genes.sn, model.sn)
fit.drg <- glmQLFit(genes.drg, model.drg)
fit.sn <- glmQLFit(genes.sn, model.sn)
test.drg <- glmQLFTest(fit.drg, coef=2:3)
test.sn <- glmQLFTest(fit.sn, coef=2:3)
test.drg$table$QValue <- p.adjust(test.drg$table$PValue, method="BH")
test.sn$table$QValue <- p.adjust(test.sn$table$PValue, method="BH")
# count the number of DEGs
sum(test.drg$table$QValue < 0.05)
sum(test.sn$table$QValue < 0.05)
# which of these genes were detected before in the "disease" factor?
sum(test.drg$table$QValue < 0.05 & disease.sel)
sum(test.sn$table$QValue < 0.05 & disease.sel)
# prepare annotations for the heatmap based on -log10 FDR
signif.drg <- -log10(test.drg$table$QValue)
signif.sn <- -log10(test.sn$table$QValue)
labels.drg <- signif.drg[!disease.sel & inter.sel]
labels.sn <- signif.sn[!disease.sel & inter.sel]

# set an alternate color scale for significance levels
col_signif <- colorRamp2(c(0, -log10(0.05), 5), c("white", "red", "yellow"))
# use the same color scale for both annotations
anova.label <- rowAnnotation(DRG=labels.drg, SN=labels.sn,
        col=list(DRG=col_signif, SN=col_signif))
Heatmap(inter_only, col=col_fun, name="Z-scored log2 CPM",
        top_annotation=col_labels, show_row_names=F, show_column_names=F,
        column_split = tissue, left_annotation = anova.label)

# -----------------------------------------------------------------------------
# Advanced PCA plots
# -----------------------------------------------------------------------------

# --- Load the edgeR and ggplot2 libraries -----------------------------------

library(edgeR)
library(ggplot2)

# --- Read in the counts table -----------------------------------------------

data <- read.delim("https://wd.cri.uic.edu/edgeR/metagenomics_counts.txt", row.names=1)
head(rownames(data), n=3)

# --- Read in the metadata table ---------------------------------------------

metadata <- read.delim("https://wd.cri.uic.edu/edgeR/metagenomics_metadata.txt",
  row.names=1)

# --- Simplify taxonomic names to just the phylum ----------------------------

rownames(data) <- gsub("^.*;p__", "", rownames(data))
head(rownames(data), n=3)

# --- Match metadata and data ------------------------------------------------

data <- data[, rownames(metadata)]

# --- Filter out low abundance taxa ------------------------------------------

data_subset <- data[rowSums(data) > 20,]

# --- Run just a few steps in edgeR to get normalized abundance --------------

taxa <- DGEList(counts = data_subset)
taxa <- calcNormFactors(taxa)
taxa.norm <- cpm(taxa, log=T)

# --- Run PCA analysis and make a data frame with metadata -------------------

pca <- prcomp(t(taxa.norm))
pca.data <- cbind(pca$x, metadata)

# --- Create plots colored by Injury and Treatment ---------------------------

ggplot(pca.data, aes(x=PC1, y=PC2, color=Injury)) + geom_point()
ggplot(pca.data, aes(x=PC1, y=PC2, color=Treatment)) + geom_point()

# --- Facet by Treatment -----------------------------------------------------

ggplot(pca.data, aes(x=PC1, y=PC2, color=Treatment)) + geom_point() +
   facet_grid(~Injury)

# --- Plot each factor, add confidence ellipse -------------------------------

ggplot(pca.data, aes(x=PC1, y=PC2, color=Treatment)) + geom_point() +
  facet_grid(~Injury) + stat_ellipse()

# --- Add counts for Bacteroidetes and Firmicutes ----------------------------

pca.data.expanded <- cbind(pca.data,
   Bacteroidetes = t(data["Bacteroidetes",]),
   Firmicutes = t(data["Firmicutes",]))

# --- Plot colored by Firmicutes/Bacteriodetes ratio -------------------------

ggplot(pca.data.expanded, aes(x=PC1, y=PC2, color=Firmicutes / Bacteroidetes)) +
   geom_point()

# --- Plot colored by the log2-ratio -----------------------------------------

ggplot(pca.data.expanded, aes(x=PC1, y=PC2, color=log2(Firmicutes / Bacteroidetes))) +
  geom_point()

# --- Make a biplot ----------------------------------------------------------

rotation <- data.frame(taxa=rownames(pca$rotation), pca$rotation)
mult.pc1 <- max(pca.data[, "PC1"]) - min(pca.data[, "PC1"]) /
           (max(rotation[, "PC1"]) - min(rotation[, "PC1"]))
mult.pc2 <- max(pca.data[, "PC2"]) - min(pca.data[, "PC2"]) /
           (max(rotation[, "PC2"]) - min(rotation[, "PC2"]))
mult <- min(mult.pc1, mult.pc2)
rotation$PC1 = mult * rotation$PC1
rotation$PC2 = mult * rotation$PC2
plot <- ggplot(pca.data, aes(x=PC1, y=PC2, color=Treatment)) + geom_point()
plot <- plot + geom_hline(yintercept=0) + geom_vline(xintercept=0)
plot <- plot + geom_text(data=rotation, aes(x=PC1, y=PC2, label=taxa),
       size=3, hjust=0, vjust=0, color="red") + coord_equal(clip='off')
plot <- plot + geom_segment(data=rotation, aes(x=0, y=0, xend=PC1, yend=PC2),
   arrow=arrow(length=unit(0.2, "cm")), color="red")
plot

