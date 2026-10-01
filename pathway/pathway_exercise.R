# =============================================================================
# Pathway Enrichment
# Research Informatics Core
# Fall 2026
#
# Generated from pathway_analysis.Rmd by purl_handouts.R -- do not edit by hand.
# =============================================================================

# =============================================================================
# SESSION 1
# =============================================================================

# -----------------------------------------------------------------------------
# Pathway Enrichment in R
# -----------------------------------------------------------------------------

# --- Enrichment test --------------------------------------------------------

# read in differential analysis results from our server
diff <- read.delim("https://wd.cri.uic.edu/pathway/RNAseq_diff.txt")
# filter to only include significant results (QValue < 0.05)
diff_sig <- subset(diff, Infection.Control...QValue < 0.05)
# Get the list of differentially expressed genes
degs <- diff_sig[,1]

# read in the antiviral gene list from RIC github repository
antiviral <- read.delim("https://wd.cri.uic.edu/pathway/antiviral_list.txt")
# remake as a vector
antiviral <- antiviral[,1]
head(antiviral)
# read in norm table from our server
norm <- read.delim("https://wd.cri.uic.edu/pathway/RNAseq_norm.txt",
   row.names=1)
all_genes <- rownames(norm)
# check how big each list is
length(degs)
length(antiviral)
length(all_genes)
# obtain true/false vectors for all genes based on intersection
# with degs or antiviral
antiviral_list <- all_genes %in% antiviral
degs_list <- all_genes %in% degs
# we can use table to confirm the number of genes in each
table(antiviral_list)
table(degs_list)
# and we can use table to get the 2x2 contingency table
fet.table <- table(data.frame(antiviral_list, degs_list))
fet.table
# run fisher's exact test
fisher.test(fet.table)

# --- More efficient processing (bonus exercise) -----------------------------

pathway_enrich <- function( degs, pathway, all_genes ){
  pathway_list <- all_genes %in% pathway
  degs_list <- all_genes %in% degs
  fet.table <- table(data.frame(pathway_list, degs_list))
  fet <- fisher.test(fet.table)
  # return the odds ratio and p-value as a vector
  result <- c(fet$p.value, fet$estimate)
  names(result)[1] = "p.value"
  return(result)
}
pathway_enrich(degs, antiviral, all_genes)

# =============================================================================
# SESSION 2
# =============================================================================

# -----------------------------------------------------------------------------
# Pathway Enrichment in GSEA
# -----------------------------------------------------------------------------

# --- Prepare input .txt and .cls files --------------------------------------

norm <- read.delim("https://wd.cri.uic.edu/pathway/RNAseq_norm.txt", row.names=1, 
                   colClasses=c(Gene.name="NULL"))
norm <- norm[rowMeans(norm) > 1,]
norm <- log2(norm + 0.1)
norm <- cbind(rownames(norm), rownames(norm), norm)
colnames(norm)[1:2] = c("NAME","DESCRIPTION")
write.table(norm, "RNAseq_for_gsea.txt", col.names=T, row.names=F, quote=F, sep="\t")
groups <- c(rep("Control",4),rep("Infection",4))
write.table(t(c(8, 2, 1)), "RNAseq_for_gsea.cls", col.names=F, row.names=F)
write.table("# Control Infection", "RNAseq_for_gsea.cls", col.names=F, row.names=F,
  append=T, quote=F)
write.table(t(groups), "RNAseq_for_gsea.cls", col.names=F, row.names=F, append=T, quote=F)

# -----------------------------------------------------------------------------
# Barplot visualization
# -----------------------------------------------------------------------------

# read in the results from DAVID
# the quote="" helps us to parse pathway names that have ' or " in them
david <- read.delim("https://wd.cri.uic.edu/pathway/david_chart.txt",
  quote="")
# check how the names are interpreted in R
colnames(david)
# sort by significance
david <- david[order(david$FDR),]
# for now, we'll focus on the KEGG pathways with FDR < 1%
david.kegg <- subset(david, Category == "KEGG_PATHWAY" & FDR < 0.01)
# log-scale the FDR and the enrichment ratio
david.kegg$logFDR <- -log10(david.kegg$FDR)
# the KEGG pathway names have IDs in them, fix it so that we just plot the name
head(david.kegg$Term)
david.kegg$Term <- gsub("mmu[0-9]*:","",david.kegg$Term)
head(david.kegg$Term)
# load the ggplot2 library
library(ggplot2)

# --- Plot $-log_{10}$ FDR ---------------------------------------------------

ggplot(david.kegg, aes(x=Term, y=logFDR)) +
  geom_col() +
  coord_flip() +
  labs(y = "-log10 FDR") +
  scale_x_discrete(limits = rev(david.kegg$Term))
library(scales)

# We only really need to do this once, but this will allow us 
# to quickly create a reverse log scale in ggplot2
revlog_trans <- trans_new("revlog",
                          function (x) { -1 * log10(x) },
                          function (x) { 10 ^ (-1 * x) }, 
                          breaks=log_breaks(base=10),
                          domain=c(1e-100, Inf))

ggplot(david.kegg, aes(x=Term, y=FDR)) + geom_col() + coord_flip() +
  scale_y_continuous("FDR corrected p-value", trans=revlog_trans)

# --- Plot the fold enrichment -----------------------------------------------

ggplot(david.kegg, aes(x=Term, y=Fold.Enrichment)) + geom_col() + coord_flip() + 
  labs(y="Fold Enrichment")

# -----------------------------------------------------------------------------
# Take home exercise: Pathway comparison plots
# -----------------------------------------------------------------------------

# --- Prepare combined data set ----------------------------------------------

kegg1 <- read.delim("https://wd.cri.uic.edu/pathway/cluster1_KEGG.txt",
  quote="")
kegg2 <- read.delim("https://wd.cri.uic.edu/pathway/cluster2_KEGG.txt",
  quote="")
kegg3 <- read.delim("https://wd.cri.uic.edu/pathway/cluster3_KEGG.txt",
  quote="")
kegg1_subset <- kegg1[,c("Term", "FDR")]
kegg2_subset <- kegg2[,c("Term", "FDR")]
kegg3_subset <- kegg3[,c("Term", "FDR")]
colnames(kegg1_subset)[2] <- "Cluster1"
colnames(kegg2_subset)[2] <- "Cluster2"
colnames(kegg3_subset)[2] <- "Cluster3"
kegg_merged <- Reduce(function(x,y) merge(x=x, y=y, by="Term", all=T),
  list(kegg1_subset, kegg2_subset, kegg3_subset))
head(kegg_merged)
kegg_merged[is.na(kegg_merged)] <- 1
kegg_subset <- kegg_merged[ apply(kegg_merged[,-1], 1, min) < 0.05, ]

kegg_subset$Term <- gsub("mmu[0-9]*:", "", kegg_subset$Term)
head(kegg_subset)

# --- Plot as side-by-side barplots ------------------------------------------

kegg_subset2 <- kegg_subset[order(apply(kegg_subset[,-1], 1, min)), ]
library(tidyr)
head(kegg_subset2)
kegg_long <- pivot_longer(kegg_subset2, !Term, names_to = "cluster")
head(kegg_long)
library(ggplot2)
library(scales)

ggplot(kegg_long, aes(x=Term, y=value, fill=cluster)) + 
  geom_col(position='dodge') + coord_flip() +
  scale_x_discrete(limits=rev(kegg_subset2$Term)) + 
  scale_y_continuous("FDR corrected p-value", trans=revlog_trans)
ggplot(kegg_long, aes(x=Term, y=value, fill=cluster)) + 
  geom_col(position='dodge') + coord_flip() +
  scale_x_discrete(limits=rev(kegg_subset2$Term)) + 
  scale_y_continuous("FDR corrected p-value", trans=revlog_trans) +
  geom_hline(yintercept = 0.05, color="red")

# --- Plot as a heatmap ------------------------------------------------------

library(circlize)
rownames(kegg_subset)<-kegg_subset$Term
kegg_subset$Term<-NULL
kegg_subset_revlog <- -log10(kegg_subset)
col_fun <- colorRamp2(c(0, -log10(0.05), max(unlist(kegg_subset_revlog))),
  c("white","yellow","red"))
library(ComplexHeatmap)
heatmap <- Heatmap(as.matrix(kegg_subset_revlog), name="-log10 FDR", col=col_fun, 
  row_names_gp = gpar(fontsize = 5))
draw(heatmap, heatmap_legend_side = "left")
# We saw the breaks were at 0, 5, 10, and 15
breaks <- c(0, 5, 10, 15)

heatmap <- Heatmap(as.matrix(kegg_subset_revlog), 
                   col=col_fun, 
                   heatmap_legend_param = list(
                     at=breaks,
                     labels=10 ^ (-1 * breaks),
                     title="FDR"),
                   row_names_gp = gpar(fontsize = 5))
draw(heatmap, heatmap_legend_side = "left")

# -----------------------------------------------------------------------------
# Heatmap of expression patterns in pathway
# -----------------------------------------------------------------------------

# --- Read in data sets ------------------------------------------------------

# read in the results from DAVID and RNA-seq normalized data
david <- read.delim("https://wd.cri.uic.edu/pathway/david_chart.txt",
  quote="")
# sort by significance
david <- david[order(david$FDR),]
david.kegg <- david[david$Category=="KEGG_PATHWAY" & david$FDR < 0.01,]
# read in normalized expression
norm <- read.delim("https://wd.cri.uic.edu/pathway/RNAseq_norm.txt",
  row.names=1, colClasses=c(Gene.name="NULL"))

# --- Make heatmap for top KEGG pathway --------------------------------------

# name of the top pathway
top.kegg.name <- david.kegg[1,"Term"]
# remove the KEGG ID
top.kegg.name <- gsub("mmu[0-9]*:","",top.kegg.name)
# get the gene list for the top pathway
# first remove the white space
top.kegg.genes <- gsub(" ", "", david.kegg[1,"Genes"])
# use string split to split the list into a vector by commas
# use unlist to turn it back into a vector
top.kegg.genes <- unlist(strsplit(top.kegg.genes,","))
length(top.kegg.genes)
# subset the norm table to these genes
top.kegg.norm <- norm[top.kegg.genes,]
dim(top.kegg.norm)
# log-scale and z-score
top.kegg.norm <- log2(top.kegg.norm + 0.1)
top.kegg.norm <- t(scale(t(top.kegg.norm)))
# plot in a heatmap
library(ComplexHeatmap)
Heatmap(top.kegg.norm, name = top.kegg.name, show_row_names = FALSE)

# -----------------------------------------------------------------------------
# Take home exercise: Pathway analysis in variant calling data
# -----------------------------------------------------------------------------

# --- Take home exercise: Filter variants in R -------------------------------

vars <- read.delim("https://wd.cri.uic.edu/pathway/example_variants.txt")
# look at vars in the R studio variable browser, or with head
head(vars)
# filtering strategy 1: damaging effects
# variants with SIFT prediction D for damaging, just the Gene column
filter1 <- as.character(vars[vars$SIFT_pred=="D","Gene"])
# some gene annotations have a comma-separated list of genes
# strsplit will split them, and unlist will give it all as a vector
filter1 <- unlist(strsplit(filter1,","))
filter1 <- unique(filter1)
head(filter1)
length(filter1)
# filtering strategy 2: mutation burden
# remove variants with synonymous or unknown change to translated sequence
filter2_all <- as.character(vars[vars$ExonicFunc!="synonymous_SNV" &
  vars$ExonicFunc!="unknown","Gene"])
# split by commas again
filter2_all <- unlist(strsplit(filter2_all,","))
# generate a count per gene using the table function
filter2_table <- table(filter2_all)
head(filter2_table)
# get the gene names from the table with counts bigger than 3
filter2 <- names(filter2_table)[filter2_table > 3]
head(filter2)
length(filter2)
# check the number of genes in common between the lists
length(intersect(filter1, filter2))
# write both lists to tables
write.table(filter1,"filter_damaging.txt",col.names=F,row.names=F,quote=F)
write.table(filter2,"filter_mutation_load.txt",col.names=F,row.names=F,quote=F)

# --- Take home exercise: Comparison visualization in R ----------------------

# read in tables; remember quote="" helps with parsing
#   pathway names with quotes in their name
gobp.damaging <- read.delim(
  "https://wd.cri.uic.edu/pathway/GOBP_damaging.txt",
  quote="")
gobp.mutation_load <- read.delim(
  "https://wd.cri.uic.edu/pathway/GOBP_mutation_load.txt",
  quote="")
# process as we did before:
# subset the columns
sub.damaging <- gobp.damaging[,c("Term","FDR")]
sub.mutation_load <- gobp.mutation_load[,c("Term","FDR")]
# name the FDR column based on the gene list
colnames(sub.damaging)[2] <- "Damaging.Variants"
colnames(sub.mutation_load)[2] <- "Mutation.Load"
# merge and replace missing values with FDR = 1
gobp.merged <- merge(x=sub.damaging, y=sub.mutation_load, by="Term", all=T)
gobp.merged[is.na(gobp.merged)] <- 1
# prepare as a data frame with term IDs in the rownames
gobp.df <- data.frame(gobp.merged[,c(2:ncol(gobp.merged))])
rownames(gobp.df) <- gobp.merged[,1]
# subset to the significant terms and log-scale FDRs
gobp.subset <- gobp.df[apply(gobp.df, 1, min) < 0.05, ]
gobp.subset <- -log10(gobp.subset)
# remove the GO IDs from the term descriptions
rownames(gobp.subset) <- gsub("GO:[0-9]*~","",rownames(gobp.subset))
# plot in a heatmap, setting a color scale first
library(circlize)
col_fun <- colorRamp2(c(0, -log10(0.05), max(unlist(gobp.subset))),
  c("white","yellow","red"))
library(ComplexHeatmap)
heatmap <- Heatmap(as.matrix(gobp.subset), name="-log10 FDR", col=col_fun,
  row_names_gp = gpar(fontsize = 5))
draw(heatmap, heatmap_legend_side = "left")

