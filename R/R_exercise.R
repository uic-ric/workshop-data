# =============================================================================
# R Essentials for Data Analysis
# Research Informatics Core
# Fall 2026
#
# Generated from R_exercise.Rmd by purl_handouts.R -- do not edit by hand.
# =============================================================================

# =============================================================================
# DAY 1
# =============================================================================

# -----------------------------------------------------------------------------
# Basic usage
# -----------------------------------------------------------------------------

# --- Creating variables -----------------------------------------------------

a <- 100
b <- "hello"
x <- TRUE
a
class(b)

# --- Multi-valued (collection) variables ------------------------------------

x <- c(6,5,4)
x
y <- c("mouse", "cat", "dog")
y
z <- list(4, "cat", TRUE, c(1,2))
z
A <- matrix(c(1:6), nrow=3)
A
df <- data.frame(x ,y)
df
length(x)
dim(A)
nrow(A)

# --- Basic operations -------------------------------------------------------

1 + 1
12/3
A
A + A 
A * A
TRUE & FALSE
T | F
y
"cat" %in% y
"kitten" %in% y
y
order(y)
y[order(y)] 
df[ order(df[,2]), ]
A
t(A)

# -----------------------------------------------------------------------------
# Combining and extracting elements/items from objects
# -----------------------------------------------------------------------------

c( c(1,2,3), 4:8 )
c( list("a","b"), list(1,c(2,3)) )
cbind(A, A)
rbind(A, A)
x[1]
x[c(1,2)]
A[1,1]
A[1:2, 1]
A[1:2, 1, drop=F]
A[1:2, ]
df[1]
df[[1]]

# --- Naming -----------------------------------------------------------------

x
names(x) <- c("first","second","third")
x
x["second"]
newx <- c("cat"=1, "dog"=2)
newx
rownames(df) <- c("r1","r2","r3")
colnames(df) <- c("size","animal")
df
df["r1","animal"]
df[c("r1","r3"), ] 
df$animal 

# --- Factors ----------------------------------------------------------------

f <- factor(c("WT","WT","KO1","KO1","KO2","KO2"))
f
levels(f)
levels(f)[3] <- "Wild-Type"
f
f[6] <- "KO3" 
f
f <- factor(c("WT","WT","KO1","KO1","KO2","KO2"), levels=c("WT","KO1","KO2"))
f
f2 <- factor(c("time1","time2","time3"), ordered=T)
f2
df
df[,2] <- factor(df[,2])
df
df[,2]

# -----------------------------------------------------------------------------
# Writing functions and using apply
# -----------------------------------------------------------------------------

# --- Write a function to find the max value of a vector. --------------------

get.max <- function(x){
	max <- x[1]
	for (i in x){
		if (i > max){
			max <- i
		}
	}
	return(max)
}
a <- c(23.3, 1, 3, 55, 6)
get.max(a)
max(a)

# --- apply statements -------------------------------------------------------

?mtcars
dim(mtcars)
head(mtcars)
apply(mtcars, 2, max) 
tapply(mtcars$mpg, mtcars$cyl, mean)
# get the cylinders
cyls <- unique(mtcars$cyl)
cyls
cyls <- cyls[order(cyls)]
cyls
lapply(cyls, function (x) sum(mtcars$cyl==x))
sapply(cyls, function (x) sum(mtcars$cyl==x))

# -----------------------------------------------------------------------------
# Read and write a file
# -----------------------------------------------------------------------------

bw_data  <- read.delim("https://wd.cri.uic.edu/R/birth_weight.txt")
head(bw_data)
data_ordered_by_age <- bw_data[order(bw_data$age), ]
head(data_ordered_by_age)
write.table(data_ordered_by_age,"birth_weight_ordered_by_age.txt",sep="\t",
   quote=F,row.names=F)

# -----------------------------------------------------------------------------
# Install CRAN and Bioconductor packages
# -----------------------------------------------------------------------------

# --- Install CRAN package: "tidyverse" --------------------------------------

install.packages("tidyverse")
library(tidyr)
library(dplyr)
library(ggplot2)

# --- Install Bioconductor package: "ComplexHeatmap" -------------------------

install.packages("BiocManager") 
BiocManager::install("ComplexHeatmap", update=F)
library(ComplexHeatmap)

# --- TAKE-HOME: How to check which packages are already installed -----------

installed.packages()
c("ggplot2", "ComplexHeatmap") %in% rownames(installed.packages())

# -----------------------------------------------------------------------------
# Cleaning and transforming a sample dataset
# -----------------------------------------------------------------------------

# --- Initial setup ----------------------------------------------------------

library(tidyverse)
head(mtcars)
data <- mtcars %>%
  rownames_to_column(var="car")
head(data)

# --- Examples ---------------------------------------------------------------

data %>%
  select(car, mpg, cyl, hp) %>% head()
data %>%
  filter(mpg > 25) %>% head()
data %>%
  mutate(power_to_weight=hp / wt) %>% head()
data %>%
  group_by(cyl) %>%
  summarise(mean_mpg=mean(mpg))
data %>%
  count(gear)
data_long <- select(data, c(car, mpg, hp, wt)) %>%
  pivot_longer(cols=c(mpg, hp, wt),
               names_to="variable",
               values_to="value")
data_long
data_long %>%
  pivot_wider(names_from=variable, values_from=value)

# --- TAKE-HOME: Additional examples -----------------------------------------

data %>%
  rename(miles_per_gallon=mpg,
         horsepower=hp) %>% head()
data %>%
  filter(cyl %in% c(4, 6)) %>% head()
data %>%
  mutate(engine_size=cut(disp, breaks=c(0, 150, 300, Inf),
                           labels=c("small", "medium", "large"),
                           right=F)) %>% head()
data_with_na <- data %>%
  mutate(mpg=replace(mpg, mpg < 20, NA))
head(data_with_na)
data_with_na %>%
  drop_na() %>% head()
data_with_na %>%
  replace_na(list(mpg=0)) %>% head()

# =============================================================================
# DAY 2
# =============================================================================

# -----------------------------------------------------------------------------
# Data Visualization
# -----------------------------------------------------------------------------

# --- Data set ---------------------------------------------------------------

sc <- read.delim("http://wd.cri.uic.edu/R/scRNA_cells.txt", row.names=1)
sc$Cluster <- factor(sc$Cluster)
head(sc)

# --- Scatterplots -----------------------------------------------------------

plot(sc$UMAP_1, sc$UMAP_2)
cell_colors <- factor(sc$Cluster)
levels(cell_colors) <- rainbow(length(levels(cell_colors)))
plot(sc$UMAP_1, sc$UMAP_2, col=cell_colors)
library(ggplot2)
ggplot(sc, aes(x=UMAP_1, y=UMAP_2, color=Cluster)) +
  geom_point()
pdf("basic_UMAP_plot.pdf")
ggplot(sc, aes(x=UMAP_1, y=UMAP_2, color=Cluster)) +
  geom_point()
dev.off()
ggplot(sc, aes(x=UMAP_1, y=UMAP_2, color=Cluster)) +
  geom_point() +
  facet_wrap( ~ Genotype)
ggplot(sc, aes(x=UMAP_1, y=UMAP_2, color=Cluster)) +
  geom_point() +
  facet_grid( Batch ~ Genotype)

# --- Box plots and violin plots ---------------------------------------------

ggplot(sc, aes(x=Cluster, y=Cxcr6, fill=Genotype)) +
  geom_boxplot()
ggplot(sc, aes(x=Cluster, y=Cxcr6, fill=Genotype)) +
  geom_violin()
baseplot <- ggplot(sc, aes(x=Cluster, y=Cxcr6, fill=Genotype))
baseplot + geom_boxplot()
baseplot + geom_violin()

# --- Bar plots --------------------------------------------------------------

ggplot(sc, aes(x=Cluster)) + geom_bar()
cluster_counts <- table(sc$Cluster)
cluster_counts
cluster_counts_df <- data.frame(cluster_counts)
cluster_counts_df
ggplot(cluster_counts_df, aes(x=Var1, y=Freq)) + geom_col()
ggplot(sc, aes(x=Cluster, fill=Genotype)) +
  geom_bar()
ggplot(sc, aes(x=Cluster, fill=Genotype)) +
  geom_bar(position="dodge")

# --- TAKE HOME: Additional plot examples ------------------------------------

ggplot(sc, aes(x=UMAP_1, y=UMAP_2, color=Cluster)) +
  geom_point() +
  theme_classic() +
  labs(x="UMAP1", y="UMAP2", color="Cell type", title="UMAP by cluster")
ggplot(sc, aes(x=UMAP_1, y=UMAP_2, color=Cxcr6)) +
  geom_point()
ggplot(sc, aes(x=UMAP_1, y=UMAP_2, color=Cxcr6)) +
  geom_point() +
  scale_color_viridis_c()
ggplot(sc, aes(x=UMAP_1, y=UMAP_2, color=Genotype)) +
  geom_point() +
  scale_color_manual(values=c("WT"="red","KO"="blue"))
gene_expr <- sc[,c(7,10:17)]
head(gene_expr)
library(tidyr)
gene_expr_long <- pivot_longer(gene_expr, !Cluster)
gene_expr_long
basic_vioplot <- ggplot(gene_expr_long, aes(x=NA,y=value,fill=Cluster)) +
  geom_violin() + facet_grid(Cluster ~ name)
basic_vioplot
basic_vioplot + 
  theme_void() + 
  coord_flip() +
  guides(fill="none") +
  theme(strip.text.x=element_text(angle=45))

# -----------------------------------------------------------------------------
# Heatmaps
# -----------------------------------------------------------------------------

library(ComplexHeatmap)
degs <- read.delim("http://wd.cri.uic.edu/R/degs.txt",row.names=1)
head(degs)
degs.z <- t(scale(t(degs)))
Heatmap(degs.z, name="Z-score")
Heatmap(degs.z, name="Z-score", show_row_names=F)

# --- TAKE HOME: Heatmap annotations -----------------------------------------

groups <- gsub("\\.[0-9]*$","",colnames(degs))
groups
group_colors <- c("Control"="blue","Model1"="orange","Model2"="purple")
column_label <- HeatmapAnnotation(df=data.frame(Group=groups),
  col=list(Group=group_colors), which="column")
Heatmap(degs.z, name="Z-score", show_row_names=F, top_annotation=column_label)

# -----------------------------------------------------------------------------
# PCA
# -----------------------------------------------------------------------------

pca <- prcomp(t(degs))
summary(pca)
pca_stats <- summary(pca)
names(pca_stats) 
importance <- pca_stats$importance
importance
pca_coords <- data.frame(pca$x)
head(pca_coords)
pca_coords$Group <- groups
ggplot(data=pca_coords, aes(x=PC1, y=PC2, color=Group)) +
    geom_point() +
    theme_classic()
screeplot(pca)

# --- TAKE HOME: Additional PCA plots ----------------------------------------

ggplot(data=pca_coords, aes(x=PC1, y=PC2, color=Group)) +
    geom_point() +
    theme_classic() +
    labs(x=paste0("PC1 (", round(100*importance[2,1], 1), "%)"),
         y=paste0("PC2 (", round(100*importance[2,2], 1), "%)"))
install.packages("ggrepel")
library(ggrepel)
pca_coords$Sample <- rownames(pca_coords)
ggplot(data=pca_coords, aes(x=PC1, y=PC2, color=Group, label=Sample)) +
    geom_point() +
    theme_classic() +
    geom_text_repel(show.legend=F) +
    labs(x=paste0("PC1 (", round(100*importance[2,1], 1), "%)"),
         y=paste0("PC2 (", round(100*importance[2,2], 1), "%)"))
import_df <- as.data.frame(t(importance)) %>% 
  rownames_to_column("axis") %>%
  mutate(axis=factor(axis, levels=colnames(importance)))

colnames(import_df)[2:4] <- c("sdev", "var", "cum_var")

ggplot(import_df, aes(x=axis, y=var)) + 
  geom_col(fill="gray", color="black") + 
  labs(y="Proportion of Variance") +
  theme_classic()

# -----------------------------------------------------------------------------
# Differences between groups (t-test, ANOVA)
# -----------------------------------------------------------------------------

# --- Data set ---------------------------------------------------------------

bw_data <- read.delim("https://wd.cri.uic.edu/R/birth_weight.txt")

# --- t-test and Wilcox test -------------------------------------------------

bw_data$smoke <- factor(bw_data$smoke)
t.test(bwt ~ smoke, data=bw_data)
ttest_result <- t.test(bwt ~ smoke, data=bw_data)
ttest_result
ttest_result$p.value
ttest_result$conf.int
wilcox.test(bwt ~ smoke, data=bw_data)
ggplot(bw_data, aes(x=smoke, y=bwt)) +
  geom_boxplot()

# --- One-way ANOVA and Kruskal-Wallis test ----------------------------------

library(dplyr)
bw_data <- bw_data %>%
  mutate(gestation_cat=cut(gestation, 
                             breaks=c(0, 37*7, 41*7, Inf), 
                             labels=c("pre term", "full term", "late term"),
                             right=F))
anova <- aov(bwt ~ gestation_cat, data=bw_data)
summary(anova)
anova_result <- summary(anova)
anova_pvalues <- anova_result[[1]][["Pr(>F)"]]
anova_pvalues
plot(anova)

# --- Kruskal-Wallis test ----------------------------------------------------

kruskal.test(bwt ~ gestation_cat, data=bw_data)
ggplot(bw_data, aes(x=gestation_cat, y=bwt)) +
  geom_boxplot()

# --- Two-way ANOVA ----------------------------------------------------------

bw_data$smoke <- factor(bw_data$smoke)
summary(aov(bwt ~ gestation_cat + smoke, data=bw_data))
summary(aov(bwt ~ gestation_cat * smoke, data=bw_data))
ggplot(bw_data, aes(x=gestation_cat, fill=smoke, y=bwt)) +
  geom_boxplot()

# --- TAKE HOME: Additional examples with t-test and ANOVA -------------------

length(bw_data$bwt)
hist(bw_data$bwt)
shapiro.test(bw_data$bwt)
qqnorm(bw_data$bwt)
qqline(bw_data$bwt, col="red")
raw.count <- read.delim("http://wd.cri.uic.edu/advanced_R/raw.count.txt",row.names=1)
count <- as.numeric(raw.count[1, ])
length(count)
hist(count)
shapiro.test(count)
qqnorm(count)
qqline(count, col="red")
# save the wilcox test object
wilcox_result <- wilcox.test(bwt ~ smoke, data=bw_data)
# plot with a subtitle
ggplot(bw_data, aes(x=smoke, y=bwt)) +
  geom_boxplot() +
  labs(title="Baby's birth weight vs. Mother's smoking status",
       subtitle=paste0("t-test, p=", signif(ttest_result$p.value, 3), 
                       ", Wilcox p=", signif(wilcox_result$p.value, 3)))

# -----------------------------------------------------------------------------
# Correlation and linear regression
# -----------------------------------------------------------------------------

# --- Correlation ------------------------------------------------------------

cor.test(bw_data$height, bw_data$weight)
cor.test(bw_data$height, bw_data$weight, method="spearman")
cor.test(bw_data$height, bw_data$weight, method="kendall")
ggplot(bw_data, aes(x=height, y=weight)) +
  geom_point()

# --- Linear regression ------------------------------------------------------

summary(lm(bwt ~ gestation, data=bw_data))
ggplot(bw_data, aes(x=gestation, y=bwt)) +
  geom_point()
ggplot(bw_data, aes(x=gestation, y=bwt)) +
  geom_point() +
  geom_smooth(method="lm", se=T)

# --- TAKE HOME: Multiple regression -----------------------------------------

summary(lm(bwt ~ gestation + smoke + weight, data=bw_data))
summary(lm(bwt ~ gestation * smoke * weight, data=bw_data))
summary(aov(bwt ~ gestation * smoke * weight, data=bw_data))
ggplot(bw_data, aes(x=gestation, y=bwt, color=smoke)) +
  geom_point() +
  geom_smooth(method="lm", se=T)
ggplot(bw_data, aes(x=weight, y=bwt, color=smoke)) +
  geom_point() +
  geom_smooth(method="lm", se=T)

# -----------------------------------------------------------------------------
# Fisher's Exact and Chi-squared tests
# -----------------------------------------------------------------------------

# --- Fisher's Exact test ----------------------------------------------------

bw_data <- bw_data %>% 
  mutate(smoke=factor(smoke), 
         parity=factor(parity))
head(bw_data[, c("smoke","parity")])
smoke_vs_parity_table <- table(bw_data[, c("smoke", "parity")])
smoke_vs_parity_table
fisher.test(smoke_vs_parity_table)
fet <- fisher.test(smoke_vs_parity_table)
log2(fet$estimate)
ggplot(bw_data, aes(x=smoke, fill=parity)) +
  geom_bar()

ggplot(bw_data, aes(x=smoke, fill=parity)) +
  geom_bar(position="fill") +
  labs(y="Fraction")

# --- TAKE HOME: Chi-square test ---------------------------------------------

bw_data <- bw_data %>%
  mutate(gestation_cat=cut(gestation, 
                             breaks=c(0, 37*7, 41*7, Inf), 
                             labels=c("pre term", "full term", "late term"),
                             right=F))
smoke_vs_gest_table <- table(bw_data[,c("smoke","gestation_cat")])
chisq.test(smoke_vs_gest_table)
ggplot(bw_data, aes(x=smoke, fill=gestation_cat)) +
  geom_bar(position="fill") +
  labs(y="Fraction")
install.packages("vcd")
library(vcd)
mosaic(smoke_vs_gest_table, shade=T, legend=T, gp=shading_binary)
# store groups to run tests over
cols <- ncol(smoke_vs_gest_table)
# start an empty data frame
term_stats <- data.frame()
# run a test for each column (term group) one at a time
for(i in 1:cols){
  # get the counts for this term group
  term.counts <- smoke_vs_gest_table[,i]
  # get the counts for all other term groups
  term.other <- rowSums(smoke_vs_gest_table[,-i])
  # we're specifically seeing if THIS term group has a different distribution
  # of smoking vs non-smoking mothers, using fisher's exact test
  term.table <- cbind(term.counts, term.other)
  fet <- fisher.test(term.table)
  # combine the estimate (odds ratio) and p-value) into the data frame
  # and use the term group as the row name
  term_stats <- rbind(term_stats, c(fet$estimate, fet$p.value) )
  rownames(term_stats)[nrow(term_stats)] <- colnames(smoke_vs_gest_table)[i]
}
# set the column names, and add log2-scaled odds ratio and FDR corrected p-value
colnames(term_stats) <- c("OddsRatio","P.Value")
term_stats$Log2OddsRatio <- log2(term_stats$OddsRatio)
term_stats$Q.Value <- p.adjust(term_stats$P.Value,method="fdr")
term_stats

# -----------------------------------------------------------------------------
# False discovery rate
# -----------------------------------------------------------------------------

# start empty vector
pval <- c()
# generate 10k random data sets from the
# SAME normal distribution
for (i in 1:10000){
        wt <- rnorm(20, mean=10, sd=3)
        ko <- rnorm(20, mean=10, sd=3)
        pval[i] <- t.test(wt, ko)$p.value
}
sum(pval<0.05)
fdr <- p.adjust(pval, method="fdr")
sum(fdr<0.05)

# =============================================================================
# TAKE-HOME EXERCISES
# =============================================================================

# -----------------------------------------------------------------------------
# Read/write data from/to Excel spreadsheets
# -----------------------------------------------------------------------------

install.packages("openxlsx")

# --- Read/write Excel spreadsheets using openxlsx ---------------------------

library(openxlsx)
sheet_1 <- read.xlsx("http://wd.cri.uic.edu/advanced_R/taxa_relative.xlsx", sheet=1)
sheet_1[1:10, 1:5]
sheet_L6 <- read.xlsx("http://wd.cri.uic.edu/advanced_R/taxa_relative.xlsx", sheet="L6")
sheet_L6[1:10, 1:5]
bw_data  <- read.delim("https://wd.cri.uic.edu/R/birth_weight.txt")
# Create a workbook with two empty sheets
wb <- createWorkbook()
addWorksheet(wb, "birth weight")
# Write the birth weight data to Excel
writeData(wb, sheet="birth weight", x=bw_data)

