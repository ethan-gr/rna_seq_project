## ============= LOAD PACKAGES =============

library("recount3")
library("SummarizedExperiment")
library("ggplot2")
library("limma")
library("edgeR")
library("pheatmap")


## ============ DATA DOWNLOAD ==============

# get project information through recount3
mouse_projects <- available_projects(organism = "mouse")

# create RangedSummarizedExperiment with the information at genes level
rse_gene_SRP181125 <- create_rse(
  subset(
    mouse_projects,
    project == "SRP181125" & project_type == "data_sources"
  )
)

# set counts per nucleotide to counts per lecture
assay(rse_gene_SRP181125, "counts") <- compute_read_counts(rse_gene_SRP181125)


## ======== EXPLORATION & CLEANING =========

# more easy use information
rse_gene_SRP181125 <- expand_sra_attributes(rse_gene_SRP181125)

# overview sra_attributes
summary(as.data.frame(colData(rse_gene_SRP181125)[, grepl("^sra_attribute", colnames(colData(rse_gene_SRP181125)))]))

# reproductive state to factor
rse_gene_SRP181125$sra_attribute.reproductive_state <- as.factor(rse_gene_SRP181125$sra_attribute.reproductive_state)
# tissue to factor
rse_gene_SRP181125$sra_attribute.tissue <- as.factor(rse_gene_SRP181125$sra_attribute.tissue)
# source name to factor
rse_gene_SRP181125$sra_attribute.source_name <- as.factor(rse_gene_SRP181125$sra_attribute.source_name)
# animal identification to factor
rse_gene_SRP181125$sra_attribute.animal_identification <- as.factor(rse_gene_SRP181125$sra_attribute.animal_identification)
# strain to factor
rse_gene_SRP181125$sra_attribute.strain <- as.factor(rse_gene_SRP181125$sra_attribute.strain)

# overview modifyed sra_attributes
summary(as.data.frame(colData(rse_gene_SRP181125)[, grepl("^sra_attribute", colnames(colData(rse_gene_SRP181125)))]))
# it seems to be equilibrium on the data

# quality exploration based on assiganation ratio
rse_gene_SRP181125$assigned_gene_prop <- rse_gene_SRP181125$recount_qc.gene_fc_count_all.assigned / rse_gene_SRP181125$recount_qc.gene_fc_count_all.total
summary(rse_gene_SRP181125$assigned_gene_prop)

# visualization
with(colData(rse_gene_SRP181125), plot(assigned_gene_prop, sra_attribute.source_name))
with(colData(rse_gene_SRP181125), tapply(assigned_gene_prop, sra_attribute.reproductive_state, summary))
with(colData(rse_gene_SRP181125), tapply(assigned_gene_prop, sra_attribute.tissue, summary))
with(colData(rse_gene_SRP181125), tapply(assigned_gene_prop, sra_attribute.source_name, summary))
# it seems to have good quality


## ============= FILTER DATA ===============

# backup original size
rse_gene_SRP181125_unfiltered <- rse_gene_SRP181125

# check histogram
hist(rse_gene_SRP181125$assigned_gene_prop)

# define cutoff quality
cutoff <- 0.6
table(rse_gene_SRP181125$assigned_gene_prop < cutoff)
rse_gene_SRP181125 <- rse_gene_SRP181125[, rse_gene_SRP181125$assigned_gene_prop > cutoff]
# no data were removed, we could refine more the data
# aiming to not devalance the cuantities on samples

# we calculate gene expression level means
gene_means <- rowMeans(assay(rse_gene_SRP181125, "counts"))
summary(gene_means)

# gene means
hist(gene_means)
hist(log(gene_means))
hist(log(gene_means[gene_means > 0.1]))

# remove low expressed genes under 0.1 expression
rse_gene_SRP181125 <- rse_gene_SRP181125[gene_means > 0.1, ]

# original and final dimentions
dim(rse_gene_SRP181125_unfiltered)
dim(rse_gene_SRP181125)

# recovery ratio
round(nrow(rse_gene_SRP181125) / nrow(rse_gene_SRP181125_unfiltered) * 100, 2)
# in order to find relevant things we had to remove 47.61 of the data


## ============ NORMALIZATION ==============

dge <- DGEList(
  counts = assay(rse_gene_SRP181125, "counts"),
  genes = rowData(rse_gene_SRP181125)
)
dge <- calcNormFactors(dge)


## ======== DIFFERENTIAL EXPRESSION ========

# define stadistic model
ggplot(as.data.frame(colData(rse_gene_SRP181125)), aes(y = assigned_gene_prop, x = sra_attribute.reproductive_state)) +
  geom_boxplot() +
  theme_bw(base_size = 20) +
  ylab("Assigned Gene Prop") +
  xlab("Reproductive State")

ggplot(as.data.frame(colData(rse_gene_SRP181125)), aes(y = assigned_gene_prop, x = sra_attribute.tissue)) +
  geom_boxplot() +
  theme_bw(base_size = 20) +
  ylab("Assigned Gene Prop") +
  xlab("Tissue")

# define stadistic model 2
mod <- model.matrix(
  ~ sra_attribute.reproductive_state + sra_attribute.tissue + assigned_gene_prop,
  data = colData(rse_gene_SRP181125)
)
colnames(mod)

# limma for differential expression analysis
vGene <- voom(dge, mod, plot = TRUE)


eb_results <- eBayes(lmFit(vGene))
de_results <- topTable(
  eb_results,
  coef = 2,
  number = nrow(rse_gene_SRP181125),
  sort.by = "none"
)

# explore
dim(de_results)
head(de_results)

# differential expressed genes between reproductive states with FDR < 5%
table(de_results$adj.P.Val < 0.05)

# visualize stadistic results
plotMA(eb_results, coef = 2)
volcanoplot(eb_results, coef = 2, highlight = 3, names = de_results$gene_name)

de_results[de_results$gene_name %in% c("Acer2", "Sgk1", "Kpna2-ps"), ]


## ========== VISUALIZE GENES DE ===========

# make data more easy to handle
df <- as.data.frame(colData(rse_gene_SRP181125)[, c("sra_attribute.reproductive_state", "sra_attribute.tissue")])
colnames(df) <- c("ReproductiveState", "Tissue")

# extract gene's values of interest (top 50)
exprs_heatmap <- vGene$E[rank(de_results$adj.P.Val) <= 30, ]

# heatmap
pheatmap(
  exprs_heatmap,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = FALSE,
  show_colnames = FALSE,
  annotation_col = df
)
