setwd("~/Local/projects/VaWRKY14 and Drought Tolerance in Arabidopsis")
source("functions.R")

##############################################################
#####################load_working_files#######################
# Get raw_counts from the files by get_raw_counts function 
raw_counts <- get_raw_counts("~/Local/projects/VaWRKY14 and Drought Tolerance in Arabidopsis/GSE108610_RAW/")
# Metadata extraction from the raw_counts data frame
# Metadata will contain sample names as row names, and genotype and replicates as columns.
# Create metadata dataframe
metadata_df <- data.frame( genotype = factor(c(rep("WT", 2), rep("VaWRKY14", 2)),
                                             levels = c("VaWRKY14", "WT")),
                           replicates = rep(c(1,2), 2),
                           row.names = colnames(raw_counts))
# Recheck the order of the samples order in both 
all(rownames(metadata_df) == (colnames(raw_counts)))

##############################################################
############ Data Cleaning and Multiple Imputation ###########
imputed_raw_counts <- performImputation(raw_counts)

##############################################################
###################### Data Exploration ######################
# Plot distribution of RNA counts for each sample
# Create a copy of the imputed_raw_counts dataframe for plotting
df_plot <- imputed_raw_counts
# Add gene_id as a column in the data frame
df_plot$gene_id <- rownames(df_plot)
# Reshape the data to long format (grouped by gene_id, samples as names, and counts as values)
melted_data <- pivot_longer(df_plot, colnames(imputed_raw_counts))
# Plot the distribution of RNA counts for each sample
melted_data %>%
  ggplot() +
  geom_histogram(aes(x = value, group = name), bins = 200) +
  facet_wrap(~name, scale = "free_x") +
  ylab("Number of genes") +
  xlab("Expression count")

# Explore mean-variance relationship
# Calculate mean and variance for each gene from the original data
mean_var_df <- calculateMeanAndVariance(raw_counts)
# Plot the dispersion between mean and variance for each gene
mean_var_df %>% ggplot() +
  geom_point(aes(x = mean_vals, y = var_vals)) +
  scale_x_log10() + scale_y_log10() +
  xlab("Mean per gene") +
  ylab("Variance per gene")

##############################################################
##### Data preparation for DE analysis (Quality Control) #####
# Creating DESeq2 object (DESeqDataSet = dds) 
dds<- DESeqDataSetFromMatrix(countData = imputed_raw_counts,
                             colData = metadata_df, design = ~ genotype)
# Counts normalization, FPKM is already used to normalized data to account for the gene length, but it doesn't cover all normalization procedures
# So that, we will normalized data using DESeq2 to account library size
# In other words, we will remove the impact of differences in sequencing depth or library size
# allowing for a more accurate comparison of gene expression levels between samples. 
# This is crucial for identifying true biological differences in gene expression.
# Normalized counts calculation 
# Estimate size factor(size factors are used for normalization)
dds <- estimateSizeFactors(dds)
size_factors <- sizeFactors(dds)
# Normalized counts extraction
normalized_raw_counts <- counts(dds, normalized = T)

# Unsupervised clustering analysis by variance stabilizing transformation (log transformation)
# Create object containing variance stabilizing data (vsd) 
vsd <- vst(dds, blind = T, fitType = "local") 
# Extract vsd matrix from the object (contain log transformation of the normalized data)
vsd_matirx <- assay(vsd)
# Computing pairwaise correlation that will be used in heatmap plotting
vsd_cor <- cor(vsd_matirx)
## Heatmap plotting to assess the similarity in gene expression between different samples
pheatmap(vsd_cor, annotation = dplyr::select(metadata_df, genotype))

# PCA (emphasizes on the variation present in a dataset)
plotPCA(vsd, intgroup = "genotype")

##############################################################
########################  DE analysis ########################
# When I ran the analysis, I found that the coefficients in resultsnames() used "VaWRKY14" as a reference
# This is because R will choose a reference level for factors based on alphabetical order.
# So that, we will reset it before DE analysis
dds$genotype <- relevel(dds$genotype, ref = "WT")
## Run DE analysis to figure out the work flow
dds <- DESeq(dds, fitType = "local") 
## Plot Mean-dispersion relationship
plotDispEsts(dds)

##############################################################
################# Explore DE analysis ########################
# Extract the results
res <- results(dds, contrast = c("genotype", "VaWRKY14", "WT") , cooksCutoff=FALSE)
res_shrinked <- lfcShrink(dds, res = res , type = "ashr")

# From the results, we will plot mean of the normalized counts VS log2 fold change for all gene tested for each comparison 
# Genes that are significantly DE colored blue
plotMA(res, ylim = c(-8,8), main = "res")
abline(h=c(-1,1), col="purple")
# Plotting after shrinkage
plotMA(res_shrinked, ylim = c(-4,4), main = "res_shrinked")
abline(h=c(-1,1), col= "purple" )

# Get summary of the results
summary(res)

# Get annotation
# Here we have not the ensemble gene_IDs, so that we will retrieve them from the biomaRt 
ensembl_arabidopsis <- useEnsemblGenomes(biomart = "plants_mart", dataset = "athaliana_eg_gene")
# But how we get the right name of the two arguments; biomart and dataset.
# We first list all the genomes in the ensembl to know the biomart name which is "plants_mart"
listEnsemblGenomes()
# Then, we need to figure out the database name
# We extract all the ensemble genomes in the "plants_mart"
ensemble_plants <- useEnsemblGenomes(biomart = "plants_mart")
# After that, we searching about the Arabidopsis
searchDatasets(ensemble_plants, pattern = "Arabidopsis")
# After getting the right names, we extract the Arabidopsis ensembl dataset
ensembl_arabidopsis <- useEnsemblGenomes(biomart = "plants_mart", dataset = "athaliana_eg_gene")
# Retrieve gene annotation
gene_annotations <- getBM(attributes = c("ensembl_gene_id","description", "chromosome_name", "uniprot_gn_symbol" ),
                          filters = "ensembl_gene_id",
                          values = c(rownames(raw_counts)),
                          mart = ensembl_arabidopsis)

# Extracting results with its annotation
DEGs <- data.frame(res_shrinked) %>%
  rownames_to_column( var = "ensembl_gene_id") %>%
  left_join(y = gene_annotations %>% dplyr::select(ensembl_gene_id,description,chromosome_name, uniprot_gn_symbol),
            by = "ensembl_gene_id")

# Extracting and arranging significant genes only
sig_DEGs <- DEGs %>%
  subset(padj < 0.01 & abs(log2FoldChange) >= 1) %>%
  arrange(padj)

# Extracting upregurated and downregulated genes
up_DEGs <- sig_DEGs %>% 
  subset(log2FoldChange > 0, 
         select = c("ensembl_gene_id", "log2FoldChange","description", "uniprot_gn_symbol"))

down_DEGs <- sig_DEGs %>% 
  subset(log2FoldChange < 0,
         select = c("ensembl_gene_id", "log2FoldChange","description", "uniprot_gn_symbol")) 

##############################################################
################ Results Visualization #######################
## Expression heatmap 
# Subsetting normalized counts of significant genes
sig_normalized_counts <- normalized_raw_counts[sig_DEGs$ensembl_gene_id,]
# Choose a color palette fromRColorBrewer and save it
heat_colors <- brewer.pal(6, "YlOrRd")
# Heatmap plotting 
pheatmap(sig_normalized_counts, color = heat_colors,
         cluster_rows = T, show_rownames = T,
         annotation = dplyr::select(metadata_df, c(genotype)), 
         scale = "row",
         main = "sig_DEGs")

## Volcano plotting
# Obtain logical vector regarding whether padj values are less than 0.01 which is return true or false
# Creating for loop to extract results
vol_plot_res <- DEGs %>% 
  mutate(threshold = padj < 0.01) 
#Creating for loop to plot volcano for each comparison 
ggplot(vol_plot_res, aes(x = log2FoldChange, y = -log10(padj), 
                         color = threshold)) +
  geom_point() +
  scale_color_manual(values = c("blue", "red"), breaks = c(FALSE, TRUE),
                     labels = c("Not Significant", "Significant")) +
  xlab("Log2 fold change") + 
  ylab("-Log10 adjusted p-value")+
  ggtitle("Volcano Plot for Differential Expression") +
  theme_minimal()

## Expression plot 
# Extracting the top 20 significant genes, then Merge the metadata, so that, we can color the plot by sample group
top_DEGs <- data.frame(sig_normalized_counts[1:20, ]) %>%
  rownames_to_column(var = "ensembl_gene_id") %>%
  gather(key = "samplename", value = "normalized_count", 2:5) %>%
  inner_join(rownames_to_column(metadata_df, var = "samplename"), by = "samplename")

# Creating Expression plotting 
top_DEGs %>% ggplot() + 
  geom_point(aes(x= ensembl_gene_id, y = normalized_count, color = genotype)) + 
  scale_y_log10() + 
  xlab("Genes") + 
  ylab("Normalized Counts") + 
  ggtitle("top 20 significant genes") + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  theme(plot.title = element_text(hjust = 0.5))


