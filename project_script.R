# Load the used libraries 
library("DESeq2")
library("RColorBrewer")
library("pheatmap")
library("tidyverse")
library("readr")
library("mice")
library("biomaRt")


#####load files#######
## In this project, we will only work on 4 files; two  for the normal plants and two drought stressed plant  
## Here we make a list of the selected data set
working_files <- list.files("~/Desktop/Project/GSE108610_RAW/",full.names = TRUE)
## Then, we will create a single data frame containing FPKM
## The data frame will consist of 8 columns for each sample, and the rows will be gene_id
## Creating an empty data frame
raw_counts <- data.frame(gene_id = character(), stringsAsFactors = FALSE)
for (each_file in working_files){
  # Reading each file into a data frame
  each_file_df <- read_delim(each_file, delim = "\t", escape_double = FALSE, trim_ws = TRUE)
  # Merge the data frames based on gene_id
  raw_counts <- merge(raw_counts, each_file_df[, c("gene_id", "FPKM")], by= "gene_id", all = T)
  # Change column name form FPKM to the sample name only 
  # by two function; basename() to remove the whole path of the file and keep the file name only
  # str_extract() to extract the part we need only which is WT_\\S{2},this means WT_ + the two unspecified letter after.
  new_col_name <- str_extract(basename(each_file), "_(\\w{5})") %>%
    paste0("GSM", .)
  # After that, we rename col by rename(), (!!new_col_name:) means to access the variable not to take it as it is
  raw_counts <- raw_counts %>% rename(!!new_col_name:= FPKM)
}
# Set the row names of the normalized counts data frame as gene_id
rownames(raw_counts) <- raw_counts$gene_id
# Removing the repeated column of gene_id
raw_counts$gene_id <- NULL

## Metadata extraction from the raw_counts dataframe 
# Metadata will contain sample names as rownames and genotype and condition as a columns
metadata_df <- data.frame( genotype = factor(c(rep("WT", 4), rep("VaWRKY14", 4)),
                                             levels = c("VaWRKY14", "WT")),
                           condition = factor(
                             c("normal", "normal", "stressed", "stressed", 
                               "normal", "normal", "stressed", "stressed"), 
                             levels = c("stressed", "normal")),
                           replicates = rep(c(1,2), 4),
                           row.names = colnames(raw_counts))

# Recheck the order of the samples order in both 
all(rownames(metadata_df) == (colnames(raw_counts)))
################################################################################
################################################################################
################################################################################
#####Data Cleaning#######
## Checking missing values and Multiple imputation
sum(is.na(raw_counts))
# Here, we have 17968 missing values. which is a large percent of our data points.
# We can not removing them as this may be lead to lose a significant data. 
# so that we will do Multiple imputation, PMM is a method that imputes missing values by matching observed values to the predicted values from a model.
# Create a mids (multiple imputed datasets) object for multiple imputation
imputed_dataset <- mice(raw_counts, method = "pmm", m = 5, seed=42)
# Complete the imputation process using the mids
imputed_raw_counts <- mice::complete(imputed_dataset)
# Recheck missing values again to make sure everything is ok
sum(is.na(imputed_raw_counts))
# Convert all the values in integers
imputed_raw_counts <- as.data.frame(lapply(imputed_raw_counts, as.integer))
#Reset rawnames again
rownames(imputed_raw_counts) <- rownames(raw_counts)
## Checking duplication in gene_id
duplicates_gene_id <- imputed_raw_counts[duplicated(rownames(imputed_raw_counts)) | duplicated(rownames(imputed_raw_counts), fromlast = TRUE), ]
# There is no duplication 
################################################################################
################################################################################
################################################################################
#####Data exploration#######
# Plot distribution of RNA counts for each sample
# Copying the last data frame to use it in the plotting.  
df_plot <- imputed_raw_counts
# Add the row names (gene_id) as a column in the data frame.
df_plot$gene_id <- rownames(df_plot)
# Start to create new data frame containing the long format (grouped by gene_id, the name will be samples, and values are counts)
melted_data <- pivot_longer(df_plot, colnames(imputed_raw_counts))
melted_data %>% ggplot() + 
  geom_histogram(aes(x=value, group=name), bins = 200) +
  facet_wrap(~name, scale="free_x")+
  ylab("Number of genes") +
  xlab("Expression count")
# Explore mean-variance relationship
# Calculate mean and variance for each gene from the original data
mean_raw_counts <- apply (raw_counts[, 1:8], 1, mean)
variance_raw_counts <- apply(raw_counts[, 1:8], 1, var)
# Creating dataframe for mean and variance of each gene
mean_var_df <- data.frame(mean_raw_counts, variance_raw_counts)
# Dispersion plotting between mean and variance for each gene
ggplot(mean_var_df) + 
  geom_point(aes(x= mean_raw_counts, y = variance_raw_counts))+
  scale_x_log10()+ scale_y_log10()+
  xlab("Mean per each gene")+
  ylab("Variance per each gene")
################################################################################
################################################################################
################################################################################
#####Data preparation for DE analysis (Quality Control)#######
# Creating DESeq2 object (DESeqDataSet = dds) 
dds<- DESeqDataSetFromMatrix(countData = imputed_raw_counts,colData = metadata_df, design = ~ genotype + condition + genotype:condition)
# Counts normalization, FPKM is already used to normalized data to account for the gene length, but it doesn't cover all normalization procedures
# So that, we will normalized data using DESeq 2 to account library size
# In other words, we will remove the impact of differences in sequencing depth or library size
#allowing for a more accurate comparison of gene expression levels between samples. 
#This is crucial for identifying true biological differences in gene expression.
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
# Computing pairwais correlation that will be used in heatmap plotting
vsd_cor <- cor(vsd_matirx)
## Heatmap plotting to assess the similarity in gene expression between different samples
pheatmap(vsd_cor, annotation = dplyr::select(metadata_df, condition, genotype))
## PCA (emphasizes on the variation present in a dataset)
plotPCA(vsd, intgroup = c("condition", "genotype"))
################################################################################
################################################################################
################################################################################
#####DE analysis #######
# When I ran the analysis, I found that the coefficients in resultsnames() used "VaWRKY14" as a reference
# This is because R will choose a reference level for factors based on alphabetical order.
# So that, we will reset it before DE analysis
dds$condition <- relevel(dds$condition, ref = "normal")
dds$genotype <- relevel(dds$genotype, ref = "WT")
## Run DE analysis to figure out the work flow
dds <- DESeq(dds, fitType = "local")
##The workflow as the results of the DE analysis
# Using pre-existing size factors(The function calculates size factors based on the library sizes of your samples)
# The estimated size factors are then stored within the DESeqDataSet object (dds) for downstream use as the following
## Estimate dispersions 
dds <- estimateDispersions(dds,fitType = "local")
# Estimate dispersion for each gene(variability of gene expression across samples)
## Gene-wise dispersion estimates
gene_wise_dispersion <- dispersions(dds)
## Plot Mean-dispersion relationship
plotDispEsts(dds)
## Final dispersion estimates
dds <- nbinomWaldTest(dds)
## Fitting model and testing
res <- results(dds, alpha = 0.05)
################################################################################
################################################################################
################################################################################
#####Exploring DE results #######
# Check the coefficients for the comparison
coefficients <- resultsNames(dds)
# Specifying the coefficients of interest
coefficients_names <- coefficients[2:3]
# Create an empty list to store the results 
all_res <- list()
# Extract results for the selected coefficients using a loop
for (each_name in coefficients_names){
  res <- results(dds, name = each_name)
  all_res[[each_name]] <- res
} 
# Here we have the results of the possible two main effects only which are 
# genotype VaWRKY14 vs WT (normal condition is the reference)= the genotype effect for normal condition
# and condition stressed vs normal (WT genotype is the reference)= the condition effect for WT genotype
# So we need to extract: (main effect + interaction)
# the genotype effect for stressed condition("genotype_VaWRKY14_vs_WT"+"genotypeVaWRKY14.conditionstressed")
# and the condition effect for VaWRKY14 genotype("condition_stressed_vs_normal"+"genotypeVaWRKY14.conditionstressed")
# Extract results for the rest using a loop
for 




# Making more estimates of the fold change which represent the expression of one sample group relative to another.
# DESeq2 contrasts
res <- results(dds, alpha = 0.05)
# From the results, we will plot mean of the normalized counts VS log2 fold change for all gene tested
plotMA(res, ylim = c(-8,8))
# Genes that are significantly DE colored blue
# LFC shrinkage(to generate more likely, lower, log2 fold change estimates, similar to what we did with dispersion.) 
res <- lfcShrink(dds,  contrast = c("condition", "stressed", "normal"), res = res, type = "normal")
# From the results, we will plot mean of the normalized counts VS log2 fold change for all gene tested
plotMA(res, ylim = c(-8,8))
# Get summary of the significant results
summary(res)







# Get annotation
# Here we have not the ensemble gene_IDs, so that we will retrieve them from the biomaRt 
# ensembl_arabidopsis <- useEnsemblGenomes(biomart = "plants_mart", dataset = "athaliana_eg_gene")
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
# Extracting results
res_all <- data.frame(res) %>%
  rownames_to_column( var = "ensembl_gene_id") %>%
  left_join(y = gene_annotations %>% dplyr::select(ensembl_gene_id,description,chromosome_name, uniprot_gn_symbol),
            by = "ensembl_gene_id")
# Extracting and arranging significant genes only
res_sig <- res_all %>%
  subset(padj < 0.05) %>% 
  arrange(padj)



################################################################################
################################################################################
################################################################################
#####Data Visualization #######
## Expression heatmap 
# Subsetting normalized counts of significant genes
sig_normalized_counts <- normalized_raw_counts[ res_sig$ensembl_gene_id,]
# Choose a color palette fromRColorBrewer and save it
heat_colors <- brewer.pal(6, "YlOrRd")
# Heatmap plotting 
pheatmap(sig_normalized_counts, color = heat_colors,
         cluster_rows = T, show_rownames = F,
         annotation = dplyr::select(metadata_df, condition), 
         scale = "row")
## Volcano plotting
# Obtain logical vector regarding whether padj values are less than 0.05 which is return true or false
res_all <- res_all %>% 
  mutate(threshold = padj < 0.05) 
# Volcano plot 
ggplot (res_all) + geom_point(aes(x = log2FoldChange, y = -log10(padj), color = threshold)) +
  xlab("Log2 fold change") + 
  ylab("-Log10 adjusted p-value") + 
  theme(legend.position = "none" , plot.title = element_text(size = rel(1.5), hjust = 0.5), axis.title = element_text(size = rel(1.25)))
## Expression plot 
# Extracting the top 20 significant genes, then Merge the metadata, so that, we can color the plot by sample group
top_sig_20 <- data.frame(sig_normalized_counts[1:20, ]) %>%
  rownames_to_column( var = "ensembl_gene_id") %>%
  gather(key = "samplename", value = "normalized_count", 2:5) %>%
  inner_join(rownames_to_column(metadata_df, var = "samplename"), by = "samplename")
# Expression plotting
ggplot(top_sig_20) + geom_point(aes(x= ensembl_gene_id, y = normalized_count, color = condition)) + 
  scale_y_log10() + 
  xlab("Genes") + 
  ylab("Normalized Counts") + 
  ggtitle("Top 20 Significant DE Genes") + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  theme(plot.title = element_text(hjust = 0.5))




