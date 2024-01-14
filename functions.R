# Load the used libraries 
library("DESeq2")
library("RColorBrewer")
library("pheatmap")
library("tidyverse")
library("readr")
library("mice")
library("biomaRt")
library("ggVennDiagram")
library("clusterProfiler") #Gene Ontology (GO) Analysis
library("ReactomePA")      #Pathway Analysis
library("igraph")          #Network Analysis

#################################################
################load_working_files###############
# Function to rad and merge gene expression data from multiple files

get_raw_counts <- function(dir_path){
  # List all the files in the directory
  working_files <- list.files(dir_path, full.names = TRUE)
  # Create an empty data frame to store gene expression data
  raw_counts <- data.frame(gene_id = character(), stringsAsFactors = FALSE)
  # Loop through each file in the list
  for(each_file in working_files){
    # Read file into a data frame
    tryCatch({
      each_file_df <- read.delim(each_file, sep = "\t",header = TRUE)
      # Merge the data frames based on gene_id
      raw_counts <- merge(raw_counts, each_file_df[, c("gene_id", "FPKM")],
                          by = "gene_id",
                          all = TRUE) 
      # Extract sample name from the file name
      sample_name <- str_extract(basename(each_file), "_(\\w{5})") %>%
        paste0("GSM", .)
      # Rename the FPKM column with the extracted sample name
      colnames(raw_counts)[colnames(raw_counts) == "FPKM"] <- sample_name
    }, error = function(e){
      cat("Error reading file:", each_file, "\n")
      cat("Error message", conditionMessage(e), "\n\n")
    })
  }
  # Set the row names of the data frame as gene_id
  rownames(raw_counts) <- raw_counts$gene_id
  # Removing the repeated column of gene_id 
  raw_counts$gene_id <- NULL
   return(raw_counts)                        
}


#################################################
##### Data Cleaning and Multiple Imputation #####
# Function to clean data and perform multiple mutation 

performImputation <- function(raw_counts, method = "pmm", m = 5, seed = 42) {
  # Check for missing values
  missing_values <- sum(is.na(raw_counts))
  cat("Number of missing values:", missing_values, "\n")
  # Perform Multiple Imputation using PMM method
  imputed_dataset <- mice::mice(raw_counts, method = method, m = m, seed = seed)
  # Complete the imputation process
  imputed_raw_counts <- mice::complete(imputed_dataset)
  # Recheck missing values after imputation
  missing_values_after_imputation <- sum(is.na(imputed_raw_counts))
  cat("Number of missing values after imputation:", missing_values_after_imputation, "\n")
  # Convert imputed values to integers
  imputed_raw_counts <- as.data.frame(lapply(imputed_raw_counts, as.integer))
  # Reset row names
  rownames(imputed_raw_counts) <- rownames(raw_counts)
  # Check for duplications in gene_id
  duplicates_gene_id <- imputed_raw_counts[duplicated(rownames(imputed_raw_counts)) | 
                                             duplicated(rownames(imputed_raw_counts), 
                                                        fromlast = TRUE), ]
  # Handle duplications if they exist
  if (nrow(duplicates_gene_id) > 0) {
    warning("Duplications found in gene_id. Duplicates are being removed.")
  } else {
    # Print a message indicating no duplications
    cat("No duplications in gene_id.\n")
  }
  return(imputed_raw_counts)
}

##############################################################
###################### Data Exploration ######################

# Function to calculate mean and variance
calculateMeanAndVariance <- function(data) {
  mean_vals <- apply(data, 1, mean)
  var_vals <- apply(data, 1, var)
  return(data.frame(mean_vals, var_vals))
}

