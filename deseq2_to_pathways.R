if (!("clusterProfiler" %in% installed.packages())) {
  # Install this package if it isn't installed yet
  BiocManager::install("clusterProfiler", update = FALSE)
}

if (!("msigdbr" %in% installed.packages())) {
  # Install this package if it isn't installed yet
  BiocManager::install("msigdbr", update = FALSE)
}

if (!("org.Hs.eg.db" %in% installed.packages())) {
  # Install this package if it isn't installed yet
  BiocManager::install("org.Hs.eg.db", update = FALSE)
}


# Attach the library
library(clusterProfiler)

# Package that contains MSigDB gene sets in tidy format
library(msigdbr)

# Human annotation package we'll use for gene identifier conversion
library(org.Hs.eg.db)

# We will need this so we can use the pipe: %>%
library(magrittr)



wd <- "~/Documents/summer_research_2023-main/bulkrna/deseq2"
setwd(wd)


folders <- c("Control-vs-PMC",
             "Control-vs-oxLDL-24h",
             "Control-vs-oxLDL-48h",
             "Control-vs-oxLDL-PMC-24h",
             "Control-vs-oxLDL-PMC-48h")

#new column for control and test group average
for (f in folders){
  csv_file <- read.csv(paste(wd, f, "Differential_expression_analysis_table.csv", sep="/"))
  csv_file$ControlMean <- rowMeans(csv_file[,5:7])
  csv_file$TestGroupMean <- rowMeans(csv_file[,8:10])
  csv_file <- csv_file[,-5:-10]
  write.csv(csv_file,file=paste(wd, f, "Differential_expression_analysis_table_new.csv", sep="/"), row.names=FALSE)
  
  
  
  
  data <- csv_file
  
  
  msigdbr_species()
  
  mm_hallmark_sets <- msigdbr(
    species = "Homo sapiens", # Replace with species name relevant to your data
    category = "C5"
  )
  
  # First let's create a mapped data frame we can join to the differential
  # expression stats
  dge_mapped_df <- data.frame(
    gene_symbol = mapIds(
      # Replace with annotation package for the organism relevant to your data
      org.Hs.eg.db,
      keys = data$ID,
      # Replace with the type of gene identifiers in your data
      keytype = "ENSEMBL",
      # Replace with the type of gene identifiers you would like to map to
      column = "SYMBOL",
      # This will keep only the first mapped value for each Ensembl ID
      multiVals = "first"
    )
  ) %>%
    # If an Ensembl gene identifier doesn't map to a gene symbol, drop that
    # from the data frame
    dplyr::filter(!is.na(gene_symbol)) %>%
    # Make an `Ensembl` column to store the rownames
    tibble::rownames_to_column("Ensembl") %>%
    # Now let's join the rest of the expression data
    dplyr::inner_join(data, by = c("Ensembl" = "ID"))
  
  
  dup_gene_symbols <- dge_mapped_df %>%
    dplyr::filter(duplicated(gene_symbol)) %>%
    dplyr::pull(gene_symbol)
  
  dge_mapped_df %>%
    dplyr::filter(gene_symbol %in% dup_gene_symbols) %>%
    dplyr::arrange(gene_symbol)
  
  filtered_dge_mapped_df <- dge_mapped_df %>%
    # Sort so that the highest absolute values of the log2 fold change are at the
    # top
    dplyr::arrange(dplyr::desc(abs(log2FoldChange))) %>%
    # Filter out the duplicated rows using `dplyr::distinct()`
    dplyr::distinct(gene_symbol, .keep_all = TRUE)
  
  print(paste("duplicated =", any(duplicated(filtered_dge_mapped_df$Gene.name))))
  
  
  
  # Let's create a named vector ranked based on the log2 fold change values
  lfc_vector <- filtered_dge_mapped_df$log2FoldChange
  names(lfc_vector) <- filtered_dge_mapped_df$gene_symbol
  
  # We need to sort the log2 fold change values in descending order here
  lfc_vector <- sort(lfc_vector, decreasing = TRUE)
  
  # Set the seed so our results are reproducible:
  set.seed(2023)
  
  gsea_results <- GSEA(
    geneList = lfc_vector, # Ordered ranked gene list
    minGSSize = 25, # Minimum gene set size
    maxGSSize = 2000, # Maximum gene set set
    pvalueCutoff = 0.01, # p-value cutoff
    eps = 0, # Boundary for calculating the p value
    seed = TRUE, # Set seed to make results reproducible
    pAdjustMethod = "BH", # Benjamini-Hochberg correction
    TERM2GENE = dplyr::select(
      mm_hallmark_sets,
      gs_name,
      gene_symbol
    )
  )
  
  # We can access the results from our `gsea_results` object using `@result`
  gsea_result_df <- data.frame(gsea_results@result)
  # Then write to file
  write.csv(gsea_result_df,file=paste(wd, f, "GSEA_result.csv", sep="/"), row.names=FALSE)
  
}


