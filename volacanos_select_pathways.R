library(msigdbr)
library(EnhancedVolcano)


wd <- "~/Documents/summer_research_2023-main/bulkrna"
setwd(wd)

folders <- c("Control-vs-PMC",
             "Control-vs-oxLDL-24h",
             "Control-vs-oxLDL-48h",
             "Control-vs-oxLDL-PMC-24h",
             "Control-vs-oxLDL-PMC-48h",
             "oxLDL-24h-vs-oxLDL-PMC-24h",
             "oxLDL-48h-vs-oxLDL-PMC-48h")


hs_hallmark_sets <- msigdbr(
  species = "Homo sapiens", # Replace with species name relevant to your data
  category = "C5"
)

#replace with relevant pathways
pathways <- c("GOBP_APOPTOTIC_SIGNALING_PATHWAY", "GOBP_LIPID_METABOLIC_PROCESS", "GOBP_INFLAMMATORY_RESPONSE") 


genes <- hs_hallmark_sets[hs_hallmark_sets$gs_name %in% pathways,][c("gs_name", "gene_symbol")]



for (f in folders){
  csv_file <- read.csv(paste(wd, "deseq2", f, "Differential_expression_analysis_table.csv", sep="/"))

  csv_file <- csv_file[csv_file$Gene.name %in% genes[genes$gs_name %in% pathways,]$gene_symbol,]
  
  pdf(file = paste(wd, "/images/volcano-select-pathways/", f, ".pdf", sep=""), width = 16, height = 16)
  
  
  
  
  plot <- EnhancedVolcano(csv_file,
                          lab = csv_file$Gene.name,
                          selectLab = genes$gene_symbol,
                          x = 'log2FoldChange',
                          y = 'pvalue',
                          pCutoff = 0.05,
                          FCcutoff = 1,
                          title = f,
                          legendPosition = "right",
                          legendLabSize = 12,
                          labSize = 4,
                          drawConnectors = TRUE,
                          max.overlaps = 50)
  
  print(plot)
  dev.off()
}
