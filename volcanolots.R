library(EnhancedVolcano)
wd <- "~/Documents/summer_research_2023-main/bulkrna/deseq2"
setwd(wd)


folders <- c("Control-vs-PMC",
             "Control-vs-oxLDL-24h",
             "Control-vs-oxLDL-48h",
             "Control-vs-oxLDL-PMC-24h",
             "Control-vs-oxLDL-PMC-48h")


for (f in folders){
  csv_file <- read.csv(paste(wd, f, "Differential_expression_analysis_table.csv", sep="/"))
  
  pdf(file = paste(wd, "/images/volcanoplots/", f, ".pdf", sep=""), width = 10)
  
  
  plot <- EnhancedVolcano(csv_file,
                          lab = csv_file$Gene.name,
                          x = 'log2FoldChange',
                          y = 'pvalue',
                          pCutoff = 0.05,
                          FCcutoff = 1,
                          title = f,
                          legendPosition = "right",
                          legendLabSize = 12,
                          labSize = 4,
                          drawConnectors = TRUE,
                          max.overlaps = 45)
  print(plot)
  dev.off()
}
