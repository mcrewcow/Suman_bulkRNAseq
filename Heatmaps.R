library(ComplexHeatmap)
library(msigdbr)



wd <- "~/Documents/summer_research_2023-main/bulkrna/deseq2"
setwd(wd)



folders <- c("Control-vs-PMC",
             "Control-vs-oxLDL-24h",
             "Control-vs-oxLDL-48h",
             "Control-vs-oxLDL-PMC-24h",
             "Control-vs-oxLDL-PMC-48h",
             "oxLDL-24h-vs-oxLDL-PMC-24h",
             "oxLDL-48h-vs-oxLDL-PMC-48h")




for(i in 1:length(folders)) {                              # Head of for-loop
  assign(paste0("csv", i),                                   # Read and store data frames
         read.csv(paste(wd, folders[i], "Differential_expression_analysis_table.csv", sep="/"))[c("ID","log2FoldChange")])
  csv <- get(paste0("csv", i))
  colnames(csv)[2] <- folders[i]
  assign(paste0("csv", i), csv)
}

csv_list <- list(csv1, csv2, csv3, csv4, csv5, csv6, csv7)

genes_full <- Reduce(function(...) merge(..., by='ID', all=TRUE), csv_list)

genes_full[is.na(genes_full)] = 0


annotations <- read.csv(paste(wd, folders[1], "Differential_expression_analysis_table.csv", sep="/"))[c("ID","Gene.name")]

genes_full <- merge(genes_full, annotations, by="ID")


nrow(genes_full[duplicated(genes_full$Gene.name),])



genes <- c("TSPAN6",
           "DPM1",
           "SCYL3")

m <- genes_full[genes_full$Gene.name %in% genes,]
rownames(m) <- m$Gene.name


Heatmap(as.matrix(m[,2:8]))


hs_hallmark_sets <- msigdbr(
  species = "Homo sapiens", # Replace with species name relevant to your data
  category = "C5"
)

#replace with relevant pathways
pathways <- c("GOBP_APOPTOTIC_SIGNALING_PATHWAY", "GOBP_LIPID_METABOLIC_PROCESS", "GOBP_INFLAMMATORY_RESPONSE") 

genes <- hs_hallmark_sets[hs_hallmark_sets$gs_name %in% pathways,][c("gs_name", "gene_symbol")]


m <- genes_full[genes_full$Gene.name %in% genes[genes$gs_name %in% pathways,]$gene_symbol,]
m$mean <- apply(m[,2:8],1,mean)
m$mean <- abs(m$mean)
m <-m[order(m$mean, decreasing = TRUE),]
rownames(m) <- m$Gene.name

Heatmap(as.matrix(m[1:50,2:8]))

