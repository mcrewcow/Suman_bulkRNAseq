library(ggplot2)
library(viridisLite)
library(viridis)

wd <- "~/Documents/summer_research_2023-main/bulkrna/deseq2"
setwd(wd)


folders <- c("Control-vs-PMC",
             "Control-vs-oxLDL-24h",
             "Control-vs-oxLDL-48h",
             "Control-vs-oxLDL-PMC-24h",
             "Control-vs-oxLDL-PMC-48h")

library(EnhancedVolcano)
for (f in folders){
  csv_file <- read.csv(paste(wd, f, "Differential_expression_analysis_table.csv", sep="/"))
  
  print(EnhancedVolcano(csv_file,
                        lab = csv_file$ID,
                        x = 'log2FoldChange',
                        y = 'pvalue',
                        pCutoff = 0.1,
                        FCcutoff = 1,
                        title = f))
  
}


# Load all results


csv1 <- read.csv(paste(wd, "gsea_res", "Control_vs_PMC.csv", sep="/"))[c("ID","NES","pvalue")]
csv2 <- read.csv(paste(wd, "gsea_res", "Control_vs_oxLDL24h.csv", sep="/"))[c("ID","NES","pvalue")]
csv3 <- read.csv(paste(wd, "gsea_res", "Control_vs_oxLDL48h.csv", sep="/"))[c("ID","NES","pvalue")]
csv4 <- read.csv(paste(wd, "gsea_res", "Control_vs_oxLDL_PMC24h.csv", sep="/"))[c("ID","NES","pvalue")]
csv5 <- read.csv(paste(wd, "gsea_res", "Control_vs_oxLDL_PMC48h.csv", sep="/"))[c("ID","NES","pvalue")]

# Check their lengths

nrow(csv1)
nrow(csv2)
nrow(csv3)
nrow(csv4)
nrow(csv5)

# make copies
csv1_1 <- csv1
csv2_1 <- csv2
csv3_1 <- csv3
csv4_1 <- csv4
csv5_1 <- csv5

# barplot


csv1_1$group <- "PMC"
csv2_1$group <- "oxLDL-24h"
csv3_1$group <- "oxLDL-48h"
csv4_1$group <- "oxLDL-PMC-24h"
csv5_1$group <- "oxLDL-PMC-48h"

df1 <- rbind(csv1_1,csv2_1,csv3_1,csv4_1,csv5_1)

df1$NES <- abs(df1$NES)

#add control
u <- unique(df1$ID) 
control <- data.frame(ID = u, NES = 1, pvalue = 1, group = "Control")
df1 <- rbind(df1, control)

#specify ID
ids <- c("GOBP_ANIMAL_ORGAN_MORPHOGENESIS", "GOBP_CELL_JUNCTION_ORGANIZATION", "HP_ABNORMAL_PALATE_MORPHOLOGY", 
         "HP_DOWNSLANTED_PALPEBRAL_FISSURES", "GOMF_ION_TRANSMEMBRANE_TRANSPORTER_ACTIVITY", "GOBP_QUINONE_METABOLIC_PROCESS") #chosen at random

x <- df1[df1$ID %in% ids, ]


row_order <- c("Control", "oxLDL-24h", "oxLDL-PMC-24h", "oxLDL-48h", "oxLDL-PMC-48h", "PMC")
x$group <- factor(x$group, levels = rev(row_order))


#plot
ggplot(data=x, aes(x=NES, y=ID, fill=group)) +
  geom_bar(stat="identity", position=position_dodge(), width = 0.5) +
  scale_fill_viridis(discrete = TRUE, breaks = row_order) +
  theme_minimal() +
  scale_x_continuous(limits = c(0.9, 2.5), 
                     oob = scales::squish) 




# radar chart

#import packages
library(fmsb)
library(tibble)

#make copy 
csv1_2 <- csv1
csv2_2 <- csv2
csv3_2 <- csv3
csv4_2 <- csv4
csv5_2 <- csv5

#rename columns

csv_list <- list(csv1_2, csv2_2, csv3_2, csv4_2, csv5_2)
new_names <- c("PMC_NES", "oxLDL_24h_NES", "oxLDL_48h_NES", "oxLDL_PMC_24h_NES", "oxLDL_PMC_48h_NES")
new_names1 <- c("PMC_pvalue", "oxLDL_24h_pvalue", "oxLDL_48h_pvalue", "oxLDL_PMC_24h_pvalue", "oxLDL_PMC_48h_pvalue")

for (i in 1:length(csv_list)) {
  names(csv_list[[i]])[names(csv_list[[i]]) == 'NES'] <- new_names[i]
}

for (i in 1:length(csv_list)) {
  names(csv_list[[i]])[names(csv_list[[i]]) == 'pvalue'] <- new_names1[i]
}


#merge into one data frame

df2 <- Reduce(function(...) merge(..., by='ID', all=TRUE), csv_list)

#replace NA with 0
df2[is.na(df2)] <- 0
#take absolute value
df2[,2:6] <- abs(df2[,2:6])

write.csv(df2,file=paste("full_gsea_results.csv"), row.names=FALSE)



#specify id
ids <- c("GOBP_ANIMAL_ORGAN_MORPHOGENESIS", "GOBP_CELL_JUNCTION_ORGANIZATION", "HP_ABNORMAL_PALATE_MORPHOLOGY", 
         "HP_DOWNSLANTED_PALPEBRAL_FISSURES", "GOMF_ION_TRANSMEMBRANE_TRANSPORTER_ACTIVITY", "GOBP_QUINONE_METABOLIC_PROCESS") #chosen at random



x <- df2[df2$ID %in% ids, ]

x <- add_column(x, Control = 1, .after = 1)
rownames(x) <- x$ID


x2 <- data.frame(t(x[-1]))



#color vectors
colors_border=c( rgb(0.2,0.5,0.5,0.9))
colors_in=c( rgb(0.2,0.5,0.5,0.2))




par(mar=c(1, 2, 2, 1)) #decrease default margin 
par(mfrow=c(2,3)) #layout

#loop over rows to draw them, add 2.5 as max and 0 as min for each var

lapply(1:6, function(i) { 
  radarchart(rbind(rep(2.5,6), rep(0,6), x2[i,]), 
             axistype = 1,
             pcol=colors_border , pfcol=colors_in ,
             plwd=2 , plty=1,
             #custom the grid
             cglcol="grey", cglty=1, axislabcol="grey", caxislabels=seq(0,2.5,0.5), cglwd=0.8,
             #custom labels
             vlcex=0.8,
             title = rownames(x2)[i],
             seg = 5
  )
}) 




# save full results
#write.csv(df2,file=paste("gsearesults.csv"), row.names=FALSE)




w <-unique(mm_hallmark_sets[, c("gs_name", "gs_description")])

# Lysosomal pathways
lysosomal_subset <- w[grep("LYSOSOM", w$gs_description, ignore.case = TRUE),] 
#View(lysosomal_subset)


# Autophagy pathway
autophagy_subset <- w[grep("AUTOPHAG", w$gs_description, ignore.case = TRUE),] 
#View(autophagy_subset)

# Lipid metabolism
lipid_subset <- w[grep("_LIPID_METABOL", w$gs_name),] 
#View(lipid_subset)

# Antioxidant pathway
antioxidant_subset <- w[grep("ANTIOXID", w$gs_description, ignore.case = TRUE),] 
#View(antioxidant_subset)

# Mitochondrial pathways
mitochondrial_subset <- w[grep("MITOCHOND", w$gs_name, ignore.case = TRUE),] 
#View(mitochondrial_subset)

# Cell death pathways
apoptosis_subset <- w[grep("apoptos", w$gs_description, ignore.case = TRUE),] 
#View(apoptosis_subset)
necrosis_subset <- w[grep("necro", w$gs_description, ignore.case = TRUE),] 
#View(necrosis_subset)

subsets <- rbind(lysosomal_subset,autophagy_subset,lipid_subset,antioxidant_subset,mitochondrial_subset,apoptosis_subset,necrosis_subset)

# !!!Check for unwanted pathways
View(subsets)
#write.csv(subsets,file=paste("subsets.csv"), row.names=FALSE)



inter <- intersect(df2$ID, subsets$gs_name)

x <- df1[df1$ID %in% inter, ]
#x <- rbind(x, df1[1,])       #pmc check


#plot
ggplot(data=x, aes(x=NES, y=ID, fill=group)) +
  geom_bar(stat="identity", position=position_dodge(), width = 0.5) +
  scale_fill_viridis(discrete = TRUE) +
  theme_minimal() +
  scale_x_continuous(limits = c(1, 2.5), 
                     oob = scales::squish)
