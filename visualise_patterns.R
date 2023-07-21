library("ggsignif")  
library("stringr")   
library("ggplot2")
library("dplyr")
library("spradarchart")



wd <- "~/Documents/summer_research_2023-main/bulkrna"
setwd(wd)



folders <- c("Control-vs-PMC",
             "Control-vs-oxLDL-24h",
             "Control-vs-oxLDL-48h",
             "Control-vs-oxLDL-PMC-24h",
             "Control-vs-oxLDL-PMC-48h")


group_names <- c("PMC", "oxLDL-24h", "oxLDL-48h", "oxLDL-PMC-24h", "oxLDL-PMC-48h", "Control")

# Load all gsea results

for(i in 1:length(folders)) { # Head of for-loop
  assign(paste0(gsub("-", "_", group_names[i]), "_csv"),                                   # Read and store data frames
         read.csv(paste("deseq2", folders[i], "GSEA_result.csv", sep="/"))[c("ID","NES","pvalue")])
}

#common dataframe
NES_df_consecutif <- data.frame()

for (i in 1:length(folders)){
  df <- paste0(gsub("-", "_", group_names[i]), "_csv")
  
  # Check their lengths
  print(paste("nrow", df, "=", nrow(get(df))))
  # make copies
  df_with_group <- get(df)
  #group column 
  df_with_group$group <- group_names[i]

  #add into one common dataframe
  NES_df_consecutif <- rbind(NES_df_consecutif, df_with_group)
  
}


#add control group
u <- unique(NES_df_consecutif$ID) 
control <- data.frame(ID = u, NES = 1, pvalue = 1, group = "Control")
NES_df_consecutif <- rbind(control,NES_df_consecutif)




# load pattern pathways 
pattern_folder <- "pattern_folder"
# same as .csv file name 
pattern_names <- c("Oxygen_Metabolism", "Metabolism", "Inflammation", "Mitochondria", "Apoptosis", "Lipid_Metabolism", "Retinol_Metabolism", "Autophagy")


for(name in pattern_names) {                              # Head of for-loop
  assign(name,                                   # Read and store data frames
         read.csv(paste(wd, pattern_folder, paste(name, ".csv", sep=""), sep="/")))["ID"]
}

#list of vectors, one for each pattern with their respective pathways
pathwayList <- lapply(pattern_names, get)
pathwayList <- lapply(pathwayList, function(df) df[, 1])



#colors for each condition (used in all plots)
colors <-  c("Control" = "#AAAAAA", "oxLDL-24h" = "#7090FF","oxLDL-48h" = "#0040AA",
               "PMC" = "#ED7014","oxLDL-PMC-24h" = "#E3242B","oxLDL-PMC-48h" = "#AA2080")


row_order <- c("Control", "oxLDL-24h", "oxLDL-PMC-24h", "oxLDL-48h", "oxLDL-PMC-48h", "PMC")


#create directories in wd given a path
createDirectories <- function(path) {
  directories <- unlist(strsplit(path, "/"))
  currentPath <- "."
  
  for (dir in directories) {

      currentPath <- file.path(currentPath, dir)
      
    if (!dir.exists(currentPath)) {
      dir.create(currentPath)
    }
  }
}



#expression over all pathways for each pattern as barplot
#saved as pdf in path
pathway_expression <- function(conditions, path, NES_df_consecutif, pathwayList, pattern_names){
  
  createDirectories(paste(path, "pathway-expression", sep="/"))
  
  for (i in 1:length(pattern_names)) {
    
    # vector of pathways for the pattern
    pathways <- as.vector(pathwayList[i][[1]])
    
    #extract NES value from common dataframe
    x <- NES_df_consecutif[NES_df_consecutif$ID %in% pathways, ]
    #keep requested conditions
    x <- x[x$group %in% conditions,]
    #reorder by order given in row_order
    x$group <- factor(x$group, levels = rev(row_order))
    
  
    pdf(file = paste(".", path, "pathway-expression", paste0(pattern_names[i], ".pdf"), sep="/"), width = 10)
    
    #plot
    plot <- ggplot(data=x, aes(x=NES, y=gsub("_", " ", ID), fill=group)) +
      geom_bar(stat="identity", position=position_dodge(), width = 0.5) +
      guides(fill = guide_legend(reverse=TRUE)) +
      scale_fill_manual(values = colors) +
      theme_minimal() +
      scale_x_continuous(limits = c(min(x$NES)-0.1, max(x$NES)+0.1), oob = scales::squish) + #set x-axis limits here
      scale_y_discrete(labels = function(x) str_wrap(x, width = 50)) +
      xlab("Normalised Enrichment Value") +
      ylab("Pathway") +
      ggtitle(pattern_names[i])
    
    
    print(plot)
    
    dev.off()
  }
}


#average expression of pathways for each pattern as barplot
#saved as pdf in path
average_expression <- function(conditions, path, NES_df_consecutif, pathwayList, pattern_names){
  
  createDirectories(paste(path, "average-expression", sep="/"))
  
  for (i in 1:length(pattern_names)) {
    
    
    # vector of pathways for the pattern
    pathways <- as.vector(pathwayList[i][[1]])
    
    #extract NES value from common dataframe
    x <- NES_df_consecutif[NES_df_consecutif$ID %in% pathways, ]
    #keep requested conditions
    x <- x[x$group %in% conditions,]
    #reorder by order given in row_order
    x$group <- factor(x$group, levels = row_order)
    
    
    
    y <- aggregate(x$NES~x$group, x, FUN=mean)
    colnames(y) <- c("group", "average")
    
    
    
    
    pdf(file = paste(".", path, "average-expression", paste0(pattern_names[i], ".pdf"), sep="/"), width = 10)
    plot <- ggplot(data=y, aes(x=average, y=gsub("_", " ", group), fill=group)) +
      geom_bar(stat="identity", position=position_dodge(), width = 0.5) +
      guides(fill = guide_legend(reverse=TRUE)) +
      scale_fill_manual(values = colors) +
      theme_minimal() +
      scale_x_continuous(limits = c(min(y$average)-0.1, max(y$average)+0.1), 
                         oob = scales::squish) +
      scale_y_discrete(labels = function(x) str_wrap(x, width = 50), limits=rev(intersect(row_order, conditions))) +
      xlab("Average NES") +
      ylab("group") +
      ggtitle(pattern_names[i])
    
    print(plot)
    
    dev.off()
  }
  
}


#average expression of pathways for each pattern as polygon plot
polygon_plot <- function(conditions, NES_df_consecutif, pathwayList, pattern_names){

  avgdf <- data.frame(matrix(ncol = length(conditions), nrow = 0))
  
  colnames(avgdf) <- conditions
  
  for (i in 1:length(pattern_names)) {
    group <- as.vector(pathwayList[i])
    
    # vector of pathways for the pattern
    pathways <- as.vector(pathwayList[i][[1]])
    
    #extract NES value from common dataframe
    x <- NES_df_consecutif[NES_df_consecutif$ID %in% pathways, ]
    #keep requested conditions
    x <- x[x$group %in% conditions,]
    #reorder by order given in row_order
    x$group <- factor(x$group, levels = rev(row_order))
  
    y <- aggregate(x$NES~x$group, x, FUN=mean)
    colnames(y) <- c("group", pattern_names[i])
    
    y2 <- data.frame(t(y[-1]))
    colnames(y2) <- y[, 1]
    y2 <- rev(y2)
  
    avgdf <- rbind(avgdf, y2)
  }
  avgdf <- cbind(pattern = rownames(avgdf), avgdf)
  
  
  
  borderDash <- list(c(0,0),c(20,3),c(20,3),c(0,0),c(3,3),c(3,3))   #for dotted line 
  
  
  
  spChartJSRadar(avgdf, scaleStartValue = min(avgdf[,2:ncol(avgdf)]), maxScale = max(avgdf[,2:ncol(avgdf)]),
                       polyAlpha = 0, BorderColor = unname(colors), borderDash = borderDash) 

}





g <- group_names[c(2,3,6)] #select conditions here

pathway_expression(g, "images/24vs48", NES_df_consecutif, pathwayList, pattern_names)
average_expression(g, "images/24vs48", NES_df_consecutif, pathwayList, pattern_names)
polygon_plot(g, NES_df_consecutif, pathwayList, pattern_names)


g <- group_names[c(1:6)] 

pathway_expression(g, "images/all-groups", NES_df_consecutif, pathwayList, pattern_names)
average_expression(g, "images/all-groups", NES_df_consecutif, pathwayList, pattern_names)
polygon_plot(g, NES_df_consecutif, pathwayList, pattern_names)









