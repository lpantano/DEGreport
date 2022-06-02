library(getopt)
library(DESeq2)
library(DEGreport)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(plotly)

# arguments
spec <- matrix(c(
    'deseq2Object', 'd', 1, "character", "DESeq2 object (required)",
    'qValue', 'q', 1, "numeric", "Q-value threshold to filter elements from the DESeq2 object (required)",
    'log2FoldChange', 'l', 1, "numeric", "Min absolute log2 fold-change threshold to filter elements from the DESeq2 object (required)",
    'minElements', 'm', 1, "numeric", "Min number of elements to form a cluster (required)",
    'explanatoryVar', 'e', 1, "character", "Explanatory variable in DESeq2 object (e.g. time). Will be used on the x-axis of the plot (required)",
    'groupVar', 'g', 1, "character", "Group variable in DESeq2 object (e.g. genotype). Determines how samples are grouped (optional)",
    'clusterMethod', 'c', 1, "character", "Cluster method. Either consensus or diana (required)",
    'removeOutliers', 'r', 1, "logical", "Remove outliers from the cluster distribution. (required)",
    'scale', 's', 1, "logical", "Scale DESeq2 normalized count matrix. (required)",
    'plotsPerColumn', 'n', 1, "numeric", "Number of plots per column in `summaryPlot´ file (required)",
    'plotsPerRow', 'w', 1, "numeric", "Number of plots per row in `summaryPlot´ file (required)",
    'degReport', 'o', 1, "character", "Output degReport table (required)",
    'elementClusterMap', 'u', 1, "character", "Output table with element - cluster relationship (required)",
    'clusterCount', 't', 1, "character", "Output table summarizing the number of elements per cluster (required)",
    'summaryPlot', 'p', 1, "character", "Output summary cluster plot (required)",
    'plotlyDirectory', 'y', 1, "character", "Output directory to store plotly cluster plots (required)",
    'help', 'h',0, "logical", "this help"),ncol=5,byrow=T)
opt = getopt(spec)

# Check parameters
if (!is.null(opt$help) || is.null(opt$deseq2Object) || is.null(opt$qValue) 
    || is.null(opt$log2FoldChange) || is.null(opt$minElements) || is.null(opt$explanatoryVar) 
    || is.null(opt$clusterMethod) || is.null(opt$removeOutliers) || is.null(opt$scale) || is.null(opt$plotsPerColumn) 
    || is.null(opt$plotsPerRow) || is.null(opt$degReport) || is.null(opt$elementClusterMap) || is.null(opt$clusterCount) 
    || is.null(opt$summaryPlot) || is.null(opt$plotlyDirectory)) {
    stop(cat(paste(getopt(spec, usage=T),"\n")));
}

# Libraries version
print(paste("DESeq2 version: ", package.version("DESeq2"), sep= ""))
print(paste("DEGreport version: ", package.version("DEGreport"), sep= ""))
print(paste("dplyr version: ", package.version("dplyr"), sep= ""))

# deseq2Object
print("Reading DESeq2 object...")
deseq2Object <- readRDS(file = opt$deseq2Object)

# cluster analysis
print("Extracting desing used to generate DESeq2 object...")
design <- as.data.frame(colData(deseq2Object))
print(design)

# check explanatory variable
print(paste("Checking explanatory variable '", opt$explanatoryVar, "' is present in the DESeq2 design...", sep=""))
if (!(opt$explanatoryVar %in% names(design))){
    stop(paste(opt$explanatoryVar, " is not present in DESeq2 object's design.", sep=""))
}

# relevel factor based on design
print("Releveling explanatory variable based on DESeq2 design...")
expVar <- dplyr::pull(design, opt$explanatoryVar)
expVar <- factor(expVar, unique(as.character(expVar)))
design[which(names(design) == opt$explanatoryVar)] <- expVar

# check group variable
if (!is.null(opt$groupVar)) {
    print(paste("Checking group variable '", opt$groupVar, "' is present in the DESeq2 design...", sep=""))
    if (!(opt$groupVar %in% names(design))){
        stop(paste(opt$design, " is not present in DESeq2 object's design.", sep=""))
    }
}

# results
print("Extracting results from DESeq2 object...")
res <- results(deseq2Object)
resDataFrame <- as.data.frame(res)
resDataFrame <- resDataFrame[complete.cases(resDataFrame),]
print(paste("Number of elements in DESeq2 object: ", nrow(resDataFrame), sep=""))

# log2 threshold
print(paste("Filtering elements in DESeq2 results table by q-value: ", opt$qValue, sep=""))
resDataFrame <- resDataFrame[(resDataFrame$padj < opt$qValue),]
print(paste("Filtering elements in DESeq2 results table by log2FoldChange: ", opt$log2FoldChange, sep=""))
resDataFrame <- resDataFrame[(abs(resDataFrame$log2FoldChange) > opt$log2FoldChange),]
print(paste("Number of elements in DESeq2 object after filtering: ", nrow(resDataFrame), sep=""))

if (nrow(resDataFrame) < opt$minElements) {
    stop(paste("Not enough elements to create clusters: ", nrow(resDataFrame), sep=""))
}

# rlog
print("Applying log2 transformation to count matrix ...")
ma <- assay(rlog(deseq2Object))[row.names(resDataFrame),]

# degPatterns
print("Running degPatterns...")
print(paste("Min number of elements to form a cluster: ", opt$minElements, sep=""))
print(paste("Remove outliers of the cluster distribution: ", opt$removeOutliers, sep=""))
print(paste("Scale normalized count matrix: ", opt$scale, sep=""))
if (opt$clusterMethod == "consensus") {
    print("Cluster Method: ConsensusClusterPlus")
    opt$clusterMethod <- TRUE
} else {
    print("Cluster Method: cluster::diana()")
    opt$clusterMethod <- FALSE
}
clusters <- degPatterns(ma, design, time = opt$explanatoryVar, consensusCluster = opt$clusterMethod, 
                            reduce = opt$removeOutliers, col = opt$groupVar, minc = opt$minElements, scale = opt$scale)

# degReport
print("Exporting degReport...")
degReport <- clusters$normalized
names(degReport) <- c("Id", names(degReport)[2:ncol(degReport)])
write.table(degReport, file = opt$degReport, sep = "\t", quote = FALSE, row.names = FALSE)

# elementClusterMap
print("Exporting element-cluster map table...")
clusteredElements <- paste(degReport$Id,degReport$cluster,sep = "_-_")
uniqueClusteredElements <- unique(clusteredElements)
uniqueClusteredElementsSplit <- strsplit(uniqueClusteredElements, "_-_")
element <- vector(mode = "character", length = length(uniqueClusteredElementsSplit))
clusterElement <- vector(mode = "character", length = length(uniqueClusteredElementsSplit))
for(elementIndex in 1:length(uniqueClusteredElementsSplit)) {
  element[elementIndex] <- uniqueClusteredElementsSplit[[elementIndex]][1]
  clusterElement[elementIndex] <- uniqueClusteredElementsSplit[[elementIndex]][2]
}
elementClusterMap <- data.frame(id = element, cluster = clusterElement)
write.table(elementClusterMap, file = opt$elementClusterMap, sep = "\t", quote = FALSE, row.names = FALSE)

# clusterCount
elementClusterMap$cluster <- as.character(elementClusterMap$cluster)
clusterCount <-  count(elementClusterMap, cluster)
names(clusterCount) <- c("cluster", "count")
write.table(clusterCount, file = opt$clusterCount, sep = "\t", quote = FALSE, row.names = FALSE)

# clusterPlot
print("Generating single cluster plots...")
numberOfCluster <- sort(unique(degReport$cluster))
clusterPlotList <- list()
plotCount <- 0
groupAndColorVariable <- opt$groupVar
if(is.null(opt$groupVar)) {
    groupAndColorVariable <- opt$explanatoryVar
}
print("Cluster identified:")
print(numberOfCluster)
for(clusterId in numberOfCluster) {
    plotCount <- plotCount + 1
    clusterSubset <- degReport %>% dplyr::filter(cluster == clusterId)
    clusterSubsetSummary <- clusterSubset %>% group_by_at(c(opt$explanatoryVar, opt$groupVar)) %>% summarize (mean = mean(value), sd = sd(value)) %>% as.data.frame()
    clusterPlot <- ggplot(data = clusterSubsetSummary, aes_string(x = opt$explanatoryVar, y = "mean", group = groupAndColorVariable, color= groupAndColorVariable)) +
                        geom_line(size=1) +
                        geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.5, position=position_dodge(0.05)) +
                        geom_point() +
                        ggtitle(paste("Cluster ", clusterId, sep = "")) +
                        theme(legend.title = element_text(size = 7)) +
                        theme(axis.title.x = element_text(size = 7)) +
                        theme(axis.title.y = element_text(size = 7)) +
                        theme(legend.text = element_text(size = 5)) +
                        theme(axis.text = element_text(size = 5))
    clusterPlotList[[plotCount]] <- clusterPlot
    clusterPlotly <- ggplotly(clusterPlot) %>% config(displaylogo=FALSE,
                                            modeBarButtonsToRemove  = list("sendDataToCloud", "hoverClosest3d", 
                                                                            "autoScale2d", "hoverClosestCartesian", "hoverCompareCartesian", "zoom2d", 
                                                                            "pan2d", "select2d", "lasso2d", "zoomIn2d", "zoomOut2d", "toggleSpikelines"))
    samplePlotlyFile <- paste(opt$plotlyDirectory, "cluster-", clusterId, ".html", sep = "")                                                     
    htmlwidgets::saveWidget(clusterPlotly, samplePlotlyFile)
}
# summaryPlot
print("Generating cluster plots summary...")
pdf(onefile = TRUE,  file = opt$summaryPlot)
if (length(clusterPlotList) > 0) {
    names(clusterPlotList) <- seq_along(clusterPlotList)
    clusterPlotList[sapply(clusterPlotList, is.null)] <- NULL
    plotIndex <- 0
    pageCount <- 1
    plotCountLimit <- (opt$plotsPerColumn * opt$plotsPerRow)
    while((plotIndex + plotCountLimit) < length(clusterPlotList)) {
        print(paste("Summary plot page: ", pageCount, sep = ""))
        if (length(clusterPlotList) > plotIndex + plotCountLimit) {
            plotsInPageList <- clusterPlotList[(plotIndex + 1):(plotIndex + plotCountLimit)]
            summaryPlot <- ggarrange(plotlist = plotsInPageList, nrow = opt$plotsPerColumn, ncol = opt$plotsPerRow)
            print(summaryPlot)
        }
       plotIndex <- plotIndex + plotCountLimit
       pageCount <- pageCount + 1
    }
    if (plotIndex < length(clusterPlotList)) {
        print(paste("Summary plot page: ", pageCount, sep = ""))
        plotsInPageList <- clusterPlotList[(plotIndex + 1):length(clusterPlotList)]
        summaryPlot <- ggarrange(plotlist = plotsInPageList, nrow = opt$plotsPerColumn, ncol = opt$plotsPerRow)
        print(summaryPlot)
    }
}
dev.off()
