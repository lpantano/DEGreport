test_degcovariates <- function(){
    library(DESeq2)
    data(humanSexDEedgeR)
    idx <- c(1:5, 75:80)
    dse <- DESeqDataSetFromMatrix(humanSexDEedgeR$counts[1:1000, idx],
                                  humanSexDEedgeR$samples[idx,],
                                  design=~group)
    res <- degCovariates(log2(counts(dse) + 0.5), colData(dse))
    checkTrue(nrow(res$cor_matrix[res$cor_matrix$fdr < 0.05,]) == 3)
}

test_genes <- function(){
    data(humanSexDEedgeR)
    library(DESeq2)
    idx <- c(1:10, 75:85)
    dse <- DESeqDataSetFromMatrix(humanSexDEedgeR$counts[1:1000, idx],
                                  humanSexDEedgeR$samples[idx,], design=~group)
    dse <- DESeq(dse)
    res <- results(dse)
    checkTrue(class(degPlotWide(dse, rownames(dse)[1:10], group="group"))[[2]] == "ggplot")
    checkTrue(class(degPlot(dse, res = res, n = 3, xs = "group", group="group"))[[2]] == "ggplot")
}

test_singleFunctions <- 
    function() 
{
    data(DEGreportSet)
    data(humanSexDEedgeR)
    
    checkTrue(class(degMean(DEGreportSet$deg[,4],
        DEGreportSet$counts))[[2]] == "ggplot")
    checkTrue(class(degVar(DEGreportSet$deg[,4],
        DEGreportSet$counts))[[2]] == "ggplot")
    checkTrue(class(degMV(humanSexDEedgeR$samples$group,DEGreportSet$deg[,4],
        DEGreportSet$counts))[[2]] == "ggplot")
    detag <- row.names(DEGreportSet$deg[1:10,])
    checkTrue(class(degMB(detag,DEGreportSet$g1,DEGreportSet$g2,
        DEGreportSet$counts))[[2]] == "ggplot")
    checkTrue(class(degVB(detag,DEGreportSet$g1,DEGreportSet$g2,
        DEGreportSet$counts))[[2]] == "ggplot")
    checkTrue(class(degVolcano(DEGreportSet$deg[,c("logFC", "PValue")]))[1] == "gtable")
    }
