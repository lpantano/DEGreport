context("Plots")

library(DESeq2)
data(humanSexDEedgeR)
idx <- c(1:5, 75:80)
dse <- DESeqDataSetFromMatrix(humanSexDEedgeR$counts[1:1000, idx],
                              humanSexDEedgeR$samples[idx,],
                              design = ~group) %>% DESeq
res <- results(dse)

test_that("degcovariates", {
    resCov <- degCovariates(log2(counts(dse) + 0.5), colData(dse), min_pc_pct = 10)
    expect_true(nrow(resCov[["corMatrix"]][resCov[["corMatrix"]][["fdr"]] < 0.05,]) == 3)
})

test_that("test_genes", {
    expect_true(degPlotWide(dse, rownames(dse)[1:10], group = "group") %>%
                    class %>% .[[2]] == "ggplot")
    expect_true(degPlot(dse, res = res, n = 3, xs = "group", group = "group") %>%
                    class %>% .[[2]] == "ggplot")
})

test_that("singleFunctions", 
    {

        expect_true(degMean(res[,4], counts(dse)) %>%
                        class %>% .[[2]] == "ggplot")
        expect_true(degVar(res[,4], counts(dse)) %>%
                              class %>% .[[2]] == "ggplot")
        expect_true(degMV(colData(dse)[["group"]], res[,4],
                          counts(dse)) %>%
                        class %>% .[[2]] == "ggplot")
        detag <- row.names(res[1:10,])
        expect_true(degMB(detag, colData(dse)[["group"]][1:5],
                          colData(dse)[["group"]][6:10],
                          counts(dse)) %>%
                        class %>% .[[2]] == "ggplot")
        expect_true(degVB(detag, colData(dse)[["group"]][1:5],
                          colData(dse)[["group"]][6:10],
                          counts(dse)) %>%
                        class %>% .[[2]] == "ggplot")
        expect_true(degVolcano(as.data.frame(res[,c("log2FoldChange", "pvalue")])) %>%
                        class %>% .[1] == "gtable")
    })
