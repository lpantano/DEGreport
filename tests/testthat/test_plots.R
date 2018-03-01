context("Plots")

data(humanGender)
idx <- c(1:5, 75:80)
counts <- assays(humanGender)[[1]]
dse <- DESeqDataSetFromMatrix(counts[1:1000, idx],
                              colData(humanGender)[idx,],
                              design = ~group) %>% DESeq

res <- results(dse)

plot_cor <- function(dse){
    ggplot(as.data.frame(counts(dse)),
           aes(x = NA20502, y = NA20504)) +
        geom_point() +
        geom_cor(method = "kendall")
}

test_that("degcovariates", {
    resCov <- degCovariates(log2(counts(dse) + 0.5), colData(dse), minPC = 10)
    expect_true(nrow(resCov[["corMatrix"]][resCov[["corMatrix"]][["fdr"]] < 0.05,]) == 3)
    expect_error(degCovariates(log2(counts(dse) + 0.5), colData(dse)[2:1,]))
})

test_that("testGenes", {
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
        expect_true(degMB(detag, colData(dse)[["group"]],
                          counts(dse)) %>%
                        class %>% .[[2]] == "ggplot")
        expect_true(degVB(detag, colData(dse)[["group"]],
                          counts(dse)) %>%
                        class %>% .[[2]] == "ggplot")
        expect_true(degVolcano(as.data.frame(res[,c("log2FoldChange", "pvalue")])) %>%
                        class %>% .[[2]] == "ggplot")
        expect_true(degPCA(counts(dse)) %>%
                        class %>% .[[2]] == "ggplot")
        expect_true(degMDS(counts(dse)) %>%
                        class %>% .[[2]] == "ggplot")
       
        expect_true(plot_cor(dse) %>%
                        class %>% .[[2]] == "ggplot")
    })
