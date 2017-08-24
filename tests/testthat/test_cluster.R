context("Clustering")

library(DESeq2)
data(humanSexDEedgeR)
idx <- c(1:5, 75:80)
dse <- DESeqDataSetFromMatrix(humanSexDEedgeR$counts[1:1000, idx],
                              humanSexDEedgeR$samples[idx,],
                              design = ~group) %>% DESeq

test_that("transform", {
    expect_gte(mean(.scale(counts(dse)[1,])), -0.5)
    expect_lte(mean(.scale(counts(dse)[1,])), 0.5)
    countsGroup <- .summarize_scale(counts(dse)[1:100,],
                                    colData(dse)[["group"]])
    expect_equal(mean(countsGroup), 0)
    hc <- .make_clusters(countsGroup)
    cluster <- .select_genes(hc, countsGroup, minc = 5)
    expect_equal(unique(cluster), c(1, 2))
    df <- data.frame(cluster = cluster, genes = names(cluster))
    expect_equal(.median_per_cluster(countsGroup, df) %>% mean, 0)
})
    