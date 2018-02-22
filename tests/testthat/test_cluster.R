context("Clustering")

data(humanGender)
idx <- c(1:5, 75:80)
counts <- assays(humanGender)[[1]]
dse <- DESeqDataSetFromMatrix(counts[1:1000, idx],
                              colData(humanGender)[idx,],
                              design = ~group) %>% DESeq

test_that("transform", {
    expect_gte(mean(.scale(counts(dse)[1,])), -0.5)
    expect_lte(mean(.scale(counts(dse)[1,])), 0.5)
    countsGroup <- .summarize_scale(counts(dse)[1:100,],
                                    colData(dse)[["group"]])
    expect_equal(mean(countsGroup), 0)
    hc <- .make_clusters(countsGroup)
    cluster0 <- .select_genes(hc, countsGroup, minc = 5, reduce = TRUE)
    cluster <- .select_genes(hc, countsGroup, minc = 5)
    expect_equal(unique(cluster), c(1, 2))
    expect_equal(unique(cluster0), c(1, 2))
    df <- data.frame(cluster = cluster, genes = names(cluster))
    expect_equal(.median_per_cluster(countsGroup, df) %>% mean, 0)
    expect_equal(.filter(df, 50) %>% nrow, 69)
    expect_equal(.group_metadata(as.data.frame(colData(dse)),
                                 "group", "group", "group") %>% nrow, 2)
})

test_that("groupDifference", {
    ma <- matrix(rnorm(50), ncol = 2)
    ma[1:20, 2] <- 1000 + ma[1:20, 2]
    expect_equal(.remove_low_difference(ma, 500) %>% nrow(), 20)
})
