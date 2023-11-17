context("Clustering")
library(DESeq2)

data(humanGender)
idx <- c(1:5, 75:80)
counts <- assays(humanGender)[[1]]
dse <- DESeqDataSetFromMatrix(counts[1:1000, idx],
                              colData(humanGender)[idx,],
                              design = ~group) %>% DESeq
ma <- assays(dse)[[1]]
des <- colData(dse)
des[["other"]] <- sample(c("a", "b"), length(idx), replace = TRUE)
res <- degPatterns(ma, des, time="group", col = "other", plot = FALSE)

test_that("cluster_plot", {
    expect_error(degPlotCluster(res, "group", "other"))
    expect_is(degPlotCluster(res$normalized, "group", "other"), "gg")
})

# test_that("missingfactor", {
#    des_missing <- des
#    levels(des_missing$group) <- c(levels(des_missing$group), "missing")
#    expect_is(degPatterns(ma, des_missing, time="group", col = NULL, plot = FALSE), "list")
# 
# })

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
    expect_equal(.remove_low_difference(ma, 500, FALSE) %>% nrow(), 20)
})

test_that("pct_variance", {
    clusters <- c(rep(1, 500), rep(2, 500))
    names(clusters) <- row.names(counts[1:1000, idx])
    pct <- .pct_var(counts[1:1000, idx], clusters)
    expect_true(pct < 100)
})

test_that("process", {
    library(dplyr)
    library(tidyr)
    library(tibble)
    table <- rownames_to_column(as.data.frame(ma), "genes") %>%
        gather("sample", "expression", -genes) %>%
        right_join(distinct(res$df[,c("genes", "cluster")]),
                   by = "genes") %>%
        left_join(rownames_to_column(as.data.frame(des), "sample"),
                  by = "sample")
    expect_true("value"  %in%  names(.process(table, "group", NULL)))
})
