context("Classes and Methods")
library(tibble)

test_that("DEGSet",{
    library(DESeq2)
    dds <- makeExampleDESeqDataSet(betaSD = 1)
    colData(dds)[["treatment"]] <- sample(colData(dds)[["condition"]], 12)
    design(dds) <-  ~ condition + treatment
    dds <- DESeq(dds)
    res <- degComps(dds, combs = c("condition"))

    expect_true(abs(deg(res[[1]], "raw")[["log2FoldChange"]][[1]]) >
                    abs(deg(res[[1]], "shrunk")[["log2FoldChange"]][[1]]))
    expect_match(degDefault(res[[1]]), "shrunk")
    expect_error(deg(res[[1]], "fake"))
    expect_identical(res[[1]][["shrunk"]], deg(res[[1]]))
    expect_true(deg(res[[1]], tidy = "tibble") %>% is_tibble)
    expect_equal(deg(res[[1]], top = 5) %>% nrow, 5)
    expect_type(significants(res[[1]]), "character")
    expect_true(plotMA(res[[1]], diff = 4) %>%
                    class %>% .[[2]] == "ggplot")
    expect_true(plotMA(res[[1]], diff = 4, raw = TRUE) %>%
                    class %>% .[[2]] == "ggplot")
})