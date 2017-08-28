context("Classes and Methods")

test_that("DEGSet",{
    library(DESeq2)
    dds <- makeExampleDESeqDataSet(betaSD = 1)
    colData(dds)[["treatment"]] <- sample(colData(dds)[["condition"]], 12)
    design(dds) <-  ~ condition + treatment
    dds <- DESeq(dds)
    res <- degComps(dds, combs = c("condition"))

    expect_true(abs(degTable(res[[1]], "raw")[["log2FoldChange"]][[1]]) >
                    abs(degTable(res[[1]], "shrunk")[["log2FoldChange"]][[1]]))
    expect_match(degDefault(res[[1]]), "shrunk")
    expect_error(degTable(res[[1]], "fake"))
    expect_identical(res[[1]][["shrunk"]], degTable(res[[1]]))
    expect_type(degSign(res[[1]]), "character")
})