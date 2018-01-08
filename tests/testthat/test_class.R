context("Classes and Methods")
library(tibble)

y <- matrix(rnbinom(80L, size = 1L, mu = 10L), nrow = 20L)
d <- DGEList(counts = y, group = rep(1L:2L, each = 2L),
             lib.size = rep(c(1000L:1001L), 2L))
rownames(d[["counts"]]) <- paste("gene", 1L:nrow(d[["counts"]]), sep=".")
d <- estimateCommonDisp(d)
de <- exactTest(d)
rawR <- topTags(de, n = Inf)

dds <- makeExampleDESeqDataSet(betaSD = 1)
colData(dds)[["treatment"]] <- sample(colData(dds)[["condition"]], 12)
design(dds) <-  ~ condition + treatment
dds <- DESeq(dds)
raw <- results(dds)
res <- degComps(dds, combs = c("condition"))


test_that("DEGSet",{

    expect_message(DEGSetFromDESeq2(raw, default = "raw"),
                "shrunken")
    expect_type(DEGSetFromDESeq2(raw, extras = list(shrunken = raw)), "list")
    expect_error(DEGSetFromDESeq2(raw, extras = list(shrunken = rawR)),
                 "DESeqResults")
    
    expect_message(DEGSetFromEdgeR(rawR, default = "raw"),
                "shrunken")
    expect_type(DEGSetFromEdgeR(rawR,  extras = list(shrunken = rawR)), "list")
    expect_error(DEGSetFromEdgeR(rawR,  extras = list(shrunken = raw)), "TopTags")
    
    expect_message(DEGSet(list(raw = raw), "raw"))
    
    expect_true(abs(deg(res[[1]], "raw")[["log2FoldChange"]][[1]]) >
                    abs(deg(res[[1]], "shrunken")[["log2FoldChange"]][[1]]))
    expect_match(degDefault(res[[1]]), "shrunken")
    expect_error(deg(res[[1]], "fake"))
    expect_identical(res[[1]][["shrunken"]], deg(res[[1]]))
    expect_true(deg(res[[1]], tidy = "tibble") %>% is_tibble)
    expect_equal(deg(res[[1]], top = 5) %>% nrow, 5)
    
    expect_type(significants(res[[1]]), "character")
    expect_type(significants(res[[1]][[1]]), "character")
    expect_type(significants(res), "character")
    expect_true(significants(res, full = TRUE) %>% is_tibble)
    
    expect_true(plotMA(res[[1]], diff = 4) %>%
                    class %>% .[[2]] == "ggplot")
    expect_true(plotMA(res[[1]], diff = 4, raw = TRUE) %>%
                    class %>% .[[2]] == "ggplot")
    expect_true(plotMA(res[[1]], diff = 4, correaltion = TRUE) %>%
                    class %>% .[[2]] == "ggplot")
    expect_error(plotMA(res[[1]], diff = 4, raw = TRUE, correlation = TRUE))
})