context("Classes and Methods")
library(tibble)
library(DESeq2)

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

    expect_type(as.DEGSet(raw, extras = list(shrunken = raw)), "list")
    expect_error(as.DEGSet(raw, extras = list(shrunken = rawR)),
                 "DESeqResults")
    
    expect_type(as.DEGSet(rawR,  extras = list(shrunken = rawR)), "list")
    expect_error(as.DEGSet(rawR,  extras = list(shrunken = raw)), "TopTags")
    
    expect_true(abs(deg(res, "raw")[["log2FoldChange"]][[1]]) >
                    abs(deg(res, "shrunken")[["log2FoldChange"]][[1]]))
    expect_match(degDefault(res), "shrunken")
    expect_error(deg(res, "fake"))
    expect_identical(res[["shrunken"]], deg(res))
    expect_true(deg(res, tidy = "tibble") %>% is_tibble)
    expect_equal(deg(res, top = 5) %>% nrow, 5)
    
    expect_type(significants(res), "character")
    expect_type(significants(res[[1]]), "character")
    expect_type(significants(res), "character")
    expect_true(significants(res, full = TRUE) %>% is_tibble)
    
    # expect_s4_class(degCorrect(res, "lfdr-stat"), "DEGSet")
    
    expect_true(degMA(res, diff = 4) %>%
                    class %>% .[[2]] == "ggplot")
    expect_true(degMA(res, diff = 4, raw = TRUE) %>%
                    class %>% .[[2]] == "ggplot")
    expect_true(degMA(res, diff = 4, correlation = TRUE) %>%
                    class %>% .[[2]] == "ggplot")
    expect_error(degMA(res, diff = 4, raw = TRUE, correlation = TRUE))
})