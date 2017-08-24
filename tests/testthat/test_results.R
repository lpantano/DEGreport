context("results")
library(DESeq2)
dds <- makeExampleDESeqDataSet(betaSD=1)
colData(dds)[["treatment"]] <- sample(colData(dds)[["condition"]], 12)
design(dds) <-  ~ condition + treatment
dds <- DESeq(dds)
res <- results(dds)

test_that("Combinations", {
    expect_named(.guessComb(dds,
               combs = c("condition"),
               contrast = list("treatment_B_vs_A", c("condition", "A", "B")),
               pairs = TRUE), c("condition_1_vs_2",
                                "treatment_B_vs_A",
                                "condition_A_vs_B") )
    expect_type(.normalizeNames(c("condition", 2), dds), "character")
    expect_output(.createComb(dds, c("treatment", "condition")) %>% str,
                  "List of 2")
})

test_that("Results",{
    expect_match(.guessResults(dds, c("condition", "A", "B"), 0.05) %>% class,
                "DESeqResults")
    expect_match(.guessResults(dds, "condition_B_vs_A", 0.05) %>% class,
                "DESeqResults")
    expect_error(.guessResults(dds, "condition_C_vs_A", 0.05))
    expect_match(.guessShrunken(dds, c("condition", "A", "B"), res) %>% class,
                 "DESeqResults")
    expect_error(.guessResults(dds, "condition_C_vs_A", res))
    
})