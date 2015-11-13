test_singleFunctions <-
    function()
{
    data(DEGreportSet)
    checkTrue(class(degMean(DEGreportSet$deg[,4],
        DEGreportSet$counts))[[2]]=="ggplot")
    checkTrue(class(degVar(DEGreportSet$deg[,4],
        DEGreportSet$counts))[[2]]=="ggplot")
    checkTrue(class(degMV(DEGreportSet$g1,DEGreportSet$g2,DEGreportSet$deg[,4],
        DEGreportSet$counts))[[2]]=="ggplot")
    detag <- row.names(DEGreportSet$deg[1:10,])
    checkTrue(class(degMB(detag,DEGreportSet$g1,DEGreportSet$g2,
        DEGreportSet$counts))[[2]]=="ggplot")
    checkTrue(class(degVB(detag,DEGreportSet$g1,DEGreportSet$g2,
        DEGreportSet$counts))[[2]]=="ggplot")
}