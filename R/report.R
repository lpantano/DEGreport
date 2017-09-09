# Wrap figure from \code{degMV} into a Nozzle object
figurepvaluebyvarexp <- 
    function(g, pvalues, counts, out)
    {
        p <- degQC(counts, g, pvalue = pvalues)
        File = "fpvaluebyvarexp.jpg"
        HFile = "fpvaluebyvaexpr.pdf"
        jpeg(file.path(out, File) ,width=600, height=400, quality=100 );
        print(p)
        dev.off();
        pdf(file.path( out, HFile));
        print(p)
        dev.off();
        
        # create a figure and make it available for exporting
        newFigure( File, fileHighRes=HFile, exportId="PVALBYVAREXP",
                   "This figure shows the distribution of pvalues according
                   the average expression and the variability of the feature. 
                   It is taking the maximum value among group1 and group2" );
    }
# Wrap figure from \code{degMB} into a Nozzle object
figurebyexp <- 
    function(tags, g, counts, out, pop=400)
    {
        p <- degMB(tags, g, counts, pop)
        File = "fexp.jpg"
        HFile = "fexp.pdf"
        jpeg(file.path(out, File), width = 600, height = 400, quality = 100 );
        print(p)
        dev.off();
        pdf(file.path(out, HFile));
        print(p)
        dev.off();
        
        # create a figure and make it available for exporting
        newFigure(File, fileHighRes=HFile, exportId="EXP",
                   "This figure shows the distribution of expression values
                   of the two groups.");
    }
# Wrap figure from \code{degVB} into a Nozzle object
figurebyvar <- 
    function(tags, g, counts, out, pop=400)
    {
        p <- degVB(tags, g, counts, pop)
        File = "fvar.jpg"
        HFile = "fvar.pdf"
        jpeg(file.path(out, File), width = 600, height = 400, quality = 100);
        print(p)
        dev.off();
        pdf(file.path(out, HFile));
        print(p)
        dev.off();
        
        # create a figure and make it available for exporting
        newFigure(File, fileHighRes=HFile, exportId="EXP",
                   "This figure shows the distribution of SD
        of the two groups.");
    }

#' Create report of RNAseq DEG anlaysis
#' 
#' @description This function get the count matrix, pvalues, and FC of a 
#' DEG analysis and create a report to help to detect possible problems 
#'   with the data.
#' @aliases createReport
#' @param g Character vector with the group the samples belong to.
#' @param counts Matrix with counts for each samples and each gene. 
#'   Should be same length than pvalues vector.
#' @param tags Genes of DEG analysis
#' @param pvalues pvalues of DEG analysis
#' @param path path to save the figure
#' @param pop random genes for background
#' @param name name of the html file
#' @return A HTML file with all figures and tables
#' @export
createReport <- 
    function(g, counts, tags, pvalues, path,
             pop = 400, name="DEGreport")
    {
        fg3 <- figurepvaluebyvarexp(g, pvalues, counts, path)
        fg4 <- figurebyexp(tags, g, counts, path, pop)
        fg5 <- figurebyvar(tags, g, counts, path, pop)
        report <- ""
        report  <- newCustomReport( "DEG Report " );
        report  <- addTo( 
            report, addTo( newSection( "Quality of DEG", class="results" ),
                           addTo(newSubSection("Pvalue vs abundance/variation"), fg3),    
                           addTo(newSubSection("Abundance distribution"), fg4),
                           addTo(newSubSection("Variation distribution"), fg5)
            ))
        
        writeReport(report, filename = file.path(path, name))
    }

