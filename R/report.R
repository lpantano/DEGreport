#' Create figure 1
#' 
#' @param pvalues  pvalues of DEG analysis
#' @param counts  matrix with counts for each samples and each gene. Should be same length than pvalues vector.
#' @param out path to save the figure
#' @return Nozzle object
figurepvaluebyexp<-function(pvalues,counts,out){
  p<-pvalueMean(pvalues,counts)
  
  File="fpvaluebyexp.jpg"
  HFile="fpvaluebyexp.pdf"
  jpeg(paste(out, File, sep="" ) ,width=600,height=400,quality=100 );
  print(p)
  dev.off();
  pdf(paste( out, HFile, sep="" )  );
  print(p)
  dev.off();
  
  # create a figure and make it available for exporting
  FR <- newFigure( File, fileHighRes=HFile, exportId="PVALBYEXP",
                      "This figure shows the distribution of pvalues according
                      the average expression of the feature." );
  return(FR)
}
#' Create figure 2
#' 
#' @param pvalues  pvalues of DEG analysis
#' @param counts  matrix with counts for each samples and each gene. Should be same length than pvalues vector.
#' @param out path to save the figure
#' @return Nozzle object
figurepvaluebyvar<-function(pvalues,counts,out){
  p<-pvalueVar(pvalues,counts)
  
  File="fpvaluebyvar.jpg"
  HFile="fpvaluebyvar.pdf"
  jpeg(paste(out, File, sep="" ) ,width=600,height=400,quality=100 );
  print(p)
  dev.off();
  pdf(paste( out, HFile, sep="" )  );
  print(p)
  dev.off();
  
  # create a figure and make it available for exporting
  FR <- newFigure( File, fileHighRes=HFile, exportId="PVALBYVAR",
                   "This figure shows the distribution of pvalues according
                   the SD of the feature." );
  return(FR)
}
#' Create figure 3
#' 
#' @param g1 list of samples in group 1
#' @param g2 list of samples in group 2
#' @param pvalues  pvalues of DEG analysis
#' @param counts  matrix with counts for each samples and each gene. Should be same length than pvalues vector.
#' @param out path to save the figure
#' @return Nozzle object
figurepvaluebyvarexp<-function(g1,g2,pvalues,counts,out){
  p<-pvalueVarMean(g1,g2,pvalues,counts)
  
  File="fpvaluebyvarexp.jpg"
  HFile="fpvaluebyvaexpr.pdf"
  jpeg(paste(out, File, sep="" ) ,width=600,height=400,quality=100 );
  print(p)
  dev.off();
  pdf(paste( out, HFile, sep="" )  );
  print(p)
  dev.off();
  
  # create a figure and make it available for exporting
  FR <- newFigure( File, fileHighRes=HFile, exportId="PVALBYVAREXP",
                   "This figure shows the distribution of pvalues according
                   the average expression and the SD of the feature." );
  return(FR)
}
#' Create figure 4
#' 
#' @param tags  genes of DEG analysis
#' @param g1 group 1
#' @param g2 group 2
#' @param counts  matrix with counts for each samples and each gene. Should be same length than pvalues vector.
#' @param out path to save the figure
#' @param pop random genes for background
#' @return Nozzle object
figurebyexp<-function(tags,g1,g2,counts,out,pop=400){
  p<-expDElist(tags,g1,g2,counts,pop)
  File="fexp.jpg"
  HFile="fexp.pdf"
  jpeg(paste(out, File, sep="" ) ,width=600,height=400,quality=100 );
  print(p)
  dev.off();
  pdf(paste( out, HFile, sep="" )  );
  print(p)
  dev.off();
  
  # create a figure and make it available for exporting
  FR <- newFigure( File, fileHighRes=HFile, exportId="EXP",
                   "This figure shows the distribution of expression values
                   of the two groups." );
  return(FR)

}
#' Create figure 5
#' 
#' @param tags  genes of DEG analysis
#' @param g1 group 1
#' @param g2 group 2
#' @param counts  matrix with counts for each samples and each gene. Should be same length than pvalues vector.
#' @param out path to save the figure
#' @param pop random genes for background
#' @return Nozzle object
figurebyvar<-function(tags,g1,g2,counts,out,pop=400){
  p<-varDElist(tags,g1,g2,counts,pop)
  File="fvar.jpg"
  HFile="fvar.pdf"
  jpeg(paste(out, File, sep="" ) ,width=600,height=400,quality=100 );
  print(p)
  dev.off();
  pdf(paste( out, HFile, sep="" )  );
  print(p)
  dev.off();
  
  # create a figure and make it available for exporting
  FR <- newFigure( File, fileHighRes=HFile, exportId="EXP",
                   "This figure shows the distribution of SD
                   of the two groups." );
  return(FR)

}
#' Create figure 6
#' 
#' @param tab  table from \link{get_rank}
#' @param out path to save the figure
#' @param colors colors for each gene
#' @return Nozzle object
figurerank<-function(tab,out,colors){
  p<-plotrank(tab,colors)
  File="fcor.jpg"
  HFile="fcor.pdf"
  jpeg(paste(out, File, sep="" ) ,width=600,height=400,quality=100 );
  print(p)
  dev.off();
  pdf(paste( out, HFile, sep="" )  );
  print(p)
  dev.off();
  
  # create a figure and make it available for exporting
  FR <- newFigure( File, fileHighRes=HFile, exportId="COR",
                   "This figure shows the correlation between
                   the rank by score and tha rank by FC" );
  return(FR)

}

#' Create table 1
#' 
#' @param tab  table from \link{get_rank}
#' @param out path to save the figure
#' @return Nozzle object
tablerank<-function(tab,out){
  countsFile<-"rank.txt"
  tab<-cbind(row.names(tab),tab)
  names(tab)<-c("Gene","mean FC","FC at 2.5%","FC at 97.5%",
                "Origial FC","score")
  write.table(tab,paste0(out,"/rank.txt"),row.names=F,quote=F,sep="\t")
  TAB <- newTable(tab , file=countsFile, exportId="TABLE_COUNTS",
                  "Top genes" );
  return(TAB)
}
#' Create report
#' 
#' @param g1 group 1
#' @param g2 group 2
#' @param counts  matrix with counts for each samples and each gene. Should be same length than pvalues vector.
#' @param tags  genes of DEG analysis
#' @param pvalues  pvalues of DEG analysis
#' @param fc FC for each gene
#' @param path path to save the figure
#' @param colors data frame with colors for each gene
#' @param pop random genes for background
#' @return nothing
createReport<-function(g1,g2,counts,tags,pvalues,fc,path,colors="",pop=400){
  fg1<-figurepvaluebyexp(pvalues,counts,path)
  fg2<-figurepvaluebyvar(pvalues,counts,path)
  fg3<-figurepvaluebyvarexp(g1,g2,pvalues,counts,path)
  fg4<-figurebyexp(tags,g1,g2,counts,path,pop)
  fg5<-figurebyvar(tags,g1,g2,counts,path,pop)
  #figurebyvarvsexp()
  #figurecor()
  tabrank<-get_rank(g1,g2,counts[tags,],fc,pop)
  tb1<-tablerank(tabrank,path)
  fg6<-figurerank(tabrank,path,colors)
  
  report<-""
  report <- newCustomReport( "DEG Report " );
  report <- addTo( 
    report, addTo( newSection( "Quality of DEG", class="results" ),
                   addTo( newSubSection( "Pvalue vs abundance" ), fg1),
                   addTo( newSubSection( "Pvalue vs variation" ), fg2),
                   addTo( newSubSection("Pvalue vs abundance/variation"), fg3),    
                   addTo( newSubSection( "Abundance distribution" ), fg4),
                   addTo( newSubSection( "Variation distribution" ), fg5),
                   addTo( newSubSection("Rank"), tb1),
                   addTo( newSubSection("FC vs rank"), fg6)
    ))
    
    writeReport( report, filename=paste(path,"/DEGReport",sep=""))
    
  return(0)
  
}

