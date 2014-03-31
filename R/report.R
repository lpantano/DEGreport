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

figurepvaluebyvarexp<-function(pvalues,counts,out){
  p<-pvalueVarMean(pvalues,counts)
  
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

figurebyvar<-function(tags,g1,g2,counts,out,pop=400){
  p<-varDElist(tags,g1,g2,counts,pop)
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
                   "This figure shows the distribution of SD
                   of the two groups." );
  return(FR)

}

plotrank<-function(tab){
  p<-expDElist(tags,g1,g2,counts,pop)
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


tablerank<-function(tab,path){
  countsFile<-"rank.txt"
  write.table(tab,paste0(path,"/rank.txt"),row.names=F,quote=F,sep="\t")
  TAB <- newTable(tab , file=countsFile, exportId="TABLE_COUNTS",
                  "Top genes" );
  return(TAB)
}

createReport<-function(g1,g2,counts,tags,pvalues,fc,path,pop=400){
  fg1<-figurepvaluebyexp(pvalues,counts,path)
  fg2<-figurepvaluebyvar(pvalues,counts,path)
  fg3<-figurepvaluebyvarexp(pvalues,counts,path)
  fg4<-figurebyexp(tags,g1,g2,counts,path,pop)
  fg5<-figurebyvar(tags,g1,g2,counts,path,pop)
  #figurebyvarvsexp()
  #figurecor()
  tabrank<-get_rank(g1,g2,counts[detags],fc)
  tb1<-tablerank(tabrank)
  
  report <- addTo( 
    report, addTo( newSection( "Quality of DEG", class="results" ),
                   addTo( newSubSection( "" ), fg1),
                   addTo( newSubSection( "" ), fg2),
                   addTo( newSubSection(""), fg3),    
                   addTo( newSubSection( "" ), fg4),
                   addTo( newSubSection( "" ), fg5),
                   addTo( newSubSection(""), tab1)
    ))
    
    writeReport( report, filename=paste(path,"/DEGReport.html",sep=""))
    
  
  
}

