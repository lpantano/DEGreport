#' Distribution of pvalues by expression range
#' 
#' @param pvalues  pvalues of DEG analysis
#' @param counts  matrix with counts for each samples and each gene. Should be same length than pvalues vector.
#' @return ggplot2 object
pvalueMean<-function(pvalues,counts){
  
  meanv<-apply(counts,1,mean)
  meanvfac<-cut(meanv,
                breaks=c(0,quantile(meanv,seq(.1,1,.1))),
                labels=paste0(seq(10,100,10),"%"),right=T)
  pvalfac<-cut(pvalues,
               breaks=c(0,seq(.1,1,.1)),
               labels=seq(0.1,1,0.1),right=T)
  d<-data.frame(pvalues=factor(pvalfac),
                meanv=factor(meanvfac))
  p<-ggplot(d,aes(pvalues,fill=meanv))+  
    geom_bar( position="fill")+
    theme_bw() + scale_fill_brewer(palette="RdYlBu")
  return(p)
}
#' Distribution of pvalues by SD range
#'
#' @param pvalues  pvalues of DEG analysis
#' @param counts  matrix with counts for each samples and each gene. Should be same length than pvalues vector.
#' @return ggplot2 object
pvalueVar<-function(pvalues,counts){
  
  sdv<-apply(counts,1,sd)
  sdvfac<-cut(sdv,
                breaks=c(0,quantile(sdv,seq(.1,1,.1))),
                labels=paste0(seq(10,100,10),"%"),right=T)
  pvalfac<-cut(pvalues,
               breaks=c(0,seq(.1,1,.1)),
               labels=seq(0.1,1,0.1),right=T)
  d<-data.frame(pvalues=factor(pvalfac),
                sdv=factor(sdvfac))
  p<-ggplot(d,aes(pvalues,fill=sdv))+  
    geom_bar()+
    theme_bw() + scale_fill_brewer(palette="RdYlBu")
  return(p)
}
#' Correlation of SD and mean of a set of genes
#'
#' @param g1 list of samples in group 1
#' @param g2 list of samples in group 2
#' @param pvalues  pvalues of DEG analysis
#' @param counts  matrix with counts for each samples and each gene. Should be same length than pvalues vector.
#' @return ggplot2 object
pvalueVarMean<-function(g1,g2,pvalues,counts){
  
  sdt1<-apply(counts[,g1],1,sd)
  sdt2<-apply(counts[,g2],1,sd)
  sdv<-apply(cbind(sdt1,sdt2),1,max)
  mt1<-apply(counts[,g1],1,mean)
  mt2<-apply(counts[,g2],1,mean)  
  meanv<-apply(cbind(mt1,mt2),1,max)
  pv<-cut(pvalues,breaks=c(-1,0.01,1.1),
          labels=c("<0.01","NoSig"))
  d<-data.frame(pvalues=pv,sdv=log2(sdv),
                meanv=log2(meanv))
  p<-ggplot(d,aes(meanv,sdv,colour=pvalues))+  
    geom_point()+
    scale_color_manual(values=c("red",rgb(0.9,0.9,0.9,0.6)))+
    theme_bw()+
    stat_quantile(aes(meanv,sdv),colour="blue",
                  quantiles = c(0.025,0.975),
                  linetype=2)
  return(p)
}
#' Distribution of DE genes expression compared the background
#'
#' @param tags  list of genes that are DE
#' @param g1 list of samples in group 1
#' @param g2 list of samples in group 2
#' @param counts  matrix with counts for each samples and each gene. Should be same length than pvalues vector.
#' @param pop number of random samples taken for background comparison
#' @return ggplot2 object
expDElist<-function(tags,g1,g2,counts,pop=400){
  delen<-length(tags)
  g<-""
  rand<-sample(row.names(counts),pop)
  g1var<-apply(counts[tags,g1],1,mean)    
  g2var<-apply(counts[tags,g2],1,mean)    
  
  rand.s1<-apply(counts[rand,g1],1,mean)
  rand.s2<-apply(counts[rand,g2],1,mean)
  res<-rbind(data.frame(g="g1",mean=g1var),data.frame(g="g2",mean=g2var))
  res<-rbind(res,data.frame(g="r1",mean=rand.s1),
             data.frame(g="r2",mean=rand.s2))
  
  p<-ggplot(res,aes(g,mean,fill=g,colour=g))+  
    geom_violin(alpha=0.2)+
    scale_y_log10()+
    theme_bw()
  return(p)
}
#' Distribution of DE genes SD compared the background
#'
#' @param tags  list of genes that are DE
#' @param g1 list of samples in group 1
#' @param g2 list of samples in group 2
#' @param counts  matrix with counts for each samples and each gene. Should be same length than pvalues vector.
#' @param pop number of random samples taken for background comparison
#' @return ggplot2 object
varDElist<-function(tags,g1,g2,counts,pop=400){
  delen<-length(tags)
  g<-""
  rand<-sample(row.names(counts),pop)
  g1var<-apply(counts[tags,g1],1,sd)    
  g2var<-apply(counts[tags,g2],1,sd)    
  
  rand.s1<-apply(counts[rand,g1],1,sd)
  rand.s2<-apply(counts[rand,g2],1,sd)
  res<-rbind(data.frame(g="g1",var=g1var),data.frame(g="g2",var=g2var))
  res<-rbind(res,data.frame(g="r1",var=rand.s1),
             data.frame(g="r2",var=rand.s2))
  
  p<-ggplot(res,aes(g,var,fill=g,colour=g))+  
    geom_violin(alpha=0.2)+
    scale_y_log10()+
    theme_bw()
  return(p)
}

#' Get number of potential combinations of two vectors
#'
#' @param g1 list of samples in group 1
#' @param g2 list of samples in group 2
#' @return maximum number of combinations of two vectors
combinations<-function(g1,g2){
  return(g1*g2)
}
#' Get random combinations of two groups
#'
#' @param g1 list of samples in group 1
#' @param g2 list of samples in group 2
#' @param pop number of combinations to be return
#' @return matrix with different combinatios of two vector 
getcomb<-function(g1,g2,pop){
  if (combinations(length(g1),length(g2))<pop){
    return(expand.grid(a=g1,b=g2))
  }else{
    #take 10% of each group, 400 times
    g1p=3
    g2p=3
    g1p<-ceiling(0.1*length(g1)+1)
    if(g1p==2){g1p=3}
    g2p<-ceiling(0.1*length(g2)+1)
    if(g2p==2){g2p=3}
    r1<-sapply(1:pop,function(x){sample(g1,g1p)})
    r2<-sapply(1:pop,function(x){sample(g2,g2p)})
    return(list(r1,r2))
  }
}
#' get the FC for each gene between two groups
#'
#' @param g1 list of samples in group 1
#' @param g2 list of samples in group 2
#' @param genes list of genes to be analized
#' @param popsize number of combinations to generate
#' @return FC for different combinations of samples in each group for each gene
get_ratio<-function(g1,g2,genes,popsize){
  pop<-getcomb(g1,g2,popsize)
  if (is.list(pop)){
  popfc<-as.data.frame(sapply(1:popsize,function(x){
    r<-rowMeans(genes[,pop[[1]][,x]])/(rowMeans(genes[,pop[[2]][,x]])+1 )
    r[is.infinite(r)]<-NaN
    return(r)
    }))
  }else{
    popfc<-as.data.frame(apply(pop,1,function(x){
      r<-genes[,x[1]]/(genes[,x[2]]+1)
      r[is.infinite(r)]<-NaN
      return(r)
    }))
  }

  return(popfc)
}

#' Get the estimates of the FC mean from a FC distribution using bayesian inference
#'
#' @param fc list of FC
#' @return matrix with values from \link{do_estimate}
get_estimate<-function(fc){
  e<-apply(fc,1,do_estimate)
  return(e)
}

#' apply bayesian inference to estimate the average FC of a distribution
#' http://www.johnmyleswhite.com/notebook/2010/08/20/using-jags-in-r-with-the-rjags-package/
#' http://public.wsu.edu/~jesse.brunner/classes/bio572/Lab7_Bayesian.html
#' @param x list of values
#' @return vector with mu and its confidence intervales (2.5% and 97.5%)
do_estimate<-function(x){
  x<-(as.numeric(x))
  #print(x[1:10])
  mx<-min(x[!is.infinite(x)],na.rm=T)
  x[is.infinite(x)]<-mx  
  sink("model.bug")
  cat(paste0("
model {
  for (i in 1:N) {
		x[i] ~ dnorm(mu, tau)
	}
	mu ~ dunif(",min(x),",",max(x),")
	tau <- pow(sigma, -2)
	sigma ~ dunif(0, 100)
}
      "),fill=TRUE)
  sink()
  #print(x)
   jags <- jags.model('model.bug',
                     data = list('x' = x,
                                 'N' = length(x)),
                     n.chains = 4,
                     n.adapt = 500)
  
   #update(jags, 10000)
   
   coda <- coda.samples(jags, variable.names = c("mu", "tau"), 
                              n.iter = 500) 
   #summary(coda)
   #plot(coda)
   g<-gelman.diag(coda)

  return(c(summary(coda)$statistics[1:2,1],
              summary(coda)[[2]][1,c(1,5)],g$psrf[1,2]))
   
   
  
}
#' Get rank data frame with best score on the top
#' 
#' @param g1 list of samples in group 1
#' @param g2 list of samples in group 2
#' @param counts count matrix for each gene and each sample
#' @param fc list of FC from the DEG analysis
#' @param popsize number of combinations to generate
#' @return data frame with the output of \link{do_estimate} for each gene
get_rank<-function(g1,g2,counts,fc,popsize){
  popfc<-get_ratio(g1,g2,counts,popsize)
  e.tab<-get_estimate(log2(popfc))
  e.tab.t<-as.data.frame(t(e.tab))
  mix<-cbind(e.tab.t[,c(1,3:4)],fc)  
  mix$sc<-apply(mix[,1:3],1,function(x){
    if (x[2]*x[3] < 0){
      d<-abs(x[2])+abs(x[3])
      s<-d/mean(abs(x[1]))
    }else{
      d<-abs(max(abs(x))-min(abs(x)))
      s<-d/abs(x[1])
    }
    return(s)  
  })
  
  
  return(mix[order(mix$sc),])
  
}
#' plot the correlation between the rank according estimator and the rank according FC
#'
#' @param mix output from get_rank function
#' @param colors colour used for each gene
#' @return ggplot2 object
plotrank<-function(mix,colors=""){
  idsc<-row.names(mix[order(abs(mix$sc)),])
  idfc<-row.names(mix[order(abs(mix[,3]),decreasing=T),])
  sc<-""
  fc<-""
  col=""
  tab.o<-data.frame(row.names=idsc,sc=1:length(idsc),fc=0.0,col="")
  tab.o[idfc,2]<-1:length(idfc)
  if (colors!="" & is.data.frame(colors)){
   tab.o[as.character(colors$genes),3]<-colors$colors
  }else{
    tab.o[,3]<-"black"
  }
  p<-ggplot(tab.o, aes(sc,fc,colour=factor(col))) +
    geom_point()+
    theme_bw()+
    scale_color_brewer(palette="Set1")+
    labs(list(y="rank by FC",x="rank by score"))
  return(p)
}


