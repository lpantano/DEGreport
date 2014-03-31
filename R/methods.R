#' distribution of pvalues by expression range
#'
#' @param pvalues  pvalues of DEG analysis
#' @param counts  matrix with counts for each samples and each gene. Should be same length than pvalues vector.
#' @examples
#' \dontrun{
#' pvalueMean(pvalues,counts)
#' }
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
#' distribution of pvalues by SD range
#'
#' @param pvalues  pvalues of DEG analysis
#' @param counts  matrix with counts for each samples and each gene. Should be same length than pvalues vector.
#' @examples
#' \dontrun{
#' pvalueVar(pvalues,counts)
#' }
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
#' @param pvalues  pvalues of DEG analysis
#' @param counts  matrix with counts for each samples and each gene. Should be same length than pvalues vector.
#' @examples
#' \dontrun{
#' pvalueVarMean(pvalues,counts)
#' }
pvalueVarMean<-function(pvalues,counts){
  
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
#' distribution of DE genes expression compared the background
#'
#' @param tags  list of genes that are DE
#' @param g1 list of samples in group 1
#' @param g1 list of samples in group 2
#' @param counts  matrix with counts for each samples and each gene. Should be same length than pvalues vector.
#' @param pop number of random samples taken for background comparison
#' @examples
#' \dontrun{
#' expDElist(tags,g1,g2,counts)
#' }
expDElist<-function(tags,g1,g2,counts,pop=400){
  delen<-length(tags)
  rand<-sample(row.names(counts),pop)
  g1var<-apply(counts[tags,g1],1,mean)    
  g2var<-apply(counts[tags,g2],1,mean)    
  
  rand.s1<-apply(counts[rand,g1],1,mean)
  rand.s2<-apply(counts[rand,g2],1,mean)
  res<-rbind(data.frame(g="g1",mean=g1var),data.frame(g="g2",mean=g2var))
  res<-rbind(res,data.frame(g="r1",mean=rand.s1),
             data.frame(g="r2",mean=rand.s2))
  
  p<-ggplot(res,aes(mean,fill=g,colour=g))+  
    geom_density(alpha=0.2)+
    xlim(c(0,300))+
    theme_bw()
  return(p)
}
#' distribution of DE genes SD compared the background
#'
#' @param tags  list of genes that are DE
#' @param g1 list of samples in group 1
#' @param g1 list of samples in group 2
#' @param counts  matrix with counts for each samples and each gene. Should be same length than pvalues vector.
#' @param pop number of random samples taken for background comparison
#' @examples
#' \dontrun{
#' varDElist(tags,g1,g2,counts)
#' }
varDElist<-function(tags,g1,g2,counts,pop=400){
  delen<-length(tags)
  rand<-sample(row.names(counts),pop)
  g1var<-apply(counts[tags,g1],1,sd)    
  g2var<-apply(counts[tags,g2],1,sd)    
  
  rand.s1<-apply(counts[rand,g1],1,sd)
  rand.s2<-apply(counts[rand,g2],1,sd)
  res<-rbind(data.frame(g="g1",var=g1var),data.frame(g="g2",var=g2var))
  res<-rbind(res,data.frame(g="r1",var=rand.s1),
             data.frame(g="r2",var=rand.s2))
  
  p<-ggplot(res,aes(var,fill=g,colour=g))+  
    geom_density(alpha=0.2)+
    xlim(c(0,300))+
    theme_bw()
  return(p)
}

#' get number of potential combinations of two vectors
#'
#' @param g1 list of samples in group 1
#' @param g1 list of samples in group 2
combinations<-function(g1,g2){
  return(g1*g2)
}
#' get random combinations of two groups
#'
#' @param g1 list of samples in group 1
#' @param g1 list of samples in group 2
getcomb<-function(g1,g2){
  if (combinations(length(g1),length(g2))<400){
    return(expand.grid(a=g1,b=g2))
  }else{
    #take 10% of each group, 400 times
    g1p=3
    g2p=3
    g1p<-ceiling(0.1*length(g1)+1)
    if(g1p==2){g1p=3}
    g2p<-ceiling(0.1*length(g2)+1)
    if(g2p==2){g2p=3}
    r1<-sapply(1:400,function(x){sample(g1,g1p)})
    r2<-sapply(1:400,function(x){sample(g2,g2p)})
    return(list(r1,r2))
  }
}
#' get the FC for each gene between two groups
#'
#' @param g1 list of samples in group 1
#' @param g1 list of samples in group 2
get_ratio<-function(g1,g2,genes){
  pop<-getcomb(g1,g2)
  if (is.list(pop)){
  popfc<-as.data.frame(sapply(1:400,function(x){
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

  #popfc$mean<-apply(popfc,1,glmMean)
  #popfc$mean<-apply(popfc,1,mean,na.rm=T)
  #popfc$var<-apply(popfc,1,sd,na.rm=T)
  #var<-sort(var)
  return(popfc)
}

#' get the estimates of the FC mean from a FC distribution using bayesian inference
#'
#' @param fc list of FC
get_estimate<-function(fc){
  e<-apply(fc,1,do_estimate)
  return(e)
}

#' apply bayesian inference to estimate the average FC of a distribution
#'
#' @param x list of values
do_estimate<-function(x){
  library('rjags')
  sink("model.bug")
  x<-(as.numeric(x))
  mx<-min(x[!is.infinite(x)],na.rm=T)
  x[is.infinite(x)]<-mx  
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
                     n.adapt = 5000)
  
   #update(jags, 10000)
   
   coda <- coda.samples(jags, variable.names = c("mu", "tau"), 
                              n.iter = 5000) 
   #summary(coda)
   #plot(coda)
   g<-gelman.diag(coda)

  return(c(summary(coda)$statistics[1:2,1],
              summary(coda)[[2]][1,c(1,5)],g$psrf[1,2]))
   
   
  
}

#' get a rank list according to the confidence of the bayesian estimator
#'
#' @param g1 list of samples in group 1
#' @param g2 list of samples in group 2
#' @param counts count matrix for each gene and each sample
#' @param fc list of FC from the DEG analysis
get_rank<-function(g1,g2,counts,fc){
  popfc<-get_ratio(g1,g2,counts)
  e.tab<-get_estimate(log2(popfc))
  e.tab.t<-as.data.frame(t(e.tab))
#   qv<-apply(e.tab.t[,1:2],1,function(x){
#     q<-quantile(rnorm(1000,x[1],x[2]^-0.5),
#                 c(.025,.975))
#     return(q)
#   })
#   e.tab.t$log2min<-qv[1,]
#   e.tab.t$log2max<-qv[2,]
#   
  mix<-cbind(e.tab.t[detags,c(1,3:4)],fc)
  #mix$state<-apply(mix[,1:2],1,function(x){
  
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
  
  
  return(mix)
  
}
#' plot the correlation between the rank according estimator and the rank according FC
#'
#' @param mix output from get_rank function
#' @param colors colour used for each gene
plotrank<-function(mix,colors){
  idsc<-row.names(mix[order(abs(mix$sc)),])
  idfc<-row.names(mix[order(abs(mix[,3]),decreasing=T),])
  
  tab.o<-data.frame(row.names=idsc,sc=1:length(idsc))
  tab.o[idfc,2]<-1:length(idfc)
  tab.o[colors$genes,3]<-colors$colors
  #tab.o[grep("[0-9]",tab.o$V3),3]<-"Aut"
  
  p<-ggplot(tab.o, aes(sc,V2,colour=factor(V3))) +
    geom_point()
  return(p)
}


#######

plotscore<-function(var,genes,g1,g2){
  var1<-var[1:500];g1<-1:3;g2<-4:6
  g1<-genes[names(var1),g1]
  v1<-apply(g1,1,sd)
  m1<-rowMeans(g1);
  g2<-genes[names(var1),g2]
  v2<-apply(g2,1,sd)
  m2<-rowMeans(g2)
  g1up<-m1+var1
  g1down<-m1-var1
  g2up<-m2+var1
  g2down<-m2-var1
  d<-data.frame(p=1:length(m1),m=m1,u=g1up,d=g1down,group=factor("g1",levels=c("g1","g2")))
  d<-rbind(d,data.frame(p=1:length(m1),m=m2,u=g2up,d=g2down,group="g2"))
  qplot(p,m, data=d, colour=factor(group)) +
   geom_ribbon(aes(ymin=log2(d), ymax=log2(u),colour=group))
  
  
}


meanE<-function(m,dat){
  var=sum(abs(dat-m))
  return(var)
}

glmMean<-function(d){  
  
  d<-d[!is.infinite(d)]
  d<-d[!is.na(d)]
  #print(summary(d))
  m=mean(d)
  m.like <- optim(par = m, fn = meanE, dat = d,method="Brent",lower=min(d),upper=max(d))
  return(m.like$par)
}


getBest<-function(i,dat){
  d<-dat[-i]
  m=max(d)
  m.like <- optim(par = m, fn = meanE, dat = d,method="Brent",lower=min(d),upper=max(d))
  return(m.like$par)
}

meanRemoving<-function(i,dat){
  d<-dat[-i]
  return(mean(d))
}

detectOutlier<-function(g,d,g1,g2){
  dat<-unlist(d[g,])
  m1<-mean(dat[g1])
  mValue1<-which.max(abs(m1-(unlist(lapply(1:length(dat[g1]),meanRemoving,dat[g1])))))
  mTrim1<-meanRemoving(mValue1,dat[g1])
  #FC=(m/minValue)
  m2<-mean(dat[g2])
  mValue2<-which.max(abs(m2-(unlist(lapply(1:length(dat[g2]),meanRemoving,dat[g2])))))
  mTrim2<-meanRemoving(mValue2,dat[g2])
  FC.c<-mTrim1/mTrim2
  FC<-m1/m2
  return(c(FC.c,FC))
}


detectVar<-function(g,d,g1,g2){
  dat<-unlist(d[g,])
  v1<-log2(var(dat[g1]))
  v2<-log2(var(dat[g2]))
  return(max(var))
}


minFCq25<-function(g,d,g1,g2){
  dat<-unlist(d[g,])
  groups<-list()
  groups[[1]]<-g1
  groups[[2]]<-g2
  maxvar<-which.max(c(var(dat[g1]),var(dat[g2])))
  minvar<-which.min(c(var(dat[g1]),var(dat[g2])))
  qmaxg<-quantile(dat[groups[[maxvar]]],c(.25,.75))
  mming<-mean(dat[groups[[minvar]]])
  minfc<-min(mming/qmaxg)
  if (minvar==2){
    minfc<-min(qmaxg/mming)
  }

  return(minfc)
}

overlap<-function(g,d,g1,g2){
  dat<-unlist(d[g,])
  q1<-quantile(dat[g1],c(0.25,0.75))
  q2<-quantile(dat[g2],c(0.25,0.75))
  m75<-min(q1[2],q2[2])
  m25<-max(q1[1],q2[1])
  a<-min(c(q1[2]-q1[1],q2[2]-q2[1])) 
  ratio<-(m75-m25)/a
  return(ratio)
}


filterData<-function(exp,gen1,gen2,gene,diff=0.3,over=0.40){
  fc.cor<-as.data.frame(t(as.data.frame((lapply(gene,detectOutlier,exp,gen1,gen2)))))
  row.names(fc.cor)<-gene
  overlap.cor<-unlist(lapply(gene,overlap,exp,gen1,gen2))
  names(overlap.cor)<-gene
  minfc.cor<-unlist(lapply(gene,minFCq25,exp,gen1,gen2))
  names(minfc.cor)<-gene
  var.cor<-unlist(lapply(gene,detectVar,exp,gen1,gen2))
  names(var.cor)<-gene
  tab.res<-cbind(fc.cor,minfc.cor,var.cor,overlap.cor)
  ##tab.res$ratio.cor<-(log2(tab.res[,1])/log2(tab.res[,2]))
  ##tab.res$ratio.cor2<-(log2(tab.res[,2])/log2(tab.res[,4]))
  ##tab.res$score<-apply(tab.res[,c(3,5,6)],1,mean)
  #tab.fil<-tab.res[tab.res[,1]>1.3 & tab.res[,4]>1.3 & tab.res[,2]>1.3 & (tab.res[,2]>0.2 | tab.res[,2]<5) ,]
  #tab.fil<-rbind(tab.fil,tab.res[tab.res[,1]<0.7 & tab.res[,4]<0.7 & tab.res[,2]<0.7 & (tab.res[,2]>0.2 | tab.res[,2]<5),])
  tab.fil<-tab.res[tab.res[,1]>(1+diff) & tab.res[,3]<over & tab.res[,2]>(1+diff) ,]
  tab.fil<-rbind(tab.fil,tab.res[tab.res[,1]<(1-diff) & tab.res[,3]<over & tab.res[,2]<(1-diff) ,])
  #tab.fil<-tab.res
  #tab.fil[order(tab.fil[,3],-tab.fil[,1]),]
  return(tab.fil)
}