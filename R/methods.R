#' Distribution of pvalues by expression range
#' @aliases degMean
#' @usage degMean(pvalues,counts)
#' @param pvalues  pvalues of DEG analysis
#' @param counts  matrix with counts for each samples and each gene.
#' row number should be the same length than pvalues vector.
#' @return ggplot2 object
#' @examples
#' data(DEGreportSet)
#' degMean(DEGreportSet$deg[,4],DEGreportSet$counts)
degMean <-
    function(pvalues,counts)
{
    meanv  <-  apply(counts,1,mean)
    meanvfac <- cut(meanv,
                breaks=c(0,quantile(meanv,seq(.1,1,.1))),
                labels=paste0(seq(10,100,10),"%"),
                right=TRUE)
    pvalfac  <-  cut(pvalues,
               breaks=c(0,seq(.1,1,.1)),
               labels=seq(0.1,1,0.1),right=TRUE)
    d  <-  data.frame(pvalues=factor(pvalfac),
                meanv=factor(meanvfac))
    suppressWarnings(ggplot(d,aes(pvalues,fill=meanv))+
        geom_bar()+
        theme_bw() +
        scale_fill_brewer("mean\nquantiles",palette="RdYlBu")) +
        labs(list(x="p-values",y="# genes"))
}

#' Distribution of pvalues by standard desviation range
#' @aliases degVar
#' @usage degVar(pvalues,counts)
#' @param pvalues  pvalues of DEG analysis
#' @param counts  matrix with counts for each samples and each gene.
#' row number should be the same length than pvalues vector.
#' @return ggplot2 object
#' @examples
#' data(DEGreportSet)
#' degVar(DEGreportSet$deg[,4],DEGreportSet$counts)
degVar <-
    function(pvalues,counts)
{
    sdv <- apply(counts,1,sd)
    sdvfac <- cut(sdv,
                breaks=c(0,quantile(sdv,seq(.1,1,.1))),
                labels=paste0(seq(10,100,10),"%"),
                right=TRUE)
    pvalfac <- cut(pvalues,
               breaks=c(0,seq(.1,1,.1)),
               labels=seq(0.1,1,0.1),right=TRUE)
    d <- data.frame(pvalues=factor(pvalfac),
                sdv=factor(sdvfac))
    suppressWarnings(ggplot(d,aes(pvalues,fill=sdv))+
        geom_bar()+
        theme_bw() +
        scale_fill_brewer("variance\nquantiles",palette="RdYlBu") +
        labs(list(x="p-values",y="# genes")))
}

#' Correlation of the standard desviation and the mean of the abundance of a
#' set of genes.
#' @aliases degMV
#' @usage degMV(g1,g2,pvalues,counts)
#' @param g1 list of samples in group 1
#' @param g2 list of samples in group 2
#' @param pvalues  pvalues of DEG analysis
#' @param counts  matrix with counts for each samples and each gene.
#' row number should be the same length than pvalues vector.
#' @return ggplot2 object
#' @examples
#' data(DEGreportSet)
#' degMV(DEGreportSet$g1,DEGreportSet$g2,DEGreportSet$deg[,4],
#'     DEGreportSet$counts)
degMV <-
    function(g1,g2,pvalues,counts)
{
    sdt1 <- apply(counts[,g1],1,sd)
    sdt2 <- apply(counts[,g2],1,sd)
    sdv <- apply(cbind(sdt1,sdt2),1,max)
    mt1 <- apply(counts[,g1],1,mean)
    mt2 <- apply(counts[,g2],1,mean)
    meanv <- apply(cbind(mt1,mt2),1,max)
    pv <- cut(pvalues,breaks=c(-1,0.01,1.1),
          labels=c("<0.01","NoSig"))
    d <- data.frame(pvalues=pv,sdv=log2(sdv),
                meanv=log2(meanv))
    suppressWarnings(ggplot(d,aes(meanv,sdv,
    colour=pvalues))+
    geom_point()+
    scale_color_manual(values=c("red",rgb(0.9,0.9,0.9,0.6)))+
    theme_bw()+
    stat_quantile(aes(meanv,sdv),colour="blue",
        quantiles = c(0.025,0.975),
        linetype=2,formula=y ~ x))
}

#' Distribution of expression of DE genes compared to the background
#'
#' @aliases degMB
#' @usage degMB(tags,g1,g2,counts,pop=400)
#' @param tags  list of genes that are DE
#' @param g1 list of samples in group 1
#' @param g2 list of samples in group 2
#' @param counts  matrix with counts for each samples and each gene.
#' Should be same length than pvalues vector.
#' @param pop number of random samples taken for background comparison
#' @return ggplot2 object
#' @examples
#' data(DEGreportSet)
#' detag <- row.names(DEGreportSet$deg[1:10,])
#' degMB(detag,DEGreportSet$g1,DEGreportSet$g2,DEGreportSet$counts)
degMB <-
    function(tags,g1,g2,counts,pop=400)
{
    delen <- length(tags)
    g <- ""
    rand <- sample(row.names(counts),pop)
    g1var <- apply(counts[tags,g1,drop=FALSE],1,mean)
    g2var <- apply(counts[tags,g2,drop=FALSE],1,mean)

    rand.s1 <- apply(counts[rand,g1,drop=FALSE],1,mean)
    rand.s2 <- apply(counts[rand,g2,drop=FALSE],1,mean)
    res <- rbind(data.frame(g="g1",mean=g1var),
             data.frame(g="g2",mean=g2var))
    res <- rbind(res,data.frame(g="r1",mean=rand.s1),
             data.frame(g="r2",mean=rand.s2))

    suppressWarnings(ggplot(res,aes(g,mean,
        fill=g,colour=g))+
        geom_violin(alpha=0.2)+
        scale_y_log10()+
        theme_bw())
}

#' Distribution of the standard desviation of
#'     DE genes compared to the background
#' @aliases degVB
#' @usage degVB(tags,g1,g2,counts,pop=400)
#' @param tags  list of genes that are DE
#' @param g1 list of samples in group 1
#' @param g2 list of samples in group 2
#' @param counts  matrix with counts for each samples and each gene.
#'     Should be same length than pvalues vector.
#' @param pop number of random samples taken for background comparison
#' @return ggplot2 object
#' @examples
#' data(DEGreportSet)
#' detag <- row.names(DEGreportSet$deg[1:10,])
#' degVB(detag,DEGreportSet$g1,DEGreportSet$g2,DEGreportSet$counts)
degVB <-
    function(tags,g1,g2,counts,pop=400)
{
    delen <- length(tags)
    g <- ""
    rand <- sample(row.names(counts),pop)
    g1var <- apply(counts[tags,g1,drop=FALSE],1,sd)
    g2var <- apply(counts[tags,g2,drop=FALSE],1,sd)

    rand.s1 <- apply(counts[rand,g1,drop=FALSE],1,sd)
    rand.s2 <- apply(counts[rand,g2,drop=FALSE],1,sd)
    res <- rbind(data.frame(g="g1",var=g1var),
             data.frame(g="g2",var=g2var))
    res <- rbind(res,data.frame(g="r1",var=rand.s1),
             data.frame(g="r2",var=rand.s2))

    suppressWarnings(ggplot(res,aes(g,var,fill=g,colour=g))+
        geom_violin(alpha=0.2)+
        scale_y_log10()+
        theme_bw())
}

#' Get number of potential combinations of two vectors
#' @aliases degNcomb
#' @usage degNcomb(g1,g2)
#' @param g1 list of samples in group 1
#' @param g2 list of samples in group 2
#' @return maximum number of combinations of two vectors
degNcomb <-
    function(g1,g2)
{
    return(g1*g2)
}

#' Get random combinations of two groups
#' @aliases degComb
#' @usage degComb(g1,g2,pop)
#' @param g1 list of samples in group 1
#' @param g2 list of samples in group 2
#' @param pop number of combinations to be return
#' @return matrix with different combinatios of two vector
degComb <-
    function(g1,g2,pop)
{
    if (degNcomb(length(g1),length(g2))<pop){
        return(as.data.frame(expand.grid(a=g1,b=g2)))
    }else{
        #take 10% of each group, 400 times
        g1p=3
        g2p=3
        g1p <- ceiling(0.1*length(g1)+1)
        if(g1p==2){g1p=3}
        g2p <- ceiling(0.1*length(g2)+1)
        if(g2p==2){g2p=3}
        r1 <- sapply(1:pop,function(x){sample(g1,g1p)})
        r2 <- sapply(1:pop,function(x){sample(g2,g2p)})
        return(list(r1,r2))
    }
}

#' get the FC for each gene between two groups
#' @aliases degFC
#' @usage degFC(g1,g2,counts,popsize)
#' @param g1 list of samples in group 1
#' @param g2 list of samples in group 2
#' @param counts count matrix of deregulated genes
#' @param popsize number of combinations to generate
#' @return FC for different combinations of samples in each group for each gene
degFC <-
    function(g1,g2,counts,popsize)
{
    pop <- degComb(g1,g2,popsize)
    if (!is.data.frame(pop)){
        popfc <- as.data.frame(lapply(1:popsize,function(x){
            r <- rowMeans(counts[,pop[[1]][,x],drop=FALSE]+1)/
                (rowMeans(counts[,pop[[2]][,x],drop=FALSE])+1)
            r[is.infinite(r)] <- NaN
            return(r)
        }))
    }else{
        popfc <- as.data.frame(apply(pop,1,function(x){
            r <- (counts[,x[1],drop=FALSE]+1)/(counts[,x[2],drop=FALSE]+1)
            r[is.infinite(r)] <- NaN
            return(r)
        }))
        if (nrow(counts) == 1){
            popfc <- t(popfc)
        }
        row.names(popfc)<-row.names(counts)
    }

    return(split(popfc,row.names(popfc)))
}

#' Get the estimates of the fold change (FC) mean from a FC distribution using
#' bayesian inference
#' @aliases degBI
#' @param fc list of FC
#' @param iter number of iteration in the mcmc model
#' @param ncores number of cores to use
#' @return matrix with values from \link{degBIcmd}
degBI <-
    function(fc,iter=1000,ncores=NULL)
{
    if (is.null(ncores)){
        e <- lapply(fc,degBIcmd,iter)
    }else{
        message("using multiple threads")
        e <- bplapply(fc, degBIcmd, BPPARAM = MulticoreParam(ncores),
                      iter = iter)
    }
    if (file.exists("model.bug")){file.remove("model.bug")}
    return(do.call(rbind.data.frame, e))
}

#' Apply bayesian inference to estimate the average fold change (FC) of
#' a distribution
#' @description code based on
#' http://www.johnmyleswhite.com/notebook/2010/08/20/
#'     using-jags-in-r-with-the-rjags-package/
#' http://public.wsu.edu/~jesse.brunner/classes/bio572/Lab7_Bayesian.html
#' @param x list of values
#' @param iter number of iteration in the mcmc model
#' @return vector with mu and its confidence intervales (2.5% and 97.5%)
degBIcmd <-
    function(x,iter=1000)
{
    x <- (as.numeric(x))
    #print(x[1:10])
    mx <- min(x[!is.infinite(x)],na.rm=TRUE)
    x[is.infinite(x)] <- mx
    if (max(x)!=min(x)){
        sink("model.bug")
        cat(paste0("
        model {
        for (i in 1:N) {
        x[i] ~ dnorm(mu, tau)
        }
        mu ~ dunif(",min(x),",",max(x),")
        tau  <-  pow(sigma, -2)
        sigma ~ dunif(0, 100)
        }
        "),fill=TRUE)
        sink()
        jags  <-  suppressMessages(jags.model('model.bug',
                                            data = list('x' = x,
                                                        'N' = length(x)),
                                            n.chains = 4,
                                            n.adapt = iter,quiet=TRUE ))

        #update(jags, 10000)
        coda  <-  suppressWarnings(coda.samples(jags,
                                              variable.names = c("mu", "tau"),
                                              n.iter = iter,quiet=TRUE ) )
        g <- gelman.diag(coda)

        return(c(summary(coda)$statistics[1:2,1],
                 summary(coda)[[2]][1,c(1,5)],g$psrf[1,2]))

    }else{
        return(c(mean(x),0,min(x),max(x),1))
    }
}

#' Get rank data frame with best score on the top
#' @aliases degRank
#' @param g1 list of samples in group 1
#' @param g2 list of samples in group 2
#' @param counts count matrix for each gene and each sample that is deregulated
#' @param fc list of FC of deregulated genes.
#'     Should be same length than counts \code{row.names}
#' @param popsize number of combinations to generate
#' @param iter number of iteration in the mcmc model
#' @param ncores number of cores to use
#' @return data frame with the output of \link{degBIcmd} for each gene
#' @examples
#' data(DEGreportSet)
#' degRank(DEGreportSet$g1,DEGreportSet$g2,
#'     DEGreportSet$counts[DEGreportSet$detag[1:5],],
#'     DEGreportSet$deg[DEGreportSet$detag[1:5],1],400,500)
degRank <-
    function(g1,g2,counts,fc,popsize=400,iter=1000,ncores=NULL)
{
    if (!is.element("rjags", installed.packages()[,1])){
        message("you need to install jags and rjags.")
        return(TRUE)
    }
    require("rjags")
    popfc <- degFC(g1,g2,counts,popsize)
    e.tab <- degBI(lapply(popfc,log2),iter,ncores)
    names(e.tab) <- c("mu","tau","q5","q95","conv")
    full <- cbind(e.tab[,c(1,3:4)],fc,row.names = row.names(counts))
    full$sc <- apply(full[,1:3],1,function(x){
        if (x[2]*x[3] < 0){
            d <- abs(x[2])+abs(x[3])
            s <- d/abs(x[1])
        }else{
            d <- abs(max(abs(x))-min(abs(x)))
            s <- d/abs(x[1])
        }
        return(s)
    })
    return(full[order(full$sc),])
}

#' plot the correlation between the rank according estimator and
#'     the rank according FC
#' @aliases degPR
#' @usage degPR(rank,colors)
#' @param rank output from \code{degRank} function
#' @param colors colour used for each gene
#' @return ggplot2 object
#' @examples
#' data(DEGreportSet)
#' degPR(DEGreportSet$rank)
degPR <-
    function(rank,colors="")
{
    idsc <- row.names(rank[order(abs(rank$sc)),])
    idfc <- row.names(rank[order(abs(rank[,3]),
                               decreasing=TRUE),])
    sc <- ""
    fc <- ""
    col=""
    tab.o <- data.frame(row.names=idsc,sc=1:length(idsc),
                      fc=0.0,col="",stringsAsFactors = FALSE)
    tab.o[idfc,2] <- 1:length(idfc)
    if (is.data.frame(colors)){
        tab.o[as.character(colors$genes),3] <- as.character(colors$colors)
    }else{
        tab.o[,3] <- "black"
    }
    suppressWarnings(ggplot(tab.o,
        aes(sc,fc,colour=factor(col))) +
        geom_point()+
        theme_bw()+
        scale_color_brewer(palette="Set1")+
        labs(list(y="rank by FC",x="rank by score")))
}

#' Create a deg object that can be used to plot expression values
#'     at shiny server:runGist(9930881)
#' @aliases degObj
#' @usage degObj(counts,design,outfile)
#' @param counts output from get_rank function
#' @param design colour used for each gene
#' @param outfile file that will contain the object
#' @return R object to be load into vizExp
degObj <-
    function(counts,design,outfile)
{
    deg <- NULL
    deg <- list(counts, design)
    save(deg, file=outfile)
    return(TRUE)
}

