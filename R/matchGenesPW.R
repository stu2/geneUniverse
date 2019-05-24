#' Get a Matched Control Gene-Set
#' 
#' If the likelihood of a gene's inclusion in a gene-set is confounded with other 
#' variables like expression magnitude or variance, then in order to perform 
#' unbiased gene-set testing, a control gene-set is needed having a comparable 
#' distribution of the variable that is confounded with gene discoverability. 
#' \code{matchGenesPW} generates a control gene-set whose distribution of some metric 
#' approximately matches that of the gene-set of interest. PW stands for 'pair-wise',
#' because the algorithm chooses the closest available gene for every 
#' gene-of-interest in each iteration.This algorithm gives the best match in terms of 
#' the distribution of \code{metric} between the genes-of-interest and the returned 
#' control gene-set, however local spikes in control-gene density will form around 
#' isolated genes-of-interest. Also note, because there is no smoothing involved, 
#' monotonically increasing selection cannot be forced (use \code{matchGenesSmth} if 
#' monotonicity is desired).
#' 
#' @param x A data-frame having minimally 2 columns: a gene ID column and the column containing the metric that is confounded with gene discovery.
#' @param goi A vector containing the list of the genes of interest. These gene-names should be present in the column of \code{x} specified in \code{geneIDCol}.
#' @param metric The name of the column of \code{x} containing the metric that is confounded with gene discovery.
#' @param geneIDCol The name of the column of \code{x} containing the gene IDs.
#' @param ngenes The minimum number of genes that \code{matchGenesPW} will attempt to return.
#' @param maxDiscord The maximum Kolmogorov-Smirnov test statistic allowable between the frequency distribution of \code{metric} between the test gene-set and the new control gene-set. If set, \code{ngenes} will be ignored.
#' @param verbose Output extra information, such as testing for the correlation between \code{metric} and inclusion in \code{goi}.
#' @param plotRes Output plots of the iterations and distributions of \code{metric} between test and control gene-sets, and discarded genes.
#' @param highRankedGenes A vector containing an expanded list of gene-names passing a less stringent significance threshold than \code{goi}, to model the distribution of \code{metric} for discoverable genes. The gene-names should be present in the column of \code{x} specified in \code{geneIDCol}. \code{highRankedGenes} can be used if there are not enough genes in \code{goi} to generate a smooth /reliable distribution of \code{metric}. 
#' @param sweepFrom Strategy for prioritizing genes in \code{goi} for matching. \code{median}: start from median ranked gene and work outwards. \code{rnd} order is random. \code{dens}: start with the most regions of the explanatory variable richest in discovered genes. \code{leftRight}:
#' @param returnAll Return a named list with control gene-names, KS-test results, plots etc. and modifed input dataframe.
#' @return A character vector containing gene IDs of a matched control gene-set (default). If \code{returnAll} is set to TRUE, a list containing all the outputs from the algorithm is returned.
#' @examples 
#' # Make biased multimodal dataset (first 300 are selected as DE and have higher expression)
#' set.seed(1)
#' df<-data.frame(Gene.ID=paste0('gene', 1:10000), 
#'                AveExpr=c(rgamma(100,shape = 3,scale = 0.5)+2,
#'                          rgamma(9200,shape = 2,scale = 1),
#'                          rgamma(700,shape = 2,scale = 0.5)+5), 
#'                diffReg=c(rep(TRUE,300),rep(FALSE,9700)))
#' controlSet <- matchGenesPW(df, df$Gene.ID[df$diffReg], metric="AveExpr", 
#'                         ngenes=5000,verbose = TRUE, plotRes = TRUE)
#' 
#' @author Stuart Archer
#' @importFrom graphics abline legend lines par plot
#' @importFrom stats density ks.test median wilcox.test approxfun
#' @export

matchGenesPW <- function(x, goi, metric='AveExpr', geneIDCol='Gene.ID',
                       ngenes=NULL, maxDiscord=NULL, verbose=F, plotRes=F, 
                       sweepFrom='median', highRankedGenes=NULL, returnAll=F){
  # brief sanity check on supplied args
  if (!sweepFrom %in% c('median', 'rnd', 'dens','leftRight', 'fullrnd')){
    stop("sweepFrom arg must be one of median, dens, rnd, leftRight")
  }
  x = .checkAndReformatX(x = x, metric = metric, geneIDCol = geneIDCol)
  goi = .vetGeneNamesArgument(gList = goi, allGeneNames = x[,geneIDCol], 
                            listNameForWarnings = 'genes-of-interest')
  if (is.null(highRankedGenes)){
    selGenes=goi
    selGenesDescr='genes-of-interest'
  } else {
    highRankedGenes=.vetGeneNamesArgument(gList = highRankedGenes, allGeneNames = x[,geneIDCol], 
                                         listNameForWarnings = 'high-ranked-genes')
    selGenes=highRankedGenes
    selGenesDescr='high-ranked genes'
  }
  x$goi     = x[,geneIDCol] %in% goi
  x$selGene = x[,geneIDCol] %in% selGenes

  # decide on runmode: reach nGenes first or MaxDiscord first?
  runMode='maxDiscord'
  numberOfSweeps=0
  if (!is.null(ngenes)){
    if(!is.null(maxDiscord)){
      runMode='maxDiscord'
      cat('Warning: both ngenes and maxDiscord were specified. Ignoring ngenes.\n')
    } else {
      runMode='ngenes'
      numberOfSweeps=ceiling(ngenes/length(selGenes)) # change this to fractional algorithm later
    }
  }
  
  if (numberOfSweeps*length(selGenes)>sum(!x$goi)){
    stop('Number of genes to be selected is more than the number available.')
  }
  # x$iter will be the order in which we consider selGenes, 
  # start from median-ranked and work outwards
  x$iter <- NA
  if (sweepFrom == 'rnd'){
    x$iter[x$selGene] = sample(1:sum(x$selGene))
  } else if (sweepFrom == 'median'){
    x$iter[x$selGene] <- rank(abs(x$rank[x$selGene] - median(x$rank[x$selGene]) ), 
                              ties.method = 'random')
  } else if (sweepFrom == 'dens'){
    stop('dens sweep not yet supported')
  } else if (sweepFrom == 'leftRight'){
    x$iter[x$selGene] <- rank(x$rank[x$selGene])  
  }
  x$controlSet = F
  
  # report on possible bias if verbose
  # maybe entropy would be a better way to do this
  if(verbose){
    cat('----------------------------------------------\n')
    cat(paste0('K-S test for difference in metric distribution between ', selGenesDescr, 
               ' and all other genes:\n'))
    kst <- ks.test(x=x[x$selGene,metric], y=x[!x$selGene,metric])
    print(kst)
    cat('----------------------------------------------\n')
  }
  
  runningDiff=0
  kstStats=list()
  if (runMode=='maxDiscord'){
    numberOfSweeps=nrow(x) # this number of iterations will never be reached
  }
  for (j in 1:numberOfSweeps){
    if (sweepFrom == 'leftRight'){
      if (j %% 2 == 1){
        x$iter[x$selGene] <- rank(x$rank[x$selGene]) 
      } else {
        x$iter[x$selGene] <- rank(0-x$rank[x$selGene]) # reverse order every 2nd sweep
      }
    } else if (sweepFrom == 'fullrnd'){
      x$iter[x$selGene] = sample(1:sum(x$selGene))
    }
    controlSet_thisSweep=rep(F,nrow(x))
    if (max(x$iter, na.rm = T) + sum(x$controlSet, na.rm=T) > nrow(x)){
      print(paste0("Maxed out the number of genes at ",j,"th iteration. Exiting."))
      break
    }
    for (i in 1:max(x$iter, na.rm = T)){
      # get one gene from selGenes (at iter==i) for which to get matching control
      rnk=x$rank[which(x$iter==i)]  
      score=x$rank-rnk+runningDiff
      score[x$controlSet]=Inf
      score[x$goi]=Inf
      score[controlSet_thisSweep]=Inf
      selected = which(abs(score)==min(abs(score)))[1]
      controlSet_thisSweep[selected]=T
      runningDiff = runningDiff - (rnk-x$rank[selected])
    }
    # test discordance
    kst <- ks.test(x=x[x$selGene,metric], y=x[(x$controlSet | controlSet_thisSweep),metric])
    kstStats[[j]]=list(statistic=kst$statistic, p.value=kst$p.value)
    if (runMode=='maxDiscord' && kst$statistic > maxDiscord){
      if (verbose){
        cat(paste0('Exceeded max K-S test discordance at iteration ',j, '.\n'))
      }
      break
    } else {
      x$controlSet=x$controlSet | controlSet_thisSweep
    }
  }
  if(verbose){
    cat (paste0(sum(x$controlSet), ' genes selected for control set.\n'))
    cat (paste0('Final difference between ', selGenesDescr, 
               ' and matched-control set (K-S test):\n'))
    ksT <- ks.test(x=x[x$selGene,metric], y=x[x$controlSet,metric])
    print(ksT)
  }
  
  if (plotRes || returnAll){
    ksDF  <-  as.data.frame(t(matrix(unlist(kstStats), nrow=2)))
    colnames(ksDF)=c('KSstat','KS-pval')
    ksDF$iter=1:nrow(ksDF)
    x$geneSet='Discarded genes'
    x$geneSet[x$controlSet]='Control gene-set'
    if (!is.null(highRankedGenes)){
      x$geneSet[x$selGene] <- 'High-ranked Genes'
    }
    x$geneSet[x$goi]='Genes of interest'
    x$geneSet <- as.factor(x$geneSet) #'Control gene-set','Discarded genes','Genes of interest','High-ranked Genes')
    x$geneSet = factor(x$geneSet,levels(x$geneSet)[c(3,4,1,2)])
    if (requireNamespace("ggplot2", quietly = TRUE) && 
        requireNamespace("gridExtra", quietly = TRUE)) {

      gr1 <- ggplot2::ggplot(ksDF) + 
        ggplot2::geom_point(ggplot2::aes_string(x='iter',y='KSstat')) +
        ggplot2::theme_bw()
      if (runMode=='maxDiscord'){
        gr1=gr1+ggplot2::geom_hline(yintercept=maxDiscord, linetype='dashed', col='red')
      }
      gr2 = ggplot2::ggplot(x) + 
        ggplot2::geom_freqpoly(ggplot2::aes_string( x='rank' ), bins=30 ) +   #x='AveExpr'
        ggplot2::facet_wrap(facets='geneSet', ncol=1, scales = "free_y") + 
        ggplot2::theme_bw()
        gridExtra::grid.arrange(gr1,gr2, ncol=2)
    } else {
      denList=list()
      mxdens=0
      mxx=0
      for (gset in unique(x$geneSet)){
        denList[[gset]] <- density(x$metric[x$geneSet==gset])
        mxdens=max(mxdens, max(denList[[gset]]$y))
        mxx=max(mxx,max(denList[[gset]]$x))
      }
      par(mfrow=c(1,2),mar=c(3.5,3.5,0,0)+1) #bottom, left, top and right margins
      plot(denList[['Selected gene-set']], col='red', 
           ylim=c(0,mxdens), xlim=c(0,mxx), main="", sub="", xlab=metric, ylab='gene density')
      lines(denList[['Control gene-set']], col='blue')
      lines(denList[['Discarded genes']], col='green')
      legend(mxx*0.5, mxdens, legend=c('selGenes', 'Control-set', 'Discarded'),
             col=c('red', 'blue', 'green'), lty=1, cex=0.8)
      plot(ksDF$iter,y=ksDF$KSstat, xlab='iteration', ylab='K-S test statistic')
      if(!is.null(maxDiscord)){
          abline(h = maxDiscord, col='red',lty=2)
      }
    }
  }
  if (returnAll){
    return(list(controlGenes=x[x$controlSet,geneIDCol], x=x, ksDF=ksDF))
  } else {
    return(x[x$controlSet,geneIDCol])
  }
}