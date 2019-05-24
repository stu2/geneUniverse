#' Get a Matched Control Gene-Set
#' 
#' If the likelihood of a gene's inclusion in a gene-set is confounded with other 
#' variables like expression magnitude or variance, then in order to perform 
#' unbiased gene-set testing, a control gene-set is needed having a comparable 
#' distribution of the variable that is confounded with gene discoverability. 
#' \code{matchGenesSmth} generates a control gene-set whose distribution of some metric 
#' approximately matches a smoothed distribution for this metric in the gene-set of 
#' interest. The smoothed distribution can be forced to be monotonically increasing
#' prior to control gene selection.   
#' 
#' @param x A data-frame having minimally 2 columns: a gene ID column and the column containing the metric that is confounded with gene discovery.
#' @param goi A vector containing the list of the genes of interest. These gene-names should be present in the column of \code{x} specified in \code{geneIDCol}.
#' @param metric The name of the column of \code{x} containing the metric that is confounded with gene discovery.
#' @param geneIDCol The name of the column of \code{x} containing the gene IDs.
#' @param ngenes The target number of genes that \code{matchGenesPW} will attempt to return. Over-rides \code{saturateAt}
#' @param saturateAt (default 0.95). This value of the observed density curve of GOIs will be adjusted to 1.0 in the target density curve for the control-gene set. Lower values return larger control gene lists. This value will be over-ridden if ngenes is specified.
#' @param mode Either 'sweep','probabilistic' or 'returnAll'. Genes are incorporated into the control set by sweeping from left to right on \code{metric}, otherwise each gene is selected purely probabilistically. 
#' @param n If \code{mode} is set to 'probabilistic', setting n>1 allows for several iterations of probabilistic selection of control gene-sets. In this case, a logical matrix is returned. 
#' @param verbose Output extra information, such as testing for the correlation between \code{metric} and inclusion in \code{selGenes}.
#' @param plotRes Output plots of the iterations and distributions of \code{metric} between test and control gene-sets, and discarded genes.
#' @param highRankedGenes A vector containing an expanded list of gene-names passing a less stringent significance threshold than \code{goi}, to model the distribution of \code{metric} for discoverable genes. The gene-names should be present in the column of \code{x} specified in \code{geneIDCol}. \code{highRankedGenes} can be used if there are not enough genes in \code{selGenes} to generate a smooth distribution of \code{metric}. 
#' @param forceMonot The density of selected genes can never decrease as \code{metric} increases.
#' @param flankWin The number of flanking GOIs to the left and right of a GOI over which local density of GOIs is aggregated. Default (10) equates to a window size of 21 GOIs (up to 10 on left + central GOI + up to 10 on right) which may include any number of non-GOIs. This is similar to nominating kernel size for kernel smoothing.
#' @return A character vector containing gene IDs of the generated control gene-set. If \code{mode} is set to 'probabilistic' and \code{n} is set to > 1, a logical matrix is returned, where rownames are GeneIDs and columns are iterations. If \code{mode} is set to 'returnAll', a named list of all outputs (including plots etc) is returned.
#' @examples 
#' # Make biased multimodal dataset (first 300 are selected as DE and have higher expression)
#' set.seed(1)
#' df<-data.frame(Gene.ID=paste0('gene', 1:10000), 
#'                AveExpr=c(rgamma(100,shape = 3,scale = 0.5)+2,
#'                          rgamma(9200,shape = 2,scale = 1),
#'                          rgamma(700,shape = 2,scale = 0.5)+5), 
#'                diffReg=c(rep(TRUE,300),rep(FALSE,9700)))
#' controlSet <- matchGenes(df, df$Gene.ID[df$diffReg], metric="AveExpr", 
#'                         ngenes=5000,verbose = TRUE, plotRes = TRUE)
#' 
#' @author Stuart Archer
#' @importFrom graphics abline legend lines par plot
#' @importFrom stats density ks.test median wilcox.test approxfun
#' @export

matchGenesSmth <- function(x, goi, metric='AveExpr', geneIDCol='Gene.ID', ngenes=NULL, saturateAt=0.95, 
                           mode='sweep', n=1, verbose=F, plotRes=F, highRankedGenes=NULL, forceMonot=F,
                           flankWin=10){
  # brief sanity check on supplied args
  if (mode == 'prob'){mode = 'probabilistic' } # thats hard to type!
  if (!mode %in% c('sweep','probabilistic', 'returnAll')){
    stop('mode must be either "sweep" or "probabilistic" or "returnAll"')
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
  x$controlSet = F
  
  # report on possible bias if verbose
  # maybe entropy would be a better way to do this
  if(verbose){
    cat('----------------------------------------------\n')
    cat(paste0('K-S test for difference in metric distribution between ', selGenesDescr, 
               ' and all other genes:\n'))
    kst <- ks.test( x = x[ x$selGene , metric ], y = x[ !x$selGene, metric ] )
    print(kst)
    cat('----------------------------------------------\n')
  }
  
  x = .makeRankDensMid(x = x, flankWin = flankWin)
  x$TargetDens = pmin(1, x$DensCentredWindow/(max(x$DensCentredWindow)*saturateAt)) # normalize DensCentredWindow
  if(forceMonot){
    for (rw in 1:nrow(x)){
      x$TargetDens[rw] = max(x$TargetDens[1:rw], na.rm = T)
    }
  } 
  
  x$pChoose=x$TargetDens
  x$pChoose[x$goi] = 0 # don't re-select GOIs
  if (!is.null(ngenes)){
    for (i in 1:100){  # 100 iterations is normally plenty to converge
      saturateAt = sum(x$pChoose) / ngenes
      x$pChoose = pmin( 1, x$pChoose / saturateAt )
    }
  }
  x$pChoose = pmin( 1, x$pChoose / saturateAt )
  
  output=list(sweep=NULL, probabilistic=NULL)
  
  if (mode %in% c("sweep","returnAll")){
    x = .selectSweep(x)  # add column controlSetBest, T/F
    output[['sweep']] = x[x$controlSetBest,geneIDCol]
  }
  if (mode %in% c("probabilistic","returnAll")){
    controlSetsProb = matrix(runif(n = n*nrow(x), min = 0, max = 1), ncol = n) <= x$pChoose
    colnames(controlSetsProb) = paste0('controlSet', 1:n)
    rownames(controlSetsProb) = x[,geneIDCol]
    if (n==1){
      output[['probabilistic']] = rownames(controlSetsProb)[controlSetsProb[,1]]
    } else {
      output[['probabilistic']] = controlSetsProb
    }
  }

  if(verbose){
    if (mode %in% c("sweep","returnAll")){
      cat (paste0(sum(x$controlSetBest), ' genes selected for control set.\n'))
      cat (paste0('Final difference between ', selGenesDescr, 
               ' and matched-control set (K-S test):\n'))
      kst <- ks.test(x=x[x$selGene,metric], y=x[x$controlSetBest,metric])
      print(kst)
      output[['KsTest_sweep']] =kst
    }
    if (mode %in% c("probabilistic","returnAll")){
      cat (paste0( '\nNumber of genes in probabilistically matched-control set, iteration 1: ',sum(controlSetsProb[,1]), '\n'))
      cat (paste0(' \nFinal difference between ', selGenesDescr, 
                    ' and probabilistically matched-control set, iteration 1 (K-S test):\n'))
      kst <- ks.test(x=x[x$selGene,metric], y=x[controlSetsProb[,1],metric])
      print(kst)
      output[['KsTest_prob']] = kst
    }
  }
  
  if (plotRes){
    if (requireNamespace("ggplot2", quietly = TRUE)) {
      x$geneSet='Discarded genes'
      if (!is.null(highRankedGenes)){
        x$geneSet[x$selGene] <- 'High-ranked Genes'
      }
      x$geneSet[x$goi]='Genes of interest'
    if (mode %in% c('sweep','returnAll')){
      x$geneSetSweep=x$geneSet
      x$geneSetSweep[x$controlSetBest]='Control gene-set (sweep)'
      x$geneSetSweep <- as.factor(x$geneSetSweep) #'Control gene-set','Discarded genes','Genes of interest','High-ranked Genes')
      x$geneSetSweep = factor(x$geneSetSweep,levels(x$geneSetSweep)[c(3,4,1,2)])
        #requireNamespace("gridExtra", quietly = TRUE)) {
        gr2 = ggplot2::ggplot(x) + 
          ggplot2::geom_freqpoly(ggplot2::aes_string( x='rank' ), bins=30 ) +   #x='AveExpr'
          ggplot2::facet_wrap(facets='geneSetSweep', ncol=1, scales = "free_y") + 
          ggplot2::theme_bw()
        gridExtra::grid.arrange(gr2, ncol=1)
        output[["sweepPlot"]] = gr2
      }
      if (mode %in% c('probabilistic','returnAll')){  # && n==1
        x$geneSetProb=x$geneSet
        x$geneSetProb[controlSetsProb[,1]]='Control gene-set (Probabilistic, iteration 1)'
        x$geneSetProb <- as.factor(x$geneSetProb) #'Control gene-set','Discarded genes','Genes of interest','High-ranked Genes')
        x$geneSetProb = factor(x$geneSetProb,levels(x$geneSetProb)[c(3,4,1,2)])
        #requireNamespace("gridExtra", quietly = TRUE)) {
        gr3 = ggplot2::ggplot(x) + 
          ggplot2::geom_freqpoly(ggplot2::aes_string( x='rank' ), bins=30 ) +   #x='AveExpr'
          ggplot2::facet_wrap(facets='geneSetProb', ncol=1, scales = "free_y") + 
          ggplot2::theme_bw()
        gridExtra::grid.arrange(gr3, ncol=1)
        output[["ProbPlot"]] = gr3
      }
    }
  }
  if (mode == 'returnAll'){
    output[['x']] <- x
    return(output)  # return the whole named list
  } else if (mode == 'sweep'){
    return(output[['sweep']])
  } else if (mode == 'probabilistic'){
    # If n=1, return a geneID vector. 
    # If n>1, return a logical matrix, rownames are GeneIDs, columns are iterations
    return(output[['probabilistic']])  
  }
}  
