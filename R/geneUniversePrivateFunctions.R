#' @importFrom stats approxfun
#' 
########################################################
#           private functions for all:
########################################################

.vetGeneNamesArgument <- function(gList, allGeneNames, listNameForWarnings='genes-of-interest'){
  lnSG=length(gList)
  gList=unique(gList)
  if( length(gList) < lnSG ){
    warning(paste0('Some ',listNameForWarnings,' are not unique. Purging duplicates.'))
    lnSG=length(gList)
  }
  gList=gList[gList %in% allGeneNames]
  if(length(gList)==0){
    stop(paste0('None of the ',listNameForWarnings,' were found in gene ID column of x'))
  } else if (length(gList) != lnSG){
    warning(paste0('Some of the ',listNameForWarnings,' were not found in gene ID column of x'))
  }
  return(gList)
}

###

.checkAndReformatX <- function(x, metric, geneIDCol, sgnf=NULL){
  if (!any(metric == colnames(x))){
    stop(paste0("supplied metric ('",metric,"') is not a column name of x"))
  } 
  if (!any(geneIDCol == colnames(x))){
    stop(paste0("'", geneIDCol, "' is not a column name of x"))
  }
  if (!is.null(sgnf)){
    if (!any(sgnf == colnames(x))){
      stop(paste0("'", sgnf, "' is not a column name of x"))
    } else {
      x$sgnf = x[,sgnf]
    }
  }
  reservedColnames=c("metric", "rank", "goi", "selGene", "controlSet", "DensCentredWindow", "TargetDens", "pChoose", "controlSetBest", "geneSet", "geneSetSweep", "geneSetProb")
  if (any (colnames(x) %in% reservedColnames)){
    warning(paste0("colnames of x will be overwritten: ", colnames(x)[colnames(x) %in% reservedColnames]), immediate. = T)
  }
  x$metric = x[,metric]
  x$rank <- rank(x$metric)
  x=x[order(x$metric),]
  #x=x[,colnames(x) %in% c(geneIDCol, 'metric', 'rank')]
  return(x)
}


########################################################
#   private functions for matchGenes():
########################################################


.forceMonotonic_OLD <- function(x, verbose=T){ # minCount is the sliding window size (number of genes)
  mxd = 0
  lastMaxima = 0
  runningScore = 0
  x$DensMonoT = 0
  for (rw in 1:nrow(x)){
    if (!is.na(x$TargetDens[rw]) && mxd < x$TargetDens[rw]) {
      mxd = x$TargetDens[rw]
      runningScore=0  # reset interval
      lastMaxima=rw-1 # reset interval
    }
    #if (x$controlSet[rw]){
    #  runningScore=runningScore+1
    #} else if (runningScore/(rw-lastMaxima) < mxd && ! x$goi[rw] ) { 
    if (runningScore/(rw-lastMaxima) < mxd && ! x$goi[rw] ) { 
      runningScore=runningScore+1
      x$controlSet[rw]=T
    } #else {
      #x$controlSet[rw]=F # reset previously chosen controlSet genes to avoid density going too high
    #}
  }
  retList <- list(x=x)
  return(retList)
}

###

# Gene density of GOIs in a rank-window of the up to 10 GOIs left, and 10 GOIs right of, the current GOI.
# flankWin*2 + 1 is the number of GOIs in the window. The proportion of intervening genes that are GOIs/(GOIs + non-GOIs) is returned

###

.makeRankDensMid <- function(x, flankWin=10){ # x must have colnames 'goi' and 'rank'
  conRanks = x$rank[ x$goi]
  if(length(conRanks) < (1+flankWin*2)){
    stop('too few control genes to estimate density for forceMonotonic')
  }
  rankDensities=NA
  for (rw in 1:length(conRanks)){
    aggRange=c(max(1, rw-flankWin), min(length(conRanks), rw+flankWin))
    rankDensities[rw] = (aggRange[2]-aggRange[1]) / (conRanks[aggRange[2]]-conRanks[aggRange[1]])
  }
  x$DensCentredWindow = NA
  x$DensCentredWindow[ x$goi ] = rankDensities # only populate DensLeft20 column for rows corresponding to controlSet genes
  x$DensCentredWindow[ !x$goi ] = approx(x = x$rank[x$goi], y=x$DensCentredWindow[x$goi], xout = x$rank[ !x$goi ], 
                                         rule=1:2, method = 'linear' )[[2]]
  x$DensCentredWindow [is.na ( x$DensCentredWindow ) ] = 0
  return(x)
}

.selectSweep <- function(x){
  runningExp=0
  runningSelected=0
  x$controlSetBest=F
  for (rw in 1:nrow(x)){
    if (!x$goi[rw]){
      runningExp = runningExp + x$pChoose[rw]
      if (runningExp > 1){
        x$controlSetBest[rw]=T
        runningExp=runningExp-1
      }
    }
  }
  return(x)
}

########################################################

# select a FDR cut-off for modelling the expr-ranks for 'nearly' significant genes
# by KS-testing of expr-rank distributions of progressive strata of FDR thresholds
# look at both cumulative (FDR < bin upper limit ) and local (bin limit < FDR < bin upper limit) 
# by KS test against expr-rank distribution of GOIs.

# For now this will be a private function but could become a main function

.selectSgnfByKS <- function(x, goiSgnf, metric='metric', sgnf = 'adj.P.Val', geneIDCol='Gene.ID', nbs=50, 
                            breakOn=0.005, binSize=NULL, doLocal=T, upperQ=0.95, lowerQ=0.05, pthresh=0.5){
  x2 = .checkAndReformatX(x, metric, geneIDCol, sgnf=sgnf)
  if (!any(x2$sgnf <= goiSgnf)){
    stop(paste0("No genes are significant (threshold: ", goiSgnf, ")"))
  }
  x2 = x2[order(x2$sgnf),]  # rank by FDR not metric
  x2$goi = (x2$sgnf <= goiSgnf)
  goiRanks = x2$rank[x2$goi] # 'rank' is the rank by metric (e.g. by Ave.Expr)
  if (length(goiRanks)<50){
    warning("Warning: fewer than 50 significant genes available for extrapolating distribution over Metric")
  }
  if (nbs>0){
    goiRanksBS=list()
    for (bs in 1:nbs){
      goiRanksBS[[bs]]=sample(goiRanks, replace = T)
    }
  }
  
  x2=x2[!x2$goi,] # remove gois from x2
  
  if (is.null(binSize)){ #assume binSize=length(goi)
    binSize=length(goiRanks)
  } 
  nbin=floor(nrow(x2)/binSize)
  ksRes <- list()
  # prepare tmp bootstrapping dataframe in advance to speed up
  if (nbs>0){
    bsKStmp=data.frame(D=rep(Inf,nbs),P=rep(-1,nbs),DCum=rep(Inf,nbs),PCum=rep(-1,nbs))
  }
  for (bn in 1:nbin){
    testRanks = x2$rank[ ((bn-1)*binSize+1) : (bn*binSize)] # KS test of just this FDR stratum's expr-ranks vs GOI expr-ranks
    ksLoc <- ks.test(goiRanks, testRanks)
    
    testRanksCum = x2$rank[ 1 : (bn*binSize)] # KS test of expr-ranks of all FDR strata up to this bin (bn) vs GOI expr-ranks
    ksCum <- ks.test(goiRanks, testRanksCum)
    if(nbs>0){
      for (bs in 1:nbs){ # bootstrap nbs times. bstrapped data can have ties due to replacement
        # local:
        if(doLocal){
          ks <- ks.test(goiRanksBS[[bs]], sample(testRanks,replace = T))
          bsKStmp$D[bs]=ks$statistic
          bsKStmp$P[bs]=ks$p.value
        }
        # cumulative:
        ks <- ks.test(goiRanksBS[[bs]], sample(testRanksCum,replace = T))
        bsKStmp$DCum[bs]=ks$statistic
        bsKStmp$PCum[bs]=ks$p.value
      }
      if (! doLocal){
        bsKStmp$D=NA
        bsKStmp$P=NA
      }
      ksRes[[bn]] = list(RankBinEnd=bn*binSize,
                         p=ksLoc$p.value, 
                         D=ksLoc$statistic, 
                         pCum = ksCum$p.value, 
                         DCum = ksCum$statistic,
                         pBS0.95 = quantile(bsKStmp$P, upperQ, na.rm=T),
                         pBS0.05 = quantile(bsKStmp$P, lowerQ, na.rm=T),
                         DBS0.95 = quantile(bsKStmp$D, upperQ, na.rm=T),
                         DBS0.05 = quantile(bsKStmp$D, lowerQ, na.rm=T),
                         pCumBS0.95 = quantile(bsKStmp$PCum, upperQ),
                         pCumBS0.05 = quantile(bsKStmp$PCum, lowerQ),
                         DCumBS0.95 = quantile(bsKStmp$DCum, upperQ),
                         DCumBS0.05 = quantile(bsKStmp$DCum, lowerQ))
      if(ksRes[[bn]]$pCumBS0.95 < breakOn){
        break
      }
    } else { # don't bootstrap
      ksRes[[bn]] = list(RankBinEnd=bn*binSize,
                         p=ksLoc$p.value, 
                         D=ksLoc$statistic, 
                         pCum = ksCum$p.value, 
                         DCum = ksCum$statistic)
      if(ksRes[[bn]]$pCum < breakOn){
        break
      }
    }
    
  }
  ksRes <- as.data.frame(do.call(rbind, ksRes))
  for (cn in 1:ncol(ksRes)){
    ksRes[,cn] <- unlist(ksRes[,cn])
  }
  gr = ggplot(ksRes, aes(x=RankBinEnd)) + 
    geom_line(aes(y=p),col='red', linetype='dashed') + 
    geom_line(aes(y=pCum)) + 
    theme_bw()
  if (nbs>0){
    gr=gr+geom_line(aes(y=pCumBS0.95),col='blue') + geom_line(aes(y=pCumBS0.05),col='blue') 
  }
  if (any(ksRes$pCum > pthresh)){
    lastBin=max(which(ksRes$pCum > pthresh)) # get the last FDR stratum with cumulative KS p-val > thresh (0.5 usually)
    apf=approxfun(y = ksRes$RankBinEnd[0:1+lastBin], x=ksRes$pCum[0:1+lastBin]) 
    RankThresh = apf(pthresh) # more precise cutoff by interpolating between last OK FDR stratum and first not-OK stratum
    sgnfThreshCO=x2$sgnf[floor(RankThresh)]
  } else { 
    sgnfThreshCO = 0
    warning('No nearly-significant gene bins showed similar distribution of "metric" as significant genes by KS test.')
  }
  return(list(ksRes=ksRes,sgnfThreshCO=sgnfThreshCO, gr=gr))
}

########################################################