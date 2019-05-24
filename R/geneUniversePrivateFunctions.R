#' @importFrom MASS kde2d
#' @importFrom stats approxfun
#' @importFrom pracma interp2
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


.use_2dSmoothing <- function(x, nbinxIn=20, nbinyIn=20, nbinyOut=500, nbinxOut=250, sgnfThresh=0.05,
                             yMarg=NULL, xMarg=NULL, kernCdfAt1xbw=0.05){
  bwAdjust = .kde2dKernelBWadjust(kernCdfAt1xbw = kernCdfAt1xbw)
  # Use 2d smoothing: AveExpr (x)  vs FDR (y)
  if (is.null(yMarg)){
    yMarg =  2 / nbinyIn
  }
  if (is.null(xMarg)){
    xMarg =  2 / nbinxIn
  }
  xFactor = 1 + 2 * xMarg
  yFactor = 1 + 2 * yMarg
  nbinxOutReal = as.integer(nbinxOut * xFactor)
  nbinyOutReal = as.integer(nbinyOut * yFactor)
  
  smth2d <- MASS::kde2d(x = x$rank, 
                  y = x$adj.P.Val, 
                  h=c(nrow(x)/nbinxIn, 1/nbinyIn)/bwAdjust,   # adjusting by 0.41079873 makes cdf = 5% at 1 bw from point (remember though, cdf only = 50% at 0 bw from point)
                  n=c(nbinxOutReal, nbinyOutReal),
                  lims = c(c(0-xMarg, 1+xMarg) * nrow(x),
                           c(0-yMarg, 1+yMarg) * 1) )
  smth2d$zYcum = t(apply(smth2d$z, MARGIN = 1, FUN = cumsum))  # convert to cumulative sum in y-axis (FDR axis)
  #lims = c(1,nrow(x),0,1))
  # note:  kernel(bw*0.7530855) = 1%;   kernel(bw*0.6064393)=5%;     kernel(bw*0.5310085)=10%;     kernel(bw*0.2891057)=50%
  # areas: cdfKernel(bw*-0.5810079)=1%; cdfKernel(bw*-0.4107987)=5%; cdfKernel(bw*-0.3200473)=10%; cdfKernel(0) = 50%
  # test this using bw=100:
  # tst <- kde2d(x=500, y=500, h=c(100,100), lims = c(0,1000,0,1000), n = c(1000,1000))
  # ap3 = approxfun(y=(200:800-500)/100, x=cumsum(tst$z[500,200:800]/sum(tst$z[500,])))
  # ap3(c(0.01,0.05,0.1,0.5))
  print('finished 2d kernel smoothing..')
  
  xp=rep((0:nbinxOutReal-xMarg*nbinxOut) * nrow(x)/nbinxOut ,   nbinyOutReal+1)
  yp=rep((0:nbinyOutReal-yMarg*nbinyOut) * 1 / nbinyOut ,  each=nbinxOutReal+1)
  #NEEDS PRACMA
  smth2dAF <- pracma::interp2(x=smth2d$x, y=smth2d$y, Z = t(smth2d$z), xp=xp, yp=yp)
  smth2dAF_df <- data.frame(x=xp,y=yp,z=smth2dAF)
  gr_tile = ggplot(smth2dAF_df) + geom_tile(aes(x=x,y=y,fill=z)) + 
    scale_fill_gradientn(colours = rev(rainbow(4))) + theme_bw() + 
    geom_hline(yintercept = sgnfThresh, linetype='dashed', col='black')
  
  print('finished data-frame of 2d smooth')
  
  vol_underThresh=data.frame(rnkbn=(0:nbinxOut) * nrow(x)/nbinxOut, underThresh=NA)
  vol_underThresh$atThresh = pracma::interp2(x=smth2d$x, y=smth2d$y, Z = t(smth2d$z), # probably not a useful metric
                                     xp=vol_underThresh$rnkbn, 
                                     yp=rep(sgnfThresh,nbinxOut+1))
  vol_underThresh$underThresh = pracma::interp2(x=smth2d$x, y=smth2d$y, Z = t(smth2d$zYcum), 
                                        xp=vol_underThresh$rnkbn, 
                                        yp=rep(sgnfThresh,nbinxOut+1))
  print('finished y-axis cumsum-based estimation')
  
  ### following old (slow) process generates very, very similar density shape; scaled differently
  # for (rw in 1:nrow(vol_underThresh)){
  #   xp = rep(vol_underThresh$rnkbn[rw],nbinyOutReal+1)  
  #   yp = (0:nbinyOutReal/nbinyOutReal) * yFactor - yMarg
  #   smth2d_Interp2 <- pracma::interp2(x=smth2d$x, y=smth2d$y, Z = t(smth2d$z), xp=xp, yp=yp)
  #   ap=approxfun(x = yp, y = smth2d_Interp2)  # x and y are now FDR and density (FDR was y before now)
  #   vol_underThresh$underThresh[rw] = integrate(ap, lower = -yMarg, upper = sgnfThresh)$value
  # }
  # print('finished integrate-based analysis')
  
  gr = ggplot(vol_underThresh) + 
    geom_line(aes(x=rnkbn, y=underThresh)) + 
    theme_bw()
  #pltly <- plot_ly(smth2dAF_df, x = ~x, y = ~y, z = ~z, type = 'scatter3d', mode = 'lines',
  #                 color = ~(as.integer(y<0.05)-1)*(-z-0.00002))
  #pltly2 <- plot_ly(smth2dAF_df, x = ~x, y = ~y, z = ~z, type = 'scatter3d', mode = "markers",
  #                  color = ~(as.integer(y<sgnfThresh)-1)*(-z-0.00002), marker = list(size=3))
  vol_underThresh$underThresh = vol_underThresh$underThresh/max(vol_underThresh$underThresh)
  vol_underThresh$atThresh = vol_underThresh$atThresh/max(vol_underThresh$atThresh)
  return(list(df = vol_underThresh, gr=gr, smth2d=smth2d, gr_tile=gr_tile, smth2dAF_df=smth2dAF_df)) #, plt3d = pltly2))
}


# get empirical bw vs cdf relationship from kde2d
.kde2dKernelBWadjust <- function(kernCdfAt1xbw){ 
  # make standard from a single point, smoothed with bw of 100
  kw <- MASS::kde2d(x=500, y=500, h=c(100,100), lims = c(0,1000,0,1000), n = c(1000,1000))
  # just get central row of matrix (cross-section)
  xSection = kw$z[500,] 
  # make approxfun to yield number of bandwidths (y) given a certain CDF value (x)
  ap = approxfun(y=(200:800-500)/100, x=cumsum(xSection[200:800]/sum(xSection)))
  return(abs(ap(kernCdfAt1xbw))) # divide the h arg of kde2d by this returned factor, and the resulting cdf of kernel will equal 'kernCdfAt1xbw' at position in 'h' 
}