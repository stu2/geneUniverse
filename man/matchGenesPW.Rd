% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/matchGenesPW.R
\name{matchGenesPW}
\alias{matchGenesPW}
\title{Get a Matched Control Gene-Set}
\usage{
matchGenesPW(x, goi, metric = "AveExpr", geneIDCol = "Gene.ID",
  ngenes = NULL, maxDiscord = NULL, verbose = F, plotRes = F,
  sweepFrom = "median", highRankedGenes = NULL, returnAll = F)
}
\arguments{
\item{x}{A data-frame having minimally 2 columns: a gene ID column and the column containing the metric that is confounded with gene discovery.}

\item{goi}{A vector containing the list of the genes of interest. These gene-names should be present in the column of \code{x} specified in \code{geneIDCol}.}

\item{metric}{The name of the column of \code{x} containing the metric that is confounded with gene discovery.}

\item{geneIDCol}{The name of the column of \code{x} containing the gene IDs.}

\item{ngenes}{The minimum number of genes that \code{matchGenesPW} will attempt to return.}

\item{maxDiscord}{The maximum Kolmogorov-Smirnov test statistic allowable between the frequency distribution of \code{metric} between the test gene-set and the new control gene-set. If set, \code{ngenes} will be ignored.}

\item{verbose}{Output extra information, such as testing for the correlation between \code{metric} and inclusion in \code{goi}.}

\item{plotRes}{Output plots of the iterations and distributions of \code{metric} between test and control gene-sets, and discarded genes.}

\item{sweepFrom}{Strategy for prioritizing genes in \code{goi} for matching. \code{median}: start from median ranked gene and work outwards. \code{rnd} order is random. \code{dens}: start with the most regions of the explanatory variable richest in discovered genes. \code{leftRight}:}

\item{highRankedGenes}{A vector containing an expanded list of gene-names passing a less stringent significance threshold than \code{goi}, to model the distribution of \code{metric} for discoverable genes. The gene-names should be present in the column of \code{x} specified in \code{geneIDCol}. \code{highRankedGenes} can be used if there are not enough genes in \code{goi} to generate a smooth /reliable distribution of \code{metric}.}

\item{returnAll}{Return a named list with control gene-names, KS-test results, plots etc. and modifed input dataframe.}
}
\value{
A character vector containing gene IDs of a matched control gene-set (default). If \code{returnAll} is set to TRUE, a list containing all the outputs from the algorithm is returned.
}
\description{
If the likelihood of a gene's inclusion in a gene-set is confounded with other 
variables like expression magnitude or variance, then in order to perform 
unbiased gene-set testing, a control gene-set is needed having a comparable 
distribution of the variable that is confounded with gene discoverability. 
\code{matchGenesPW} generates a control gene-set whose distribution of some metric 
approximately matches that of the gene-set of interest. PW stands for 'pair-wise',
because the algorithm chooses the closest available gene for every 
gene-of-interest in each iteration.This algorithm gives the best match in terms of 
the distribution of \code{metric} between the genes-of-interest and the returned 
control gene-set, however local spikes in control-gene density will form around 
isolated genes-of-interest. Also note, because there is no smoothing involved, 
monotonically increasing selection cannot be forced (use \code{matchGenesSmth} if 
monotonicity is desired).
}
\examples{
# Make biased multimodal dataset (first 300 are selected as DE and have higher expression)
set.seed(1)
df<-data.frame(Gene.ID=paste0('gene', 1:10000), 
               AveExpr=c(rgamma(100,shape = 3,scale = 0.5)+2,
                         rgamma(9200,shape = 2,scale = 1),
                         rgamma(700,shape = 2,scale = 0.5)+5), 
               diffReg=c(rep(TRUE,300),rep(FALSE,9700)))
controlSet <- matchGenesPW(df, df$Gene.ID[df$diffReg], metric="AveExpr", 
                        ngenes=5000,verbose = TRUE, plotRes = TRUE)

}
\author{
Stuart Archer
}
