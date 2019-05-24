\name{geneUniverse-package}
\alias{geneUniverse-package}
\alias{geneUniverse}
\docType{package}
\title{
\packageTitle{geneUniverse}
}
\description{
\packageDescription{geneUniverse}
}
\details{
If the likelihood of a gene's inclusion in a gene-set is confounded with other variables like expression magnitude or variance, then in order to perform unbiased gene-set testing, a control gene-set is needed having a comparable distribution of the variable that is confounded with gene discoverability. 

\code{geneUniverse} provides functions \code{matchGenesPW()} and \code{matchGenesSmth()}, which can select a subset of control genes whose distribution of some variable matches that of the genes-of-interest. The selection of control genes is non-parametric (rank-based). If the number of significant genes is small, an expanded set can also be given as an argument. For example, if the test-set has only about 10 genes with FDR < 0.01 (not enough to model a reliable distribution of \code{metric}), but there are hundreds with FDR < 0.05, the latter set can be used for generating the matched control genes using the \code{highRankedGenes} argument. 

Note that if both highly- and lowly- expressed genes are under-represented in the significant gene-set, there should be good reason to believe that this is due to false-negatives in both cases if using \code{matchGenesPW()} to generate a control set. If the under-representation of highly-expressed genes is rather due to there being genuinely less expression change in these genes, then removal of these from the matched-control set may not be appropriate, and a monotonically increasing model for selection bias (e.g. as produced by \code{matchGenesSmth()} with \code{forceMonot=TRUE}, or, if only performing Gene Ontology enrichment analysis, using the package \code{goseq}) may be preferred.
  
}
\author{
Stuart K. Archer
}
\references{
}
\keyword{ package 
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
controlSetMonotonic <- matchGenesSmth(x=df, goi=df$Gene.ID[df$diffReg], metric="AveExpr", 
                                      geneIDCol = "Gene.ID", ngenes=5000, verbose = TRUE, 
                                      plotRes = TRUE, forceMonot = T)
}
