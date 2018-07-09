\name{sidier-package}
\alias{sidier-package}
\alias{sidier}
\docType{package}
\title{SIDIER: Substitution and Indel Distances to Infer Evolutionary Relationships}
\description{Package for evolutionary reconstruction and bar code analyses based on substitutions and insertion-deletion (indels) in a distance-based framework.}
\details{

The DESCRIPTION file:
\packageDESCRIPTION{sidier}
\packageIndices{sidier}
Functions can be classified according to the following groups:

- \emph{Barcode analysis}:
\code{\link{barcode.gap}}; 
\code{\link{barcode.quality}}; 
\code{\link{barcode.summary}}; 
\code{\link{compare.dist}}; 
\code{\link{genbank.sp.names}}; 
\code{\link{inter.intra.plot}}; 
\code{\link{rule}}. 

- \emph{Example files}:
\code{\link{alignExample}}; 
\code{\link{Example_Spatial.plot_Alignment}}; 
\code{\link{ex_alignment1}}; 
\code{\link{ex_Coords}}.

- \emph{Indel coding methods}:
\code{\link{BARRIEL}}; 
\code{\link{FIFTH}}; 
\code{\link{MCIC}}; 
\code{\link{SIC}}.

- \emph{Matrix/network manipulation}:
\code{\link{distance.comb}}; 
\code{\link{mergeNodes}}; 
\code{\link{nt.gap.comb}}; 
\code{\link{simplify.network}}.

- \emph{Network-from-distance methods}:
\code{\link{NINA.thr}}; 
\code{\link{perc.thr}}; 
\code{\link{zero.thr}}; 
\code{\link{single.network}}; 
\code{\link{single.network.module}}.

- \emph{Network visualization}:
\code{\link{mutation.network}}; 
\code{\link{pie.network}}; 
\code{\link{double.plot}}; 
\code{\link{colour.scheme}}; 
\code{\link{spatial.plot}}.

- \emph{Sequence/haplotype analysis}:
\code{\link{FilterHaplo}}; 
\code{\link{FindHaplo}}; 
\code{\link{GetHaplo}}; 
\code{\link{HapPerPop}}; 
\code{\link{mutationSummary}}; 
\code{\link{pop.dist}}; 
\code{\link{simuEvolution}}.

}
\author{
\packageAuthor{sidier}

Maintainer: \packageMaintainer{sidier}
}
\references{
Muñoz-Pajares, A. J. (2013). SIDIER: substitution and indel distances to infer evolutionary relationships. Methods in Ecology and Evolution 4, 1195-1200. doi: 10.1111/2041-210X.12118
}
%\keyword{ distance matrices, networks, evolution }
%\seealso{
%~~ Optional links to other man pages, e.g. ~~
%~~ \code{\link[<pkg>:<pkg>-package]{<pkg>}} ~~
%}
%\examples{
%~~ simple examples of the most important functions ~~
%}
