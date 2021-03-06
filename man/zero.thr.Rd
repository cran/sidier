\name{zero.thr}
\alias{zero.thr}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Zero distance networks}
\description{Given a distance matrix, this function computes a network connecting nodes showing distances equal to zero.}
\usage{
zero.thr(dis,ptPDF=TRUE,ptPDFname="zero_Network.pdf",cex.label=1,cex.vertex=1,
bgcol="white",label.col="black",label=colnames(dis),modules=FALSE,moduleCol=NA,
modFileName="Modules_summary.txt",ncs=4,na.rm.row.col=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dis}{the input distance matrix}
  \item{ptPDF}{a logical, must the resulting network be saved as a pdf file?}
  \item{ptPDFname}{if ptPDF=TRUE, the name of the pdf file containing the resulting network to be saved ("zero_Network.pdf", by default)}
  \item{cex.label}{a numeric; the size of the node labels.}
  \item{cex.vertex}{a numeric; the size of the nodes.}
  \item{bgcol}{the background colour for each node in the network. Can be equal for all nodes (if only one colour is defined), customized (if several colours are defined), or can represent different modules (see "modules" option).
}
  \item{label.col}{vector of strings defining the colour of labels for each node in the network. Can be equal for all nodes (if only one colour is defined) or customized (if several colours are defined).
}
  \item{label}{vector of strings, labels for each node. By default are the column names of the distance matrix (dis). (See the 'substr' function in base package to automatically reduce name lengths).
}
  \item{modules}{a logical, must nodes belonging to different modules be represented with different colours? If TRUE, a text file containing information on modules for each node is also produced.
}
  \item{moduleCol}{(if modules=TRUE) a vector of strings, defining the colour of nodes belonging to different modules in the network. If 'NA' or less colours than modules are defined, colours are automatically defined.}
  \item{modFileName}{(if modules=TRUE) the name of the text file containing a summary of module results
}
  \item{ncs}{a numeric; number of decimal places to display threshold in plot title.}
  \item{na.rm.row.col}{a logical; if TRUE, missing values are removed before the computation proceeds.}
}

\details{
In some circumstances you may get distance matrices showing off-diagonal zeros. In such cases you may consider that the
existence of these off-diagonal zeros suggests that some of the groups you defined (e.g., populations) are not genetically different.
Thus, you must re-define groups to get a matrix composed only by different  groups using the 'mergeNodes' function and 
estimate a percolation network using the 'perc.thr' function. On the other hand, you may consider that, despite the off-
diagonal zeros, the groups you defined are actually different. In that case you may not be able to estimate a 
percolation threshold, but you can represent the original distance matrix using the 'NINA.thr' or the 'zero.thr' functions.

'mergeNodes' select all rows (and columns) showing a distance equal to zero and generates a new row (and column). The distance between the new merged and the remaining rows (or columns) in the matrix is estimated as the arithmetic mean of the selected elements. The biological interpretation of the new matrix could be hard if the original matrix shows a large number of off-diagonal zeros.

'perc.thr' estimates a threshold to represent a distance matrix as a network. To estimate this threshold, the algorithm represents as a link all distances lower than a range of thresholds (by default, select 101 values from 0 to 1), defined as the percentage of the maximum distance in the input matrix. For each threshold a network is built and the number of clusters (that is, the number of isolated groups of nodes) in the network is also estimated. Finally, the algorithm selects the lower threshold connecting a higher number of nodes. Note that the resulting network may show isolated nodes if it is necessary to represent a large number of links to connect a low number of nodes.

'NINA.thr' is identical to 'perc.thr', but, in the last step, the algorithm selects the lower threshold connecting all nodes in a single cluster. The information provided by this function may be limited if the original distance matrix shows high variation.

\strong{'zero.thr'} represents as a link only distances equal to zero. The information provided by this function may be limited if the original matrix shows few off-diagonal zeros.
}
\value{A network connecting nodes showing a distance equal to zero.}
\author{A. J. Muñoz-Pajares}
\seealso{\code{\link{NINA.thr}}, \code{\link{perc.thr}}, \code{\link{mergeNodes}}}
\examples{

#EXAMPLE 1: FEW OFF-DIAGONAL ZEROS
#Generating a distance matrix:
Dis1<-matrix(c(
0.00,0.77,0.28,0.94,0.17,0.14,0.08,0.49,0.64,0.01,
0.77,0.00,0.12,0.78,0.97,0.02,0.58,0.09,0.36,0.33,
0.28,0.12,0.00,0.70,0.73,0.06,0.50,0.79,0.80,0.94,
0.94,0.78,0.70,0.00,0.00,0.78,0.04,0.42,0.25,0.85,
0.17,0.97,0.73,0.00,0.00,0.30,0.55,0.12,0.68,0.99,
0.14,0.02,0.06,0.78,0.30,0.00,0.71,1.00,0.64,0.88,
0.08,0.58,0.50,0.04,0.55,0.71,0.00,0.35,0.84,0.76,
0.49,0.09,0.79,0.42,0.12,1.00,0.35,0.00,0.56,0.81,
0.64,0.36,0.80,0.25,0.68,0.64,0.84,0.56,0.00,0.62,
0.01,0.33,0.94,0.85,0.99,0.88,0.76,0.81,0.62,0.00),ncol=10)
colnames(Dis1)<-c(paste("Pop",c(1:10),sep=""))
row.names(Dis1)<-colnames(Dis1)

# No percolation threshold can be found.
#perc.thr(Dis1)

#Check Dis1 and merge populations showing distances equal to zero:
Dis1
Dis1_Merged<-mergeNodes(dis=Dis1)
#Check the merged matrix. A new "population" has been defined merging populations 4 and 5.
#Distances between the merged and the remaining populations are estimated as the arithmetic mean.
Dis1_Merged
# It is now possible to estimate a percolation threshold
perc.thr(dis=Dis1_Merged,ptPDF=FALSE, estimPDF=FALSE, estimOutfile=FALSE) 

# EXAMPLE 2: TOO MANY OFF-DIAGONAL ZEROS
#Generating a distance matrix:
Dis2<-matrix(c(
0.00,0.77,0.28,0.00,0.17,0.14,0.00,0.49,0.64,0.01,
0.77,0.00,0.12,0.00,0.97,0.02,0.00,0.09,0.36,0.33,
0.28,0.12,0.00,0.70,0.73,0.06,0.50,0.79,0.00,0.94,
0.00,0.00,0.70,0.00,0.00,0.78,0.04,0.00,0.00,0.00,
0.17,0.97,0.73,0.00,0.00,0.30,0.55,0.12,0.00,0.00,
0.14,0.02,0.06,0.78,0.30,0.00,0.71,1.00,0.64,0.00,
0.00,0.00,0.50,0.04,0.55,0.71,0.00,0.35,0.84,0.00,
0.49,0.09,0.79,0.00,0.12,1.00,0.35,0.00,0.56,0.81,
0.64,0.36,0.00,0.00,0.00,0.64,0.84,0.56,0.00,0.62,
0.01,0.33,0.94,0.00,0.00,0.00,0.00,0.81,0.62,0.00),ncol=10)
colnames(Dis2)<-c(paste("Pop",c(1:10),sep=""))
row.names(Dis2)<-colnames(Dis2)

# # No percolation threshold can be found
# #perc.thr(Dis2)
# 
# #Check Dis2 and merge populations showing distances equal to zero:
# Dis2
# Dis2_Merged<-mergeNodes(dis=Dis2)
# 
# #Check the merged matrix. Many new "populations" have been defined and both the new
# #matrix and the resulting network are difficult to interpret:
# Dis2_Merged
# perc.thr(dis=Dis2_Merged,ptPDF=FALSE, estimPDF=FALSE, estimOutfile=FALSE) 
# 
# #Instead of percolation network, representing zeros as the lowest values may be informative:
# zero.thr(dis=Dis2,ptPDF=FALSE)
# # Adjusting sizes and showing modules:
# zero.thr(dis=Dis2,ptPDF=FALSE,cex.label=0.8,cex.vertex=1.2,modules=TRUE)
# 
# #In the previous example, the 'zero.thr' method is unuseful: 
# zero.thr(dis=Dis1,ptPDF=FALSE)
# 
# #In both cases, the 'No Isolation Nodes Allowed' method yields an informative matrix:
# NINA.thr(dis=Dis1)
# NINA.thr(dis=Dis2)
}
