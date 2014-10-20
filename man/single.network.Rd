\name{single.network}
\alias{single.network}
\title{ Plot a network given a threshold
}
\description{This function plots a network connecting nodes showing distances equal or lower than the defined threshold value.
}
\usage{
single.network(dis, threshold = NA, ptPDF = TRUE, ptPDFname = "Network.pdf",
bgcol = "white", label.col = "black", label = colnames(dis), modules = FALSE,
moduleCol = NA, modFileName = "Modules_summary.txt", na.rm.row.col = FALSE)}
\arguments{
  \item{dis}{the distance matrix to be represented}
  \item{threshold}{a numeric between 0 and 1, is the value of the maximum distance to be considered as a link. This value is referred to the maximum distance in the input matrix (e.g., a value of 0.32 will represent a link between nodes showing distances equal or lower than 32\% of the maximum distance found in the distance matrix).}
  \item{ptPDF}{a logical, must the percolation network be saved as a pdf file?}
  \item{ptPDFname}{if ptPDF=TRUE, the name of the pdf file containing the percolation network to be saved ("percolationNetwork.pdf", by default)}
  \item{bgcol}{string defining the colour of the background for each node in the network. Can be equal for all nodes (if only one colour is defined), customized (if several colours are defined), or can represent different modules (see modules option).}
  \item{label.col}{vector of strings defining the colour of labels for each node in the network. Can be equal for all nodes (if only one colour is defined) or customized (if several colours are defined).}
  \item{label}{vector of strings, labels for each node. By default are the column names of the distance matrix (dis). (See substr function in base package to automatically reduce name lengths).}
  \item{modules}{a logical, must nodes belonging to different modules be represented as different colours?}
  \item{moduleCol}{(if modules=TRUE) vector of strings, defining the colour of nodes belonging to different modules in the network.}
  \item{modFileName}{(if modules=TRUE) the name of a generated file containing a summary of module results}
  \item{na.rm.row.col}{a logical; if TRUE, missing values are removed before the computation proceeds.}
}
\author{A. J. Muñoz-Pajares}
\seealso{\code{\link{perc.thr}}, \code{\link{NINA.thr}}}
\examples{
#generating distance matrix:
dis<-matrix(nrow=12,c(0.0000,0.5000,0.1875,0.5000,0.6250,0.5000,0.2500,0.6250,
0.3750,0.3750,0.3750,0.3750,0.5000,0.0000,0.7500,0.0000,0.6250,0.0000,0.8750,
0.6250,0.3750,0.3750,0.3750,0.3750,0.1875,0.7500,0.0000,0.7500,0.8750,0.7500,
0.2500,0.8750,0.6250,0.6250,0.6250,0.6250,0.5000,0.0000,0.7500,0.0000,0.6250,
0.0000,0.8750,0.6250,0.3750,0.3750,0.3750,0.3750,0.6250,0.6250,0.8750,0.6250,
0.0000,0.6250,0.5000,0.0000,0.2500,0.2500,0.2500,0.2500,0.5000,0.0000,0.7500,
0.0000,0.6250,0.0000,0.8750,0.6250,0.3750,0.3750,0.3750,0.3750,0.2500,0.8750,
0.2500,0.8750,0.5000,0.8750,0.0000,0.5000,0.5000,0.5000,0.5000,0.5000,0.6250,
0.6250,0.8750,0.6250,0.0000,0.6250,0.5000,0.0000,0.2500,0.2500,0.2500,0.2500,
0.3750,0.3750,0.6250,0.3750,0.2500,0.3750,0.5000,0.2500,0.0000,0.0000,0.0000,
0.0000,0.3750,0.3750,0.6250,0.3750,0.2500,0.3750,0.5000,0.2500,0.0000,0.0000,
0.0000,0.0000,0.3750,0.3750,0.6250,0.3750,0.2500,0.3750,0.5000,0.2500,0.0000,
0.0000,0.0000,0.0000,0.3750,0.3750,0.6250,0.3750,0.2500,0.3750,0.5000,0.2500,
0.0000,0.0000,0.0000,0.0000))
row.names(dis)<-c("Population1_sequence1","Population1_sequence2",
"Population1_sequence3","Population1_sequence4","Population2_sequence1",
"Population2_sequence2","Population2_sequence3","Population2_sequence4",
"Population3_sequence1","Population3_sequence2","Population3_sequence3",
"Population3_sequence4")
colnames(dis)<-row.names(dis)

#Representing distances equal or lower than 37% of the maximum distance:
single.network(dis=dis,threshold=0.37,label=paste(substr(row.names(dis),11,11),
substr(row.names(dis),21,21),sep="-"))

}

