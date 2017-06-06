\name{single.network.module}
\alias{single.network.module}
\title{Get modules and network given a threshold
}
\description{Gets details on modules and connections in the network connecting nodes showing distances equal or lower than the defined threshold value.
}
\usage{
single.network.module(dis,threshold=NA,refer2max=TRUE,out="module",
save.file=FALSE,modFileName="Modules_summary.txt")}

\arguments{
  \item{dis}{the distance matrix to be represented}
  \item{threshold}{a numeric between 0 and 1, is the value of the maximum distance to be considered as a link. This value is referred to the maximum distance in the input matrix (e.g., a value of 0.32 will represent a link between nodes showing distances equal or lower than 32\% of the maximum distance found in the distance matrix).}
 \item{refer2max}{
a logic, "TRUE" to refer the threshold value to the maximum distance in the input matrix (e.g., a value of 0.32 will represent a link between nodes showing distances equal or lower than 32\% of the maximum distance found in the distance matrix). "FALSE" to refer the threshold to a specific value (e.g., a value of 0.32 will represent a link between nodes showing distances equal or lower than 0.32, regardless the maximum distance found in the distance matrix).
}
 \item{out}{a string, the type of output, "module" to get a matrix with two columns giving each sequence name and the module it belongs to, and "network" to get a square matrix representing connection (1) or lack of connection (0) between sequences in the network.}
  \item{save.file}{
a logic, "TRUE" to save the summary of network modules, attributing every individual to a module.
}

 \item{modFileName}{(if modules=TRUE) the name of a generated file containing a summary of module results}
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
single.network.module(dis=dis,threshold=0.37)
single.network.module(dis=dis,threshold=0.37,out="network")

# Compare these outputs with:
single.network(dis=dis,threshold=0.37,label=paste(substr(row.names(dis),11,11),
substr(row.names(dis),21,21),sep="-"))
}

