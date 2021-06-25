\name{barcode.quality}
\alias{barcode.quality}
\title{Estimates of barcode quality
}
\description{
Provides several estimates of the quality of a barcode classification, comparing network modules with attributed species names  
}
\usage{
barcode.quality(dismat=NA,threshold=NA,refer2max=FALSE,save.file=FALSE,
modFileName="Modules_summary.txt",verbose=FALSE,output="list")

}
\arguments{
  \item{dismat}{
a matrix containing the pairwise genetic distances between individual sequences
} 
 \item{threshold}{
a numeric between 0 and 1, is the value of the maximum distance to be represented as a link in the network
}
 \item{refer2max}{
a logic, "TRUE" to refer the threshold value to the maximum distance in the input matrix (e.g., a value of 0.32 will represent a link between nodes showing distances equal or lower than 32\% of the maximum distance found in the distance matrix). "FALSE" to refer the threshold to a specific value (e.g., a value of 0.32 will represent a link between nodes showing distances equal or lower than 0.32, regardless the maximum distance found in the distance matrix).
}
  \item{save.file}{
a logic, "TRUE" to save the summary of network modules, attributing every individual to a module.
}
  \item{modFileName}{
if save.file=TRUE, a string: the name of the file containing the summary of network modules.
}
  \item{verbose}{
a logic, "TRUE" to obtain a complete report of the quality estimation (see details).
}
  \item{output}{
if verbose=TRUE, a string controlling the type of object produced for the output, being either "matrix" or "list".
}

}
\details{
This function assumes that the species names reflect the "real" taxonomic status and compare these names with the modules obtained in the network analysis. The quality is evaluated using different estimators: 
 \deqn{Accuracy= \frac{T_{+} + T_{-}}{T_{+} + T_{-} + F_{+} + F_{-}}}
 \deqn{Precision= \frac{T_{+}}{T_{+} + F_{+}}}
 \deqn{Fscore= \frac{T_{+}}{T_{+} + F_{+} + F_{-}}}
 \deqn{Qvalue = \frac{1}{N}\sum_{1}^{N}\frac{S_{link}}{S_{all}+S_{unlink}}}

where T+ is the number of true positives (number of sequences with the same species name and classified in the same module); T- is the number of true negatives (number of sequences with different species name and classified in different modules); F+ represents false positive (number of sequences with different species name classified in the same module); F- is the number of false negative (number of sequences with the same species name classified in different modules); N is the number of nodes in the network, Slink is the number of nodes of the same species connected to the node i; Sunlink is the number of nodes of the same species belonging to a different module; and Sall is the number of all possible connections to other nodes of the same species.

}
\value{
If verbose is set to "FALSE", a matrix with the estimators of the barcode quality. If verbose is set to "TRUE", either a matrix or a list (depending on the output option selected) containing the following elements:
\item{Number.of.modules}{Number of modules found in the network analysis.}
\item{Number.of.species.per.module}{A matrix containing: The number of species classified in only one module (N.sp.mod.1); the maximum number of species found in a module (N.sp.mod.MAX); and the mean number of species found per module (N.sp.mod.MED).}\item{Number.of.species}{The number of species defined for the analysis.}
\item{Number.of.modules.per.species}{A matrix containing: The number of modules composed of only one species (N.mod.sp.1); the maximum number of modules containing the same species; the mean number of modules containing the same species.}
\item{Number.of.modules.fitting.defined.species}{The number of modules containing only one species but all the individuals of this species.}
\item{Quality.estimates}{A matrix containing the Qvalue, Accuracy, Precision and Fscore of the barcode classification.}

}
\author{
A.J. Muñoz-Pajares
}

\examples{
# my.dist<-matrix(abs(rnorm(100)),ncol=10,
# dimnames=list(paste("sp",rep(1:5,2),sep=""),
# paste("sp",rep(1:5,2),sep="")))
# my.dist<-as.matrix(as.dist(my.dist))
# 
# barcode.quality(dismat=my.dist,threshold=0.2,refer2max=FALSE,save.file=TRUE,
# modFileName="Modules_summary.txt",verbose=FALSE,output="list")
}
