\name{barcode.summary}
\alias{barcode.summary}
\title{Summary of the inter- and intraspecific distances
}
\description{For every species, provides the minimum, maximum, mean and median values of inter- and intraspecific distances.
}
\usage{
barcode.summary(dismat=NULL,save.distances=FALSE,folder.name="distance_matrices")
}
\arguments{
  \item{dismat}{
a symmetric matrix containing the pairwise genetic distances between individual sequences
}
  \item{save.distances}{
a logic, "TRUE" to save the pairwise distances estimated per species (one file per species)
}
  \item{folder.name}{
a string, if save.distance=TRUE, the name of the folder to save distances
}

}

\value{
A list with two elements:
\item{Intraspecific}{a matrix containing information about the intraspecific distances.
}
\item{Interspecific}{a matrix containing information about the interspecific distances.
}
In both cases, the information provided is the minimum, maximum, median, mean, first and third quartile values.
}
\author{
A.J. Muñoz-Pajares
}

\seealso{
\code{\link{barcode.gap}}
}
\examples{

my.dist<-matrix(abs(rnorm(100)),ncol=10,
dimnames=list(paste("sp",rep(1:5,2),sep=""),
paste("sp",rep(1:5,2),sep="")))
my.dist<-as.matrix(as.dist(my.dist))
barcode.summary(my.dist)

}
