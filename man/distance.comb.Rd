\name{distance.comb}
\alias{distance.comb}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Distance matrices combiantion
}
\description{
This function allows combining distance matrices. The weight of each matrix must be defined by user.
}
\usage{
distance.comb(matrices = NA, alphas = NA, method = "Corrected",
saveFile = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{matrices}{
a vector of strings containing the names of the matrices to be combined}
  \item{alphas}{
a vector of numerics containing the weight of each matrix}
  \item{method}{
a string defining whether each distance matrix must be divided by its maximum value before the combination ("Corrected") or not ("Uncorrected"). Consequently, if the "Corrected" method is chosen, both matrices will range between 0 and 1 before to be combined.
}
  \item{saveFile}{
a logical; if TRUE (default), each ouput matrix is saved in a different text file.
}
}

\value{
A matrix containing the weighted combination of the original matrices}
\author{
A.J. Munoz-Pajares
}

\examples{
mat1<-matrix(rep(1,16),ncol=4)
mat2<-matrix(rep(2,16),ncol=4)
mat3<-matrix(rep(3,16),ncol=4)
mat4<-matrix(rep(6,16),ncol=4)

distance.comb(matrices=c("mat1","mat2","mat3","mat4"),alphas=rep(0.25,4),
saveFile=FALSE,method="Uncorrected")
distance.comb(matrices=c("mat1","mat2","mat3","mat4"),alphas=rep(0.25,4),
saveFile=FALSE,method="Corrected")
distance.comb(matrices=c("mat1","mat2","mat3","mat4"),alphas=c(0.66,0.33,0,0),
saveFile=FALSE,method="Uncorrected")
}
