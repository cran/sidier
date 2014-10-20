\name{distance.comb}
\alias{distance.comb}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Distance matrices combination
}
\description{
This function allows combining distance matrices. The weight of each matrix must be defined by user.
}
\usage{
distance.comb(matrices = NA, alphas = NA, method = "Corrected",
saveFile = TRUE, na.rm = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{matrices}{
a vector of strings containing the names of the matrices to be combined}
  \item{alphas}{
a vector of numerics containing the weight of each matrix. By default, combination is done giving the same same weight to all matrices.}
  \item{method}{
a string defining whether each distance matrix must be divided by its maximum value before the combination ("Corrected") or not ("Uncorrected"). Consequently, if the "Corrected" method is chosen, both matrices will range between 0 and 1 before being combined.
}
  \item{saveFile}{
a logical; if TRUE (default), the output matrix is saved a text file.
}
  \item{na.rm}{
a logical; if TRUE, missing values are removed before the computation proceeds.
}
}

\details{Matrices do not require to contain rows in the same order. The algorithm will search and combine rows of matrices with identical names. If row names of any input matrix are not defined, the algorithm will combine rows in order.}

\value{
A matrix containing the weighted combination of the original matrices}

\author{
A. J. Muñoz-Pajares
}

\seealso{
nt.gap.comb
}
\examples{
mat1<-matrix(rep(1,16),ncol=4)
mat2<-matrix(rep(2,16),ncol=4)
mat3<-matrix(rep(3,16),ncol=4)
mat4<-matrix(rep(6,16),ncol=4)
mat5<-matrix(c(rep(1,4),rep(3,4),rep(1,4),rep(1,4)),ncol=4)
mat6<-matrix(c(rep(1,4),rep(4,4),rep(4,4),rep(4,4)),ncol=4)
mat7<-mat6
colnames(mat5)<-c("a","b","c","d")
colnames(mat6)<-c("b","a","c","d")
row.names(mat5)<-c("a","b","c","d")
row.names(mat6)<-c("b","a","c","d")

# Matrices with information about the same elements (a-d) but in different order, are automatically
# sorted before combination....

distance.comb(matrices=c("mat5","mat6"),saveFile=FALSE,method="Uncorrected")

# ... but this is not possible if any of the matrices lack row names:

distance.comb(matrices=c("mat5","mat7"),saveFile=FALSE,method="Uncorrected")

# More examples:
distance.comb(matrices=c("mat1","mat2","mat3","mat4"),alphas=rep(0.25,4),
saveFile=FALSE,method="Uncorrected")

distance.comb(matrices=c("mat1","mat2","mat3","mat4"),alphas=rep(0.25,4),
saveFile=FALSE,method="Corrected")

distance.comb(matrices=c("mat1","mat2","mat3","mat4"),alphas=c(0.66,0.33,0,0),
saveFile=FALSE,method="Uncorrected")

distance.comb(matrices=c("mat1","mat2","mat3","mat6"),alphas=c(0.66,0.33,0,0),
saveFile=FALSE,method="Uncorrected")
}
