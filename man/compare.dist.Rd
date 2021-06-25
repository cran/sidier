\name{compare.dist}
\alias{compare.dist}
\title{Threshold to discriminate species comparing intra- and interspecific distance distributions
}
\description{
This function implements the Lefebure's method to quantify the overlap between two distributions and to determine the best threshold value to discriminate them.
}
\usage{
compare.dist(distr1=NULL,distr2=NULL,N=50,
normalize=TRUE,main=NA,col1="gray",col2="black",
col.border1="gray",col.border2="black",
col.line1="gray",col.line2="black",
Ylab=c("Abundance","Abundance","Success"),
Xlab=c("data1","data2","Threshold"))

}
\arguments{
  \item{distr1}{
a matrix containing the pairwise genetic distances between individual sequences
}
  \item{distr2}{
a vector containing the minimum and maximum value in the x-axis
}
  \item{N}{
a numeric, the number of categories for the x-axis.
}
  \item{normalize}{
a logic, "TRUE" to display percentage and "FALSE" for number of occurrences in the Y axis.
}
  \item{main}{
a vector with two elements containing the main titles of both plots.
}
  \item{col1}{a string, the color to fill the histogram for distribution 1}
  \item{col2}{a string, the color to fill the histogram for distribution 2}
  \item{col.border1}{a string, the color for the border around the histogram for distribution 1}
  \item{col.border2}{a string, the color for the border around the histogram for distribution 2}
  \item{col.line1}{a string, the color for the line representing distribution 1}
  \item{col.line2}{a string, the color for the line representing distribution 2}

  \item{Ylab}{
a three strings vector containing the labels of Y axes for the three plots to be represented.
}
  \item{Xlab}{
a three strings vector containing the labels of X axes for the three plots to be represented
}

}
\value{
The estimated threshold and its success of identification.
}
\references{
Lefébure T, Douady CJ, Gouy M, Gibert J (2006). Relationship between morphological taxonomy and molecular divergence within Crustacea: Proposal of a molecular threshold to help species delimitation. Mol Phylogenet Evol 40: 435–447. 
}
\author{
A.J. Muñoz-Pajares
}

\examples{

# ## Weak overlap
# intra<-rnorm(mean=0.08,sd=0.04,n=100)
# inter<-rnorm(mean=0.38,sd=0.10,n=1000)
# intra[intra<0]<-0
# inter[inter<0]<-0
# compare.dist(distr1=intra,distr2=inter,N=50)
# 
# # Strong overlap
# distr1<-rnorm(5000,mean=0.25,sd=0.070)
# distr2<-rnorm(5000,mean=0.31,sd=0.075)
# N<-50
# compare.dist(distr1,distr2,N)
# 
}
