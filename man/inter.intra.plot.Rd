\name{inter.intra.plot}
\alias{inter.intra.plot}
\title{Histogram of the intra- and interspecific distances
}
\description{Plot histogram for inter and intra-specific distances together}
\usage{
inter.intra.plot(dismat=NA, xlim=NULL,ylim=NULL,
intra.col="gray",intra.density=0,intra.n=30,plot="N",
inter.col="black",inter.density=0,inter.n=30,legend=TRUE,
main="",xlab="Genetic distances",ylab=NULL)
}
\arguments{
  \item{dismat}{
a symmetric matrix containing the pairwise genetic distances between individual sequences.
}
  \item{xlim}{
a vector containing the minimum and maximum value in the x-axis
}
  \item{ylim}{
a vector containing the minimum and maximum value in the y-axis
}
  \item{intra.col}{
the colour for the intraspecific distance distribution
}
  \item{intra.density}{
a numeric, the density of shading lines for the intraspecific distance distribution
}
  \item{intra.n}{
a numeric, the number of categories to represent in the intraspecific distance distribution
}
  \item{plot}{
a string, "freq" to represent frequency values in the y-axis and "N" for number of occurrences
}
  \item{inter.col}{
the colour for the interspecific distance distribution
}
  \item{inter.density}{
a numeric, the density of shading lines for the interspecific distance distribution
}
  \item{inter.n}{
a numeric, the number of categories to represent in the interspecific distance distribution
}
  \item{legend}{
a logic, "TRUE" to show plot legend
}
  \item{main}{
a string containing the title of the plot
}
  \item{xlab}{
a string with the label of the x-axis
}
  \item{ylab}{
a string with the label of the x-axis
}

}
\value{
A list with two elements:
\item{Intraspecific}{a vector containing all the intraspecific distances.
}
\item{Interspecific}{a vector containing all the interspecific distances.
}
}
\author{
A.J. Muñoz-Pajares
}

\examples{
# # Generating a distance matrix:
# 
# my.mat<-matrix(nrow=100,ncol=100,
# dimnames=list(paste("sp",rep(1:2,50),
# sep=""),paste("sp",rep(1:2,50),sep="")))
# L<-my.mat[seq(1,nrow(my.mat),2),seq(1,ncol(my.mat),2)]
# my.mat[seq(1,nrow(my.mat),2),seq(1,ncol(my.mat),2)]<-rnorm(0.15,n=L,sd=0.01)
# my.mat[seq(2,nrow(my.mat),2),seq(2,ncol(my.mat),2)]<-rnorm(0.15,n=L,sd=0.01)
# my.mat[seq(1,nrow(my.mat),2),seq(2,ncol(my.mat),2)]<-rnorm(0.3,n=L,sd=0.04)
# my.mat[seq(2,nrow(my.mat),2),seq(1,ncol(my.mat),2)]<-rnorm(0.3,n=L,sd=0.04)
# #Converting to symmetric
# my.mat<-as.matrix(as.dist(my.mat))  
# inter.intra.plot(dismat=my.mat)
# inter.intra.plot(dismat=my.mat,intra.n=10)
# inter.intra.plot(dismat=my.mat,plot="Freq",intra.n=10)
}
