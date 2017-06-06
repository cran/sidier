\name{barcode.gap}
\alias{barcode.gap}
\title{Barcode gap identification 
}
\description{
Identifies barcode gaps based on representing intra- and interspecific distances. Species above the 1:1 line are considered to show a barcode gap.
}
\usage{
barcode.gap(summary=NULL,stat.intra="max",stat.inter="min",
xlab=NULL, ylab=NULL, legend=TRUE, lab.nodes="nogap")
}
\arguments{
  \item{summary}{
a list produced by \code{\link{barcode.summary}}. From this list, the maximum intraspecific and the minimum interspecific distances per species are represented. To use any other intra- and interspecific distance, use the "inter" and "intra" options.
}
 \item{stat.intra}{a string, the inter-specific statistic used to estimate the quotient interspecific/intraspecific. Accepted values are "max", "min", "median", and "mean"
}
 \item{stat.inter}{a string, the inter-specific statistic used to estimate the quotient interspecific/intraspecific. Accepted values are "max", "min", "median", and "mean".
}
  \item{xlab}{
a string, the x-axis label
}
  \item{ylab}{
a string, the y-axis label
}
  \item{legend}{
a logic, to show information about species showing and lacking barcode within the plot
}
  \item{lab.nodes}{
a string to select the name of species to be represented in the plot: "gap" to represent species showing barcode gap; "nogap" to represent species lacking barcode gap; "all" for representing all species names. Other value will represent no names.
}

}

\value{
A list with two elements:
\item{no.barcode.gap}{a matrix containing the name of the species lacking barcode gap and their mean intra- and interspecific distances.
}
\item{barcode.gap}{a matrix containing the name of the species showing barcode gap and their mean intra- and interspecific distances.
}
}

\author{
A.J. Muñoz-Pajares
}

\seealso{
\code{\link{barcode.summary}}
}
\examples{

my.dist<-matrix(abs(rnorm(100)),ncol=10,
dimnames=list(paste("sp",rep(1:5,2),sep=""),
paste("sp",rep(1:5,2),sep="")))
my.dist<-as.matrix(as.dist(my.dist))
sum<-barcode.summary(my.dist)
barcode.gap(sum)

}
