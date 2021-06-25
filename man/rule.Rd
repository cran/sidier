\name{rule}
\alias{rule}
\title{Threshold to discriminate species.
}
\description{Threshold to discriminate species showing a ratio interspecific/intraspecific distances higher than a given value.
}
\usage{
rule(summary=NULL,rule=NULL,stat.intra="max",
stat.inter="min",pch.intra=16, pch.inter=16,
pch.out=21,col.intra="gray",col.inter="black",
col.out="black",label=F)

}
\arguments{
  \item{summary}{
a list produced by \code{\link{barcode.summary}}. From this list, the maximum intraspecific and the minimu interspecific distances per species are represented. To use any other intra- and interspecific distance, use the "inter" and "intra" options.
}
 \item{rule}{
a numeric. Only species showing interspecific distances higher than 'rule' times the intraspecific distances will be considered for threshold estimation.
}
 \item{stat.intra}{a string, the inter-specific statistic used to estimate the quotient interspecific/intraspecific. Accepted values are "max", "min", "median", and "mean"
}
 \item{stat.inter}{a string, the inter-specific statistic used to estimate the quotient interspecific/intraspecific. Accepted values are "max", "min", "median", and "mean".
}
\item{pch.intra}{Either an integer or single character defining the symbol to be used for intraspecific distances. Only species showing a ratio interspecific/intraspecific higher than the value defined by 'rule' are affected by 'pch.intra'.
}
\item{pch.inter}{Either an integer or single character defining the symbol to be used for interspecific distances Only species showing a ratio interspecific/intraspecific higher than the value defined by 'rule' are affected by 'pch.inter'.
}
\item{pch.out}{Either an integer or single character defining the symbol to be used for species showing a ratio interspecific/intraspecific lower than the value defined by 'rule'.
}
\item{col.intra}{Either an integer or string defining the colour for intraspecific distances showing a ratio interspecific/intraspecific higher than the value defined by 'rule'.
}
\item{col.inter}{Either an integer or string defining the colour for interspecific distances showing a ratio interspecific/intraspecific higher than the value defined by 'rule'.
}
\item{col.out}{Either an integer or string defining the colour forspecies showing a ratio interspecific/intraspecific lower than the value defined by 'rule'.
}
\item{label}{a string to set node labels on those species showing a ratio interspecific/intraspecific higher than the value defined ("rule"), on species showing a ratio interspecific/intraspecific higher than the value defined ("norule"), on all species ("all"). Any other value will produce no lable 
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

\examples{

my.dist<-matrix(c(0,0.3,0.24,0.45,0.23,0.01,0.11,0.34,0.64,0.34,
0.3,0,0.32,0.75,0.65,0.53,0.012,0.52,0.15,0.52,0.24,0.32,
0,0.92,0.36,0.62,0.85,0.008,0.82,0.65,0.45,0.75,0.92,0,
0.22,0.56,0.74,0.46,0.005,0.73,0.23,0.65,0.36,0.22,0,
0.34,0.24,0.42,0.35,0.009,0.01,0.53,0.62,0.56,0.34,0,
0.23,0.73,0.23,0.63,0.11,0.012,0.85,0.74,0.24,0.23,0,
0.25,0.63,0.54,0.34,0.52,0.008,0.46,0.42,0.73,0.25,0,
0.32,0.41,0.64,0.15,0.82,0.005,0.35,0.23,0.63,0.32,0,
0.23,0.34,0.52,0.65,0.73,0.009,0.63,0.54,0.41,0.23,0),
ncol=10,dimnames=list(paste("sp",rep(1:5,2),sep=""),
paste("sp",rep(1:5,2),sep="")))

# rule(barcode.summary(my.dist),rule=10)
}
