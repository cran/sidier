\name{colour.scheme}
\alias{colour.scheme}
\title{internal function for node colour scheme}
\description{This function is called during network representations to set node colours. If the number of colours defined by user do not match with the number of elements, the algorithm provide a default set of colours.}
\usage{
colour.scheme(def=NA, N=NA, colors=c("green2","red","yellow","blue","DarkOrchid1",
"gray51","chocolate","cyan4","saddle brown","aquamarine","chartreuse","chocolate1",
"DarkOrchid3","gray18","gold","DarkOrchid4","green4","gray29", "sienna3","tan1","blue4",
"limegreen","gray73","bisque3","deeppink","red4","OliveDrab4","gray95", "salmon",
"DeepPink4","green yellow","gray4","hot pink","pink2","dark orange","gold3"))
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{def}{a vector containing the set of colours defined by user}
  \item{N}{a numeric representing the number of elements to be coloured}
  \item{colors}{a vector with default colours to be used if 'def' is different from 'N'}
}

\details{
If the number of elements is higher than the number of colours (35 by default), colours are randomly selected.
}

\value{
a vector of strings representing 'N' colours 
}

\author{A. J. Muñoz-Pajares}

\examples{
colour.scheme(def=c("blue","red"),N=4)
Colors<-colour.scheme(def=c("blue","red"),N=4,colors=c("black","gray33","gray66","orange","red"))
plot(c(1:4),col=Colors,pch=16)

#Given 10 individuals classified into three groups,
#this will provide the colour for each individual:
group<-c(1,1,1,2,2,2,1,2,3,3) # defining groups
colour.scheme(N=length(unique(group)))[group]

}
