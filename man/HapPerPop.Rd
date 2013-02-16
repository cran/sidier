\name{HapPerPop}
\alias{HapPerPop}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Returns the number of haplotypes per population.
%%  function to do ... 
}
\description{Given a two column matrix, this function returns the number of haplotypes per population. The input matrix must contain one row per individual. The first column must contain the population name, while the second must contain the name of the haplotype. The desired matrix can be obtained using "FindHaplo".
%%   A concise (1-5 lines) description of what the function does. 
}
\usage{
HapPerPop(readfile = T, sep = " ", header = F, inputFile = NA, input = NA, saveFile = T, Wname = NA, Iname = NA)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{readfile}{a logical; if TRUE (default) the input matrix is provided in a text file. If FALSE, the matrix is provided as an R object.
%%     Describe \code{readfile} here
}
  \item{sep}{the character separating columns in the input matrix (space, by default).
%%     Describe \code{sep} here
}
  \item{header}{a logical value indicating whether the input matrix contains the names of the variables as its first line. (Default=FALSE).
%%     Describe \code{header} here
}
  \item{inputFile}{(if readfile=TRUE) the name of the file containing the input matrix.(if readfile=TRUE) the name of the file containing the input matrix.
%%     Describe \code{inputFile} here
}
  \item{input}{a logical; if TRUE (default), the two ouput matrices computed are saved as two different text files.
%%     Describe \code{input} here
}
  \item{saveFile}{a logical; if TRUE (default), the two ouput matrices computed are saved as two different text files.
%%     Describe \code{saveFile} here
}
  \item{Wname}{the name given to the output weighted matrix file.
%%     Describe \code{Wname} here
}
  \item{Iname}{the name given to the output interaction matrix file
%%     Describe \code{Iname} here
}
}
\details{Two output matrices are estimated, one giving the abundance of each haplotype per population (named weighted matrix) and the other representing presence/absence of each haplotype per population by 1/0 (named interaction matrix).
%%   If necessary, more details than the description above 
}
\value{A list containing two matrices. 
%%  Describe the value returned
%%  If it is a LIST, use
  \item{comp1 }{The first matrix contains the weighted matrix, that is, the number of haplotypes (columns) found per population (rows). }
  \item{comp2 }{The second is the interaction matrix, containing information about the presence or absence of each haplotype (columns) per population (rows).}
%% ...
}
\author{A.J. Munoz-Pajares
%%  who you are
}


\seealso{FindHaplo
%% objects to See Also as \code{\link{help}}, 
}
\examples{
cat("Sequence.Name Haplotype.Name",
"Population1 H1",
"Population1 H2",
"Population1 H3",
"Population1 H2",
"Population2 H4",
"Population2 H5",
"Population2 H6",
"Population2 H4",
"Population3 H7",
"Population3 H7",
"Population3 H7",
"Population3 H7",
     file = "3_FindHaplo_Example2_modified.txt", sep = "\n")
     example2_2 <- read.table("3_FindHaplo_Example2_modified.txt",header=TRUE)
	 FH<-read.table("3_FindHaplo_Example2_modified.txt",header=TRUE)

# Reading the alignment from an object and saving the two computed distance matrices:
HapPerPop(readfile=FALSE,input=FH,header=TRUE,saveFile=FALSE)

# Reading the alignment directly from file, displaying only the weighted matrix:
HapPerPop(readfile=TRUE,inputFile="3_FindHaplo_Example2_modified.txt",header=TRUE,saveFile=FALSE)[[1]]
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ kwd1 }
%\keyword{ kwd2 }% __ONLY ONE__ keyword per line
