\name{HapPerPop}
\alias{HapPerPop}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Returns the number of haplotypes per population.
}
\description{Given a two column matrix, this function returns the number of haplotypes per population. The input matrix must contain one row per individual. The first column must contain the population name, while the second must contain the name of the haplotype. The desired matrix can be obtained using "FindHaplo".
}
\usage{
HapPerPop(inputFile = NA, sep = " ", header = F, input = NA,
NameIniPopulations = NA, NameEndPopulations = NA, saveFile = T,
Wname = NA, Iname = NA)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{inputFile}{the name of the file containing the two columns input matrix.
}
  \item{sep}{the character separating columns in the input matrix (space, by default).
}
  \item{header}{a logical value indicating whether the input matrix contains the names of the variables as its first line. (Default=FALSE).
}
  \item{input}{the two coumns input matrix stored in memory
}
  \item{NameIniPopulations}{Position within the input matrix rownames of the initial character referring population origin. 
}
  \item{NameEndPopulations}{Position within the input matrix rownames of the last character referring population origin.
}
  \item{saveFile}{a logical; if TRUE (default), the two ouput matrices computed are saved as two different text files.}
  \item{Wname}{the name given to the output weighted matrix file.
}
  \item{Iname}{
the name given to the output interaction matrix file
}
}
\details{
If both NameIniPopulations and NameEndPopulations are not defined, complete row names are used as population identifiers.}
\value{A list containing two matrices:
  \item{Weighted }{The first matrix (named weighted matrix) contains the abundance of each haplotype per population, represented by the number of haplotypes (columns) found per population (rows). }
  \item{Interaction }{The second matrix (named interaction matrix) contains information about the presence or absence of each haplotype (columns) per population (rows) represented by 1 or 0, respectively.}
}
\author{
A.J. Munoz-Pajares
}

\seealso{FindHaplo
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

# Reading the alignment directly from file, displaying only the 
# weighted matrix:
	HapPerPop(inputFile="3_FindHaplo_Example2_modified.txt",header=TRUE,
saveFile=FALSE)[[1]]

# Reading the alignment from an object and saving the two computed
# distance matrices:
	FH<-read.table("3_FindHaplo_Example2_modified.txt",header=TRUE)
	HapPerPop(input=FH,header=TRUE,saveFile=FALSE)

}


