\name{HapPerPop}
\alias{HapPerPop}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Returns the number of haplotypes per population.
}
\description{Given a two column matrix, this function returns the number of haplotypes per population (weighted matrix) and the presence/absence of haplotypes per population (interaction matrix). The input matrix must contain one row per individual. The first column must contain the population name, while the second must contain the name of the haplotypes. This input matrix can be obtained using the "FindHaplo" function.
}
\usage{
HapPerPop(inputFile = NA, sep = " ", input = NA, header = FALSE,
NameIniPopulations = NA, NameEndPopulations = NA, saveFile = TRUE,
Wname = NA, Iname = NA)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{inputFile}{the name of the file containing the two columns input matrix. Alternatively you can provide the name of a matrix stored in memory using the "input" option.
}
  \item{sep}{the character separating columns in the input matrix (space, by default).
}
  \item{input}{the two columns input matrix stored in memory. Alternatively you can provide the name of the file containing the input matrix using the "inputFile" option.
}
  \item{header}{a logical value indicating whether the input matrix contains the names of the variables as its first line. (Default=FALSE).
}
  \item{NameIniPopulations}{Position within the input matrix rownames of the initial character referring population name. This option is useful if names contained in the first column includes more information than the population name (e.g., marker name, individual details...).
}
  \item{NameEndPopulations}{Position within the input matrix rownames of the last character referring population name. This option is useful if names contained in the first column includes more information than the population name (e.g., marker name, individual details...).
}
  \item{saveFile}{a logical; if TRUE (default), the two ouput matrices computed are saved as two different text files.}
  \item{Wname}{the name given to the output weighted matrix file.
}
  \item{Iname}{
the name given to the output interaction matrix file
}
}
\details{
If both NameIniPopulations and NameEndPopulations are not defined, the names contained in the input matrix first column are used as population identifiers.}
\value{A list containing two matrices:
  \item{Weighted }{The first matrix (named weighted matrix) contains the abundance of each haplotype per population, represented by the number of haplotypes (columns) found per population (rows). }
  \item{Interaction }{The second matrix (named interaction matrix) contains information about the presence or absence of each haplotype (columns) per population (rows) represented by 1 or 0, respectively.}
}
\author{
A. J. Muñoz-Pajares
}

\seealso{\code{\link{FindHaplo}}
}
\examples{
\dontrun{
# cat("Sequence.Name Haplotype.Name",
# "Population1 H1",
# "Population1 H2",
# "Population1 H3",
# "Population1 H2",
# "Population2 H4",
# "Population2 H5",
# "Population2 H6",
# "Population2 H4",
# "Population3 H7",
# "Population3 H7",
# "Population3 H7",
# "Population3 H7",
#      file = "3_FindHaplo_Example2_modified.txt", sep = "\n")
# 
# # Reading the alignment directly from file:
# 	HapPerPop(inputFile="3_FindHaplo_Example2_modified.txt",header=TRUE,
# saveFile=FALSE)
# 
# cat("Sequence.Name Haplotype.Name",
# "Population1id1 H1",
# "Population1id2 H2",
# "Population1id3 H3",
# "Population1id4 H2",
# "Population2id1 H4",
# "Population2id2 H5",
# "Population2id3 H6",
# "Population2id4 H4",
# "Population3id1 H7",
# "Population3id2 H7",
# "Population3id3 H7",
# "Population3id4 H7",
#      file = "4_FindHaplo_Example2_modified.txt", sep = "\n")
# 
# # Reading the alignment directly from file. First column includes population
# # and individual names. Consequently, 12 populations are considered:
# 	HapPerPop(inputFile="4_FindHaplo_Example2_modified.txt",header=TRUE,
# saveFile=FALSE)
# 
# # Population names within the input matrix first column goes from 
# # character 1 to 11. Now 3 populations are considered:
# 	HapPerPop(inputFile="4_FindHaplo_Example2_modified.txt",header=TRUE,
# saveFile=FALSE,NameIniPopulations=1, NameEndPopulations=11)
# 
# # If population names are set from character 1 to 3, all samples would
# # be treated as a single population
# 	HapPerPop(inputFile="4_FindHaplo_Example2_modified.txt",header=TRUE,
# saveFile=FALSE,NameIniPopulations=1, NameEndPopulations=3)
# 
# # Reading the alignment directly from file, displaying only the 
# # weighted matrix:
# 	HapPerPop(inputFile="4_FindHaplo_Example2_modified.txt",header=TRUE,
# saveFile=FALSE,NameIniPopulations=1, NameEndPopulations=11)[[1]]
# 
# # Reading the alignment from an object and saving the two computed
# # distance matrices:
# 	FH<-read.table("3_FindHaplo_Example2_modified.txt",header=TRUE)
# 	HapPerPop(input=FH,header=TRUE,saveFile=FALSE)
# }
}


