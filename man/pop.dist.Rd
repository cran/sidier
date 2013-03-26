\name{pop.dist}
\alias{pop.dist}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Distances among populations
}
\description{This function computes the among population distance matrix based on the frequency of haplotypes per population and the among haplotypes distance matrix. It is mandatory to define haplotype and population names in the input file. See example for details.
}
\usage{
pop.dist(DistFile = NA, distances = NA, HaploFile = NA, Haplos = NA,
 outType = "O", logfile = TRUE, saveFile = TRUE, NameIniPopulations
 = NA, NameEndPopulations = NA, NameIniHaplotypes = NA,
 NameEndHaplotypes = NA)
}

\arguments{
  \item{DistFile}{the name of the file containing the distance matrix among haplotypes to be analysed.
}
  \item{distances}{the distance matrix among haplotypes (stored in memory) to be analysed.
}
  \item{HaploFile}{the name of the file containing the matrix with the number of haplotypes found per population.
}
  \item{Haplos}{the name of the matrix (stored in memory) containing the number of haplotypes found per population.
}
  \item{outType}{a strig; the format of output matrix. "L" for lower diagonal hemi-matrix; "7" for upper diagonal hemi-matrix; "O" for both hemi-matrices (default).
}
  \item{logfile}{a logical; if TRUE (default), it saves a file containing matrix names used (inputDist and HaploFile)
}
  \item{saveFile}{a logical; if TRUE (default), function output is saved as a text file.
}
  \item{NameIniPopulations}{a numeric indicating the position of the initial character of population names within the individual name in the matrix containing the number of haplotypes found per population (see exaple for details).
}
  \item{NameEndPopulations}{a numeric indicating the position of the last character of population names within the individual name in the matrix containing the number of haplotypes found per population (see exaple for details).
}
  \item{NameIniHaplotypes}{a numeric indicating the position of the initial character of haplotype names within the individual name in the distance matrix (see exaple for details).
}
  \item{NameEndHaplotypes}{a numeric indicating the position of the last character of haplotype names within the individual name in the distance matrix (see exaple for details).
}
}
\details{Each element in the population distance matrix is calculated as the arithmetic mean of the distances among all the sequences sampled in the two compared populations, as follows:
\deqn{dist(i,j) = \frac{\sum_{k=1}^{m}{\sum_{l=1}^{n}{dist(H_{ki},H_{lj})}}}{m*n}}
where \emph{dist(i,j)} represents the distance between populations \emph{i} and \emph{j}, \emph{m} and \emph{n} are the number of sequences in populations \emph{i} and \emph{j}, respectively, and \emph{dist(H_{ki},H_{lj})} is the distance between the \emph{k-th} sequence found in population \emph{i} and the \emph{l-th} sequence found in population \emph{j}. 
}
\value{A matrix containing the genetic distances among populations, based on the haplotype distances and their frequencies per populations.

}
\author{A.J. Munoz-Pajares
}
\examples{
cat(" H1 H2 H3 H4 H5",
"Population1 1 2 1 0 0",
"Population2 0 0 0 4 1",
"Population3 0 1 0 0 3",
     file = "4_Example3_HapPerPop_Weighted.txt", sep = "\n")


cat("H1 H2 H3 H4 H5",
"H1 0 1 2 3 1",
"H2 1 0 3 4 2",
"H3 2 3 0 1 1",
"H4 3 4 1 0 2",
"H5 1 2 1 2 0",
     file = "4_Example3_IndelDistanceMatrixMullerMod.txt", sep = "\n")
     example3_2 <- read.table("4_Example3_IndelDistanceMatrixMullerMod.txt"
,header=TRUE)

# Checking row names to estimate NameIniHaplotypes,NameEndHaplotypes:
 row.names(read.table(file="4_Example3_IndelDistanceMatrixMullerMod.txt"))
## [1] "H1" "H2" "H3" "H4" "H5" NameIniHaplotypes=1. NameEndHaplotypes=2
# Checking row names to estimate NameIniPopulations, and NameEndPopulations
 row.names(read.table(file="4_Example3_HapPerPop_Weighted.txt"))
## [1] "Population1" "Population2" "Population3"
## NameIniPopulations=1 NameEndPopulations =11


# Reading files. Distance matrix must contain haplotype names. Abundance
# matrix must contain both, haplotype and population names:

pop.dist (DistFile="4_Example3_IndelDistanceMatrixMullerMod.txt", 
HaploFile="4_Example3_HapPerPop_Weighted.txt", outType="O",
NameIniHaplotypes=1,NameEndHaplotypes=2,NameIniPopulations=1,
NameEndPopulations=11)

}
