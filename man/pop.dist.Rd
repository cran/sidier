\name{pop.dist}
\alias{pop.dist}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Distances among populations
%%  function to do ... 
}
\description{This function computes the among population distance matrix based on the frequency of haplotypes per population and the among haplotypes distance matrix. It is mandatory to define haplotype and population names in the input file. See example for details
%%   A concise (1-5 lines) description of what the function does. 
}
\usage{
pop.dist(DistFile = T, inputDist = NA, distances = NA, HaploFile = T, inputHaplo = NA, Haplos = NA, outType = "O", logfile = TRUE, saveFile = TRUE, NameIni = NA, NameEnd = NA)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{DistFile}{a logical; if TRUE (default) input distance matrix among haplotypes is provided as a matrix in a text file. If FALSE, the matrix must be provided as an R object.

%%     Describe \code{DistFile} here
}
  \item{inputDist}{the name of the file containing the distance matrix among haplotypes.

%%     Describe \code{inputDist} here
}
  \item{distances}{the name of the distance matrix among haplotypes to be analysed (if "DistFile" is set to FALSE,).
%%     Describe \code{distances} here
}
  \item{HaploFile}{a logical; if TRUE (default) the input matrix containing the number of haplotypes found per population is provided as a matrix in a text file. If FALSE, the matrix must be provided as an R object. See HapPerPop for details on how to estimate such matrix.
%%     Describe \code{HaploFile} here
}
  \item{inputHaplo}{the name of the file containing the matrix with the number of haplotypes found per population.
%%     Describe \code{inputHaplo} here
}
  \item{Haplos}{the name of the matrix containing the number of haplotypes found per population (if "DistFile" is set to FALSE,).
%%     Describe \code{Haplos} here
}
  \item{outType}{a strig; the format of output matrix. "L" for lower diagonal hemi-matrix; "7" for upper diagonal hemi-matrix; "O" for both hemi-matrices (default).
%%     Describe \code{outType} here
}
  \item{logfile}{a logical; if TRUE (default), it saves a file containing matrix names used (inputDist and HaploFile)
%%     Describe \code{logfile} here
}
  \item{saveFile}{a logical; if TRUE (default), function output is saved as a text file.
%%     Describe \code{saveFile} here
}
  \item{NameIni}{a numeric indicating the position of the initial character of population name within the individual name in the distance matrix.
%%     Describe \code{NameIni} here
}
  \item{NameEnd}{a numeric indicating the position of the last character of population name within the individual name in the distance matrix.
%%     Describe \code{NameEnd} here
}
}
\value{A matrix containing the genetic distances among populations, based on the haplotype distances and their frequencies per populations.
%%  Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\author{A.J. Munoz-Pajares
%%  who you are
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
     example3_2 <- read.table("4_Example3_IndelDistanceMatrixMullerMod.txt",header=TRUE)

cat(">Population1_sequence1",
"GGGGAGCTGTCGGGCTAGTAGCTGTATCAGTCGTACGTAGTAGTCGTGTCGATCGATGGCGCGGCGCATC--------------------TAGCGCTAGCTGATGCTAGTAGCGTAGAGTATG",
">Population1_sequence2",
"GGGGAGCTGTCGGGCTAGTAGCTGTATCAGTCGTACGTAGTAGTCGTGTCGATCGATGGCGCGGCGCATC--------------------TAGCGCTAGCTGATGCTAGTAGCGTAGAGTATG",
">Population1_sequence3",
"GGGGAGCTGTCGGGCTAGTAGCTGTATCAGTCGTACGTAGTAGTCGTGTCGATCGATGGCGCGGCGCATC--------------------TAGCGCTAGCTGATGCTAGTAGCGTAGAGTATG",
">Population1_sequence4",
"GGGGAGCTGTCGGGCTAGTAGCTGTATCAGTCGTACGTAGTAGTCGTGTCGATCGATGGCGCGGCGCATC--------------------TAGCGCTAGCTGATGCTAGTAGCGTAGAGTATG",
">Population2_sequence1",
"TTATAGCTGTCGGGCTAGTAGCTGTATCAGTC--------------------TCGATGGCGCGGCGCATCAATATTATATCGGCGATGCGTAGCGCTAGCTGATGCTAGTAGCGTAGAGTATG",
">Population2_sequence2",
"TTATAGCTGTCGGGCTAGTAGCTGTATCAGTC--------------------TCGATGGCGCGGCGCATCAATATTATATCGGCGATGCGTAGCGCTAGCTGA----------GTAGAGTATG",
">Population2_sequence3",
"TTATAGCTGTCGGGCTAGTAGCTGTATCAGTC--------------------TCGATGGCGCGGCGCATCAATATTATATCGGCGATGCGTAGCGCTAGCTGATGCTAGTAGCGTAGAAAAAA",
">Population2_sequence4",
"TTATAGCTGTCGGGCTAGTAGCTGTATCAGTC--------------------TCGATGGCGCGGCGCATCAATATTATATCGGCGATGCGTAGCGCTAGCTGATGCTAGTAGCGTAGAGTATG",
">Population3_sequence1",
"GGGGAGCTGTCGGGCTAGTAGCTGTATCAGTCGTACGTAGTAGTCGTGTCGATCGATGGCGCGGCGCATC--------------------TAGCGCTAGCTGATGCTAGTAGCGTAGAGTATG",
">Population3_sequence2",
"GGGGAGCTGTCGGGCTAGTAGCTGTATCAGTCGTACGTAGTAGTCGTGTCGATCGATGGCGCGGCGCATC--------------------TAGCGCTAGCTGATGCTAGTAGCGTAGAGTATG",
">Population3_sequence3",
"GGGGAGCTGTCGGGCTAGTAGCTGTATCAGTCGTACGTAGTAGTCGTGTCGATCGATGGCGCGGCGCATC--------------------TAGCGCTAGCTGATGCTAGTAGCGTAGAGTATG",
">Population3_sequence4",
"GGGGAGCTGTCGGGCTAGTAGCTGTATCAGTCGTACGTAGTAGTCGTGTCGATCGATGGCGCGGCGCATC--------------------TAGCGCTAGCTGATGCTAGTAGCGTAGAGTATG",
     file = "5_Example4.fas", sep = "\n")
     example4 <- read.table("5_Example4.fas",header=TRUE)

# Reading files. Distance matrix must contain haplotype names. Abundance matrix must contain both, haplotype and population names:

pop.dist (DistFile=TRUE, inputDist="4_Example3_IndelDistanceMatrixMullerMod.txt", HaploFile=TRUE, inputHaplo="4_Example3_HapPerPop_Weighted.txt", outType="O", logfile=FALSE, saveFile=FALSE)

# It may be convenient to manually modify files to get the appropriate names. However, an automated example from fasta sequence names (using the substr function) is shown below:

# Estimating distances between unique haplotypes
uniqueHaplo<-GetHaplo(input="5_Example4.fas",saveFile=FALSE)
distGap<-MCIC(readfile=FALSE,align=uniqueHaplo,saveFile=FALSE)
dist.nt <-dist.dna(uniqueHaplo,model="raw",pairwise.deletion=TRUE)
DISTnt<-as.matrix(dist.nt)

# Replacing sequence names by haplotype names in both distance matrices
for (Hi in 1:length(colnames(distGap)))
colnames(distGap)[Hi]<-FindHaplo(input="5_Example4.fas",saveFile=FALSE)[which(colnames(distGap)[Hi]==FindHaplo(input="5_Example4.fas",saveFile=FALSE)[,1]),2]
row.names(distGap)<-colnames(distGap)
for (Hi in 1:length(colnames(DISTnt)))
colnames(DISTnt)[Hi]<-FindHaplo(input="5_Example4.fas",saveFile=FALSE)[which(colnames(DISTnt)[Hi]==FindHaplo(input="5_Example4.fas",saveFile=FALSE)[,1]),2]
row.names(DISTnt)<-colnames(DISTnt)

#Combining distance matrices and setting haplotype names
CombinedDistance<-as.data.frame(nt.gap.comb(DISTgap=distGap, range=0.5, method="Corrected", saveFile=FALSE, DISTnuc=DISTnt)[[2]])
colnames(CombinedDistance)<-row.names(CombinedDistance)

# Estimating haplotype abundance per population and setting population names:
Haplotypes<-FindHaplo(input="5_Example4.fas",saveFile=FALSE)
Haplotypes[,1]<-substr(Haplotypes[,1],1,11)
Weighted<-as.data.frame(HapPerPop(readfile=FALSE,header=TRUE,input=Haplotypes)[1])
colnames(Weighted)<-substr(colnames(Weighted),10,11)

# Estimating population distances
pop.dist (DistFile=FALSE, distances=CombinedDistance, HaploFile=FALSE, Haplos=Weighted, outType="O", logfile=FALSE, saveFile=FALSE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ kwd1 }
%\keyword{ kwd2 }% __ONLY ONE__ keyword per line
