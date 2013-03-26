\name{nt.gap.comb}
\alias{nt.gap.comb}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{substitution and indel distance combinations
}
\description{This function obtains a lineal combination from two original matrices. The weight of each matrix in the combination must be defined. If it is a range of values, several matrices are computed.
}
\usage{
nt.gap.comb(DISTnuc = NA, DISTgap = NA, alpha = seq(0, 1, 0.1),
method = "Corrected", saveFile = TRUE)
}
\arguments{
  \item{DISTnuc}{a matrix containing substitution genetic distances. See "dist.dna" in "ape" package.
}
  \item{DISTgap}{a matrix containing indel genetic distances. See MCMC function in this package.
}
  \item{alpha}{a numeric between 0 and 1, is the weight given to the indel genetic distance matrix in the combination. By definition, the weight of the substitution genetic matrix is the complementary value (i.e., 1-alpha). It is also possible to define multiple weights to estimate different combinations.}
  \item{method}{a string defining whether each distance matrix must be divided by its maximum value before the combination ("Corrected") or not ("Uncorrected"). Consequently, if the "Corrected" method is chosen, both matrices will range between 0 and 1 before to be combined.
}
  \item{saveFile}{a logical; if TRUE (default), each ouput matrix is saved in a different text file.
}}
\value{If "alpha" is a single value, this function generates a data frame containing the estimated combination of substitution and indel distance matrices. If "alpha" is a vector of values, this function generates a list of data frames.
}
\author{A.J. Munoz-Pajares
}

\seealso{MCIC
}
\examples{
cat(">Population1_sequence1",
"TTATAAAATCTA----TAGC",
">Population1_sequence2",
"TAAT----TCTA----TAAC",
">Population1_sequence3",
"TTATAAAAATTA----TAGC",
">Population1_sequence4",
"TAAT----TCTA----TAAC",
">Population2_sequence1",
"TTAT----TCGAGGGGTAGC",
">Population2_sequence2",
"TAAT----TCTA----TAAC",
">Population2_sequence3",
"TTATAAAA--------TAGC",
">Population2_sequence4",
"TTAT----TCGAGGGGTAGC",
">Population3_sequence1",
"TTAT----TCGA----TAGC",
">Population3_sequence2",
"TTAT----TCGA----TAGC",
">Population3_sequence3",
"TTAT----TCGA----TAGC",
">Population3_sequence4",
"TTAT----TCGA----TAGC",
     file = "ex2.fas", sep = "\n")

 # Estimating indel distances after reading the alignment from file:
distGap<-MCIC(input="ex2.fas",saveFile=FALSE)
 # Estimating substitution distances after reading the alignment from file:
library(ape)
align<-read.dna(file="ex2.fas",format="fasta")
dist.nt<-dist.dna(align,model="raw",pairwise.deletion=TRUE)
DISTnt<-as.matrix(dist.nt)
 # Obtaining 11 corrected combined matrices using a range of alpha values:
nt.gap.comb(DISTgap=distGap, alpha=seq(0,1,0.1), method="Corrected", 
saveFile=FALSE, DISTnuc=DISTnt)
 # Obtaining the arithmetic mean of both matrices using both the corrected
 # and the uncorrected methods.
nt.gap.comb(DISTgap=distGap, alpha=0.5, method="Uncorrected", saveFile=FALSE,
 DISTnuc=DISTnt)
 # Obtaining a range of combinations...
Range01<-nt.gap.comb(DISTgap=distGap, alpha=seq(0,1,0.1), method="Uncorrected",
 saveFile=FALSE, DISTnuc=DISTnt)
 # ...and displaying the arithmetic mean (alpha=0.5 is the element number 6
 # in the resulting data frame):
Range01[[6]]
}
