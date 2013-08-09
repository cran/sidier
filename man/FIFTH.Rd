\name{FIFTH}
\alias{FIFTH}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Indel distances following the fifth state rationale
}
\description{
This function computes an indel distance matrix following the rationale of the fifth state. For that, each gap within the alignment is treated as an independent mutation event.
}
\usage{
FIFTH(inputFile = NA, align = NA, saveFile = T,
outname = paste(inputFile,"IndelDistanceFifthState.txt",
 sep = "_"), addExtremes = F)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{inputFile}{
the name of the fasta file to be analysed.
}
  \item{align}{
the name of the alignment to be analysed. See "read.dna" in ape  package for details about reading alignments.
}
  \item{saveFile}{
a logical; if TRUE (default), function output is saved as a text file.
}
  \item{outname}{
if "saveFile" is set to TRUE (default), contains the name of the output file.
}
  \item{addExtremes}{
a logical; if TRUE, additional nucleotide sites are included in both extremes of the alignment. This will allow estimating distances for alignments showing gaps in terminal positions.
}

}
\details{
It is recommended to estimate this distance matrix after removing repeated sequences from the alignment. Repeated sequences increase computation time but do not provide additional information (because they produce duplicated rows and columns in the final distance matrix).
}
\value{
A matrix containing the genetic distances estimated as indels pairwise differences.
}
\author{
A.J. Munoz-Pajares
}

\seealso{
BARRIEL, MCIC, SIC
}
\examples{

cat(">Population1_sequence1",
"A-AGGGTC-CT---G",
">Population1_sequence2",
"TAA---TCGCT---G",
">Population1_sequence3",
"TAAGGGTCGCT---G",
">Population1_sequence4",
"TAA---TCGCT---G",
">Population2_sequence1",
"TTACGGTCG---TTG",
">Population2_sequence2",
"TAA---TCG---TTG",
">Population2_sequence3",
"TAA---TCGCTATTG",
">Population2_sequence4",
"TTACGGTCG---TTG",
">Population3_sequence1",
"TTA---TCG---TAG",
">Population3_sequence2",
"TTA---TCG---TAG",
">Population3_sequence3",
"TTA---TCG---TAG",
">Population3_sequence4",
"TTA---TCG---TAG",
     file = "ex3.fas", sep = "\n")

# Reading the alignment directly from file and saving no output file:
library(ape)
FIFTH (align=read.dna("ex3.fas",format="fasta"), saveFile = FALSE)

# Analysing the same dataset, but using only unique sequences:
uni<-GetHaplo(inputFile="ex3.fas",saveFile=FALSE)
FIFTH (align=uni, saveFile = FALSE)
}
