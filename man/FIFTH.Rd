\name{FIFTH}
\alias{FIFTH}
\title{Indel distances following the fifth state rationale
}
\description{
This function computes an indel distance matrix following the rationale of the fifth state. For that, each gap within the alignment is treated as an independent mutation event.
}
\usage{
FIFTH(inputFile = NA, align = NA, saveFile = TRUE,
outname = paste(inputFile,"IndelDistanceFifthState.txt",
 sep = "_"), addExtremes = FALSE)
}
\arguments{
  \item{inputFile}{
the name of the fasta file to be analysed. Alternatively you can provide the name of a "DNAbin" class alignment stored in memory using the "align" option.
}
  \item{align}{
the name of the alignment to be analysed. See "read.dna" in ape  package for details about reading alignments. Alternatively you can provide the name of the file containing the alignment in fasta format using the "inputFile" option.
}
  \item{saveFile}{
a logical; if TRUE (default), it produces an output text file containing the resulting distance matrix.
}
  \item{outname}{
if "saveFile" is set to TRUE (default), contains the name of the output file.
}
  \item{addExtremes}{
a logical; if TRUE, additional nucleotide sites are included in both extremes of the alignment. This will allow estimating distances for alignments showing gaps in terminal positions.
}

}
\details{
It is recommended to estimate this distance matrix using only the unique sequences in the
alignment. Repeated sequences increase computation time but do not provide additional 
information (because they produce duplicated rows and columns in the final distance matrix).

It is of critical importance to correctly identify indels homology in the provided alignment.
For this reason, addExtremes is set to false by default, and computation may not be done
unless flanking regions were homologous.
}
\value{
A matrix containing the genetic distances estimated as indels pairwise differences.
}
\author{
A. J. Muñoz-Pajares
}

\seealso{
\code{\link{BARRIEL}}, \code{\link{MCIC}}, \code{\link{SIC}}
}
\examples{

# This will generate an example file in your working directory:
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
