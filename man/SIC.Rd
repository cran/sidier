\name{SIC}
\alias{SIC}
\title{Indel distances following the Simple Index Coding method
}
\description{
This function codifies gapped positions in a sequence alignment following the rationale of the method described by Simmons and Ochoterrena (2000). Based on the yielded indel coding matrix, this function also computes a pairwise indel distance matrix.
}
\usage{
SIC(inputFile = NA, align = NA, saveFile = T,
outnameDist=paste(inputFile,"IndelDistanceSIC.txt",
sep = "_"), outnameCode = paste(inputFile,
"SIC_coding.txt", sep = "_"), addExtremes = F)
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
  \item{outnameDist}{
if "saveFile" is set to TRUE (default), contains the name of the distance output file.
}
  \item{outnameCode}{
if "saveFile" is set to TRUE (default), contains the name of the indel coding output file.

}
  \item{addExtremes}{
a logical; if TRUE, additional nucleotide sites are included in both extremes of the alignment. This will allow estimating distances for alignments showing gaps in terminal positions.
}
}
\details{
It is recommended to estimate this distance matrix after removing repeated sequences from the alignment. Repeated sequences increase computation time but do not provide additional information (because they produce duplicated rows and columns in the final distance matrix).
}
\value{
A list with two elements:
\item{indel conding matrix }{Describes the initial and final site of each gap and its presence or absence per sequence.
}
\item{distance matrix }{Contains genetic distances based on comparing indel presence/ausence between sequences.
}
}
\references{Simmons, M., Ochoterena, H. & Carr, T. (2001). Incorporation, relative homoplasy, and effect of gap characters in sequence-based phylogenetic analyses. Systematic Biology, 50, 454-62.
}
\author{
A.J. Munoz-Pajares
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
BARRIEL, MCIC, FIFTH
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
library(ape)
SIC (align=read.dna("ex3.fas",format="fasta"), saveFile = FALSE)

# Analysing the same dataset, but using only unique sequences:
uni<-GetHaplo(inputFile="ex3.fas",saveFile=FALSE)
SIC (align=uni, saveFile = FALSE)
}
