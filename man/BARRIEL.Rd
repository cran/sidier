\name{BARRIEL}
\alias{BARRIEL}
\title{Indel distances following Barriel method
}
\description{
This function codifies gapped positions in a sequence alignment following the rationale of the method described by Barriel (1994). Based on the yielded indel coding matrix, this function also computes a pairwise indel distance matrix.
}
\usage{
BARRIEL(inputFile = NA, align = NA, saveFile = T,
outnameDist = paste(inputFile, "IndelDistanceBarriel.txt",
sep = "_"), outnameCode = paste(inputFile, "Barriel_coding.txt",
 sep = "_"), addExtremes = F)
}
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
\references{
Barriel, V., 1994. Molecular phylogenies and how to code insertion/
deletion events. Life Sci. 317, 693-701.
}
\author{
A.J. Munoz-Pajares
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
MCIC, SIC, FIFTH
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

BARRIEL (align=read.dna("ex3.fas",format="fasta"), saveFile = FALSE)

# Analysing the same dataset, but using only unique sequences:
uni<-GetHaplo(inputFile="ex3.fas",saveFile=FALSE)
BARRIEL (align=uni, saveFile = FALSE)
}
