\name{mutationSummary}
\alias{mutationSummary}
\title{
Summary of observed mutations
}
\description{
This function computes the number of substitutions and indels observed in a given alignment.
}
\usage{
mutationSummary(align, addExtremes = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{align}{
the name of the alignment to be analysed. See "read.dna" in ape  package for details about reading alignments.
}
  \item{addExtremes}{
a logical; if TRUE, additional nucleotide sites are included in both extremes of the alignment. This will allow estimating distances for alignments showing gaps in terminal positions.
}
}
\details{
It is recommended to estimate this distance matrix after removing repeated sequences from the alignment. Repeated sequences increase computation time but do not provide additional information (because they produce duplicated rows and columns in the final distance matrix).
}
\value{
\item{Sites }{A matrix containing: the number of sites per sequence (Length); the number of constant and variable sites; the number of singletons and informative sites.}
\item{Events }{A matrix containing: the number of substitution (singletons, informative, and total) and indel (singletons, informative, and total) events}
}
\author{
A.J. Munoz-Pajares
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

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
mutationSummary (align=read.dna("ex3.fas",format="fasta"))

}
