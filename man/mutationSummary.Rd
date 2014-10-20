\name{mutationSummary}
\alias{mutationSummary}
\title{
Summary of observed mutations
}
\description{
This function computes the number of substitutions and indels observed in a given alignment.
}
\usage{
mutationSummary(align, addExtremes = FALSE, output = "brief")
}
\arguments{
  \item{align}{
the name of the "DNAbin" alignment to be analysed. See "?read.dna" in the ape  package for details about reading alignments.}
  \item{addExtremes}{
a logical; if TRUE, additional nucleotide sites are included in both extremes of the alignment. This will allow estimating distances for alignments showing gaps in terminal positions.}
  \item{output}{a string; defines the kind of output. Two values are accepted:
- "brief" (default) produces an output showing the number of mutations (sites and events).
- "detailed" produces an output showing the number of mutations (sites and events), the position of each mutation, and the state of these sites per sequence (A, T, C, G or - for substitutions and 1 or 0 for indels).}
}
\value{A list containing:
	\item{Sites }{A matrix containing: the number of sites per sequence (Length); the number of constant and variable sites; the number of singletons and informative sites.}
	\item{Events }{A matrix containing: the number of substitution (singletons, informative, and total) and indel (singletons, informative, and total) events}
	\item{Constants.Alignment}{A matrix showing constant sites in the alignment (Only shown if output=="detailed").}
	\item{Variables.Alignment}{A matrix showing variable sites in the alignment (Only shown if output=="detailed").}
	\item{Singletons.Alignment}{A matrix showing singleton sites in the alignment (Only shown if output=="detailed").}
	\item{Inforatives.Alignment}{A matrix showing informative sites in the alignment (Only shown if output=="detailed").}
	\item{Substitutions}{A matrix showing substitution sites in the alignment (Only shown if output=="detailed").}
	\item{Subst.Single}{A matrix showing singleton substitution sites in the alignment (Only shown if output=="detailed").}
	\item{Subst.Info}{A matrix showing informative substitution sites in the alignment (Only shown if output=="detailed").}
	\item{Gaps}{A matrix showing gap sites in the alignment (Only shown if output=="detailed").}
	\item{Gaps.Single}{A matrix showing singleton gap sites in the alignment (Only shown if output=="detailed").}
	\item{Gaps.Info}{A matrix showing informative gap sites in the alignment (Only shown if output=="detailed").}
}
\author{A. J. Muñoz-Pajares}

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
mutationSummary (align=read.dna("ex3.fas",format="fasta"),output="detailed")

#A more complex alignment
data(ex_alignment1) # this will read a fasta file with the name 'alignExample'
mutationSummary(align=alignExample,addExtremes=TRUE)

}
