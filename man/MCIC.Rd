\name{MCIC}
\alias{MCIC}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Modified Complex Indel Coding as distance matrix
%%  function to do ... 
}
\description{This function computes the insertion-deletion (indel) distance matrix following the rationale of the Modified Complex Indel Coding (Muller, 2006) to estimate transition matrices, as described in Munoz-Pajares (2013). 
%%   A concise (1-5 lines) description of what the function does. 
}
\usage{
MCIC(inputFile = NA, align = NA, saveFile = T, outname =
                 paste(inputFile, "IndelDistanceMatrixMullerMod.txt"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{inputFile}{the name of the fasta file to be analysed.
}
  \item{align}{the name of the alignment to be analysed. See "read.dna" in ape  package for details about reading alignments.
}
  \item{saveFile}{a logical; if TRUE (default), function output is saved as a text file.
}
  \item{outname}{if "saveFile" is set to TRUE (default), contains the name of the output file.
}
}
\value{A matrix containing the genetic distances estimated as indels pairwise differences.
}
\details{It is recommended to estimate this distance matrix after removing repeated sequences from the alignment. Repeated sequences increase computation time but do not provide additional information (because they produce duplicated rows and columns in the final distance matrix).
}
\references{
Muller K. (2006). Incorporating information from length-mutational events into phylogenetic analysis. Molecular Phylogenetics and Evolution, 38, 667-676.

Munoz-Pajares, AJ. (2013) Erysimum mediohispanicum at the evolutionary crossroad: Phylogeography, phenotype and pollinators. PhD Thesis.
}
\author{A.J. Munoz-Pajares
}

\examples{

cat(">Population1_sequence1",
"AAAGGGTCGCT---G",
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
MCIC (input="ex3.fas", saveFile = FALSE)

# Analysing the same dataset, but using only unique sequences:
uni<-GetHaplo(inputFile="ex3.fas",saveFile=FALSE)
MCIC (align=uni, saveFile = FALSE)

}
