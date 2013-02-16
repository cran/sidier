\name{MCIC}
\alias{MCIC}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Modified Complex Indel Coding as distance matrix
%%  function to do ... 
}
\description{This function computes the insertion-deletion (indel) distance matrix following the rationale of the Modified Complex Indel Coding (Muller, 2006) to estimate transition matrices, as described in Munoz-Pajares. 
%%   A concise (1-5 lines) description of what the function does. 
}
\usage{
MCIC(readfile = T, input = NA, align = NA, saveFile = T, outname = paste(input, "IndelDistanceMatrixMullerMod.txt"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{readfile}{a logical; if TRUE (default) input alignment is provided as a fasta format in a text file. If FALSE, the alignment is provided as an R object.
%%     Describe \code{readfile} here
}
  \item{input}{the name of the fasta file to be analysed.
%%     Describe \code{input} here
}
  \item{align}{the name of the alignment to be analysed (if "readfile" is set to FALSE,). See "read.dna" in ape  package for details about reading alignments.
%%     Describe \code{align} here
}
  \item{saveFile}{a logical; if TRUE (default), function output is saved as a text file.
%%     Describe \code{saveFile} here
}
  \item{outname}{if "SaveFile" is set to TRUE (default), contains the name of the output file.
%%     Describe \code{outname} here
}
}
\value{A matrix containing the genetic distances estimated as indels pairwise differences.
%%  Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
Muller K. (2006). Incorporating information from length-mutational events into phylogenetic analysis. Molecular Phylogenetics and Evolution, 38, 667-676.

Munoz-Pajares, AJ. SIDIER: Substitution and Indel Distances to Infer Evolutionary Relationships.
%% put references to the literature/web site here 
}
\author{A.J. Munoz-Pajares
%%  who you are
}

\examples{

cat(">Population1_sequence1",
"TTATAGCTGTCGGGCTAGTAGCTGTATCAGTCGTACGTAGTAGTCGTGTCGATCGATGGCGCGGCGCATC--------------------TAGCGCTAGCTGATGCTAGTAGCGTAGAGTATG",
">Population1_sequence2",
"TTATAGCTGTCGGGCTA------GTATCAGTCGTACGTAGTAGTCGTGTCGATCGATGGCGCGGCGCATC--------------------TAGCGCTAGCTGATGCTAGTAGCGTAGAGTATG",
">Population1_sequence3",
"GGGGAGCTGTCGGGCTAGTAGCTGTATCAGTCGTACGTAGTAGTCGTGTCGATCGATGGCGCGGCGCATC--------------------TAGCGCTAGCTGATGCTAGTAGCGTAGAGTATG",
">Population1_sequence4",
"TTATAGCTGTCGGGCTA------GTATCAGTCGTACGTAGTAGTCGTGTCGATCGATGGCGCGGCGCATC--------------------TAGCGCTAGCTGATGCTAGTAGCGTAGAGTATG",
">Population2_sequence1",
"TTATAGCTGTCGGGCTAGTAGCTGTATCAGTC--------------------TCGATGGCGCGGCGCATCAATATTATATCGGCGATGCGTAGCGCTAGCTGATGCTAGTAGCGTAGAGTATG",
">Population2_sequence2",
"TTATAGCTGTCGGGCTAGTAGCTGTATCAGTC--------------------TCGATGGCGCGGCGCATCAATATTATATCGGCGATGCGTAGCGCTAGCTGA----------GTAGAGTATG",
">Population2_sequence3",
"TTATAGCTGTCGGGCTAGTAGCTGTATCAGTC--------------------TCGATGGCGCGGCGCATCAATATTATATCGGCGATGCGTAGCGCTAGCTGATGCTAGTAGCGTAGAAAAAA",
">Population2_sequence4",
"TTATAGCTGTCGGGCTAGTAGCTGTATCAGTC--------------------TCGATGGCGCGGCGCATCAATATTATATCGGCGATGCGTAGCGCTAGCTGATGCTAGTAGCGTAGAGTATG",
">Population3_sequence1",
"TTATAGCTGTCGGGCTAGTAGCTGTATCAGTC--------------------TCGATGGCGCGGCGCATC--------------------TAGCGCTAGCTGATGCTAGTAGCGTAGAGTATG",
">Population3_sequence2",
"TTATAGCTGTCGGGCTAGTAGCTGTATCAGTC--------------------TCGATGGCGCGGCGCATC--------------------TAGCGCTAGCTGATGCTAGTAGCGTAGAGTATG",
">Population3_sequence3",
"TTATAGCTGTCGGGCTAGTAGCTGTATCAGTC--------------------TCGATGGCGCGGCGCATC--------------------TAGCGCTAGCTGATGCTAGTAGCGTAGAGTATG",
">Population3_sequence4",
"TTATAGCTGTCGGGCTAGTAGCTGTATCAGTC--------------------TCGATGGCGCGGCGCATC--------------------TAGCGCTAGCTGATGCTAGTAGCGTAGAGTATG",
     file = "ex2.fas", sep = "\n")

# Reading the alignment directly from file and saving no output file:
MCIC (input="ex2.fas", saveFile = FALSE)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ kwd1 }
%\keyword{ kwd2 }% __ONLY ONE__ keyword per line
