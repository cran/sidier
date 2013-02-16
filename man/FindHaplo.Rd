\name{FindHaplo}
\alias{FindHaplo}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Find equal haplotypesFind equal haplotypes
%%  function to do ... 
}
\description{This function assigns the same name to equal haplotypes in a sequence alignment.
%%   A concise (1-5 lines) description of what the function does. 
}
\usage{
FindHaplo(readfile = T, input = NA, align = NA, saveFile = T, outname = "FindHaplo.txt")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{readfile}{a logical; if TRUE (default), the input alignment is provided as a fasta format in a text file. If FALSE, the alignment is provided as an R object.
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
  \item{outname}{if "SaveFile" is set to TRUE (default), contains the name of the output file ("FindHaplo.txt" by default).
%%     Describe \code{outname} here
}
}
\value{A matrix showing the assigned haplotype name to each sequence in the alignment.
%%  Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\author{A.J. Munoz-Pajares
%%  who you are
}

%% Make other sections like Warning with \section{Warning }{....} 

\seealso{HapPerPop
%% objects to See Also as \code{\link{help}}, 
}
\examples{

cat(">Population1_sequence1",
"TTATAGCTGTCGGGCTAGTAGCTGTATCAGTCGTACGTAGTAGTCGTGTCGATCGATGGCGCGGCGCATC--------------------TAGCGCTAGCTGATGCTAGTAGCGTAGAGTATG",
">Population2_sequence1",
"TTATAGCTGTCGGGCTAGTAGCTGTATCAGTC--------------------TCGATGGCGCGGCGCATCAATATTATATCGGCGATGCGTAGCGCTAGCTGATGCTAGTAGCGTAGAGTATG",
">Population3_sequence1",
"TTATAGCTGTCGGGCTAGTAGCTGTATCAGTC-------GTA----------TCGATGGCGCGGCGCATC--------------------TAGCGCTAGCTGATGCTAGTAGCGTAGAGTATG",
">Population4_sequence1",
"TTATAGCTGTCGGGCTAGTAGCTGTATCAGTC----------------------------------------------------------TAGCGCTAGCTGATGCTAGTAGCGTAGAGTATG",
">Population5_sequence1",
"TTATAGCTGTCGGGCTAGTAGCTGTATCAGTCGTACGTA--------------------------------------TATCGGCGATGCGTAGCGCTAGCTGATGCTAGTAGCGTAGAGTATG",
     file = "ex1.fas", sep = "\n")
	 example1 <- read.dna("ex1.fas", format = "fasta")
     # Reading the alignment from an object:
	 FindHaplo(readfile=FALSE,align=example1)
	 # Reading the alignment directly from file:
	 FindHaplo(input="ex1.fas") 

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ kwd1 }
%\keyword{ kwd2 }% __ONLY ONE__ keyword per line
