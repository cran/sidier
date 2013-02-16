\name{nt.gap.comb}
\alias{nt.gap.comb}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{substitution and indel distance combinations
%%  function to do ... 
}
\description{This function obtains a lineal combination from two original matrices. The weight of each matrix in the combination must be defined. If it is a range of values, several matrices are computed.
%%   A concise (1-5 lines) description of what the function does. 
}
\usage{
nt.gap.comb(DISTnuc = NA, DISTgap = NA, range = seq(0, 1, 0.1), method = "Corrected", saveFile = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{DISTnuc}{a matrix containing substitution genetic distances. See "dist.dna" in "ape" package.
%%     Describe \code{DISTnuc} here
}
  \item{DISTgap}{a matrix containing indel genetic distances. See MCMC function in this package.
%%     Describe \code{DISTgap} here
}
  \item{range}{a numeric between 0 and 1, is the weights given to the indel genetic distance matrix in the combination. By definition, the weight of the substitution genetic matrix is the complementary value.
%%     Describe \code{range} here
}
  \item{method}{a string defining whether each distance matrix must be divided by its maximum value before the combination ("Corrected") or not ("Uncorrected"). Consequently, if the "Corrected" method is chosen, both matrices will range between 0 and 1 before to be combined.
%%     Describe \code{method} here
}
  \item{saveFile}{a logical; if TRUE (default), each ouput matrix is saved in a different text file.

%%     Describe \code{saveFile} here
}
}
\value{A list containing the estimated combination of substitution and indel distance matrices.
%%  Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\author{A.J. Munoz-Pajares
%%  who you are
}

\seealso{MCIC
%% objects to See Also as \code{\link{help}}, 
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
     example2 <- read.dna("ex2.fas", format = "fasta")

 # Estimating indel distances after reading the alignment from file:
distGap<-MCIC(input="ex2.fas",saveFile=FALSE)
 # Estimating substitution distances after reading the alignment from file:
align<-read.dna(file="ex2.fas",format="fasta")
dist.nt<-dist.dna(align,model="raw",pairwise.deletion=TRUE)
DISTnt<-as.matrix(dist.nt)
# Obtaining 11 corrected combined matrices using a range of alpha values:
nt.gap.comb(DISTgap=distGap, range=seq(0,1,0.1), method="Corrected", saveFile=FALSE, DISTnuc=DISTnt)
 # Obtaining the arithmetic mean of both matrices using both the corrected and the uncorrected methods.
nt.gap.comb(DISTgap=distGap, range=0.5, method="Both", saveFile=FALSE, DISTnuc=DISTnt)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ kwd1 }
%\keyword{ kwd2 }% __ONLY ONE__ keyword per line
