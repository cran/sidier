\name{GetHaplo}
\alias{GetHaplo}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Get sequences of unique haplotypes
%%  function to do ... 
}
\description{This function returns the subset of unique sequences composing a given alignment. 
%%   A concise (1-5 lines) description of what the function does. 
}
\usage{
GetHaplo(readfile = T, input = NA, align = NA, saveFile = T, outname = "Haplotypes.txt", format = "fasta", seqsNames = NA)
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
  \item{outname}{if "SaveFile" is set to TRUE (default), contains the name of the output file ("Haplotypes.txt" by default).
%%     Describe \code{outname} here
}
  \item{format}{format of the DNA sequences to be saved: "interleaved", "sequential", or "fasta" (default). See "write.dna" in ape  package for details.
%%     Describe \code{format} here
}
  \item{seqsNames}{names for each DNA sequence saved: Three choices are possible: if  n unique sequences are found, "Inf.Hap" assign names from H1 to Hn (according to input order). The second option is to define a vector containing n names. By default, input sequence names are used.
%%     Describe \code{seqsNames} here
}
}
\details{If two equal sequences are not identically aligned, they will be considered as different haplotypes. To avoid misleading results in uncertain alignments it is recommended to use as input the original unaligned sequences, including gaps after the last nucleotide of short sequences to make all sequence lengths equal.
%%   If necessary, more details than the description above 
}
\value{A file containing unique sequences from the input file.
%%  Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
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
     example2 <- read.dna("ex2.fas", format = "fasta")
	 alin<-read.dna(file="ex2.fas",format="fasta")

 # Reading the alignment from an object and saving haplotypes names as sequential numbers:
GetHaplo(readfile=FALSE,align=alin,outname="Haplotypes_sequentialNames.txt",seqsNames="Inf.Hap")

# Reading the alignment directly from file and saving using sequence input names:
GetHaplo(input="ex2.fas")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ kwd1 }
%\keyword{ kwd2 }% __ONLY ONE__ keyword per line
