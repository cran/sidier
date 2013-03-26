\name{FindHaplo}
\alias{FindHaplo}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Find equal haplotypes
}
\description{This function assigns the same name to equal haplotypes in a sequence alignment.
}
\usage{
FindHaplo(inputFile = NA, align = NA, saveFile = T,
outname = "FindHaplo.txt")
}
\arguments{
  \item{inputFile}{the name of a sequence alignment file in fasta format to be analysed.
}
  \item{align}{the name of a sequence alignment stored in memory to be analysed. See "read.dna" in ape  package for details about reading alignments.
}
  \item{saveFile}{a logical; if TRUE (default), function output is saved as a text file.
}
  \item{outname}{if "saveFile" is set to TRUE (default), contains the name of the output file ("FindHaplo.txt" by default).
}
}
\value{A two columns matrix containing the original sequence name and the haplotype name assigned to each sequence in the input alignment.
}
\author{A.J. Munoz-Pajares
}
\seealso{HapPerPop
}
\examples{

	#generating an alignment file:
cat(">Population1_sequence1",
"TTATAGGTAGCTTCGATATTG",
">Population2_sequence1",
"TTA---GTAGCTTCGAAATTG",
">Population3_sequence1",
"TTA---GTA---TCG---TAG",
">Population4_sequence1",
"TTATAGGTA---TCG---TTG",
">Population5_sequence1",
"TTA------------AAATTG",
file = "ex1.fas", sep = "\n")

	 # Reading the alignment directly from file:
	 FindHaplo(input="ex1.fas") 

	 # Reading the alignment from an object:
library(ape)
example1<-read.dna(file="ex1.fas",format="fasta")
FindHaplo(align=example1)


}

