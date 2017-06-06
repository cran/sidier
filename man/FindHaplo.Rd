\name{FindHaplo}
\alias{FindHaplo}
\title{Find equal haplotypes
}
\description{This function assigns the same name to equal haplotypes in a sequence alignment.
}
\usage{
FindHaplo(inputFile = NA, align = NA,
saveFile = TRUE, outname = "FindHaplo.txt")
}
\arguments{
  \item{inputFile}{the name of a sequence alignment file in fasta format to be analysed. Alternatively you can provide the name of a "DNAbin" class alignment stored in memory using the "align" option.
}
  \item{align}{the name of the "DNAbin" alignment to be analysed. See "?read.dna" in the ape  package for details about reading alignments. Alternatively you can provide the name of the file containing the alignment in fasta format using the "inputFile" option.
}
  \item{saveFile}{a logical; if TRUE (default), function output is saved as a text file.
}
  \item{outname}{if "saveFile" is set to TRUE (default), contains the name of the output file ("FindHaplo.txt" by default).
}
}
\details{The algorithm identifies identical sequences even if they are wrongly aligned (see example).}
\value{A two columns matrix containing the original sequence name and the haplotype name assigned to each sequence in the input alignment.
}
\author{A. J. Muñoz-Pajares
}
\seealso{\code{\link{GetHaplo}}, \code{\link{HapPerPop}}
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
	 FindHaplo(inputFile="ex1.fas") 

	 # Reading the alignment from an object:
library(ape)
example1<-read.dna(file="ex1.fas",format="fasta")
FindHaplo(align=example1)

	#generating a new alignment file with identical sequences wrongly aligned:
cat(">Pop1_seq1",
"TTATTCTA--------TAGC",
">Pop1_seq2",
"TTAT----TCTA----TAGC",
">Pop1_seq3",
"TAAT----TCTA------AC",
     file = "ex2.2.fas", sep = "\n")

# Note that seq1 and seq2 are equal, but the alignment is different.
# However, this function identifies seq1 and seq2 as identical:
FindHaplo(inputFile="ex2.2.fas")

}

