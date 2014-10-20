\name{GetHaplo}
\alias{GetHaplo}
\title{Get sequences of unique haplotypes
}
\description{This function returns the subset of unique sequences composing a given alignment. 
}
\usage{
GetHaplo(inputFile = NA, align = NA, saveFile = TRUE,
outname = "Haplotypes.txt", format = "fasta", seqsNames = NA, silent = FALSE)
}
\arguments{
  \item{inputFile}{the name of a sequence alignment file in fasta format to be analysed. Alternatively you can provide the name of a "DNAbin" class alignment stored in memory using the "align" option.
}
  \item{align}{the name of the "DNAbin" alignment to be analysed. See "?read.dna" in the ape  package for details about reading alignments. Alternatively you can provide the name of the file containing the alignment in fasta format using the "inputFile" option.
}
  \item{saveFile}{a logical; if TRUE (default), function output is saved as a text file
}
  \item{outname}{if "saveFile" is set to TRUE (default), contains the name of the output file ("Haplotypes.txt" by default).
}
  \item{format}{format of the DNA sequences to be saved: "interleaved", "sequential", or "fasta" (default). See "write.dna" in ape  package for details.
}
  \item{seqsNames}{names for each DNA sequence saved: Three choices are possible: if n unique sequences are found, "Inf.Hap" assigns names from H1 to Hn (according to input order). The second option is to define a vector containing n names. By default, input sequence names are used.}
  \item{silent}{a logical; if TRUE (default), it prints the number of unique sequences found and the name of the output file.
}
}
\details{The algorithm identifies identical sequences even if they are wrongly aligned (see example).
}
\value{A file containing unique sequences from the input file.
}
\author{A. J. Muñoz-Pajares
}

\seealso{\code{\link{FindHaplo}}
}

\examples{
	#generating an alignment file:
cat(">Population1_sequence1",
"TTATAAAATCTA----TAGC",
">Population1_sequence2",
"TAAT----TCTA----TAAC",
">Population1_sequence3",
"TTATAAAAATTA----TAGC",
">Population1_sequence4",
"TAAT----TCTA----TAAC",
">Population2_sequence1",
"TTAT----TCGAGGGGTAGC",
">Population2_sequence2",
"TAAT----TCTA----TAAC",
">Population2_sequence3",
"TTATAAAA--------TAGC",
">Population2_sequence4",
"TTAT----TCGAGGGGTAGC",
">Population3_sequence1",
"TTAT----TCGA----TAGC",
">Population3_sequence2",
"TTAT----TCGA----TAGC",
">Population3_sequence3",
"TTAT----TCGA----TAGC",
">Population3_sequence4",
"TTAT----TCGA----TAGC",
     file = "ex2.fas", sep = "\n")

# Getting unique haplotypes reading the alignment from a file and setting
#haplotype names:
	GetHaplo(inputFile="ex2.fas",outname="ex2_unique.fas",seqsNames=
c("HaploK001","HaploK002","HaploS001","HaploR001","HaploR002","HaploR003"))
# Reading the alignment from an object and using original sequence names:
    library(ape)
    example2 <- read.dna("ex2.fas", format = "fasta")
	GetHaplo(align=example2,outname="Haplotypes_DefaultNames.txt")
# Reading the alignment from an object and using new haplotype names:
	GetHaplo(align=example2,outname="Haplotypes_sequentialNames.txt",
seqsNames="Inf.Hap")


	#generating a new alignment file with identical sequences wrongly alined:
cat(">Pop1_seq1",
"TTATTCTA--------TAGC",
">Pop1_seq2",
"TTAT----TCTA----TAGC",
">Pop1_seq3",
"TAAT----TCTA------AC",
     file = "ex2.2.fas", sep = "\n")

# Note that seq1 and seq2 are equal, but the alignment is different.
# However, this function identifies seq1 and seq2 as identical:
	a<-GetHaplo(inputFile="ex2.2.fas",saveFile=FALSE)

}

