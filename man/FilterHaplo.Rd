\name{FilterHaplo}
\alias{FilterHaplo}
\title{Filter haplotypes by occurrence}
\description{Provides a subset of a original alignment composed only of haplotypes showing the range of occurrences provided}
\usage{
FilterHaplo(inputFile=NA, align=NA, Nmin=0, Nmax=NULL,
 saveFile=FALSE, outname="FilterHaplo.txt")
}
\arguments{
  \item{inputFile}{the name of a sequence alignment file in fasta format to be analysed. Alternatively you can provide the name of a "DNAbin" class alignment stored in memory using the "align" option.
}
  \item{align}{the name of the "DNAbin" alignment to be analysed. See "?read.dna" in the ape  package for details about reading alignments. Alternatively you can provide the name of the file containing the alignment in fasta format using the "inputFile" option.
}
 \item{Nmin}{Minimum occurrence of an haplotype to be included in the subset}
 \item{Nmax}{Maximum occurrence of an haplotype to be included in the subset}
  \item{saveFile}{a logical; if TRUE (default), function output is saved as a text file in fasta format
}
  \item{outname}{if "saveFile" is set to TRUE (default), contains the name of the output file.}
}
\author{A. J. Muñoz-Pajares}
\examples{
# cat(">Population1_sequence1",
# "TTATAAAATCTA----TAGC",
# ">Population1_sequence2",
# "TAAT----TCTA----TAAC",
# ">Population1_sequence3",
# "TTATAAAAATTA----TAGC",
# ">Population1_sequence4",
# "TAAT----TCTA----TAAC",
# ">Population2_sequence1",
# "TTAT----TCGA----TAGC",
# ">Population2_sequence2",
# "TTAT----TCGA----TAGC",
# ">Population2_sequence3",
# "TTAT----TCGA----TAGC",
# ">Population2_sequence4",
# "TTAT----TCGA----TAGC",
# ">Population3_sequence1",
# "TTAT----TCGAGGGGTAGC",
# ">Population3_sequence2",
# "TAAT----TCTA----TAAC",
# ">Population3_sequence3",
# "TTATAAAA--------TAGC",
# ">Population3_sequence4",
# "TTAT----TCGAGGGGTAGC",
#      file = "ex2.fas", sep = "\n")
# library(ape)
# example<-read.dna(file="ex2.fas",format="fasta")
# 
# # Exclude unique haplotypes
# FilterHaplo(align=example,Nmin=2)
# 
# # Include only unique haplotypes
# FilterHaplo(align=example,Nmax=1)
# 
# # Filter haplotypes appearing between 2 and 4 times
# FilterHaplo(align=example,Nmax=4,Nmin=2)
}
