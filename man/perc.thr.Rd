\name{perc.thr}
\alias{perc.thr}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Percolation threshold network
}
\description{This function computes the percolation network following Rozenfeld et al. (2008), as described in Munoz-Pajares (2013). 
}
\usage{
perc.thr(dis, range = seq(0, 1, 0.01), ptPDF = TRUE,
ptPDFname = "PercolatedNetwork.pdf", estimPDF = TRUE,
estimPDFname = "PercThr Estimation.pdf", estimOutfile = TRUE,
estimOutName = "PercThresholdEstimation.txt", appendOutfile 
= TRUE, plotALL = FALSE, bgcol = "white", label.col = "black",
label = colnames(dis), modules = FALSE, moduleCol = NA,
modFileName = "Modules_summary.txt")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dis}{the distance matrix to be represented
}
  \item{range}{a numeric vector between 0 and 1, is the range of thresholds (referred to the maximum distance in a matrix) to be screened (by default, 101 values from 0 to 1).
}
  \item{ptPDF}{a logical, must the percolated network be saved as a pdf file?
}
  \item{ptPDFname}{if ptPDF=TRUE, the name of the pdf file containing the percolation network to be saved ("PercolatedNetwork.pdf", by default)
}
  \item{estimPDF}{a logical, must the percolation threshold estimation be saved as a pdf file?
}
  \item{estimPDFname}{if estimPDF=TRUE (default), defines the name of the pdf file containing the percolation threshold estimation ("PercThr Estimation.pdf" by default).
}
  \item{estimOutfile}{a logical, must the value of <s> for each threshold be saved as a text file?
}
  \item{estimOutName}{if estimOutfile=TRUE (default), contains the name of the text file containing the percolation threshold estimation ("PercThr Estimation.txt" by default).
}
  \item{appendOutfile}{a logical, if estimOutfile=TRUE, it defines whether results must be appended to an existing file with the same name (TRUE) or the existing file must be replaced (FALSE).
}
  \item{plotALL}{a logical, must all the networks calculated during the percolation threshold estimation (defined by "range" option) be saved as different pdf files? (FALSE, by default). If TRUE, for each value defined in threshold, one file is generated. If TRUE, do not set modules=TRUE.
}
  \item{bgcol}{the colour of the background for each node in the network. Can be equal for all nodes (if only one colour is defined), customized (if several colours are defined), or can represent different modules (see "modules" option).
}
  \item{label.col}{the colour of labels for each node in the network. Can be equal for all nodes (if only one colour is defined) or customized (if several colours are defined).
}
  \item{label}{vector of strings, labels for each node. By default are the column names of the distance matrix (dis). (See "substr" function in base package to automatically reduce name lengths).
}
  \item{modules}{a logical, must nodes belonging to different modules be represented as different colours?
}
  \item{moduleCol}{(if modules=TRUE) vector of strings, defining the colour of nodes belonging to different modules in the network. If TRUE, do not set plotALL=TRUE
}
  \item{modFileName}{(if modules=TRUE) the name of the file to be generated containing a summary of module results (sequence name, module, and colour in network)
}
}
\details{By default, percolation threshold is estimated with an accuracy of 0.01, but it may be increased by setting the decimal places in threshold function (e.g., range=seq(0,1,0.0001)). However, it may strongly increase computation times (in this example, it is required to estimate 100,001 instead of 101 networks). It is also possible to increase accuracy with a low increase in computation time by repeating the process and increasing decimal places only in a range close to a previously estimated percolation threshold. For example, if the estimated percolation threshold is 0.48, it is possible to define a second round using range=seq(0.47,0.49,0.0001), which provide an accurary of 0.0001 estimating only 201 networks.
Because no modules can be defined when all nodes are isolated, if both, "modules" and "plotALL" options are set to TRUE, a segmentation fault is produced. 
}
\references{
Rozenfeld AF, Arnaud-Haond S, Hernandez-Garcia E, Eguiluz VM, Serrao EA, Duarte CM. (2008). Network analysis identifies weak and strong links in a metapopulation system. Proceedings of the National Academy of Sciences,105, 18824 -18829.

Munoz-Pajares, AJ. (2013) Erysimum mediohispanicum at the evolutionary crossroad: Phylogeography, phenotype and pollinators.
}
\author{A.J. Munoz-Pajares
}
\seealso{single.network
}
\examples{
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

 # Estimating indel distances after reading the alignment from file:
distGap<-MCIC(input="ex2.fas",saveFile=FALSE)
 # Estimating substitution distances after reading the alignment from file:
library(ape)
align<-read.dna(file="ex2.fas",format="fasta")
dist.nt <-dist.dna(align,model="raw",pairwise.deletion=TRUE)
DISTnt<-as.matrix(dist.nt)


 # Obtaining the arithmetic mean of both matrices using the corrected method:
CombinedDistance<-nt.gap.comb(DISTgap=distGap, alpha=0.5, method="Corrected",
saveFile=FALSE, DISTnuc=DISTnt)
 # Estimating the percolation threshold of the combined distance, modifying
 # labels:
perc.thr(dis=CombinedDistance,label=paste(substr(row.names(
CombinedDistance),11,11),substr(row.names(CombinedDistance),21,21),sep="-"))
 # The same network showing different modules as different colours
 # (randomly selected):
perc.thr(dis=as.data.frame(CombinedDistance),label=paste(substr(row.names(
as.data.frame(CombinedDistance)),11,11),substr(row.names(as.data.frame(
CombinedDistance)),21,21),sep="-"), modules=TRUE)

 # The same network showing different modules as different colours
 # (defined by user):
perc.thr(dis=as.data.frame(CombinedDistance),label=paste(substr(row.names(
as.data.frame(CombinedDistance)),11,11),substr(row.names(as.data.frame(
CombinedDistance)),21,21),sep="-"), modules=TRUE,moduleCol=c("pink",
"lightblue","lightgreen"))
 

}
