\name{perc.thr}
\alias{perc.thr}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Percolation threshold network
%%  function to do ... 
}
\description{This function computes the percolation network following Rozenfeld et al. (2008), as described in Munoz-Pajares. 
%%   A concise (1-5 lines) description of what the function does. 
}
\usage{
perc.thr(dis, threshold = seq(0, 1, 0.01), ptPDF = TRUE, ptPDFname = "PercolatedNetwork.pdf", estimPDF = TRUE, estimPDFname = "PercThr Estimation.pdf", estimOutfile = TRUE, estimOutName = "PercThresholdEstimation.txt", appendOutfile = TRUE, plotALL = FALSE, bgcol = "white", label.col = "black", label = colnames(dis), modules = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dis}{the distance matrix to be represented
%%     Describe \code{dis} here
}
  \item{threshold}{a numeric vector between 0 and 1, is the range of thresholds (referred to the maximum distance in a matrix) to be screened (by default, 101 values from 0 to 1).
%%     Describe \code{threshold} here
}
  \item{ptPDF}{a logical, must the percolated network be saved as a pdf file?
%%     Describe \code{ptPDF} here
}
  \item{ptPDFname}{if ptPDF=TRUE, the name of the pdf file containing the percolation network to be saved ("PercolatedNetwork.pdf", by default)
%%     Describe \code{ptPDFname} here
}
  \item{estimPDF}{a logical, must the percolation threshold estimation be saved as a pdf file? If estimPDF=TRUE (default) the value of <s> for each threshold is also saved

%%     Describe \code{estimPDF} here
}
  \item{estimPDFname}{if estimPDF=TRUE (default), defines the name of the pdf file containing the percolation threshold estimation ("PercThr Estimation.pdf" by default).
%%     Describe \code{estimPDFname} here
}
  \item{estimOutfile}{a logical, must the matrix containing percolation threshold estimation variables be saved as a pdf file?
%%     Describe \code{estimOutfile} here
}
  \item{estimOutName}{if estimOutfile=TRUE (default), contains the name of the text file containing the percolation threshold estimation ("PercThr Estimation.txt" by default).

%%     Describe \code{estimOutName} here
}
  \item{appendOutfile}{a logical, if estimOutfile=TRUE, it defines whether results must be appended to an existing file with the same name (TRUE) or the existing file must be replaced (FALSE).
%%     Describe \code{appendOutfile} here
}
  \item{plotALL}{a logical, must all the networks calculated during the  percolation threshold estimation be saved as different pdf files? (FALSE, by default). If TRUE, for each value defined in threshold, one file is generated.
%%     Describe \code{plotALL} here
}
  \item{bgcol}{string, defining the colour of the background for each node in the network. Can be equal for all nodes (if only one colour is defined), customized (if several colours are defined), or can represent different modules (see modules option).
%%     Describe \code{bgcol} here
}
  \item{label.col}{string, defining the colour of labels for each node in the network. Can be equal for all nodes (if only one colour is defined) or customized (if several colours are defined),
%%     Describe \code{label.col} here
}
  \item{label}{string, labels for each node. By default are the column names of the distance matrix (dis). (See substr function in base package to automatically reduce name lengths).
%%     Describe \code{label} here
}
  \item{modules}{a logical, must nodes belonging to different modules be represented as different colours?
%%     Describe \code{modules} here
}
}
\details{By default, percolation threshold is estimated with an accuracy of 0.01, but it may be increased by setting the decimal places in threshold function (e.g., seq(0,1,0.0001)). However, it may strongly increase computation times (in this example, it is required to estimate 100 001 instead of 101 networks). It is also possible to increase accuracy with a low increase in computation time by repeating the process and increasing decimal places only in a range close to a previously estimated percolation threshold.
%%   If necessary, more details than the description above 
}

\references{Rozenfeld AF, Arnaud-Haond S, Hernandez-Garcia E, Eguiluz VM, Serrao EA, Duarte CM. (2008). Network analysis identifies weak and strong links in a metapopulation system. Proceedings of the National Academy of Sciences,105, 18824 -18829.

Munoz-Pajares, AJ. SIDIER: Substitution and Indel Distances to Infer Evolutionary Relationships.
%% put references to the literature/web site here 
}
\author{A.J. Munoz-Pajares
%%  who you are
}

%% Make other sections like Warning with \section{Warning }{....} 

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
dist.nt <-dist.dna(align,model="raw",pairwise.deletion=TRUE)
DISTnt<-as.matrix(dist.nt)


 # Obtaining the arithmetic mean of both matrices using the corrected method:
CombinedDistance<-nt.gap.comb(DISTgap=distGap, range=0.5, method="Corrected", saveFile=FALSE, DISTnuc=DISTnt)
 # Estimating the percolation threshold of the combined distance, modifying labels:
#perc.thr(dis=as.data.frame(CombinedDistance$Corrected),label=paste(substr(row.names(as.data.frame(CombinedDistance$Corrected)),11,11),substr(row.names(as.data.frame(CombinedDistance$Corrected)),21,21),sep="-"),ptPDF=FALSE,estimPDF=FALSE)
 # Estimating the percolation threshold of the combined distance, modifying labels:
#perc.thr(dis=as.data.frame(CombinedDistance$Corrected),label=paste(substr(row.names(as.data.frame(CombinedDistance$Corrected)),11,11),substr(row.names(as.data.frame(CombinedDistance$Corrected)),21,21),sep="-"))
 # The same network showing different modules as different colours (randomly selected):
#perc.thr(dis=as.data.frame(CombinedDistance$Corrected),label=paste(substr(row.names(as.data.frame(CombinedDistance$Corrected)),11,11),substr(row.names(as.data.frame(CombinedDistance$Corrected)),21,21),sep="-"), module=TRUE)



}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ kwd1 }
%\keyword{ kwd2 }% __ONLY ONE__ keyword per line
