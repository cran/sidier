\name{double.plot}
\alias{double.plot}
\title{Haplotype and population networks including mutations and haplotype frequencies.}
\description{This function makes a double plot by dividing the active device into two parts. The left part is used to represent the input alignment as a haplotypic network displaying mutations. The right part is used to represent the same input alignment as a population network displaying nodes as pie charts.}
\usage{
double.plot(align = NA, indel.method = "MCIC",
substitution.model = "raw",pairwise.deletion = TRUE,
network.method.mut = "percolation", network.method.pie = "percolation",
range = seq(0, 1, 0.01), addExtremes = FALSE, alpha.mut = "info",
alpha.pie = "info", combination.method.mut = "Corrected",
combination.method.pie = "Corrected", na.rm.row.col.mut = FALSE,
na.rm.row.col.pie = FALSE, save.distance.mut = FALSE,
save.distance.name.mut = "DistanceMatrix_threshold_Mutations.txt",
save.distance.pie = FALSE, save.distance.name.pie =
"DistanceMatrix_threshold_Pies.txt", modules=FALSE,
modules.col=NA, bgcol = NA, label.col.mut = "black",
label.col.pie = "black", label.mut = NA,label.pie = NA,
label.sub.str.mut = NA,label.sub.str.pie = NA, colInd = "red",
colSust = "black", lwd.mut = 1, InScale=1, SuScale=1, lwd.edge = 1.5,
cex.mut = 1, cex.label.mut = 1, cex.label.pie = 1, cex.vertex = 1,
main=c("Haplotypes","Populations"),
NameIniPopulations = NA, NameEndPopulations = NA, NameIniHaplotypes = NA,
NameEndHaplotypes = NA, cex.pie = 1, HaplosNames = NA, offset.label = 1.5)
}
\arguments{
  \item{align}{a 'DNAbin' object; the alignment to be analysed. See "read.dna" in the ape package for details about reading alignments.}
  \item{indel.method}{a sting; the method to define indel events in your alignment. The available methods are:

	  -"MCIC": (Default) Estimates indel events following the rationale of the Modified Complex Indel Coding (Muller, 2006).

	  -"SIC": Estimates indel events following the rationale of Simmons and Ochoterrena (2001).

	  -"FIFTH": Estimates indel events following the rationale of the fifth state: each gap within the alignment is treated
as an independent mutation event.

	  -"BARRIEL": Estimates indel events following the rationale of Barriel (1994): singleton gaps are not taken into account.
}
  \item{substitution.model}{a string; the substitution evolutionary model to estimate the distance matrix. By default is set to "raw" and
estimates the pairwise proportion of variant sites. See the evolutionary models available using ?dist.dna from the ape package.
}
  \item{pairwise.deletion}{a logical; if TRUE (default) substitutions found in regions being a gap in other sequences will account
for the distance matrix. If FALSE, sites being a gap in at least one sequence will be removed before distance estimation.
}
  \item{network.method.mut}{a string; the method to build the haplotypic network. The available methods are:

		-"percolation": computes a network using the percolation network method following Rozenfeld et al. (2008).
 See ?perc.thr for details

		-"NINA": computes a network using the No Isolation Nodes Allowed method. See ?NINA.thr for details.

		-"zero": computes a network connecting all nodes showing distances equal to zero. See ?zero.thr for details.}

  \item{network.method.pie}{a string; the method to build the population network. The available methods are the same than for 'network.method.mut'}
  \item{range}{a numeric vector between 0 and 1, is the range of thresholds (referred to the maximum distance in the input matrix) to be screened (by default, 101 values from 0 to 1). This option is used for "percolation" and "NINA" network methods and ignored for "zero" method.}
  \item{addExtremes}{a logical; if TRUE, additional nucleotide sites are included in both extremes of the alignment. This will allow estimating distances for alignments showing gaps in terminal positions. This option is used for "SIC", "FIFTH" and "BARRIEL" indel methods and ignored for "MCIC" method.}
  \item{alpha.mut}{a numeric between 0 and 1, is the weight given to the indel genetic distance matrix in the combination to represent the haplotypic network. By definition, the weight of the substitution genetic matrix is the complementary value (i.e., 1-alpha). The value "info" (default) will use the proportion of informative substitutions per informative indel event as weight. It is also possible to define multiple weights to estimate different combinations.}
  \item{alpha.pie}{a numeric between 0 and 1, is the weight given to the indel genetic distance matrix in the combination to represent the population network. By definition, the weight of the substitution genetic matrix is the complementary value (i.e., 1-alpha). The value "info" (default) will use the proportion of informative substitutions per informative indel event as weight. It is also possible to define multiple weights to estimate different combinations.}
  \item{combination.method.mut}{a string defining whether each distance matrix must be divided by its maximum value before the combination ("Corrected") or not ("Uncorrected"). Consequently, if the "Corrected" method is chosen, both matrices will range between 0 and 1 before being combined. This option affects the haplotype network depiction.}
  \item{combination.method.pie}{a string defining whether each distance matrix must be divided by its maximum value before the combination ("Corrected") or not ("Uncorrected"). Consequently, if the "Corrected" method is chosen, both matrices will range between 0 and 1 before being combined. This option affects the population network depiction.}
  \item{na.rm.row.col.mut}{a logical; if TRUE, distance matrix missing values are removed.}
  \item{na.rm.row.col.pie}{a logical; if TRUE, distance matrix missing values are removed.}
  \item{save.distance.mut}{a logical; if TRUE, the distance matrix used to build the haplotypic network will be saved as a file.}
  \item{save.distance.name.mut}{a string; if save.distance.mut=TRUE, it defines the name of the file to be saved.}
  \item{save.distance.pie}{a logical; if TRUE, the distance matrix used to build the population network will be saved as a file.}
  \item{save.distance.name.pie}{a string; if save.distance.pie=TRUE, it defines the name of the file to be saved.}
  \item{modules}{a logical. If TRUE, nodes colours are set according to modules in the network of haplotypes.}
  \item{modules.col}{(if modules=TRUE) a vector of strings defining the colour of nodes belonging to different modules in the network. If 'NA' (or there are less colours than haplotypes), colours are automatically selected.}
  \item{bgcol}{a vector of strings; the colour of the background for each node in the haplotypic network. The same colours will be used to represent haplotypes in the population network. If set to 'NA' (default), colours are automatically defined.}
  \item{label.col.mut}{a vector of strings; the colour of labels for each node in the haplotypic network. Can be equal for all nodes (if only one colour is defined) or customized (if several colours are defined).}
  \item{label.col.pie}{a vector of strings; the colour of labels for each node in the population network. Can be equal for all nodes (if only one colour is defined) or customized (if several colours are defined).}
  \item{label.mut}{a vector of strings; labels for each node in the haplotypic network. By default are the sequence names. (See "substr" function in base package to automatically reduce name lengths)}
  \item{label.pie}{a vector of strings; labels for each node in the population network. By default are the sequence names. (See "substr" function in base package to automatically reduce name lengths)}
  \item{label.sub.str.mut}{a vector of two numerics; if node labels are a substring of sequence names, these two numbers represent the initial and final character of the string to be represented in the haplotypic network. See Example for details.}
  \item{label.sub.str.pie}{a vector of two numerics; if node labels are a substring of sequence names, these two numbers represent the initial and final character of the string to be represented in the population network. See Example for details.}
  \item{colInd}{a strings; the color used to represent indels.}
  \item{colSust}{a strings; the color used to represent substitutions.}
  \item{lwd.mut}{a numeric; the width of the line (marks perpendicular to the edge line) representing mutations (1 by default).}
  \item{InScale}{a numeric; the number of indels each mark represents. By default is set to 1 (that is, 1 mark= 1 indel). If set to 10, then 1 mark=10 indels. In that case, if there are 25 indels in a given edge it is represented by three marks (being one of them half width than the other two).}
  \item{SuScale}{a numeric; the number of substitutions each mark represents. By default is set to 1 (that is, 1 mark= 1 substitution). If set to 10, then 1 mark=10 substitutions In that case, if there are 25 substitutions in a given edge it is represented by three marks (being one of them half width than the other two).}
  \item{lwd.edge}{a numeric; the width of the edge linking nodes (1.5 by default).}
  \item{cex.mut}{a numeric; the length of the line (perpendicular to the edge line) representing mutations (1 by default).}
  \item{cex.label.mut}{a numeric; the size of the node labels in the haplotypic network.}
  \item{cex.label.pie}{a numeric; the size of the node labels in the population network.}
  \item{cex.vertex}{a numeric; the size of the nodes in the haplotypic network.}
  \item{main}{a vector with two elements defining the title for the left and right plots, respectively. Alternatively, may be set to "summary" to display the main options selected for representing the networks. Finally, if set to "", the algorithm will show no title for any network.}
  \item{NameIniPopulations}{a numeric; Position of the initial character of population names within sequence names. If 'NA' (default), it is set to 1.}
  \item{NameEndPopulations}{a numeric; Position of the last character of population names within sequence names. If 'NA' (default), it is set to the first "_" character in the sequences name.}
  \item{NameIniHaplotypes}{a numeric; Position of the initial character of haplotype names within sequence names. If 'NA' (default), haplotype names are defined by the algorithm and the value is set accordingly.}
  \item{NameEndHaplotypes}{a numeric; Position of the last character of haplotype names within sequence names. If 'NA' (default), haplotype names are defined by the algorithm and the value is set accordingly.}
  \item{cex.pie}{a numeric; the size of the nodes in the population network.}
  \item{HaplosNames}{a sting; the name of the haplotypes (if different from default: H1...Hn)}
  \item{offset.label}{a numeric, the separation between node and label.}
}
\references{
Barriel, V., 1994. Molecular phylogenies and how to code insertion/
deletion events. Life Sci. 317, 693-701, cited and described by
Simmons, M.P., Müller, K. & Norton, A.P. (2007)
The relative performance of indel-coding methods in simulations.
Molecular Phylogenetics and Evolution, 44, 724–740.

Muller K. (2006). Incorporating information from length-mutational events into phylogenetic analysis. Molecular Phylogenetics and Evolution, 38, 667-676.

Rozenfeld AF, Arnaud-Haond S, Hernandez-Garcia E, Eguiluz VM, Serrao EA, Duarte CM. (2008). Network analysis identifies weak and strong links in a metapopulation system. Proceedings of the National Academy of Sciences,105, 18824-18829.

Simmons, M.P. & Ochoterena, H. (2000). Gaps as Characters in Sequence-Based
Phylogenetic Analyses. Systematic Biology, 49, 369-381.
}
\author{A. J. Muñoz-Pajares}

\seealso{\code{\link{mutation.network}}, \code{\link{pie.network}}}
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
library(ape)

# Double plot (computation time may exceed 5 seconds)
#double.plot(align=read.dna(file="ex2.fas",format="fasta"))

# Other options
#data(ex_alignment1) # this will read a fasta file with the name 'alignExample'
#double.plot(alignExample,label.sub.str.mut=c(7,9))
#double.plot(alignExample,label.sub.str.mut=c(7,9),InScale=3,SuScale=2,lwd.mut=1.5)
}
