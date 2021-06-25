\name{mutation.network}
\alias{mutation.network}
\title{Haplotype network depiction including mutations}
\description{This function represents an alignment as a network and displays mutations (substitutions and indels) as marks in edges.}
\usage{
mutation.network(align = NA, indel.method = "MCIC",
substitution.model = "raw", pairwise.deletion = TRUE,
network.method = "percolation", range = seq(0, 1, 0.01),
merge.dist.zero=TRUE, addExtremes = FALSE, alpha = "info",
combination.method = "Corrected", na.rm.row.col = FALSE,
modules = FALSE, moduleCol = NA,
modFileName = "Modules_summary.txt", save.distance = FALSE,
save.distance.name = "DistanceMatrix_threshold.txt",
silent = FALSE, bgcol = "white", label.col = "black",
label = NA, label.sub.str = NA, colInd = "red",
colSust = "black", lwd.mut = 1, lwd.edge = 1.5,
cex.mut = 1, cex.label = 1, cex.vertex = 1, main = "",
InScale = 1, SuScale = 1, legend = NA, legend.bty = "o",
legend.pos="bottomright", large.range = FALSE, pies = FALSE,
NameIniPopulations = NA, NameEndPopulations = NA,
NameIniHaplotypes = NA, NameEndHaplotypes = NA, 
HaplosNames = NA,verbose = TRUE)
}
\arguments{
  \item{align}{a 'DNAbin' object; the alignment to be analysed. See "read.dna" in the ape package for details about reading alignments.}
  \item{indel.method}{a sting; the method to define indel events in your alignments. The available methods are:

	  -"MCIC": (Default) Estimates indel events following the rationale of the Modified Complex Indel Coding (Muller, 2006).

	  -"SIC": Estimates indel events following the rationale of Simmons and Ochoterrena (2000).

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
  \item{network.method}{a string; the method to build the network. The available methods are:

		-"percolation": computes a network using the percolation network method following Rozenfeld et al. (2008).
 See ?perc.thr for details

		-"NINA": computes a network using the No Isolation Nodes Allowed method. See ?NINA.thr for details.

		-"zero": computes a network connecting all nodes showing distances equal to zero. See ?zero.thr for details.
}
  \item{range}{a numeric vector between 0 and 1, is the range of thresholds (referred to the maximum distance in the input matrix) to be screened (by default, 101 values from 0 to 1). This option is used for "percolation" and "NINA" network methods and ignored for "zero" method.}
  \item{merge.dist.zero}{a logical; if TRUE, nodes showing a distance equal to zero are merged using the mergeNodes() function}
  \item{addExtremes}{a logical; if TRUE, additional nucleotide sites are included in both extremes of the alignment. This will allow estimating distances for alignments showing gaps in terminal positions. This option is used for "SIC", "FIFTH" and "BARRIEL" indel methods and ignored for "MCIC" method.}
  \item{alpha}{a numeric between 0 and 1, is the weight given to the indel genetic distance matrix for the combination. By definition, the weight of the substitution genetic matrix is the complementary value (i.e., 1-alpha). The value "info" (default) will use the proportion of informative substitutions per informative indel event as weight. It is also possible to define multiple weights to estimate different combinations.}
  \item{combination.method}{a string defining whether each distance matrix must be divided by its maximum value before the combination ("Corrected") or not ("Uncorrected"). Consequently, if the "Corrected" method is chosen, both matrices will range between 0 and 1 before being combined.}
  \item{na.rm.row.col}{a logical; if TRUE, distance matrix missing values are removed.}
  \item{modules}{a logical, If TRUE, nodes belonging to different modules are represented as different colours (defined by 'moduleCol'). Modules (defined as subsets of nodes that conform densely connected subgraphs) are estimated by means of random walks (see 'igraph' package for details).}
  \item{moduleCol}{(if modules=TRUE) a vector of strings defining the colour of nodes belonging to different modules in the network. If 'NA' (or there are less colours than haplotypes), colours are automatically selected}
  \item{modFileName}{(if modules=TRUE) a string, the name of the file to be generated containing a summary of module results (sequence name, module, and colour in network)}
  \item{save.distance}{a logical; if FALSE (default), the distance matrix used for computation is saved in a file}
  \item{save.distance.name}{if save.distace=TRUE, a string defining the file name}
  \item{silent}{a logical; if FALSE (default), displays a list containing the number of indels and substitutions represented in the network.}
  \item{bgcol}{a vector of strings; the colour of the background for each node in the network. Can be equal for all nodes (if only one colour is defined), customized (if several colours are defined), or can represent different modules (see "modules" option). If set to 'NA' (default) or if less colours than haplotypes are defined, colours are automatically selected.}
  \item{label.col}{a vector of strings; the colour of labels for each node in the network. Can be equal for all nodes (if only one colour is defined) or customized (if several colours are defined).}
  \item{label}{a vector of strings; labels for each node. By default are the sequence names. (See 'label.sub.str' to automatically reduce name lengths)}
  \item{label.sub.str}{a vector of two numerics; if node labels are a substring of sequence names, these two numbers represent the initial and final character of the string to be represented. See Example for details.}
  \item{lwd.edge}{a numeric; the width of the edge linking nodes (1.5 by default).}
  \item{colInd}{a strings; the colour used to represent indels (red by default).}
  \item{colSust}{a strings; the colour used to represent substitutions (black by default).}
  \item{lwd.mut}{a numeric; the width of the line (perpendicular to the edge line) representing mutations (1 by default).}
  \item{cex.mut}{a numeric; the length of the line (perpendicular to the edge line) representing mutations (1 by default).}
  \item{cex.vertex}{a numeric; the size of the nodes.}
  \item{cex.label}{a numeric; the size of the node labels.}
  \item{main}{if set to "summary" the main options selected for representing the network are displayed in title. The default value ("") shows no title for the network.}
  \item{InScale}{a numeric; the number of indels each mark represents. By default is set to 1 (that is, 1 mark= 1 indel). If set to 10, then 1 mark=10 indels. In that case, if there are 25 indels in a given edge it is represented by three marks (being one of them half width than the other two).}
  \item{SuScale}{a numeric; the number of substitutions each mark represents. By default is set to 1 (that is, 1 mark= 1 substitution). If set to 10, then 1 mark=10 substitutions In that case, if there are 25 substitutions in a given edge it is represented by three marks (being one of them half width than the other two).}
  \item{legend}{a logic; if TRUE, plots a legend representing the scale of marks (that is, the number of mutations represented by a mark).}
  \item{legend.bty}{a letter; the type of box to be drawn around the legend. The allowed values are 'o' (default, producing a four-sides box) and 'n' (producing no box).}
  \item{legend.pos}{a string, defines legend position ("bottomright" by default). Other possible values are: ‘"bottomright"’, ‘"bottom"’, ‘"bottomleft"’, ‘"left"’, ‘"topleft"’, ‘"top"’, ‘"topright"’, ‘"right"’ and ‘"center"’.}
  \item{large.range}{a logic, TRUE for representing node size according to three categories: haplotypes appearing less than 10 times, between 10 and 100 times and more than 100 times}
  \item{pies}{a logic, TRUE for representing nodes as pies and FALSE for representing nodes as points}
  \item{NameIniPopulations}{a numeric; Position of the initial character of population names within sequence names. If not provided, it is set to 1. It is used only if NameEndPopulations is also defined.}
  \item{NameEndPopulations}{a numeric; Position of the last character of population names within sequence names. If not provided, it is set to the first "_" character in the sequences name. It is used only if NameIniPopulations is also defined.}
  \item{NameIniHaplotypes}{a numeric; Position of the initial character of haplotype names within sequence names. If not provided, haplotype names are given and the value is set accordingly. It is used only if NameEndHaplotypes is also defined.}
  \item{NameEndHaplotypes}{a numeric; Position of the last character of haplotype names within sequence names. If not provided, haplotype names are given and the value is set accordingly. It is used only if NameIniHaplotypes is also defined.}
  \item{HaplosNames}{a sting; the name of the haplotypes (if different from default: H1...Hn)}
  \item{verbose}{a logical, if TRUE details on the calculation are shown.}
}

\details{Despite the large list of options, the only mandatory option for this function is the alignment to be represented. The remaining options can be classified into four groups:

	1- options defining the computation of both indel and substitution distances (indel.method, substitution.model, pairwise.deletion).

	2- options defining the combination of these two distance matrices (alpha, combination.method, na.rm.row.col, addExtremes, save.distance, save.distance.name).

	3- options defining the computation of the network (network.method, range).

	4- options customizing the resulting network (modules, moduleCol, modFileName, bgcol, label.col, label, label.sub.str, colInd, colSust, lwd.mut, lwd.edge, cex.mut, cex.label, cex.vertex, main).

Although the 'indel.method' option affect both the distance estimation and the number of mutations represented in the network, the 'substitution.model' and 'pairwise.deletion' options only affect the distance matrix computation. The number of substitutuions represented in the network are always estimated using the "N" model and the pairwise deletion of indels.}

\references{
Barriel, V., 1994. Molecular phylogenies and how to code insertion/
deletion events. Life Sci. 317, 693-701, cited and described by
Simmons, M.P., Müller, K. & Norton, A.P. (2007)
The relative performance of indel-coding methods in simulations.
Molecular Phylogenetics and Evolution, 44, 724–740.

Muller K. (2006). Incorporating information from length-mutational events
into phylogenetic analysis. Molecular Phylogenetics and Evolution, 38, 667-676.

Paradis, E., Claude, J. & Strimmer, K. (2004). APE: analyses of phylogenetics
and evolution in R language. Bioinformatics, 20, 289-290.

Rozenfeld AF, Arnaud-Haond S, Hernandez-Garcia E, Eguiluz VM, Serrao EA, Duarte CM. (2008). Network analysis identifies weak and strong links in a metapopulation system. Proceedings of the National Academy of Sciences,105, 18824-18829.

Simmons, M.P. & Ochoterena, H. (2000). Gaps as Characters in Sequence-Based
Phylogenetic Analyses. Systematic Biology, 49, 369-381.

}
\author{A. J. Muñoz-Pajares}

\seealso{\code{\link{plot.network}}, \code{\link{double.plot}}}
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
# "TTAT----TCGAGGGGTAGC",
# ">Population2_sequence2",
# "TAAT----TCTA----TAAC",
# ">Population2_sequence3",
# "TTATAAAA--------TAGC",
# ">Population2_sequence4",
# "TTAT----TCGAGGGGTAGC",
# ">Population3_sequence1",
# "TTAT----TCGA----TAGC",
# ">Population3_sequence2",
# "TTAT----TCGA----TAGC",
# ">Population3_sequence3",
# "TTAT----TCGA----TAGC",
# ">Population3_sequence4",
# "TTAT----TCGA----TAGC",
#      file = "ex2.fas", sep = "\n")
# 
# library(ape)
# 
# #Network with default options
# mutation.network (align=read.dna(file="ex2.fas",format="fasta"))
# 
# #Using more options:
# mutation.network (align=read.dna(file="ex2.fas",format="fasta"),modules=TRUE)
# 
# #A more complex alignment
# data(ex_alignment1) # this will read a fasta file with the name 'alignExample'
# mutation.network (align=alignExample,modules=TRUE,
# InScale=2, SuScale=2,legend=TRUE,lwd.mut=1.8)
# 
}
