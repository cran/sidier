\name{pie.network}
\alias{pie.network}
\title{Population network depiction including haplotype frequencies}
\description{This function represents an alignment as a population network and displays nodes as pie charts where haplotype frequencies are proportional to the area depicted in different colours.}
\usage{
pie.network(align = NA, indel.method = "MCIC", substitution.model = "raw",
pairwise.deletion = TRUE, network.method = "percolation",
range = seq(0, 1, 0.01), addExtremes = FALSE, alpha = "info",
combination.method = "Corrected", na.rm.row.col = FALSE,
NameIniPopulations = NA, NameEndPopulations = NA, NameIniHaplotypes = NA,
NameEndHaplotypes = NA, save.distance = FALSE,
save.distance.name = "DistanceMatrix_threshold.txt",
pop.distance.matrix = NULL, Haplos = NULL, HaplosPerPop = NULL,
col.pie = NA, label.col = "black", label = NA, label.sub.str = NA,
cex.label = 1, cex.pie = 1, main = "", HaplosNames = NA,
offset.label = 1.5, pie.size = "equal", coord = NULL, get.coord = TRUE)
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{align}{a 'DNAbin' object; the alignment to be analysed. See "read.dna" in the ape package for details about reading alignments. Other inputs are available: Use a distance matrix instead an alignment using the 'align' option or provide a list of haplotypes and frequencies per population using 'Haplos' and 'HapPerPop' options
}
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
  \item{addExtremes}{a logical; if TRUE, additional nucleotide sites are included in both extremes of the alignment. This will allow estimating distances for alignments showing gaps in terminal positions. This option is used for "SIC", "FIFTH" and "BARRIEL" indel methods and ignored for "MCIC" method.}
  \item{alpha}{a numeric between 0 and 1, is the weight given to the indel genetic distance matrix in the combination. By definition, the weight of the substitution genetic matrix is the complementary value (i.e., 1-alpha). The value "info" will use the proportion of informative substitutions per informative indel event as weight. It is also possible to define multiple weights to estimate different combinations.}
  \item{combination.method}{a string defining whether each distance matrix must be divided by its maximum value before the combination ("Corrected") or not ("Uncorrected"). Consequently, if the "Corrected" method is chosen (default option), both matrices are corrected to range between 0 and 1 before being combined.}
  \item{na.rm.row.col}{a logical; if TRUE, removes rows and columns showing missing values within the distance matrix.}
  \item{NameIniPopulations}{a numeric; Position of the initial character of population names within sequence names. If not provided, it is set to 1. It is used only if NameEndPopulations is also defined.}
  \item{NameEndPopulations}{a numeric; Position of the last character of population names within sequence names. If not provided, it is set to the first "_" character in the sequences name. It is used only if NameIniPopulations is also defined.}
  \item{NameIniHaplotypes}{a numeric; Position of the initial character of haplotype names within sequence names. If not provided, haplotype names are given and the value is set accordingly. It is used only if NameEndHaplotypes is also defined.}
  \item{NameEndHaplotypes}{a numeric; Position of the last character of haplotype names within sequence names. If not provided, haplotype names are given and the value is set accordingly. It is used only if NameIniHaplotypes is also defined.}
  \item{save.distance}{a logical; if TRUE, the distance matrix used to build the network will be saved as a file.}
  \item{save.distance.name}{a string; if save.distance=TRUE, the name of the file to be saved.}
  \item{pop.distance.matrix}{a matrix containing the population distances. Alternatively, it can be estimated from a given sequence alignment using 'align'. Alternatively, you can provide a list of haplotypes and frequencies using 'Haplos' and 'HapPerPop'}
  \item{Haplos}{a two columns matrix containing sequence names and haplotype names as reported by  \code{\link{FindHaplo}}. Alternatively, you can define an input alignment using 'align' or a distance matrix using 'pop.distance.matrix'.}
  \item{HaplosPerPop}{a matrix containing the number of haplotypes found per population, as reported by \code{\link{HapPerPop}} (Weighted matrix).Alternatively, you can define an input alignment using 'align' or a distance matrix using 'pop.distance.matrix'.}
  \item{col.pie}{a vector of strings; the colour to represent each haplotype. If 'NA' (or there are less colours than haplotyes), colours are automatically selected.}
  \item{label.col}{a vector of strings; the colour of labels for each node in the network. Can be equal for all nodes (if only one colour is defined) or customized (if several colours are defined).}
  \item{label}{a vector of strings; labels for each node. By default are the sequence names. (See "substr" function in base package to automatically reduce name lengths)}
  \item{label.sub.str}{a vector of two numerics; if node labels are a substring of sequence names, these two numbers represent the initial and final character of the string to be represented. See Example for details.}
  \item{cex.label}{a numeric; the size of the node labels.}
  \item{cex.pie}{a numeric; the size of the nodes (pie charts).}
  \item{main}{a sting; if set to "summary" the main options selected for representing the network are displayed in title. The default value ("") shows no title for the network.}
  \item{HaplosNames}{a sting; the name of the haplotypes (if different from default: H1...Hn)}
  \item{offset.label}{a numeric; the separation between node and label.}
  \item{pie.size}{a string to define the ratio of pies representing populations. Possible values are: "equal" (default) to give the same size to all pies; "radius" to make the pie radius proportional to the population sample size; "area" to make the pie area proportional to the population sample size; or "points" to display simple vertices instead of pies representing haplotypes per population.}

  \item{coord}{a two columns matrix containing coordinates where each haplotypes must be represented.}

  \item{get.coord}{a logical, TRUE to obtain coordinates of nodes within the network}
}
\details{
It is recommended to use equal length names with population and individual names separated by '_' (e.g., Pop01_id001...Pop23_id107) and set population (both, NameIniPopulations and NameEndPopulations,) and haplotype (both, NameIniHaplotypes and NameIniHaplotypes) identifiers accordingly. If any of these identifiers is not provided, the algorithm will behave as follows:

   -If only haplotype name identifiers are defined, population names are assumed between character 1 and the first symbol '_' in sequences name.

   -If only population name identifiers are defined, haplotype are automatically found and named using the 'HapPerPop' function.

   -If both are not defined, population names are assumed between character 1 and the first symbol '_' in sequences name and haplotypes are automatically found and named using the 'HapPerPop' function.

}

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
\seealso{mutation.network, double.plot}
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
# # The input format is recognized, and names identifiers can be omitted:
# 	 pie.network(align=example)
# 
# # Is identical to:
# 	 pie.network(align=example, NameIniPopulations=1,NameEndPopulations=11)
# 
# # Using different colours:
# 	 pie.network(align=example, NameIniPopulations=1,NameEndPopulations=11,
# 	 col.pie=c("red","blue","pink","orange","black","grey"))
# 
# # col.pie is omitted if less colours than haplotypes are defined:
# 	 pie.network(align=example, NameIniPopulations=1,NameEndPopulations=11,
# 	 col.pie=c("red","blue","pink"))
# 
# # and also if more colours than haplotypes are defined:
# 	  pie.network(align=example, NameIniPopulations=1,NameEndPopulations=11,
# 	 col.pie=c("red","blue","green","purple","pink","orange","gray"))
# 
}
