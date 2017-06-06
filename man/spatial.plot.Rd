\name{spatial.plot}
\alias{spatial.plot}
\title{spatial plot of populations}
\description{This function estimates the phylogeographic relationships among populations, displaying nodes according to geographic coordinates on maps.}
\usage{
spatial.plot(dis=NULL, align=NA, X=NULL, Y=NULL, indel.method="MCIC",
substitution.model="raw", pairwise.deletion=TRUE, alpha="info",
combination.method="Corrected", na.rm.row.col=FALSE, addExtremes=FALSE,
NameIniPopulations=NA, NameEndPopulations=NA, NameIniHaplotypes=NA,
NameEndHaplotypes=NA, HaplosNames=NA, save.distance=FALSE,
save.distance.name="DistanceMatrix_threshold.txt",
network.method="percolation", range=seq(0,1,0.01), modules=FALSE,
moduleCol=NA, modFileName="Modules_summary.txt", bgcol="white",
label.col="black", label=NA, label.sub.str=NA, label.pos= "b",
cex.label=1,cex.vertex=1,vertex.size="equal", plot.edges=TRUE,
lwd.edge=1,to.ggmap=FALSE, plot.ggmap=FALSE, zoom.ggmap=6,
maptype.ggmap="satellite", label.size.ggmap=3)
}

\arguments{
  \item{dis}{a matrix; the distance matrix to be analysed. Alternatively, you can define an alignment using 'align' option.}
  \item{align}{a 'DNAbin' object; the alignment to be analysed. See "read.dna" in the ape package for details about reading alignments. Alternatively, you can define a distance matrix using the 'dis' option.}
  \item{X}{a vector; longitude for each population}
  \item{Y}{a vector; latitude for each population}
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

		-"zero": computes a network connecting all nodes showing distances equal to zero. See ?NINA.thr for details.
}
  \item{range}{a numeric vector between 0 and 1, is the range of thresholds (referred to the maximum distance in the input matrix) to be screened (by default, 101 values from 0 to 1). This option is used for "percolation" and "NINA" network methods and ignored for "zero" method.}
  \item{addExtremes}{a logical; if TRUE, additional nucleotide sites are included in both extremes of the alignment. This will allow estimating distances for alignments showing gaps in terminal positions. This option is used for "SIC", "FIFTH" and "BARRIEL" indel methods and ignored for "MCIC" method.}
  \item{alpha}{a numeric between 0 and 1, is the weight given to the indel genetic distance matrix in the combination. By definition, the weight of the substitution genetic matrix is the complementary value (i.e., 1-alpha). The value "info" will use the proportion of informative substitutions per informative indel event as weight. It is also possible to define multiple weights to estimate different combinations.}
  \item{combination.method}{a string defining whether each distance matrix must be divided by its maximum value before the combination ("Corrected") or not ("Uncorrected"). Consequently, if the "Corrected" method is chosen, both matrices will range between 0 and 1 before being combined.}
  \item{na.rm.row.col}{a logical; if TRUE, distance matrix missing values are removed.}
  \item{modules}{a logical, If TRUE, nodes belonging to different modules are represented as different colours (defined by 'moduleCol').}
  \item{moduleCol}{(if modules=TRUE) a vector of strings defining the colour of nodes belonging to different modules in the network. If 'NA' (or there are less colours than haplotyes), colours are automatically selected}
  \item{modFileName}{(if modules=TRUE) a string, the name of the file to be generated containing a summary of module results (sequence name, module, and colour in network)}
  \item{save.distance}{a logical; if TRUE, the distance matrix used to build the network will be saved as a file.}
  \item{save.distance.name}{a string; if save.distance=TRUE, it defines the name of the file to be saved.}
  \item{bgcol}{a vector of strings; the colour of the background for each node in the network. Can be equal for all nodes (if only one colour is defined), customized (if several colours are defined), or can represent different modules (see "modules" option). If set to 'NA' (default) or if less colours than haplotyes are defined, colours are automatically selected.}
  \item{label.col}{a vector of strings; the colour of labels for each node in the network. Can be equal for all nodes (if only one colour is defined) or customized (if several colours are defined).}
  \item{label}{a vector of strings; labels for each node. By default are the sequence names. (See "substr" function in base package to automatically reduce name lengths)}
  \item{label.sub.str}{a vector of two numerics; if node labels are a sub-string of sequence names, these two numbers represent the initial and final character of the string to be represented. See Example for details.}
  \item{label.pos}{a sting; position for vertex labels regarding vertex position (do not affect the ggmap output). Possible values are: "b" or "below" (default), "a" or "above"; "l" or "left"; "r" or "right" and "c" or "centre"}
  \item{lwd.edge}{a numeric; the width of the edge linking nodes (1.5 by default).}
  \item{cex.label}{a numeric; the size of the node labels.}
  \item{cex.vertex}{a numeric; the size of the nodes.}
  \item{NameIniPopulations}{a numeric; Position of the initial character of population names within sequence names. If not provided, it is set to 1. It is used only if NameEndPopulations is also defined.}
  \item{NameEndPopulations}{a numeric; Position of the last character of population names within sequence names. If not provided, it is set to the first "_" character in the sequences name. It is used only if NameIniPopulations is also defined.}
  \item{NameIniHaplotypes}{a numeric; Position of the initial character of haplotype names within sequence names. If not provided, haplotype names are given and the value is set accordingly. It is used only if NameEndHaplotypes is also defined.}
  \item{NameEndHaplotypes}{a numeric; Position of the last character of haplotype names within sequence names. If not provided, haplotype names are given and the value is set accordingly. It is used only if NameIniHaplotypes is also defined.}
  \item{plot.edges}{a logical; must the edges connecting nodes be potted?}
  \item{HaplosNames}{a sting; the name of the haplotypes (if different from default: H1...Hn)}
  \item{to.ggmap}{a logical; if TRUE, the algorithm generates a list with information required to represent the resulting network using ggmap (see details).}
  \item{plot.ggmap}{a logical; if TRUE, populations (and edges producing a network if 'plot.edges' is set to TRUE) are represented within a map automatically downloaded according the population coordinates.}
  \item{zoom.ggmap}{a numeric; sets the zoom of the map (higher values mean deeper zoom)}
  \item{maptype.ggmap}{a string; types of maps implemented by 'ggplot' are: "terrain", "satellite", "roadmap", "hybrid", "toner", and "watercolor")}
  \item{label.size.ggmap}{a numeric; controls the labe size in the ggplot}
  \item{vertex.size}{a string to define the ratio of vertices representing populations. Possible values are: "equal" (default) to give the same size to all vertices; or "area" to make the vertex area proportional to the population sample size.}

}

\details{

Despite the large list of options, the only mandatory options for this function are the geographic coordinates ('X' and 'Y' options) of the studied populations and either the alignment or the distance matrix ('align' or 'dis', respectively). The remaining options can be classified into five groups:

	1- options defining the computation of both indel and substitution distances (indel.method, substitution.model, pairwise.deletion).

	2- options defining the combination of these two distance matrices (alpha, combination.method, na.rm.row.col, addExtremes, NameIniPopulations, NameEndPopulations, NameIniHaplotypes, NameEndHaplotypes, HaplosNames, save.distance, save.distance.name).

	3- options defining the computation of the network (network.method, range).

	4- options customizing the resulting network (modules, moduleCol, modFileName, bgcol, label.col, label, label.sub.str, cex.label, cex.vertex, vertex.size, plot.edges, lwd.edge).

	5- options dealing with map representation (to.ggmap, plot.ggmap, zoom.ggmap, maptype.ggmap, label.size.ggmap).

Although the 'indel.method' option affects both the distance estimation and the number of mutations represented in the network, the 'substitution.model' and 'pairwise.deletion' options only affect the distance matrix computation. 

This function provides limited options for representing of the resulting population network within a map using the 'ggmap' package. To take advantage of the additional options implemented in ggmap, the 'to.ggmap' option generates a list with the following information:

	1- location: centroid of the population coordinates (required to center the map)

	2- colours: the colour to represent each population (useful, for example to represent modules)

	3- coordinates: the geographic coordinates of the studied populations

	4- network: the resulting population network, represented as a 1/0 matrix

	5- links: a two column matrix representing the edges within the resulting network. Each row provides information on the two elements that are connected by a link.
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

\examples{

library(ggplot2)
data(ex_Coords)
data(ex_alignment1) # this will read a fasta file with the name 'alignExample'

# A simple plot of the population network using geographic coordinates:
spatial.plot (align=alignExample,X=ex_Coords[,2],Y=ex_Coords[,3])

# Changing vertex names and location:
spatial.plot (align=alignExample,X=ex_Coords[,2],Y=ex_Coords[,3],
cex.vertex=2,label=c(1:8),label.pos="c",modules=TRUE)

# Plotting network on a map:
# Uncomment the lines below. It would take more than 5 seconds to run
# spatial.plot (align=alignExample,X=ex_Coords[,2],Y=ex_Coords[,3],
# cex.vertex=2,label=c(1:8),modules=TRUE, plot.ggmap=TRUE)

# Displaying only population coordinates (sampling desing).
# Uncomment the lines below. It would take more than 5 seconds to run
# spatial.plot (align=alignExample,X=ex_Coords[,2],Y=ex_Coords[,3],
# cex.vertex=2,label=c(1:8), plot.ggmap=TRUE,plot.edges=FALSE,
# bgcol=c("red","orange","green4","green1","yellow","brown","blue","purple"))



}
