\name{simuEvolution}
\alias{simuEvolution}
\title{
Simulate sequences evolution
}
\description{
This function simulates the evolution of a set of sequences. It is necessary to define evolution topology, substitution rate, indel rate and insertion/deletion rate in a matrix (see details).
}
\usage{
simuEvolution(input, seqL, iLength, nReplicates)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{input}{
Matrix defining evolution topology and mutation rates.
}
  \item{seqL}{
Length of the simulated sequences.
}
  \item{iLength}{
Length of indel mutations.
}
  \item{nReplicates}{
Number of independent sequence sets to be simulated.
}
}
\details{
Evolution details must be provided in a file consisting in five columns separated by spaces. The first two columns define topology by indicating the ancestor and the derived sequence, respectively. The remaining columns provide rates for substitutions and indels as well as the ratio between insertions and deletions.
The simulation is performed over the complete alignment. To test the effect of alignment method over the simulated sequences it will be necessary to degap the yielded sequences.
}
\value{
For each replicate, two files are generated: one containing the alignment with all the generated sequences and the other containing only tips sequences (i.e., sequences that are not the ancestor of any other sequence).
}
\author{
A. J. Muñoz-Pajares
}

\examples{

#Generating matrix defining evolution:
Input<-matrix(c(1,rep(2:8,2),2:16,rep(0.03,15),rep(0.008,15),rep(0.5,15)),ncol=5)
#Simulating 2 replicates of the evolutionary process:
simuEvolution(input=Input, seqL=1000, iLength=20, nReplicates=2)
}
