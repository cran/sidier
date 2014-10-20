\name{alignExample}
\alias{alignExample}
\title{example alignment #1 ('DNAbin' class)}
\description{object of class 'DNAbin' to test some function within this package}
\usage{data(ex_alignment1)}
\details{
Because fasta file examples are not automatically loaded into R environment, 'ex_alignment1' function generates a fasta file (named Example_Spatial.plot_Alignment) that is stored as a 'DNAbin' object named alignExample
}

\author{A. J. Muñoz-Pajares}

\seealso{
\code{\link{ex_alignment1}}, \code{\link{Example_Spatial.plot_Alignment}}
}
\examples{

data(ex_alignment1) # this will read a fasta file with the name 'alignExample'
alignExample
}
