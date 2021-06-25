\name{filter.whole.taxo}
\alias{filter.whole.taxo}
\title{Get consensus taxonomy}
\description{Given the taxonomy of multiple BLAST hits for a query sequence, provides the most likely tanoxomy for the query sequence taking BLAST percentage of identity values into account}
\usage{filter.whole.taxo(whole.taxo)}

\arguments{
  \item{whole.taxo}{data.frame containing BLAST results and taxonomy information. Can be produced by \code{\link{assign.whole.taxo}}}
}

\details{
The expected input data.frame must contain information about BLAST hits (particularly, a "pident" column with the percentage of identity) and seven additional columns containing the name of kingdom, phylum, class, order, family, genus, and species for every subject sequence.

Depending on the "pident" value, taxonomy for the subject sequence will be retained until species (if pident>=97), genus (97>pident>=90), family  (90>pident>=85), order (85>pident>=80), family (80>pident>=75), or class (75>pident>=0). For taxonomic levels showing pident lower that these thresholds, "low_pident" is returned.
}

\value{a data.frame contining all the information provided in the input data.frame and seven additional columns containing the name of kingdom, phylum, class, order, family, genus, and species for this sequence after filtering by BLAST percentage of identity.
}


\author{A. J. Muñoz-Pajares}

\seealso{
\code{\link{get.majority.taxo},\link{assign.whole.taxo}}
}

\examples{
#data(ex_BLAST)
#TAXO <- assign.whole.taxo(ex_BLAST)
#FILT_TAXO <- filter.whole.taxo(TAXO)
}

