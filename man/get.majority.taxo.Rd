\name{get.majority.taxo}
\alias{get.majority.taxo}
\title{Get majority taxonomy for a sequence}
\description{Given filtered taxonomy of multiple BLAST hits, provides the consensus taxonomy for every query sequence}

\usage{get.majority.taxo(filtered.taxo,verbose=TRUE)}

\arguments{
  \item{filtered.taxo}{a data.frame containing (at least) the output of a BLAST analysis and the filtered taxonomy for every subject sequence. This format is produced by \code{\link{filter.whole.taxo}}}
  \item{verbose}{a logical, if TRUE details on the calculation are shown.}
}

\details{
The expected input data.frame must contain information about filtered taxonomy for every subject sequence and a unique code for every query sequence. 

Nonmeaningful names (including "unidentified", "sp", "low_pident") are coerced to "Uninformative".

The output dataframe includes a column with a summary of the alternative classifications found for taxonomy at every level. For example: "10|10|10|10|10|10|9+1" means that all 10 matches have identical taxonomy up to genus, but two species has been identified, being the taxonomy of 9 of the subject sequences identical and different from the remaining subject.
}

\value{a data.frame with the followin columns: "qseqid", the unique identifier of every original query sequence; seven columns containing the filtered taxonomy ("kingdom.final", "phylum.final", "class.final", "order.final", "family.final", "genus.final", and "species.final"); "values": The frequency of the different taxonomical names provides for every level, separated by "|" (see details).  
}


\author{A. J. Muñoz-Pajares}

\seealso{
\code{\link{filter.whole.taxo},\link{assign.whole.taxo}}
}
\examples{
# data(ex_BLAST)
# TAXO <- assign.whole.taxo(ex_BLAST)
# FILT_TAXO <- filter.whole.taxo(TAXO)
# MAJ_TAXO <- get.majority.taxo(TAXO)
}


