\name{assign.whole.taxo}
\alias{assign.whole.taxo}
\title{Get taxonomy described in sequence names}
\description{Assign taxonomy to every line in a BLAST output using the information provided in the name of the subject sequences (stitle)}
\usage{assign.whole.taxo(BLAST)}

\arguments{
  \item{BLAST}{data.frame containing the output of a BLAST analysis. The first column must be the name of the sequences matching the queries and must contain information on the taxonomy of the subject sequences. See details.}
}

\details{
The expected input data.frame must contain information about
taxonomy in the first column. Additional information is accepted if separated
by "|", but taxonomy 
must be the last bit of information. Taxonomical information must be provided for kingdom,
phylum, class, order, family, genus, and species, each separated
by ";" and identified by a letter as follows: 

 optionalTEXT | optionalTEXT | k__kingdomName;p__phylumName;c__className; o__orderName;f__familyName;g__genusName;s__speciesName

This is the typical format of sequence names in several databases. Thus a BLAST output using any of these databases will automatically
produce the desired format.
}

\value{a data.frame contining all the information provided in the input data.frame and seven additional columns containing the name of kingdom, phylum, class, order, family, genus, and species for this sequence
}


\author{A. J. Muñoz-Pajares}

\seealso{
\code{\link{filter.whole.taxo}, and \link{get.majority.taxo}}
}
\examples{
# data(ex_BLAST)
# TAXO <- assign.whole.taxo(ex_BLAST)
}


