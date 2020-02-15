\name{BlastN}
\alias{BlastN}
\title{
Takes a database genome and a query genome and blasts the latter against the former.
}
\description{
Takes a query genome and a BLAST database  and tries to find hits for the query genome by slicing it up by a user defined quantity. 
}
\usage{
BlastN(query, db, slicesize=100, leaps=900)
}
\arguments{
  \item{query}{Path to the query genome file
}
  \item{db}{Path to the database genome file (must already have been turned into a BLAST database via MakeDB)}
  
  \item{slicesize}{Length of each nucleotide sequence to BLAST against the database}
  
  \item{leaps}{Number of nucleotides to skip between each BLAST query}
}

\details{
Scans through the query genome, reading in a nucleotide sequence of length slicesize, then stores that sequence in an output file. Skips leaps number of nucleotides and then repeats until the end of the query genome is reached. The output file containing all the nucleotide sequences is then queried all at the same time against the database.
}
\value{Returns a number of data files (text files). These should not be altered by the user.
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Johnathan Lo and Heath Blackmon
}

\keyword{ BLAST }

