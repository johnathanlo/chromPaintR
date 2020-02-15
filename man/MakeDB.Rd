\name{MakeDB}
\alias{MakeDB}
\title{
Takes a fasta file and turns it into a BLAST database.
}
\description{
Takes a fasta file and turns it into a BLAST database.
}
\usage{
FindXmers(dbtype = "nucl", fsa_file)
}
\arguments{
  \item{dbtype}{The type of database to create. Currently only nucleotide databases are supported.
}
  \item{fsa_file}{path to fasta file to construct BLAST database from}
}

\details{
Creates a BLAST database.
}
\value{Creates a BLAST database. No return value.
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Johnathan Lo and Heath Blackmon
}

\keyword{ BLAST }

