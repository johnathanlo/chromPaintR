\name{InfHom}
\alias{InfHom}
\title{
Inference of chromosome homology
}
\description{
Takes the output from AnalyzeBlastR and conducts a Monte Carlo simulation to determine homology.
}
\usage{
InfHom(chromdat, scaff.lim = F, dbchrom.lim = F, max.scaffs = 10, max.dbchroms = 10, resolution =1000)
}
\arguments{
  \item{chromdat}{Output produced by AnalyzeBlastR
}
  \item{scaff.lim}{Boolean option whether or not to limit the number of scaffolds considered; if true, the limit is set by max.scaffs}
  \item{dbchrom.lim}{Same as above except for the number of database chromosomes}
  \item{max.scaffs}{Maximum number of scaffolds to infer homology for. Automatically takes largest.}
  \item{max.dbchroms}{Same as above except for database chromosomes.}
  \item{resolution}{Number of iterations to conduct random sampling. Higher will lead to more statistical power and less variance, but will take more time.}
}

\details{
Takes the data provided by AnalyzeBlastR about the distribution of hits from the BLAST database to determine whether a query scaffold has more hits than one would expect from a purely random sample of the database chromosomes. 
}
\value{Returns a dataframe where columns are the database chromosomes, in decreasing size, and rows are the query scaffolds, also in decreasing size. Each cell contains the number of hits that the query scaffold contains from the corresponding database chromosome, with significance notation (*,**,***) based on the percentage of times that the Monte Carlo simulation produced a scaffold with more hits than what is observed. 
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Johnathan Lo and Heath Blackmon
}

\keyword{ BLAST}
