\name{GraphResults}
\alias{GraphResults}
\title{
Graphical depiction of results from AnalyzeBlastR
}
\usage{
GraphChroms(chromdat,
            scaff.lim = T,
            dbchrom.lim = T,
            sample = 1,
            max.scaffs = 10,
            max.dbchroms = 10)
}
\arguments{
  \item{chromdat}{output from AnalyzeBlastR}
  \item{scaff.lim}{Boolean option whether or not to limit the number of scaffolds considered; if true, the limit is set by max.scaffs}
  \item{dbchrom.lim}{Same as above except for the number of database chromosomes}
  \item{sample}{Indicates the proportion of data to graph. Can save time and also reduce noisiness.}
  \item{max.scaffs}{Maximum number of scaffolds to infer homology for. Automatically takes largest.}
  \item{max.dbchroms}{Same as above except for database chromosomes.}
}
\details{
Produces two side by side plots. The left plot depicts color coded database chromosomes, while the right plot depicts the distribution of database hits on the query scaffolds. Note that due to the limited number of colors that can feasibly be distinguished at a single time, it is recommended to plot no more than 10 database chromosomes. 
}
\examples{

}
\keyword{BLAST}
