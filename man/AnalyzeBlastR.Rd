\name{AnalyzeBlastR}
\alias{AnalyzeBlast}
\title{
Analyzes output files from BlastN.
}
\description{
This function takes output files from BlastN and organizes the hits by their database scaffold, query scaffold, location on each scaffold, and filters them based on a user input threshold of consecutive number of hits and whether or not all hits should be shown in the final output.
}
\usage{
AnalyzeBlastR(data = "results_v03.out", chromlengthfile = "chromlength.out", min.chain = 5, allhits = F)
}
\arguments{
  \item{data}{path to an output file produced by BlastN, default is "results_v03.out" in working directory.
}
  \item{chromlengthfile}{path to an output file produced by BlastN, default is "chromlength.out" in working directory.}
  \item{min.chain}{The minimum number of consecutive hits to register as a homologous sequence.}
}

\details{
Sorts through the dataframe produced by ReadFasta() to find all instances of a specific monomer-length microsatellite, e.g. 2mer, 3mer, and provides some accessory statistics.
}
\value{Returns a list containing 2 dataframes. The first dataframe is called "summary" and is organized by scaffold on the query genome. It contains information for each scaffold, in parallel order, on hits from the database genome ("Database hits"), locations of those hits on the query scaffold ("Locations on query scaffold"), locations of the hits on the database ("Locations on database"), and length of each scaffold ("Scaffold Lengths"). 
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Johnathan Lo and Heath Blackmon
}

\keyword{ microsatellites }
\keyword{ fasta file analysis }
