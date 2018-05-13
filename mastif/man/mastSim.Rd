\name{mastSim}
\alias{mastSim}
\title{Data simulation for mast model}
\description{
Simulates data for analysis by \code{\link{mastif}} in package \code{mastif}.
}
\usage{
  mastSim(sim)
}
\arguments{
  \item{sim}{\code{list} of inputs for simulation, see Details.}
}

\details{

The \code{list} \code{sim} contains the following:

\code{specNames}: \code{character vector} of species names.

\code{seedNames}: \code{character vector} of seed names.

\code{nyr = 5}:  average number of years for a plot

\code{ntree = 10}:  average number of trees in \code{specNames} on a plot

\code{plotWide = 100}:  diameter of plot

\code{nplot = 3}:  number of plots

\code{ntrap = 20}:  average number of seed traps on a plot

\code{meanDist = 15}:  mean dispersal (meters)

}

\value{

Returns an object of \code{\link{class} "mast"}, a list containing the following components:

\item{distall}{\code{matrix} of distances between all seed traps (rows) and trees (columns).  Block diagonal structure by plot, padded with the value \code{10000} for tree/trap combinations on different plots.}

\item{formulaFec}{\code{formula} for the fecundity model.}

\item{formulaRep}{\code{formula} for the maturation model.}

\item{treeData}{\code{data.frame} of tree data needed for \code{mast}.}

\item{seedData}{\code{data.frame} of seed data needed for \code{mast}.}

\item{xytree}{\code{data.frame} of tree locations needed for \code{mast}.}

\item{xytrap}{\code{data.frame} of seed locations needed for \code{mast}.}


\item{trueValues}{\code{list} holding true parameter values includes:

\code{fec}: \code{vector} of log fecundity states. 

\code{repr}: vector of maturation states.

\code{betaFec}: fecundity parameter values.

\code{betaRep}: maturation parameter values.

\code{upar}: dispersal parameter, m^2.

\code{R}: species to seed type matrix.
}

\item{xfec}{design \code{matrix} for fecundity.}

\item{xrep}{design \code{matrix} for maturation.}

\item{years}{\code{vector} of years.}

\item{sample}{\code{list} with data summary for trees and seeds.}

}


\examples{
\dontrun{
# simulate data 
seedNames  <- specNames  <- 'acerRubr'
sim <- list(nyr=10, ntree=30, nplot=5, ntrap=40, meanDist = 15,
            specNames = specNames, seedNames = seedNames)
inputs <- mastSim(sim)
output <- mastif( inputs = inputs, ng = 2000, burnin = 1000 )

plotPars <- list(trueValues = inputs$trueValues)
mastPlot(output, plotPars)
}
}
\author{
James S Clark, \email{jimclark@duke.edu}
}
\references{
Clark, J.S., S. LaDeau, and I. Ibanez. 2004. Fecundity of trees and the colonization-competition hypothesis, Ecological Monographs, 74, 415-442.

Clark, J.S., D. M Bell, M. Kwit, A. Powell, And K. Zhu. 2013. Dynamic inverse prediction and sensitivity analysis with high-dimensional responses: application to climate-change vulnerability of biodiversity.  Journal of Biological, Environmental, and Agricultural Statistics, 18, 376-404.
}
\seealso{
\code{\link{mastSim}} simulates data

A more detailed vignette is can be obtained with:

\code{browseVignettes('mastif')}

website 'http://sites.nicholas.duke.edu/clarklab/code/'.
}

