\name{mastMap}
\alias{mastMap}
\title{Map dispersal data with predictions for mastif model}
\description{
Maps dispersal data (trees and seed traps) with predictions.
}
\usage{
                        
  mastMap(mapList, mapPlot, mapYears, treeSymbol = NULL, PREDICT = F, 
          treeScale = 1, trapScale = 1, xlim = NULL, ylim = NULL, 
          MAPTRAPS = T, seedMax = NULL, fecMax = NULL, mfrow = NULL, 
          LEGEND = F, plotScale = 1, SCALEBAR = F, scaleValue = NULL, 
          COLORSCALE = T)
}
\arguments{
   \item{mapList}{\code{list} that includes \code{data.frame}s described in the \code{mastif} help page: \code{treeData}, \code{seedData}, \code{xytree}, \code{xytrap}. Alternatively, \code{mapList} can be an object of \code{class mast} simulated by \code{mastSim} or fitted by \code{mast}.}
   
      \item{mapPlot}{\code{character} for the plot name, which appears in \code{treeData$plot} and \code{seedData$plot}.}

   \item{mapYears}{\code{numeric vector} containing years to plot, which appear in \code{treeData$year} and \code{seedData$year}.}
  
   \item{treeSymbol}{\code{character} indicating a column name in \code{treeData} holding numeric values that scale the sizes of trees in the map, e.g., \code{diam} for diameter.}
   
   \item{PREDICT}{if \code{TRUE}, then predicted seed will be mapped.  Must be predicted in \code{mastif}.}
   
  \item{treeScale}{\code{numeric} value to scale sizes of tree symbols.}
  
  \item{trapScale}{\code{numeric} value to scale sizes of trap symbols.}
  
  \item{xlim}{\code{vector} with upper and lower plot limits.}
  
  \item{ylim}{\code{vector} with upper and lower plot limits.}
  
  \item{MAPTRAPS}{\code{logical} show trap symbols.}
    
  \item{seedMax}{\code{numeric} value can be used to set maximum seed symbol/prediction contour.}
  
    \item{fecMax}{\code{numeric} value can be used to set maximum fecundity symbol size.}
  
    \item{mfrow}{\code{numeric vector} of length-2 to set arrangement of maps on plot, as in \code{pars()}.}
    
  \item{LEGEND}{add \code{legend} with species names to map.}
  
  \item{plotScale}{some  level of control on size of plot is possible with values > 1.}
  
  \item{SCALEBAR}{add a scale bar with distance to bottom of map.}
  
  \item{scaleValue}{\code{numeric} value in meters used if \code{SCALEBAR = T}.}
  
  \item{COLORSCALE}{if \code{TRUE}, then a color scale is added for predicted seed.}
}

\details{

Generates of map of seed traps and trees, with symbols scaled to the sizes relative to seed counts in \code{sdata$seedNames} and \code{treeSymbol}.  Sizes are adjusted with \code{scaleTree} and \code{scaleTrap}.

If \code{PREDICT = T}, then predictions come in the object fitted in \code{mast} with \code{predictList} used to specify prediction plots and years.  See the help page for \code{mastif}. 

More detailed vignettes can be obtained with:

\code{browseVignettes('mastif')}
}

\value{

Only graphical outputs.
}
\examples{
\dontrun{
# simulate data (see \link{\code{mastSim}})
seedNames  <- specNames  <- 'acerRubr'
sim    <- list(nyr=10, ntree=30, nplot=5, ntrap=40, meanDist = 15,
               specNames = specNames, seedNames = seedNames)
inputs <- mastSim(sim)

mastMap(inputs, treeSymbol=treeData$diam, mapPlot = 'p1', 
        mapYears = inputs$years[1])

# for Pinus
d <- "https://github.com/jimclarkatduke/mast/blob/master/mast_pinu.Rdata?raw=True"
repmis::source_data(d)

mapList <- list( treeData = treeData, seedData = seedData, 
                 specNames = specNames, seedNames = seedNames, 
                 xytree = xytree, xytrap = xytrap)
mastMap(mapList, mapPlot = 'DF_HW', mapYears = c(2004:2007), 
        treeSymbol = treeData$diam, scaleTree = .5, scaleTrap=.5,
        scalePlot = 1.2, LEGEND=T)
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

