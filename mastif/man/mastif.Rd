\name{mastif}
\alias{mastif}
\alias{print.mastif}
\alias{summary.mastif}
\title{Gibbs sampler for mast data}
\description{
Estimates productivity and dispersion of seeds observed at seed traps, using information on locations, and covariates that could explain source strength.  Data can be simulated with \code{\link{mastSim}}.
}
\usage{



  mastif(inputs, formulaFec = NULL, formulaRep = NULL, predList = NULL, 
         yearEffect = NULL, randomEffect = NULL,
         modelYears = NULL, ng = NULL, burnin = NULL)
  
  \method{print}{mastif}(x, ...)
  
  \method{summary}{mastif}(object, ...)
}
\arguments{
  \item{inputs}{\code{list} that must include two \code{character vector}s, (\code{specNames}, \code{seedNames}) and four \code{data.frame}s (\code{treeData}, \code{seedData}, \code{xytree}, \code{xytrap}).  Alternatively, \code{inputs} can be output from \code{mastif}, in which case latent states are initialized from previous fit.  See details.}
  
  \item{formulaFec}{R \code{formula} for fecundity model, e.g., \code{~ x1 + x2}.}
  
  \item{formulaRep}{R \code{formula} for maturation model, e.g., \code{~ x1 + x2}.}
  
  \item{predList}{\code{list} holding \code{plots} and \code{years} for seed prediction.}
    
  \item{yearEffect}{\code{list} holding \code{specGroups} and \code{plotGroups} for year effects or AR(p) model.}
    
    \item{randomEffect}{\code{list} holding \code{randGroups} and \code{formulaRan} for random effects in fecundity.}
    
    \item{modelYears}{\code{numeric vector} can be a subset of years in data set.}

  \item{ng}{\code{numeric} number of Gibbs steps.}
  
  \item{burnin}{\code{numeric} number of burnin steps, before predictions are saved.}
  
  \item{object}{currently, also an object of \code{class mastif}.}
  
  \item{x}{object of \code{class mastif}.}
  
  \item{...}{further arguments not used here.}
}

\details{

\code{inputs} is either a fitted object from \code{mastif} or it includes the following: 

  \code{specNames} is a \code{character vector} containing names of species, \code{specNames}, that appear in the \code{treeData$species} column.

  \code{seedNames} is a \code{character vector} of seed types that appear as column names in \code{seedData}.

  \code{treeData} is a \code{data.frame} holding tree information, including predictors and tree-year identification.  Required columns are \code{plot}, \code{tree}, \code{species}, \code{year}, \code{diam}, and any other predictors for fecundity or maturation. Predictors can be covariates or factors.  A column \code{repr} can be included, with values 0 - observed immature, 1 - observed mature, NA - uncertain or not observed.
  
  \code{seedData} is a \code{data.frame} holding seed counts with seed trap and year identification.  Required columns are \code{plot}, \code{trap}, \code{year}, and \code{seedNames}, the latter holding seed counts.
  
  \code{xytree} is a \code{data.frame} holding tree locations.  Required columns are \code{plot}, \code{tree}, \code{x}, and \code{y}.
  
  \code{xytrap} is a \code{data.frame} holding seed trap locations.  Required columns are \code{plot}, \code{trap}, \code{x}, and \code{y}.

\code{formulaFec} and \code{formulaRep} specify the models for fecundity and maturation.  Variables listed in formulas appear as column headings in \code{treeData}.  Note that \code{formulaFec} and \code{formulaRep} begin with \code{~}, not \code{y ~}.  The response \code{matrix} is constructed from seed types in \code{seedData}. 

The \code{treeData$tree} column has values that are unique for a tree within a \code{plot}.  These reference the same unique identifiers in \code{xytree$tree}.  In addition to these identifiers, the \code{data.frame xytree} holds columns \code{x} and \code{y} for map locations.
 
The \code{character vector seedNames} holds the names of columns in \code{seedData} for seed counts.  The elements of \code{seedNames} are seed types produced by one or more of the species in \code{specNames}.  \code{seedData} must also include columns for \code{trap}, \code{plot}, and \code{year}, which link with  columns in \code{xytrap}, which additionally includes columns \code{x} and \code{y}.
  
\code{predList} includes the names of \code{plots} and \code{years} to be predicted.  It can include a \code{numeric} value \code{mapMeters} for the distance between lattice points in the prediction grid.  See examples.

\code{yearEffect} is a \code{list} indicating the column names in \code{treeData} for random groups in year effects or AR(p) models.  See examples.

\code{randomEffect} is a \code{list} indicating the column names in \code{treeData} for random groups in fecundity estimates, the \code{character randGroups}, and the \code{formula} for random effects.  The \code{formula} must be a subset of predictors from \code{formulaFec}.  See examples.

\code{modelYears} is a \code{numeric vector} of years to include in the analysis.

\code{ng} is the number of Gibbs steps.  \code{burnin} is the number of initial steps, must be less than \code{ng}.  

Additional arguments to \code{inputs} can include prior parameter values and limits, with default values:
    
\code{priorDist = 20} is a prior mean dispersal distance in meters.

\code{priorVDist = 10} is the prior variance on mean dispersal distance in meters. 

\code{minDist = 2} and \code{maxDist = 50} are the minimum and maximum values for the mean dispersal kernel in meters.

\code{minDiam = 2} is the minimum diameter below which a tree of unknown maturation status is assumed to be immature, in cm.  This assignment is overridden if there is a column \code{treeData$repr} and the value is not \code{NA}.

\code{maxDiam = 80} is the maximum diameter above which a tree of unknown maturation status is assumed to be mature, in cm.  This assignment is overridden if there is a column \code{treeData$repr} and the value is not \code{NA}.

\code{sigmaMu = 1}  and \code{sigmaWt = nrow(inputs$treeData)} are the prior mean and the prior weight on log fecundity variance.

\code{maxF = 1e+8}, maximum fecundity, helps stabilize analysis of especially noisy data.

More detailed vignettes can be obtained with:

\code{browseVignettes('mastif')}
}

\value{

Returns an object of \code{\link{class} "mastif"}, which is a list containing the following components:

\item{\code{inputs}}{\code{list} includes inputs to the model, including transformations not resulting from posterior simulation.}

\item{chains}{\code{list} of MCMC matrices, each with \code{ng} rows, includes:

\code{bfec}: fecundity coefficients 

\code{brep}: maturation coefficients

\code{ugibbs}: dispersal parameters

\code{sgibbs}: residual variance and, after \code{burnin}, rmspe (root mean square prediction error)

\code{rgibbs}: if multiple seed types, then rows are the \code{specNames} to \code{seedNames matrix R}.  

\code{bygibbsF, bygibbsR}: if \code{yearEffects} are included in the model, these are fixed year/lag and random year/lag effects by random group.

\code{agibbs}: if \code{randomEffects} on individuals, rows are covariance matrix.

}

\item{parameters}{\code{list} of parameter estimates summarized from \code{chains}.

\code{acfMat}: autocorrelation on fecundity by random group.

\code{betaFec}: fecundity regression coefficients (log scale).

\code{betaRep}: maturation regression coefficients (probit scale).

\code{dpars}: dispersal kernel coefficient, by random group, on meter scale.

\code{omegaList}: fecundity covariance between trees for the same plot-year.

\code{pacfMat}: partial autocorrelation on fecundity.

\code{pacfSe}: standard errors for \code{pacfMat}.

\code{pacsMat}: partial autocorrelation on seed counts.

\code{sigma}: estimate of residual log fecundity variance.

\code{sigmaList}: fecundity covariance over years.

\code{upars}: dispersal kernel coefficient, by random group, on meter^2 scale.

\code{rMu, rSe}: if more than one \code{seedName}, posterior mean and standard error on \code{R} matrix

If \code{yearEffects}, then the following will be included in \code{$parameters}:

\code{betaYrMu, betaYrSe}: posterior mean and standard errors on fixed year effects.

\code{betaYrRand, betaYrRandSE}: posterior mean and standard errors on random year effects.

\code{alphaMu, alphaSe}: posterior mean and standard errors on random tree effects.

\code{aMu, aSe}: posterior mean and standard errors on random effects covariance.

}

\item{prediction}{\code{list} of latent variable estimates and prediction:

\code{fecPred} includes \code{matrEst, fecEstMu, fecEstSe} maturation and fecundity estimates and \code{matrPred, fecPred, fecEst} maturation and fecundity predictions.

\code{seedPred} includes seed predictions by species, from estimates of latent fecundity and maturation (\code{_estMean, _estSe}) and from the fully generative model (\code{_predMean, _predSe}).

If \code{predictList} is passed to \code{mast}, then predictions are returned as seeds per m^2 (not per trap) for the regular prediction grid spaced \code{mapMeters} apart in the \code{data.frame seedPredGrid}.  Tree predictions are returned in\code{treePredGrid} are included. If the AR(p) model is used (\code{yearEffect$p} is supplied), then both \code{data.frames} include \code{p}-yr hind casts and \code{p}-yr forecasts.

}

}
\examples{
\dontrun{
# simulate data (see \link{\code{mastSim}})
seedNames  <- specNames  <- 'acerRubr'
sim <- list(nyr=10, ntree=30, nplot=5, ntrap=40, meanDist = 15,
               specNames = specNames, seedNames = seedNames)
               
inputs   <- mastSim(sim)        # simulate data
predList <- list( mapMeters = 3, plots = inputs$plots[1], 
                  years = inputs$years ) 
output   <- mastif( inputs = inputs, ng = 2000, burnin = 1000, 
                    predList = predList)
mastPlot(output)

# for Pinus data

d <- "https://github.com/jimclarkatduke/mast/blob/master/mast_pinu.Rdata?raw=True"
repmis::source_data(d)

formulaFec <- as.formula( ~ I(log(diam)) )   # fecundity model
formulaRep <- as.formula( ~ I(log(diam)) )            # maturation model

yearEffect   <- list(specGroups = 'species', plotGroups = 'region')   
randomEffect <- list(randGroups = 'treeID', 
                     formulaRan = as.formula( ~ I(log(diam)) ) )

inputs   <- list( specNames = specNames, seedNames = seedNames, 
                  treeData = treeData, seedData = seedData, 
                  xytree = xytree, xytrap = xytrap, priorDist = 8,
                  priorVDist = 1, minDist = 4, maxDist = 25, minDiam = 15,
                  maxDiam = 70, maxF = 1e+8)
output <- mastif(inputs = inputs, formulaFec, formulaRep, ng = 2500, 
                 burnin = 1000, yearEffect = yearEffect, 
                 randomEffect = randomEffect)
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

