## ----outform, echo=F-----------------------------------------------------
insertPlot <- function(file, caption){
#    outputFormat = knitr::opts_knit$get("rmarkdown.pandoc.to")
#  if(outputFormat == 'latex')
#    paste("![ ", caption, " ](", file, ")",sep="")
}
bigskip <- function(){
#  outputFormat = knitr::opts_knit$get("rmarkdown.pandoc.to")
#  if(outputFormat == 'latex')
#    "\\bigskip"
#  else
    "<br>"
}

## ----getFiles, eval = F, echo=F------------------------------------------
#  Rcpp::sourceCpp('../RcppFunctions/cppFns.cpp')
#  source('../RFunctions/mastifFunctions.r')

## ----simSetup0, eval = F-------------------------------------------------
#  seedNames  <- specNames  <- 'acerRubr'
#  sim <- list(nyr=10, ntree=30, nplot=5, ntrap=40, meanDist = 15,
#                 specNames = specNames, seedNames = seedNames)

## ----sim0, eval = F------------------------------------------------------
#  inputs     <- mastSim(sim)        # simulate dispersal data
#  seedData   <- inputs$seedData     # year, plot, trap, seed counts
#  trueValues <- inputs$trueValues   # true states and parameter values
#  treeData   <- inputs$treeData     # year, plot, tree data
#  xytree     <- inputs$xytree       # tree locations
#  xytrap     <- inputs$xytrap       # trap locations
#  formulaFec <- inputs$formulaFec   # fecundity model
#  formulaRep <- inputs$formulaRep   # maturation model

## ----formulaFec, eval = F------------------------------------------------
#  formulaFec

## ----treeData0, eval = F-------------------------------------------------
#  head(treeData)

## ----xytree, eval = F----------------------------------------------------
#  head(xytree, 5)

## ----treeData1, eval = F-------------------------------------------------
#  head(seedData)

## ----map1a, eval = F-----------------------------------------------------
#  dataTab <- table(treeData$plot,treeData$year)
#  
#  w <- which(dataTab > 0,arr.ind=T) # a plot-year with observations
#  w <- w[sample(nrow(w),1),]
#  
#  mapYears <- as.numeric( colnames(dataTab)[w[2]] )
#  mapPlot  <- rownames(dataTab)[w[1]]
#  
#  mastMap(inputs, treeSymbol=treeData$diam, mapPlot = mapPlot,
#          mapYears = mapYears, SCALEBAR = T)

## ----map3, eval=F--------------------------------------------------------
#  trueValues <- inputs$trueValues
#  mastMap(inputs, treeSymbol=trueValues$fec, mapPlot = mapPlot,
#          mapYears = mapYears, scaleTree=1, scaleTrap = 1, SCALEBAR = T)

## ----hist0, eval = F-----------------------------------------------------
#  par( mfrow=c(2,1),bty='n', mar=c(4,4,1,1) )
#  seedData <- inputs$seedData
#  seedNames <- inputs$seedNames
#  
#  hist( as.matrix(seedData[,seedNames]) ,nclass=100,
#        xlab = 'seed count', ylab='per trap', main='' )
#  hist( trueValues$fec,nclass=100, xlab = 'seeds produced', ylab = 'per tree',
#        main = '')

## ----mast2, eval = F-----------------------------------------------------
#  output   <- mastif( inputs = inputs, ng = 1500, burnin = 500 )

## ----tabPars0, eval = F--------------------------------------------------
#  summary( output )

## ----pars, eval = F------------------------------------------------------
#  plotPars <- list(trueValues = trueValues)
#  mastPlot(output, plotPars)

## ----restart, eval=F-----------------------------------------------------
#  predList <- list( mapMeters = 3, plots = mapPlot, years = mapYears )
#  output   <- mastif( inputs = output, ng = 2000, burnin = 1000,
#                    predList = predList)

## ----mapout, eval = F----------------------------------------------------
#  mastMap(output, mapPlot = predList$plots, mapYears = predList$years,
#          scaleTree = 2, scaleTrap = .5, PREDICT = T, scaleValue = 20)

## ----simSetup, eval = F--------------------------------------------------
#  specNames <- c('pinuTaeda','pinuEchi','pinuVirg')
#  seedNames <- c('pinuTaeda','pinuEchi','pinuVirg','pinuUNKN')
#  sim    <- list(nyr=5, ntree=25, nplot=8, ntrap=50, specNames = specNames,
#                    seedNames = seedNames)

## ----sim, eval = F-------------------------------------------------------
#  inputs <- mastSim(sim)        # simulate dispersal data
#  R      <- inputs$trueValues$R # species to seedNames probability matrix
#  R

## ----mast3, eval = F-----------------------------------------------------
#  output <- mastif( inputs = inputs, ng = 2000, burnin = 1000)

## ----tabPars, eval = F---------------------------------------------------
#  summary( output )

## ----pars0, eval = F-----------------------------------------------------
#  plotPars <- list(trueValues = inputs$trueValues)
#  mastPlot(output, plotPars)

## ----again, eval=F-------------------------------------------------------
#  tab   <- with( inputs$seedData, table(plot, year) )
#  years <- as.numeric( colnames(tab)[tab[1,] > 0] ) # years for 1st plot
#  predList <- list( plots = 'p1', years = years )
#  output <- mastif( inputs = output, ng = 5000, burnin = 2000,
#                    predList = predList)
#  mastPlot(output, plotPars)

## ----map21, eval=F-------------------------------------------------------
#  library(repmis)
#  d <- "https://github.com/jimclarkatduke/mast/blob/master/mast_liri.Rdata?raw=True"
#  source_data(d)
#  mapList <- list( treeData = treeData, seedData = seedData,
#                   specNames = specNames, seedNames = seedNames,
#                   xytree = xytree, xytrap = xytrap)
#  mastMap(mapList, mapPlot = 'DF_HW', mapYear = 2011:2014,
#          treeSymbol = treeData$diam, scaleTree = 1, scaleTrap=.7,
#          SCALEBAR=T, scaleValue=50)

## ----litu1, eval=F-------------------------------------------------------
#  head(treeData, 3)

## ----litu2, eval=F-------------------------------------------------------
#  head(seedData, 3)

## ----fit, eval=F---------------------------------------------------------
#  formulaFec <- as.formula( ~ canopy + I(log(diam)) )   # fecundity model
#  formulaRep <- as.formula( ~ I(log(diam)) )            # maturation model
#  
#  inputs   <- list(specNames = specNames, seedNames = seedNames,
#                   treeData = treeData, seedData = seedData, xytree = xytree,
#                   xytrap = xytrap, priorDist = 25, priorVDist = 5, minDist = 8,
#                   maxDist = 40, sigmaMu = 1, minDiam = 10, maxDiam = 60)
#  output <- mastif( formulaFec, formulaRep, inputs = inputs,  ng = 1500,
#                    burnin = 500 )

## ----more, eval=F--------------------------------------------------------
#  predList <- list( mapMeters = 10, plots = 'DF_HW', years = 2010:2015 )
#  output <- mastif( inputs = output, ng = 4000, burnin = 1500, predList = predList )

## ----plotmydata1, eval=F-------------------------------------------------
#  mastPlot(output)

## ----outpars, eval=F-----------------------------------------------------
#  summary( output )

## ----fitSum, eval=F------------------------------------------------------
#  output$fit

## ----ranEff0, eval=F-----------------------------------------------------
#  randomEffect <- list(randGroups = 'tree',
#                       formulaRan = as.formula( ~ I(log(diam)) ) )

## ----reinit0, eval=F-----------------------------------------------------
#  treeData$lastFec <- output$inputs$treeData$lastFec
#  treeData$lastRep <- output$inputs$treeData$lastRep
#  inputs$treeData  <- treeData

## ----ranEff, eval=F------------------------------------------------------
#  randomEffect <- list(randGroups = 'tree',
#                       formulaRan = as.formula( ~ I(log(diam)) ) )
#  output <- mastif( formulaFec, formulaRep, inputs = inputs,
#                    ng = 2000, burnin = 1000,
#                    randomEffect = randomEffect )

## ----ranEff2, eval=F-----------------------------------------------------
#  output <- mastif( inputs = output, ng = 4000, burnin = 1000,
#                    randomEffect = randomEffect)
#  mastPlot(output)

## ----fitSum2, eval=F-----------------------------------------------------
#  output$fit

## ----regtab, eval=F------------------------------------------------------
#  with(treeData, colSums( table(treeID, region)) )

## ----newgr, eval=F-------------------------------------------------------
#  province <- rep('mtn',nrow(treeData))
#  province[treeData$region == 'DF'] <- 'piedmont'
#  treeData$province <- as.factor(province)

## ----yrpl, eval=F--------------------------------------------------------
#  yearEffect <- list(plotGroups = 'province')

## ----reinit, eval=F------------------------------------------------------
#  treeData$lastFec <- output$inputs$treeData$lastFec
#  treeData$lastRep <- output$inputs$treeData$lastRep

## ----fit3, eval=F--------------------------------------------------------
#  inputs$treeData  <- treeData
#  output <- mastif(formulaFec, formulaRep,  inputs = inputs, ng = 1500, burnin = 500,
#                 randomEffect = randomEffect, yearEffect = yearEffect)

## ----moreYR, eval=F------------------------------------------------------
#  predList <- list( mapMeters = 10, plots = 'DF_HW', years = 1998:2014 )
#  output <- mastif(inputs = output, predList = predList,
#                 randomEffect = randomEffect, yearEffect = yearEffect,
#                 ng = 3000, burnin = 1000)
#  mastPlot(output)

## ----map2, eval=F--------------------------------------------------------
#  d <- "https://github.com/jimclarkatduke/mast/blob/master/mast_pinu.Rdata?raw=True"
#  repmis::source_data(d)
#  
#  mapList <- list( treeData = treeData, seedData = seedData,
#                   specNames = specNames, seedNames = seedNames,
#                   xytree = xytree, xytrap = xytrap)
#  mastMap(mapList, mapPlot = 'DF_HW', mapYears = c(2004:2007),
#          treeSymbol = treeData$diam, scaleTree = .5, scaleTrap=.5,
#          scalePlot = 1.2, LEGEND=T)

## ----fit0, eval=F--------------------------------------------------------
#  formulaFec <- as.formula( ~ canopy + I(log(diam)) )   # fecundity model
#  formulaRep <- as.formula( ~ I(log(diam)) )            # maturation model
#  
#  yearEffect   <- list(specGroups = 'species', plotGroups = 'region', p = 5)
#  randomEffect <- list(randGroups = 'tree',
#                       formulaRan = as.formula( ~ I(log(diam)) ) )
#  
#  seedMass <- matrix( c(0.0170,0.0270,0.0167,0.0070,0.0080), ncol=1)
#  rownames(seedMass) <- c('pinuEchi','pinuRigi','pinuStro','pinuTaed','pinuVirg')
#  colnames(seedMass) <- 'gmPerSeed'
#  
#  
#  inputs   <- list( specNames = specNames, seedNames = seedNames,
#                    treeData = treeData, seedData = seedData,
#                    xytree = xytree, xytrap = xytrap, priorDist = 20,
#                    priorVDist = 5, minDist = 15, maxDist = 30,
#                    minDiam = 15, maxDiam = 40,
#                    maxF = 1e+8, seedMass = seedMass)
#  output <- mastif(formulaFec, formulaRep, inputs = inputs, ng = 2000,
#                 burnin = 500, yearEffect = yearEffect,
#                 randomEffect = randomEffect)

## ----moreAR, eval=F------------------------------------------------------
#   output <- mastif(inputs = output, ng = 2000,
#                 burnin = 1000, yearEffect = yearEffect,
#                 randomEffect = randomEffect)
#  plotPars <- list(MAPS = F, SPACETIME=T, SAVEPLOTS=F)
#  mastPlot(output, plotPars = plotPars)

## ----moreYR2, eval=F-----------------------------------------------------
#  plots <- c('DF_EW','DF_BW','DF_HW','HF_ST')
#  years <- 1980:2025
#  predList <- list( mapMeters = 8, plots = plots, years = years )
#  output <- mastif(inputs = output, predList = predList, yearEffect = yearEffect,
#                 randomEffect = randomEffect, ng = 2000, burnin = 1000)

## ----yrPlot, eval=F------------------------------------------------------
#  mastPlot( output, plotPars = list(MAPS=F) )

## ----onemap, eval=F------------------------------------------------------
#  mastMap(output, mapPlot = 'DF_EW', mapYears = c(2011:2012),
#          PREDICT=T, scaleTree = 1, scaleTrap=.3, LEGEND=T, scaleValue=50,
#          scalePlot = 1.5, COLORSCALE = T, mfrow=c(2,1))

## ----onemap1, eval=F-----------------------------------------------------
#  mastMap(output, mapPlot = 'DF_EW', mapYears = 2014, PREDICT=T,
#          scaleTree = 1, scaleTrap=.3, LEGEND=T, scalePlot = 10,
#          SCALEBAR = T, COLORSCALE = T)

## ----outpars0, eval=F----------------------------------------------------
#  summary( output )

## ----experiments, echo=F, eval=F-----------------------------------------
#  rCons <- rMu[,1:2]
#  
#  rCons[,1] <- .9
#  rCons[,2] <- .1
#  rCons[4,1] <- .1
#  rCons[4,2] <- .1
#  
#  tmp <- getCvol(output, plot='DF_HW', rCons, npoints = 5, nyears = 5,
#                      sampleYears = 2003:2015, nrep = 50)
#  rCons
#  tmp$entropy
#  C <- tmp$C
#  index <- tmp$index
#  
#  library(corrplot)
#  
#  #C[upper.tri(C)] <- NA
#  ccor <- cov2cor(C)
#  
#  diag(ccor) <- min(ccor)
#  
#  clim <- range(ccor[lower.tri(ccor)])
#  clim[1] <- clim[1] - .0001
#  clim[2] <- clim[2] + .0001
#  
#  col1 <- colorRampPalette(c("white","white","white","white","white",
#                             "white","white",
#                             "white", "yellow", "red"))
#  
#  corrplot(ccor, method='square',  cl.lim = clim, is.cor=F, col = col1(20),
#           cl.length=21,
#           diag = F, type = 'lower', addgrid.col=NA, tl.pos = 'n')

## ----nopred, eval=F------------------------------------------------------
#  d <- "https://github.com/jimclarkatduke/mast/blob/master/mast_liri.Rdata?raw=True"
#  repmis::source_data(d)
#  
#  formulaFec   <- as.formula(~ 1)
#  formulaRep   <- as.formula( ~ I(log(diam)) )
#  yearEffect   <- list(specGroups = 'species', plotGroups = 'region')
#  randomEffect <- list(randGroups = 'treeID',
#                       formulaRan = as.formula( ~ 1 ) )
#  inputs <- list(specNames = specNames, seedNames = seedNames,
#                 treeData = treeData, seedData = seedData,
#                 xytree = xytree, xytrap = xytrap,
#                 priorDist = 10, priorVDist = 5, maxDist = 50, minDist = 5,
#                 minDiam = 25, maxF = 1e+6)
#  output <- mastif(inputs, formulaFec, formulaRep, ng = 2000, burnin = 1000,
#                 randomEffect = randomEffect, yearEffect = yearEffect )
#  mastPlot(output)

## ----land, eval=F, echo=F------------------------------------------------
#  
#  # rank distances from random locations in a plot:
#  
#  wj <- which( as.character(xytree$plot) == 'DF_EW' )
#  xy <- xytree[wj, c('x','y')]
#  lims <- apply(xy,2,range)
#  
#  nsite <- 50
#  x <- runif(nsite, lims[1,1],lims[2,1])
#  y <- runif(nsite, lims[1,2],lims[2,2])
#  xysite <- cbind(x,y)
#  colnames(xysite) <- c('x','y')
#  
#  tmp <- rankNeighbors(xy, xysite, dr = 1, nrank=10)
#  rmat <- tmp$rmat
#  mh   <- tmp$mh
#  breaks <- tmp$breaks
#  
#  par(mfrow=c(1,1), bty='n')
#  
#  plot(breaks, rmat[1,], type='s', lwd=2)
#  
#  for(k in 1:nrow(rmat)){
#    lines(breaks, rmat[k,], type='s', col=.getColor('black',1/k))
#  }
#  plot(breaks, mh, type='s')
#  
#  upar <- 15
#  f    <- 10000
#  kern <- t(upar/pi/(upar + breaks^2)^2)
#  elam <- f*sum(kern*mh)
#  vlam <- f^2*sum(mh^2*kern) - elam^2

