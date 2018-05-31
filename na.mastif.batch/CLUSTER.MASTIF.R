# NASA MASTIF
#Cluster code

library(stringr)
library(mastif, lib.loc ="r_libs/")
args   <- commandArgs(trailingOnly = FALSE)
myargs <- args[length(args)]
myargs <- sub("-","",myargs)
myargs <- eval(parse(text=myargs))

t=myargs[1]

library(Rcpp, lib.loc ="r_libs/")
load("all.treeData.rdata")
load("otherTreeData.rdata")

priorVals  <- read.table('priorParameters.txt',header=T)
traitTable <- read.table('traitTable.txt', header=T)

#run and save all species output

genFull <- c('abies','acer','ailanthus','amelanchier','betula','carpinus',
             'carya','celtis','cercis','chionanthus','cornus','fagus','fraxinus',
             'ilex','juniperus','liquidambar','liriodendron','magnolia',
             'morus','nyssa','ostrya','oxydendron','picea','pinus','prunus',
             'quercus','robinia','sassafras','sorbus','tilia','tsuga','ulmus')

#####
covTreeData <- all.treeData[[t]]
covTreeData$province <- sub('\\_.*', '', covTreeData$plot)
covTreeData$province <- str_replace(covTreeData$province, "DF", "piedmont")
covTreeData$province <- str_replace(covTreeData$province, "MH", "sApps")
covTreeData$province <- str_replace(covTreeData$province, "CW", "sApps")
covTreeData$province <- str_replace(covTreeData$province,"GSNP", "sApps")
covTreeData$province <- str_replace(covTreeData$province, 'HF', 'NE')
covTreeData$province <- str_replace(covTreeData$province, 'B', 'NE')

names(covTreeData) <- gsub("_", ".", names(covTreeData))

covTreeData$f.PETO <- (covTreeData$flowering.covs.pr.data - covTreeData$flowering.covs.pet.data)
covTreeData$s.PETO <- (covTreeData$seed.set.covs.pr.data - covTreeData$seed.set.covs.pet.data)
otherTreeData.now <- otherTreeData[[t]]

#dealing with UNKN fruits
  wf <- grep('UNKN_fruit',otherTreeData.now$seedNames)
  if(length(wf) > 0){
    otherTreeData.now$seedNames[wf] <- 'unknfruit'
    wf <- grep('UNKN_fruit',colnames(otherTreeData.now$seedData))
    colnames(otherTreeData.now$seedData)[wf] <- 'unknfruit'
  }
  
  wf <- grep('UNKN', otherTreeData.now$seedNames)
  if(length(wf) > 1)stop('multiple unknowns')
  
  wspec <- match(otherTreeData.now$specNames,traitTable$code4)
  seedMass <- as.matrix( traitTable[wspec,'gmPerSeed',drop=F] )*1000
  rownames(seedMass) <- specNames
  
  wna <- which(is.na(seedMass))
  if(length(wna) > 0){
    seedMass[wna] <- mean(seedMass,na.rm=T)
  }
  
  print(seedMass)

#fixing col names to be compatible with MASTIF
 # names(covTreeData) <- str_replace_all(names(covTreeData), "[[.]]", "")
  
setwd('output')
  
  formulaFec <- as.formula( ~ I(log(diam)) + yearlyPETO.tmin1 + yearlyPETO.tmin2 + flowering.covs.pr.data + flowering.covs.tmin.data + s.PETO)

  formulaRep <- as.formula( ~ I(log(diam)) )
  
  randomEffect <- list(randGroups = 'treeID',
                       formulaRan = as.formula( ~ I(log(diam)) ) )
  yearEffect <- list(specGroups = 'species', plotGroups = 'province')
  
  priors <- priorVals[priorVals$genus == genFull[t],]
  plots <- sort(unique(as.character(otherTreeData.now$plot)))
  years <- sort(unique(c(otherTreeData.now$year,otherTreeData.now$seedData$year)))
  covTreeData$province <- as.factor(covTreeData$province)
  
  inputs   <- list( specNames = otherTreeData.now$specNames, 
                    seedNames = otherTreeData.now$seedNames, 
                    treeData = covTreeData, seedData = otherTreeData.now$seedData,
                    xytree = otherTreeData.now$xytree, xytrap = otherTreeData.now$xytrap, priorDist = priors$priorDist, 
                    priorVDist = priors$priorVDist, minDist = priors$minDist, 
                    maxDist = priors$maxDist, minDiam = priors$minDiam, 
                    maxDiam = priors$maxDiam, maxF = priors$maxF)
predPlots <- plots
predList <- list(mapMeters = 10, plots = predPlots, years = years )

output <- mastif(formulaFec, formulaRep,  inputs = inputs, ng =20000, burnin = 1000, yearEffect = yearEffect,randomEffect = randomEffect,predList = predList)


save(output, file = paste0(spp.names[t],"output.rdata"))
            