
vec2mat <- function(xx, ROW=F){
  
  #if(ROW) make row vector
  
  if(is.matrix(xx))return(xx)
  
  cc <- names(xx)
  xx <- matrix(xx)
  rownames(xx) <- cc
  if(!ROW)xx <- t(xx)
  xx
}


.cleanRows <- function(xmat, xcol, STOP=F){
  
  ww <- which(duplicated(xmat[,xcol]))
  if(length(ww) == 0)return(xmat)
  
  if(STOP)stop(paste('duplicates in', xcol))
  warning(paste('removed duplicates in', xcol))
  xmat[-ww,]
}

gen4code <- function(xx, nn=4){
  
  #shorten genus name in genusSpecies string to nn characters
  
  FAC  <- F
  if(is.factor(xx))FAC <- T
  xx   <- as.character(xx)
  fc   <- sapply(gregexpr("[A-Z]",xx),'[',1)
  wf   <- which(fc > 5)
  if(length(wf) > 0){
    gen <- substr(xx[wf],1,4)
    spp <- substr(xx[wf],fc[wf],1000)
    xx[wf] <- columnPaste(gen,spp,'')
  }
  if(FAC)xx <- as.factor(xx)
  xx
}

combineMastData <- function(inputs, newlist=NULL, regions=NULL){
  
  # inputs and newlist both have inputs
  
  if(is.null(newlist)){
    
    inputs$treeData  <- inputs$treeData[inputs$treeData$region %in% regions,]
    inputs$specNames <- inputs$specNames[inputs$specNames %in% inputs$treeData$species]
    plots            <- sort(unique(as.character(inputs$treeData$plot)))
    inputs$seedData  <- inputs$seedData[inputs$seedData$plot %in% plots,]
    inputs$xytree    <- inputs$xytree[inputs$xytree$plot %in% plots,]
    inputs$xytrap    <- inputs$xytrap[inputs$xytrap$plot %in% plots,]
    
    ssum <- colSums( as.matrix(inputs$seedData[,inputs$seedNames]) )
    stmp <- inputs$seedData[,!colnames(inputs$seedData) %in% inputs$seedNames]
    
    inputs$seedData <- data.frame( cbind(stmp), inputs$seedData[,names(ssum)[ssum > 0]] )
    inputs$seedNames <- names(ssum)[ssum > 0]
    
    inputs$treeData <- cleanFactors(inputs$treeData)
    inputs$seedData <- cleanFactors(inputs$seedData)
    inputs$xytree <- cleanFactors(inputs$xytree)
    inputs$xytrap <- cleanFactors(inputs$xytrap)
    
    return(inputs)
    
  }
  
  newlist$treeData  <- newlist$treeData[newlist$treeData$region %in% regions,]
  plots             <- sort(unique(as.character(newlist$treeData$plot)))
  newlist$seedData  <- newlist$seedData[newlist$seedData$plot %in% plots,]
  newlist$xytree    <- newlist$xytree[newlist$xytree$plot %in% plots,]
  newlist$xytrap    <- newlist$xytrap[newlist$xytrap$plot %in% plots,]
  newlist$specNames <- newlist$specNames[newlist$specNames %in% newlist$treeData$species]
  
  treeData <- merge(inputs$treeData, newlist$treeData, by=colnames(inputs$treeData), all=T)
  treeData <- treeData[!duplicated(treeData),]
  xytree   <- merge(inputs$xytree, newlist$xytree, by=colnames(inputs$xytree), all=T)
  xytrap   <- merge(inputs$xytrap, newlist$xytrap, by=colnames(inputs$xytrap), all=T)
  
  
  sd1 <- inputs$seedData[,!names(inputs$seedData) %in% inputs$seedNames]
  sd2 <- newlist$seedData[,!names(newlist$seedData) %in% newlist$seedNames]
  
  
  w1 <- rownames(inputs$seedData)
  w2 <- rownames(newlist$seedData)
  w0 <- sort(unique(c(w1,w2)))
  
  sall <- sort(unique( c(inputs$seedNames, newlist$seedNames) ))
  
  scount <- matrix(0,length(w0),length(sall))
  rownames(scount) <- w0
  colnames(scount) <- sall
  
  for(k in 1:length(sall)){
    wi <- which(colnames(inputs$seedData) == sall[k])
    wn <- which(colnames(newlist$seedData) == sall[k])
    
    if(length(wi) == 1)scount[w1,sall[k]] <- scount[w1,sall[k]] + inputs$seedData[,wi]
    if(length(wn) == 1)scount[w1,sall[k]] <- scount[w2,sall[k]] + newlist$seedData[,wn]
  }
  scount <- scount[,colSums(scount) > 0]
  stmp   <- merge(sd1,sd2, by=colnames(sd1), all=T)
  
  plotTrapYr <- columnPaste(stmp$trapID,stmp$year)
  stmp <- stmp[!duplicated(plotTrapYr),]
  plotTrapYr <- columnPaste(stmp$trapID,stmp$year)
  
  rownames(stmp) <- plotTrapYr
  sc2 <- scount[rownames(stmp),]
  stmp <- cbind(stmp,sc2)
 # inputs$seedData <- stmp
  
  inputs$treeData  <- cleanFactors(treeData)
  inputs$seedData  <- cleanFactors(stmp)
  inputs$xytree    <- cleanFactors(xytree)
  inputs$xytrap    <- cleanFactors(xytrap)
  inputs$seedNames <- colnames(scount)
  inputs$specNames <- levels(treeData$species)
  
  inputs
}

rankNeighbors <- function(xy, xysite, dr = 5, nrank=10){
  
  # rmat is number of trees at distance (columns)
  # rows are rank: closest tree, 2nd closest, ...
  
  nsite  <- nrow(xysite)
  nrank  <- min( c(nsite, nrank) )
  breaks <- seq(0, 1000, by=dr)
  
  dj <- .distmat(xy[,'x'],xy[,'y'],xysite[,'x'], xysite[,'y']) 
  
  rr <- apply(dj, 2, order, decreasing=F)

  rmat <- matrix(NA, nrank, length(breaks)-1)
  rownames(rmat) <- paste('rank',c(1:nrank),sep='-')
  
  for(k in 1:nrank){
    
    loc <- cbind(rr[k,], 1:ncol(rr))
    tmp <- hist(dj[ loc ], breaks = breaks, plot=F)
    rmat[k,] <- tmp$density
  }
  rtot <- colSums(rmat,na.rm=T)
  rseq <- c(1:ncol(rmat))
  rseq[rtot == 0] <- 0
  rmat <- rmat[,1:max(rseq)]
  mh <- colSums(rmat)
  list(rmat = rmat, mh = mh, breaks = breaks[1:ncol(rmat)])
}


HMC <- function (ff, fMin, fMax, ep, L, tdat, sdat, ug, 
                 mu, sg, zz, R, SAMPR, distance, 
                 obsRows, obsYr, seedNames){
  
  #Hamiltonian Markov chain steps
  
  getU <- function(q, U = T){   # yq = log(fec)
    
    # for Hamiltonian
    # if U == T, return U, if U == F, return dU/df
    nseed <- ncol(R)
    fq <- exp(q)
    fq[fq > fMax] <- fMax[fq > fMax]
    fq[fq < fMin] <- fMin[fq < fMin]
    
    if(SAMPR){
      fq <- matrix(fq,length(fq),ncol=ncol(R))*R[drop=F,as.character(tdat$specPlot),]
    }else{
      fq <- matrix(fq,ncol=1)
    }
    
    uvec <- ug[ attr(distance,'group') ]
    dmat <- t(uvec/pi/(uvec + t(distance)^2)^2)
    dmat[dmat < 1e-8] <- 0
    
    lambda <- kernYrRcpp(dmat, fq*zz, years = obsYr, seedyear = sdat$year,
                         treeyear = tdat$year, seedrow = sdat$drow,
                         treecol = tdat$dcol)
    ss <- as.matrix(sdat[,seedNames])
    lambda[lambda < 1e-6] <- 1e-6
    
    if(U){
      mmat  <- matrix(0,max(tdat$plotyr),1)
      sprob <- -ss*log(lambda) + activeArea*lambda
      ii    <- rep(sdat$plotyr, nseed)
      tmp   <- .myBy(as.vector(sprob), ii, ii*0+1, summat=mmat, fun='sum')
      tprob <- 1/sg*(q - mu)^2 + tmp[tdat$plotyr]
      
      return( tprob ) 
    }
    
    kmat <- dmat[sdat$drow,tdat$dcol]
    smat <- -ss/lambda + activeArea
    
    if(nseed == 1){
      svec <- ff*colSums( kmat*as.vector(smat) )
    }else{
      svec <- rep(0,length(q))
      for(m in 1:nseed){
        sv <-  colSums(kmat*as.vector(smat[,m]) )*fq[,m]
        svec <- svec + sv
      }
    }
    svec + (q - mu)/sg
  }
  
  
#  logmax <- log(maxF)
  
  q <- log( ff )
  p <- currentP <- rnorm(length(q)) 
  
  activeArea <- sdat$area*sdat$active
  
  # half step for momentum 
  p <- p - ep*getU(q, U = F)/2
  
  # Alternate full steps for position and momentum
  wall <- 1:length(q)
  
  for (i in 1:L){
    q[wall] <- q[wall] + ep[wall]*p[wall]  
    if(i < L)p[wall] <- p[wall] - ep[wall]*getU(q, U = F)[wall] 
    wall <- which(q < log(fMax))
  }
  
  # Make a half step for momentum at the end.
  p <- p - ep*getU(q, U = F)/2
  
  # Negate momentum at end of trajectory to make the proposal symmetric
  p <- -p
  
  # Evaluate potential and kinetic energies at start and end of trajectory
  currentU  <- getU(log( ff ), U = T) 
  proposedU <- getU(q, U = T)
  currentK  <- currentP^2/2 
  proposedK <- p^2/2
  
  # Accept or reject the state at end of trajectory, returning either 
  # the position at the end of the trajectory or the initial position
  sp <- currentU - proposedU + currentK - proposedK
  ii <- tdat$plotyr
  
  pnow <- .myBy(sp*zz, ii, ii*0 + 1, fun='sum')
  
  a <- exp(pnow)
  wa <- which( runif(length(a)) < a)
  
  if(length(wa) > 0){
    wt <- which(tdat$plotyr %in% wa & zz == 1)
    ff[wt] <- exp( q[wt] )
    ep[wt] <- ep[wt]*1.1
    ep[-wt] <- ep[-wt]*.9
  }else{
    ep <- ep*.9
  }
  
  list(fg = ff, epsilon = ep, accept = length(wa))
}

msarLagTemplate <- function(plag, data, icol, jcol, gcol, ocol){
  
  # index for sampling betaYr
  # data - data.frame created by msarSetup
  # icol - individual column (integer or name)
  # jcol - time column (integer)
  # gcol - group column (integer or name)
  # ocol - indicator for observed (1) or interpolated (0)
  
  ifac <- gfac <- F
  idata <- data[,icol]
  jdata <- data[,jcol]
  gdata <- data[,gcol]
  odata <- data[,ocol]
  
  if(is.character(idata))idata <- as.factor(idata)
  if(is.character(gdata))gdata <- as.factor(gdata)
  
  if(is.factor(idata)){
    ifac   <- T
    indivs <- attr(idata,'levels')
    idata  <- match( as.character(idata), indivs )
  }else{
    indivs <- sort(unique(idata))
  }
  if(is.factor(gdata)){
    gfac   <- T
    groups <- attr(gdata,'levels')
    gdata  <- match( as.character(gdata), groups )
  }else{
    groups <- sort(unique(gdata))
  }
                   
  times  <- sort(unique(jdata))
  ngroup <- length(groups)
  nyr    <- length(times)
  lagGroup <- vector('list',ngroup)
  names(lagGroup) <- groups
  
  lagMatrix <-  numeric(0)
  lagGroup  <- numeric(0)
  
  #  OBSCOL <- F
  #  if('obs' %in% colnames(tdata))OBSCOL <- T
  
  nall <- 0
  
  for(m in 1:ngroup){
    
    wm   <- which(gdata == m & odata == 1)  #obs states required
    lagm <- numeric(0)
    
    im <- idata[wm]
    tall <- sort(unique(im))
    
    orm <- order(im,jdata[wm])
    
    imk <- match(im,tall)
    tmat <- matrix(NA,length(tall),nyr)
    tmat[ cbind(imk,jdata[wm]) ] <- wm
    rownames(tmat) <- tall
    
    for(j in (plag+1):nyr){
      
      jname <- paste(tall,times[j],sep='-')
      tj    <- tmat[,j:(j-plag),drop=F]
      rownames(tj) <- jname
      lagm <- rbind(lagm,tj)
    }
    wm <- unique( which(is.na(lagm),arr.ind=T)[,1] )
    if(length(wm) > 0)lagm <- lagm[-wm,]
  
    nall <- nall + nrow(lagm)
    
    colnames(lagm) <- paste('lag',c(0:plag),sep='-')
    lagm <- lagm[order(lagm[,1]),]
    lagMatrix <- rbind(lagMatrix,lagm)
    lagGroup  <- c(lagGroup,rep(m,nrow(lagm)))
    #   lagGroup[[m]] <- lagm
  }
  
  cat('\nNumber of full observations with AR model is: ')
  print( nall ) 
  
  if(nall < 20)stop('not enough observations for AR(p), try reducing p')
  
  #  lagGroup
  list(matrix = lagMatrix, group = lagGroup)
}

getVarType <- function(vnames, data, i, j){
  
  # 'i'  - individual variable
  # 'j'  - time variable
  # 'ij' - individual/time
  
  if(is.factor(i))i <- as.character(i)
  
  id <- sort(unique(i))
  ni <- length(id)
  yr <- sort(unique(j))
  ny <- length(yr)
  i <- match(i, id)
  j <- match(j, yr)
  
  vnew <- vector('list',length(vnames))
  names(vnew) <- vnames
  
  
  for(k in 1:length(vnames)){
    cj <- data[,vnames[k]]
    if(is.factor(cj))cj <- as.character(cj)
    mj <- matrix(NA, ni,ny)
    mj[ cbind(i, j) ] <- cj
    
 #   print(k)
 #   print(dim(mj))
    
    
    
    rr <- suppressWarnings( apply(mj,1,range,na.rm=T) )
    rc <- suppressWarnings( apply(mj,2,range,na.rm=T) )
    rowSame <- all(rr[1,] == rr[2,])
    colSame <- all(rc[1,] == rc[2,])
    if(is.na(rowSame))rowSame <- F
    if(is.na(colSame))colSame <- F
    if(!rowSame & !colSame)vnew[[k]] <- 'ij'
    if( rowSame &  colSame)vnew[[k]] <- 'i'
    if(!rowSame &  colSame)vnew[[k]] <- 'j'
    if( rowSame & !colSame)vnew[[k]] <- 'i'
  }
  
  vnew
}

msarSetup <- function(data, variables, plag, icol, jcol, gcol=NULL,
                      minGroup=10){
  
  # icol - column names in data for individual index; integer or factor
  # jcol - column name in data for time index; integer
  # gcol - column name for group index; integer or factor
  # gmat - matrix of individual columns to retain in output
  # plag - AR lag 
  # minGroup - minimum group size
  
  if(!is.data.frame(data))stop('data must be a data.frame')
  
  huge <- 1e+10
  
  data$obs <- 1
  
  if(is.factor(data[,icol]))data[,icol] <- droplevels(data[,icol])
  if(is.factor(data[,jcol]))data[,jcol] <- droplevels(data[,jcol])
  if(is.null(gcol)){
    gcol <- 'group'
    data$group <- rep(1,nrow(data))
  }else{
    if(is.factor(data[,gcol]))data[,gcol] <- droplevels(data[,gcol])
  }
  
  i <- as.character(data[,icol])
  j <- data[,jcol]
  g <- data[,gcol]
  
  io <- i
  if(is.factor(io))io <- as.character(io)
  oo <- apply(cbind(io,j),1,paste0,collapse='-')
  
  iall <- sort(unique(as.character(i)))
  jall <- c((min(j) - plag):(max(j) + plag))
  
  groups <- unique(g)
  
  baseMat <- base <- with(data, table(i, j) )
  base[base == 0] <- NA
  
  if( max(base,na.rm=T) > 1 )
    stop('must be a unique individual-time pair in data')
  
  ii <- match(i,iall)
  jj <- match(j,jall)  # original times
  ni <- length(iall)
  nj <- length(jall)
  
  #groups with sufficient times
  tmp <- table(g,jj)
  ntmp <- as.numeric(colnames(tmp))
  tmp[,ntmp < plag] <- 0
  tsum <- apply(tmp,1,max)
  

  cat('\nno. trees with > plag years, by group:\n')
  print(tsum)
  
  wlow <- which(tsum < minGroup)
  if(length(wlow) > 0)cat( paste('\nsmall group: ',names(tsum)[wlow],'\n',sep='') )
  
  back  <- matrix(0,nrow(base),plag)
  
  expMat <- zmat <- cbind(back, baseMat, back) #i by j, with plag
  colnames(expMat) <- colnames(zmat) <- jall
  
  imat <- jmat <- matrix(NA,ni,nj)
  imat[ cbind(ii, jj) ] <- c(1:nrow(data))
  jmat[ cbind(ii, jj) ] <- j
  first <- apply(jmat,1,which.min)
  last  <- apply(jmat,1,which.max)
  
  emat <- imat*0
  
  for(k in 1:ni)emat[k, (first[k]-plag):(last[k]+plag) ] <- 1
  ijFull <- which(emat == 1, arr.ind=T)
  
  vtypes    <- getVarType(colnames(data), data, ii, jj) 
  vtypes$obs <- 'ij'
  
  tpy <- as.character(columnPaste(data$treeID,data$year))
  
  tmp <- fillMissing(vtypes, data, icol, jcol, 
                      jtimes = jall, ijIndex = cbind(ii,jj), ijFull = ijFull)
  data <- tmp$data
  naVars <- tmp$naVars
  tpn    <- as.character(columnPaste(data$treeID,data$year))
  data$obs[!tpn %in% tpy] <- 0
  rm(tpy,tpn)
  
  i <- match(data[,icol],iall)
  j <- match(data[,jcol],jall)
  g <- match(data[,gcol],groups)
  
  io <- data[,icol]
  if(is.factor(io))io <- as.character(io)
 # odd <- apply(cbind(io,data[,jcol]),1,paste0,collapse='-')
 # obs[odd %in% oo] <- 1
  
  data <- cbind(data,i,j,g)
  data <- data[order(data$i, data$j),]
  
  wideGroup <- apply(with(data, table(i,g)),1,which.max) 
  
  ngroup <- length(groups)
  betaYr <- round( matrix(rnorm(ngroup*plag,0,.1),ngroup,plag), 3)
  rownames(betaYr) <- groups
  colnames(betaYr) <- paste('lag',c(1:plag),sep='-')
  
  list(xdata = data, times = jall,
       groupByInd = wideGroup, betaYr = betaYr, plag = plag)
}

fillMissing <- function(variables, data, icol, jcol, jtimes,
                        ijIndex = NULL, ijFull=NULL){
  
  # ijIndex is location in i by j matrix
  # ijFull is location in i by j matrix that includes added obs
  # jtimes - integer
  
  if(!is.data.frame(data))data <- as.data.frame(data)
  
  icc <- which(sapply(data,is.character))
  if(length(icc) > 0){
    for(k in icc)data[[k]] <- as.factor(data[[k]])
  }
  
  if(is.null(ijIndex)){
    ijIndex <- cbind(ii, jj)
    ii <- data[,icol]
    jj <- data[,jcol]
    if(is.factor(ii))ii <- as.character(ii)
    iall <- sort(unique(ii))
    ii   <- match(ii, iall)
    jj   <- match(jj, jtimes)
  }else{
    ii <- ijIndex[,1]
    jj <- ijIndex[,2]
    iall <- sort(unique(ii))
  }
  if(is.null(ijFull))ijFull <- ijIndex
  
  ny <- length(jtimes)
  ni <- max(iall)
  
  emat <- matrix(NA,ni, ny)
  newData <- vector('list',ncol(data))
  names(newData) <- names(data)
  ffact <- which( sapply(data, is.factor) )
  naVars <- character(0)
  
  for(k in 1:ncol(data)){           #expand to pre- and post-data
    
    jmat  <- matrix(NA,ni,ny,byrow=T)
    vtype <- variables[names(data)[k]]
    
    kvar <- data[,k]
    if(k %in% ffact)kvar <- as.character(kvar)
    jmat[ cbind(ii,jj) ]  <- kvar
    w0 <- which( is.na(jmat),arr.ind=T)
    w1 <- which(!is.na(jmat),arr.ind=T )
    
    if(vtype == 'i'){
      w1 <- w1[!duplicated(w1[,1]),]
      w2 <- w0
      w2[,2] <- w1[match(w2[,1],w1[,1]),2]
      jmat[ w0 ] <- jmat[ w2 ]
    }
    if(vtype == 'j'){
      w1 <- w1[!duplicated(w1[,2]),]
      w2 <- w0
      w2[,1] <- w1[match(w2[,2],w1[,2]),1]
      jmat[ w0 ] <- jmat[ w2 ]
      if(colnames(data)[k] == jcol)jmat <- matrix(jtimes,ni,ny, byrow=T)
    }
    if(vtype == 'ij'){            # use trend
      if( is.numeric(kvar) ){
        rr   <- range(jmat,na.rm=T)
        if(rr[1] < .1)rr[1] <- .1
        if(rr[2] <= rr[1])rr[2] <- rr[1] + .01
        jmat <- .interpRows(jmat, INCREASING=F,minVal=rr[1],maxVal=rr[2],
                            defaultValue=NULL,tinySlope=.00001) 
      }else{
        naVars <- c(naVars, colnames(data)[k])
      }
    }
    
    ktmp <- jmat[ijFull]
    if( k %in% ffact )ktmp <- as.factor(ktmp)
    newData[[k]] <- ktmp
  }
  
  xdata <- data.frame(newData)
  rownames(xdata) <- NULL
  xdata <- xdata[order(xdata[,icol],xdata[,jcol]),]
  list(data = xdata, naVars = naVars)
}

.propZ <- function(znow, last0first1, matYr){
  
  # repr - known repr from tdata
  # random walk proposal
  
  new <- matYr + sample( c(-1:1), nrow(znow), replace=T)
  
  new[last0first1[,'all0'] == 1] <- 0
  new[last0first1[,'all1'] == 1] <- ncol(znow) + 1
  
  ww  <- which(new < last0first1[,1])
  new[ww] <- last0first1[ww,1]
  
  ww  <- which(new > last0first1[,2])
  new[ww] <- last0first1[ww,2]
  
  down <- which(new < matYr & new > 0)
  znow[ cbind(down,new[down]) ] <- 1   # advance 1 year
  
  up <- which(new > matYr & new < ncol(znow))
  znow[ cbind(up,matYr[up]) ] <- 0     # delay 1 year
  
  znow[last0first1[,'all0'] == 1, ] <- 0
  znow[last0first1[,'all1'] == 1, ] <- 1
  
  list(zmat = znow, matYr = new)
}

.setupYear <- function(ylist, data, jtimes){
  
  specNames <- lagGroup <- NULL
  nyr <- length(jtimes)
  
  yrnames <- names(ylist)
  mcol    <- unlist(ylist)
  mcol    <- mcol[names(mcol) %in% c('specGroups','plotGroups')]
  yeGroup <- data[,mcol]
  if(length(mcol) > 1){
    if(length(mcol) > 2)stop('only two groups available for yearEffects')
    cy <- cbind(as.character(yeGroup[,1]),as.character(yeGroup[,2]))
    yeGroup <- apply( cy, 1, paste0, collapse='-')
  }
  yeGr   <- sort(as.character(unique(yeGroup)))      
  ygr    <- match(yeGroup, yeGr)
  yyr    <- match(data$year,jtimes)
  ngroup <- length(yeGr)
  data$group <- yeGroup
    
  yrIndex <- cbind(ygr, yyr)
  colnames(yrIndex) <- c('group','year')
  gy     <- apply( yrIndex, 1, paste0,collapse='-')
  gyall  <- sort(unique(gy))
  groupYr <- match(gy,gyall)
  yrIndex <- cbind(yrIndex,groupYr)
  rownames(yrIndex) <- gy
  
  betaYr <- .myBy(ygr*0+1, yrIndex[,'group'], yrIndex[,'year'], fun='sum')
  rownames(betaYr) <- yeGr
  colnames(betaYr) <- jtimes
  
  
  rmax <- apply(betaYr,1,max)
  wr   <- which(rmax < 10)
  if(length(wr) > 0)cat('\nat least one year group has < 10 trees\n')
  
  # reference class removed if specGroups are species
  yeRef <- yeGr
  dash <- grep('-', yeGr)
  if(length(dash) > 0){
    yeRef <- matrix( unlist(strsplit(yeGr,'-')), ncol=2,byrow=T) 
    wr    <- match(yeRef[,1],specNames)
    yeRef <- yeGr[ !duplicated(wr)  ]
  }
  betaYr <- betaYr*0
  
 # lagGroup <- numeric(0)
  
#  if(AR){
#    lagGroup <- .lagSetup(ylist$p, ngroup, yrIndex, data, years, yeGr)
#    betaYr <- betaYr[,1:ylist$p,drop=F]*0
#    colnames(betaYr) <- colnames(lagm)[-1]
#  }
    
  list(betaYr = betaYr, yrRef = yeGr, yrIndex = yrIndex, yeGr = yeGr,
       tdata = data)
}

.lagSetup <- function(plag, ngroup, yrInd, tdata, yr, yeGr){
  
  # tdata must have treeID, year           
  
  nyr <- length(yr)
  lagGroup <- vector('list',ngroup)
  
  OBSCOL <- F
  if('obs' %in% colnames(tdata))OBSCOL <- T
  
  for(m in 1:ngroup){
    
    if(!OBSCOL){
      toc <- rep(1,nrow(yrInd))
    }else{
      toc <- tdata$obs
    }
    
    wm   <- which(yrInd[,'group'] == m & toc == 1)  #all states required
    lagm <- numeric(0)
    
    tall <- sort(unique(tdata$treeID[wm]))
    tmat <- matrix(NA,length(tall),nyr)
    ti   <- match(tdata$treeID[wm],tall)
    yi   <- match(tdata$year[wm],yr)
    tmat[ cbind(ti,yi) ] <- wm
    rownames(tmat) <- tall
    
    for(j in (plag+1):nyr){
      
      jname <- paste(tall,yr[j],sep='-')
      tj    <- tmat[,j:(j-plag)]
      rownames(tj) <- jname
      lagm <- rbind(lagm,tj)
    }
    wm <- unique( which(is.na(lagm),arr.ind=T)[,1] )
    if(length(wm) > 0)lagm <- lagm[-wm,]
    
    colnames(lagm) <- paste('lag',c(0:plag),sep='-')
    lagGroup[[m]] <- lagm
  }
  names(lagGroup) <- yeGr

  lagGroup
}

.boxCoeffsLabs <- function( boxPars, labels, colLabs = NULL, cex=1, 
                            xadj = 0){
  
  ncols <- length(labels)
  if(is.null(colLabs))colLabs <- rep('black',ncols)
  
  xaxp <- par('xaxp')
  da <- diff(xaxp[1:2])
 # at <- xadj + (da/ncols/2 + xaxp[1] + c(0:ncols)*da/(xaxp[3] + 1))[1:ncols]
  
  at <- (xaxp[1] + c(0:ncols)*da/ncols)[-1]
  at <- at - diff(at)[1]/2
  
  yends <- boxPars$stats[c(1,nrow(boxPars$stats)),]
  yfig  <- par('usr')[3:4]
  dy    <- diff(yfig)
  dends <- rbind(yfig[1] - yends[1,], yfig[2] - yends[2,])
  sides <- apply( abs(dends), 2, which.max)
  wt  <- which(sides == 1)
  if(length(wt) > 0)text(at[wt],yends[1,wt] - dy/20,labels[wt], 
                         offset = -.1,
                         col=colLabs[wt], pos=2, srt=90, cex = cex)
  wt  <- which(sides == 2)
  if(length(wt) > 0)text(at[wt],yends[2,wt] + dy/20,labels[wt], 
                         offset = -.1,
                         col=colLabs[wt], pos=4, srt=90, cex = cex)
}

summaryData <- function(tdata, sdata, seedNames, xytree, xytrap, 
                        plotDims=NULL, plotArea=NULL){
  
  plots <- as.character(sort(unique(tdata$plot)))
  nplot <- length(plots)
  years <- sort(unique(tdata$year))
  nyr   <- length(years)
  
  if(is.null(plotArea)){
    plotArea <- matrix( plotDims[,'area']/10000, nplot, nyr)
    colnames(plotArea) <- years
    rownames(plotArea) <- plots
  }
    
  
 # cy <- as.character(years)
  
  ba    <- pi*(tdata$diam/2)^2
  
  summat <- matrix(0,nplot,nyr)
  
  ii <- match(as.character(tdata$plot),plots)
  jj <- match(tdata$year,years)
  tmp <- .myBy(ba,ii,jj,summat=summat,fun='sum')/10000  # m2 per plot
  colnames(tmp) <- years
  rownames(tmp) <- plots
  
  yr <- as.character( years[years %in% colnames(plotArea)] )
  
  baPlotYr <- round(tmp[,yr]/plotArea[drop=F,plots,yr],6)
  
  seedCount <- as.matrix(sdata[,seedNames,drop=F])
  
  ii <- match(as.character(sdata$plot),plots)
  jj <- match(sdata$year,years)
  ii <- rep(ii,length(seedNames))
  jj <- rep(jj,length(seedNames))
  aa <- rep(sdata$active*sdata$area,length(seedNames))
  tmp <- .myBy(as.vector(seedCount)/aa,ii,jj,summat = summat,fun='sum')
  colnames(tmp) <- years
  rownames(tmp) <- plots
  
  seedPerBA <- tmp[drop=F,plots,yr]/baPlotYr[drop=F,plots,yr]
  
  keep <- which(colSums(baPlotYr,na.rm=T) > 0)
  
  list(BA = baPlotYr[,keep,drop=F], seedPerBA = seedPerBA[,keep,drop=F])
}

getBins <- function(xx, nbin=15){
  
  bins <- unique( quantile(as.vector(xx[xx > 0]),seq(0,1,length=nbin)) )
  bins <- c(0,bins)
  bins[-1] <- bins[-1] - diff(bins)/2
  bins <- c(bins,max(bins)+1)
  db   <- diff(bins)
  w0   <- which(db < 1)
  if(length(w0) > 0)bins <- bins[-(w0+1)]
  sort(unique(bins))
}

.cov2Cor <- function(covmat, covInv = NULL){  
  
  # covariance matrix to correlation matrix
  # if covInv provided, return inverse correlation matrix
  
  d    <- nrow(covmat)
  di   <- diag(covmat)
  s    <- matrix(di,d,d)
  cc   <- covmat/sqrt(s*t(s))
  
  if(!is.null(covInv)){
    dc <- diag(sqrt(di))
    ci <- dc%*%covInv%*%dc
    return(ci)
  }
  cc
}

mastPlot <- function(output, plotPars = NULL){
  
  SAVEPLOTS <- F
  
  data <- treeData <- priorUgibbs <- fecPred <- plotDims <- 
    plotHaByYr <- plotArea <- NULL
  
  seedNames <- specNames <- R <- formulaFec <- formulaRep <- xfec <- xrep <-
    tdata <- seedData <- xytree <- xytrap <- distall <- ng <- burnin <- nplot <-
    ntree <- ntrap <- nyr <- maxF <- bfec <- brep <- upar <- rgibbs <- 
    betaFec <- betaRep <- rMu <- rSe <- usigma <- fecMu <- fecSe <- matrMu <- 
    seedPred <- inputs <- chains <- parameters <- predictions <- 
    upars <- dpars <- trueValues <- betaYrMu <- betaYrSe <-  
    sgibbs <- ugibbs <- omegaE <- predPlots <- betaYrRand <- 
    betaYrRandSE <- prediction <- eigenMu <- facLevels <- specPlots <- NULL
  randGroups <- formulaRan <- rnGroups <- reIndex <- xrandCols <- NULL  
  specGroups <- plotGroups <- yrIndex <- randomEffect <- yearEffect <- NULL
  YR <- AR <- RANDOM <- TV <- SAMPR <- PREDICT <- SPACETIME <- F
  pacfMat <- pacfSe <- pacsMat <- pacsSe <- obsRows <- NULL
  modelYears <- seedPredGrid <- treePredGrid <- NULL
  
  outFolder <- 'mastPlots'
  yeGr <- NULL
  plotsPerPage <- 4
  MAPS <- T
  
  for(k in 1:length(output))assign( names(output)[k], output[[k]] )
  for(k in 1:length(inputs))assign( names(inputs)[k], inputs[[k]] )
  for(k in 1:length(chains))assign( names(chains)[k], chains[[k]] )
  for(k in 1:length(parameters))assign( names(parameters)[k], parameters[[k]] )
  for(k in 1:length(prediction))assign( names(prediction)[k], prediction[[k]] )
  if('arList' %in% names(data)){
    for(k in 1:length(data$arList))
      assign( names(data$arList)[k], data$arList[[k]] )
    AR <- T
    plag <- ncol(data$arList$betaYr)
  }
  if(!is.null(plotPars)){
    for(k in 1:length(plotPars))assign( names(plotPars)[k], plotPars[[k]] )
    if( 'trueValues' %in% names(plotPars)  ){
      TV <- T
      for(k in 1:length(trueValues))
        assign( names(trueValues)[k], trueValues[[k]] )
    }
  }
  if('trueValues' %in% names(inputs)){
    TV <- T
    for(k in 1:length(inputs$trueValues))
      assign( names(inputs$trueValues)[k], inputs$trueValues[[k]] )
  }
  
  maxF      <- output$inputs$maxF
  specNames <- output$data$setupData$specNames
  seedNames <- output$data$setupData$seedNames
  if('rgibbs' %in% names(output$chains))SAMPR <- T
  bfec   <- .orderChain(bfec, specNames)
  brep   <- .orderChain(brep, specNames)
  if('randomEffect' %in% names(output$inputs)){
    RANDOM <- T
    for(k in 1:length(randomEffect))
      assign( names(randomEffect)[k], randomEffect[[k]] )
    agibbs <- .orderChain(agibbs, specNames)
  }
  if('yearEffect' %in% names(output$inputs)){
    yrIndex <- output$data$setupYear$yrIndex
    if(is.null(yrIndex))yrIndex <- output$inputs$yrIndex
    if('p' %in% names(output$inputs$yearEffect)){
      AR <- T
    }else{
      YR <- T
      yeGr <- output$data$setupYear$yeGr
    }
    ngroup <- length(yeGr)
    bygibbsR <- .orderChain(bygibbsR, yeGr)
    ugibbs   <- .orderChain(ugibbs, yeGr)
  }
  
  if(SAVEPLOTS){
    tt <- paste('\nPlots saved to ',outFolder,'/\n',sep='')
    cat(tt)
  }
  
  if(AR)YR <- F
  
  xmean <- output$data$setupData$xmean  # to unstandardize xfec, xrep
  xsd   <- output$data$setupData$xsd
  
  tdata <- data$setupData$tdata
  sdata <- data$setupData$sdata
  
  
  ###############
  rm(treeData)
  ##############
  
  
  if(is.null(yeGr))yeGr <- 'all'
  
  nspec  <- length(specNames)
  ntype  <- length(seedNames)
  years  <- sort(unique(tdata$year))
  nyr    <- length(years)
  ngroup <- length(yeGr)
  plots  <- sort(unique(as.character(tdata$plot)))
  nplot  <- length(plots)
  
  tmp1 <- as.character(tdata$plot)
  tmp2 <- tdata$year
  if(!is.null(seedPredGrid)){
    PREDICT <- T
    tmp1 <- c(tmp1, as.character(seedPredGrid$plot))
    tmp2 <- c(tmp2, seedPredGrid$year)
  }
  
  plotYrTable <- table(plot = tmp1, year = tmp2)
  rm(tmp1)
  rm(tmp2)

  
  cfun <- colorRampPalette( c('#66c2a5','#fc8d62','#8da0cb') )
  specCol <- cfun(nspec)
  names(specCol) <- specNames
  cols <- specCol
  
  
  gfun <- colorRampPalette( c("#8DD3C7", "#BEBADA", "#FB8072",
                              "#80B1D3", "#FDB462") )
  groupCol <- gfun(ngroup)
  names(groupCol) <- yeGr
  
  if(ngroup == 1 & ncol(ugibbs) > 1){
    groupCol <- gfun(ncol(ugibbs))
    names(groupCol) <- colnames(ugibbs)
  }
                     
  gfun <- colorRampPalette( c("forestgreen","#8DD3C7", "#BEBADA", "#FB8072",
                              "#80B1D3", "#FDB462", "brown") )
  plotCol <- gfun(nplot)
  names(plotCol) <- plots
  
  if(SAVEPLOTS){
    ff <- file.exists(outFolder)
    if(!ff)dir.create(outFolder)
  }
  
  # empirical summary
  
 # pdim <- getPlotDims(xytree,xytrap)
  
  tmp <- summaryData(tdata, sdata, seedNames, xytree, xytrap, plotDims, 
                     plotArea)
  
  BA <- as.matrix( tmp$BA[drop=F,plots,] )
  seedPerBA  <- as.matrix( tmp$seedPerBA[drop=F,plots,] )
  BA[!is.finite(BA)] <- NA
  seedPerBA[!is.finite(seedPerBA)] <- NA
  
  pyr <- as.character(years[years %in% colnames(BA) & 
                              years %in% colnames(seedPerBA)])
  
  graphics.off()
  
  if(SAVEPLOTS)pdf( file=.outFile(outFolder,'dataSummary.pdf') )
  
  par(bty='n', mfrow=c(1,1),mar=c(5,5,3,3),oma=c(2,2,1,2), cex=1.2)
  
  xlim <- range(BA[,pyr],na.rm=T)
  xlim[2] <- 1.1*xlim[2]
  ylim <- range(seedPerBA[,pyr],na.rm=T)
  ylim[2] <- ylim[2]*1.2
  
  xlab <- expression( paste('BA (',plain(m)^2, plain(ha)^-1,')' ) )
  ylab <- expression( paste('Seeds (', plain(m)^-2, 'BA ', plain(yr)^-1,')') )
  zlab <- expression( paste(bar(y) %+-% plain(sd), ' (for non-zero values)') )
  
  plot(NULL,xlim=xlim,ylim=sqrt(ylim),yaxt='n', xlab=xlab, ylab=ylab, yaxt='n')
  
  tt   <- sqrtSeq(1.2*sqrt(ylim[2]))
  at   <- tt$at
  labs <- tt$labs
  axis(2, at = at, labels = labs)
  
  syr <- seedPerBA[,pyr,drop=F]
  qs <- quantile(syr[syr > 0],probs = pnorm(c(-1,0,1)),na.rm=T)
  abline(h=sqrt(qs),lty=2,lwd=1)
  abline(h=sqrt(qs[2]),lwd=1)
  
  for(j in 1:nrow(BA)){
    points(BA[j,pyr],sqrt(syr[j,]),
           col=.getColor(plotCol[rownames(BA)[j]],.5),pch=16)
    b1 <- BA[j,pyr]
    xj <- .1*diff(xlim) + quantile(b1[b1 > 0], .7, na.rm=T)
    yj <- 1.2*max(sqrt(syr[j,]), na.rm=T)
    if(xj > xlim[2])xj <- xlim[2]*.95
    if(yj > ylim[2])yj <- ylim[2]*.95
    if(yj < diff(sqrt(ylim))/5)yj <- yj + diff(sqrt(ylim))/5
    text(xj,yj, rownames(BA)[j],srt=45, col=plotCol[rownames(BA)[j]])
  }
  mtext('Quantial for non-zeros',outer=T,side=4, line=-2)
  mtext('mean (solid), 68% (dashed)',outer=T,side=4, line=-1)
  
  
  if(!SAVEPLOTS){
    readline('seeds vs BA by yr-- return to continue ')
  } else {
    dev.off( )
  }
  
  ########### MCMC chains
    
  graphics.off()
  
  refVals <- NULL
  
  if(TV)refVals <- betaRep
  .chainPlot(brep, burnin, 'maturation', 
             refVals = refVals, SAVEPLOTS, outFolder)
  
  if(TV)refVals <- betaFec
  .chainPlot(bfec, burnin, 'fecundity', 
             refVals = refVals, SAVEPLOTS, outFolder)
 
  if(TV)refVals <- upar
  .chainPlot(ugibbs, burnin, 'dispersal parameter u', 
             refVals = refVals, SAVEPLOTS, outFolder)
  
  if(ngroup > 1){
  .chainPlot(priorUgibbs, burnin, 'dispersal mean and variance', 
             refVals = NULL, SAVEPLOTS, outFolder)
  }
 
  .chainPlot(sgibbs, burnin, 'variance sigma', refVals = NULL, 
             SAVEPLOTS, outFolder)
 
  if(SAMPR){
    
    mg   <- rgibbs
    posR <- attr(rMu,'posR')
    #  wc   <- which(apply(mg, 2, FUN=var) > 0)
    rff  <- NULL
    
    #  if(length(wc) > 0){
    
    #   mg  <- mg[,wc,drop=F]
    
    tmp <- columnSplit(colnames(mg),'_')
    seedCols <- tmp[,1]
    
    wdash <- grep('-',tmp[1,2])
    if(length(wdash) > 0){
      tmp  <- columnSplit(tmp[,2],'-')
      specCols <- tmp[,1]
      plotCols <- tmp[,2]
      np <- length(plots)
    }else{
      specCols <- tmp[,2]
      plotCols <- NULL
      np <- 1
    }
    
    for(j in 1:np){
      if(!is.null(plotCols)){
        wj <- which(plotCols == plots[j])
        if(length(wj) == 0)next
        rj <- range(mg[,wj])
        if(rj[1] == 1 | rj[2] == 0)next
      }else{
        wj <- 1:ncol(mg)
      }
      
      label <- paste('M matrix', plots[j])
      
      if(TV){
        tt <- columnSplit(colnames(mg)[wj])
        rff <- trueValues$R[ cbind(tt[,2],tt[,1]) ]
      }
      
      .chainPlot(mg[,wj,drop=F], burnin, label, ylim=c(0,1),
                 refVals = rff, SAVEPLOTS, outFolder)
    }
    
    graphics.off()
    
    if(SAVEPLOTS)pdf( file=.outFile(outFolder,'Mpars.pdf') )
    
    espec <- sort(unique(specCols))
    
    par(mfrow=c(length(espec),1), bty='n', mar=c(.5,4,.5,3), oma=c(4,3,2,4))
    tlab <- ''
  
    plotColors <- plotCol[plotCols]

    for(j in 1:nspec){
      
      cj <- which(specCols  == specNames[j])
      if(length(cj) == 0)next
      
      unkn <- grep('UNKN',colnames(mg)[cj])
      if(length(unkn) > 0){
        cj <- c(cj[-unkn],cj[unkn])
      }
      
      cols <- plotColors[cj]
      
      boxPars <- .boxCoeffs(mg[burnin:ng,cj,drop=F], specNames[j], xlab = tlab,
                 ylab=specNames[j], addSpec='', ylim=c(0,1),
                 cols = cols, yaxt='n')
      wt <- c(match(seedCols[cj], seedNames), 1000)
      wt[length(wt)] <- wt[length(wt)] + 1
      wt <- which(diff(wt) != 0)
      abline(v=wt+.5,col='grey')
      text(wt-.5,.5,seedCols[cj[wt]], srt=90, cex=.8)
      
      abline(h=1,col='grey', lwd=1)
      axis(2, at=c(0,1), las=2)
      tlab <- ''
    }
    
    mtext('Seed types', side=3, line=.4, outer=T, cex=1.2)
    mtext('Species', side=2, line=.4, outer=T, cex=1.2)
    
    ncol <- round(length(plots)/4) + 1
    
    cornerLegend('bottomright', plots, text.col = plotCol[plots],
                 cex=.9, bty='n', ncol=ncol)
      
    if(!SAVEPLOTS){
      readline('species -> seed type -- return to continue ')
    } else {
      dev.off( )
    }
    
    
  graphics.off()
  
  ########################
  
  
  #inverse probability species h| seedtype m
  
  fec <- output$prediction$fecPred
  
  nsim <- 40
  ksim <- sample(c(burnin:ng), nsim, replace=T)
  
  npairs <- columnSplit(colnames(rgibbs),'_')[,c(2,1)]
  
  rmat <- sprob <- sprob2 <- parameters$rMu*0
  
  for(m in 1:length(plots)){
    
    wj <- grep( paste('-',plots[m],sep=''), colnames(rgibbs))
    
    if(length(wj) <= 1)next
    mm <- rgibbs[,wj,drop=F]
    mp <- npairs[wj,]
    
    wm <- which(fec$plot == plots[m])
    ff <- fec$fecEstMu[wm]
    ss <- fec$fecEstSe[wm] + .00001
    mt <- fec$matrEst[wm]
    
    rrow <- grep(plots[m],rownames(rmat))
    
    for(k in 1:nsim){
      
      rmat <- rmat*0
      rmat[npairs[wj,]] <- mm[ksim[k],]
      rmm  <- rmat[drop=F,rrow,]
      
      fk <- .tnorm(length(ff), 0, maxF, ff, ss)*mt
      ii <- match(fec$species[wm],specNames)
      tf <- tapply(ff, ii, FUN=sum)
      tf <- tf/sum(tf)
      names(tf) <- paste(specNames[as.numeric(names(tf))],plots[m],sep='-')
      tmat <- rmm*0
      tmat[names(tf),] <- rep(tf, ncol(tmat))
      
      sf <- rmm*tmat/matrix( colSums(rmm*tmat),nrow(rmm), ncol(rmm), byrow=T )
      sf[is.na(sf)] <- 0
      sprob[rrow,] <- sprob[rrow,] + sf
      sprob2[rrow,] <- sprob2[rrow,] + sf^2
    }
  }
  seed2SpecMu <- sprob/nsim
  sse <- sprob2/nsim - seed2SpecMu^2
  sse[sse < 1e-30] <- 0
  seed2SpecSe <- sqrt( sse )
  
  if(max(seed2SpecSe) > 1e-20){
    
    if(SAVEPLOTS)pdf( file=.outFile(outFolder,'undiffSeed.pdf') )
    
    kplots <- character(0)
    
    for(m in 1:length(plots)){
      rrow <- grep(plots[m],rownames(seed2SpecMu))
      smm  <- seed2SpecMu[rrow,]
      wk   <-  which(smm > 1e-20 & smm < 1, arr.ind=T)
      if(length(wk) == 0)next
      kplots <- c(kplots,plots[m])
    }
    
    tt <- .getPlotLayout(length(kplots))
    par(mfrow=tt$mfrow, bty='n', mar=c(3,3,1,1), oma=c(3,3,1,1))
    aspec <- character(0)
    
    for(m in 1:length(kplots)){
      
      rrow <- grep(kplots[m],rownames(seed2SpecMu))
      smm  <- seed2SpecMu[rrow,]
      emm  <- seed2SpecSe[rrow,]
      wk   <-  which(smm > 1e-20 & smm < 1, arr.ind=T)
      if(length(wk) == 0)next
      
      cspec <- columnSplit(rownames(wk),'-')[,1]
      aspec <- c(aspec,cspec)
      colm <- specCol[cspec]
      tmp  <- barplot(smm[wk], beside=T, col=colm, border=colm, ylim=c(0,1),
                      yaxt='n')
      labels <- F
      if(m %in% tt$left)labels=c(0,1)
      axis(2, c(0,1), labels=labels)
      segments(tmp, smm[wk], tmp, smm[wk] + 1.96*emm[wk], lwd = 1.5, col=colm)
      segments(tmp-.1, smm[wk] + 1.96*emm[wk], tmp+.1, smm[wk] + 1.96*emm[wk],
               lwd=1.5, col=colm)
      .plotLabel(kplots[m], 'topright', above=T)
    }
    
    mtext('Species contribution', side=1, line=.4, outer=T, cex=1.2)
    mtext('Fraction of unknown seed type', side=2, line=.4, outer=T, cex=1.2)
    
    aspec <- sort(unique(aspec))
    
    cornerLegend('bottomright', aspec, text.col = specCol[aspec],
                 cex=1.1, bty='n', ncol=1)
    
    if(!SAVEPLOTS){
      readline('Species to undiff seed -- return to continue ')
    } else {
      dev.off( )
    }
  }
  }
  
  ##############################
  if(RANDOM){
    
    aMu <- parameters$aMu
    
    att <- aMu*0
    wc <- 1
    if(length(aMu) > 1){
      diag(att) <- 1
      wc <- which(att == 1)
    }
    
    vaa <- agibbs[,wc,drop=F]
    vrr <- apply(vaa,2,sd)
    if( max(vrr) > 1e-5 ){
      
      .chainPlot(agibbs[,wc,drop=F], burnin, 'random effects variance', 
                 refVals = NULL, SAVEPLOTS, outFolder)
    }
    
    graphics.off()
  }
  
  ########### coefficients
  
  if(SAVEPLOTS)pdf( file=.outFile(outFolder,'fecundityCoeffs.pdf') )
  
  betaFec <- parameters$betaFec
  betaRep <- parameters$betaRep
  fnames  <- .coeffNames(rownames(betaFec))
  rnames  <- .coeffNames(rownames(betaRep))
  
  xfec <- output$data$setupData$xfec
  xrep <- output$data$setupData$xrep
  Qf   <- ncol(xfec)/nspec
  Qr   <- ncol(xrep)/nspec
  
  if(nspec == 1){
    
    par(mfrow=c(1,2), mar=c(5,4,2,1), bty='n')
    tmp <- .boxplotQuant( brep[drop=F,burnin:ng,], add = F, xaxt = 'n',
                          xlim = NULL, ylab='', 
                          outline=F, col=.getColor('black', .2), 
                          border='black', lty=1, boxfill=NULL )
    
    axis(1, at = c(1:nrow(betaRep)), labels=rnames, las = 2)
    abline(h=0, col='grey', lwd=2, cex.axis=.6)
    title('b) Maturation')
    if(TV)abline(h=trueValues$betaRep, lwd=2, lty=2, col='grey')
    
    tmp <- .boxplotQuant( bfec[drop=F,burnin:ng,], add = F, xaxt = 'n',
                          xlim = NULL, ylab='Parameter value', 
                          outline=F, col=.getColor('black', .2), 
                          border='black', lty=1, boxfill=NULL )
    axis(1, at = c(1:nrow(betaFec)), labels=fnames, las = 2)
    abline(h=0, col='grey', lty=2)
    title('a) Fecundity')
    if(TV){
      abline(h=trueValues$betaFec, lwd=2, lty=2, col='grey')
      text(2,trueValues$betaFec[1] + 2,'true values')
    }
    
  }else{
    
    par(mfrow=c(2,2), mar=c(1,3,1,.2), bty='n', oma=c(2,2,1,1))
    nint <- grep(':',colnames(chains$bfec))
    
    intt <- c(1:ncol(bfec))
    if(length(nint) > 0)intt <- intt[-nint]
    
    bc <- bfec[burnin:ng,intt,drop=F]
    colnames(bc) <- .replaceString(colnames(bc),'species','')
    
    boxPars <- .boxCoeffs(bc, specNames, xaxt='n',
                          xlab = 'Fecundity', ylab='', 
                          cols = specCol[specNames], addSpec = '')
    .boxCoeffsLabs( boxPars, specNames, specCol[specNames], cex=.7)
    
    if(length(nint) > 0){
      bc <- bfec[burnin:ng,nint,drop=F]
      colnames(bc) <- .replaceString(colnames(bc),'species','')
      scol <- rep(character(0), ncol(bc))
      
      bcc <- .boxCoeffs(bc, specNames, xlab = '', xaxt='n',
                        ylab=' ', cols = specCol[specNames], addSpec = '' )
      gg <- grep(':',fnames)
      cnames <- fnames[gg]
      cnames <- columnSplit(cnames,':', LASTONLY=T)
      cnames <- cnames[seq(1,length(cnames),by=2)]
      stats  <- bcc$stats
      snew   <- numeric(0)
      if(length(cnames) == 1){
        .plotLabel(cnames,'bottomleft', cex=.8)
      }else{
        for(k in 1:length(cnames)){
          skk <- stats[,grep(cnames[k], colnames(stats), fixed=T), drop=F]
          mins <- apply(skk,1,min)
          maxs <- apply(skk,1,max)
          sk  <- skk[,1]*0
          sk[1] <- mins[1]
          sk[length(sk)] <- maxs[length(sk)]
          snew <- cbind(snew,sk)
        }
        colnames(snew) <- cnames
        bcc$stats <- snew
        .boxCoeffsLabs( bcc, cnames, cex=.7)
      }
    }
    
    nint <- grep(':',colnames(chains$brep))
    intt <- c(1:ncol(chains$brep))
    intt <- c(1:ncol(chains$brep))
    if(length(nint) > 0)intt <- intt[-nint]
    
    bc <- brep[burnin:ng,intt,drop=F]
    colnames(bc) <- .replaceString(colnames(bc),'species','')
    
    boxPars <- .boxCoeffs(bc, specNames, xaxt='n',
                          xlab = 'Maturation', ylab='', 
                          cols = specCol[specNames], addSpec = '')
    mtext('intercepts', 1, cex=.9)
    
    if(length(nint) > 0){
      bc <- brep[burnin:ng,nint,drop=F]
      colnames(bc) <- .replaceString(colnames(bc),'species','')
      for(k in 1:nspec){
        wk <- grep(specNames[k], colnames(bc))
        colnames(bc)[wk] <- specNames[k]
      }
      bcc <- .boxCoeffs(bc, specNames, xlab = '', xaxt='n',
                        ylab=' ', cols = specCol[specNames], addSpec = '')
      gg <- grep(':',rnames)
      cnames <- rnames[gg]
      cnames <- columnSplit(cnames,':', LASTONLY=T)
      cnames <- cnames[seq(1,length(cnames),by=2)]
      stats  <- bcc$stats
      snew   <- numeric(0)
      if(length(cnames) == 1){
        mtext(cnames, 1, cex=.9)
      }else{
        for(k in 1:length(cnames)){
          gk  <- grep(cnames[k], colnames(stats), fixed=T)
          if(length(gk) == 0)next
          skk <- stats[,, drop=F]
          mins <- apply(skk,1,min)
          maxs <- apply(skk,1,max)
          sk  <- skk[,1]*0
          sk[1] <- mins[1]
          sk[length(sk)] <- maxs[length(sk)]
          snew <- cbind(snew,sk)
        }
        if(length(snew) > 0){
          colnames(snew) <- cnames
          bcc$stats <- snew
        }
        .boxCoeffsLabs( bcc, cnames, cex=.7)
      }
    }
    mtext('Posterior estimate',side=2, outer=T)
  }
  
  if(!SAVEPLOTS){
    readline('fecundity, maturation -- return to continue ')
  } else {
    dev.off( )
  }
  
  ############ maturation
  
  
  if(SAVEPLOTS)pdf( file=.outFile(outFolder,'maturation.pdf') )
  
  mfrow <- .getPlotLayout(nspec)$mfrow
  par(mfrow=mfrow, bty='n', mar=c(3,3,1,1),  oma=c(1,1,1,1.1))
  
  nsim <- 1000
  lseq <- seq(0, 1, length=100)

  
  for(j in 1:nspec){
    
    fspec <- fecPred[obsRows,]
    tspec <- tdata[obsRows,]
    xspec <- xrep[obsRows,]
    pspec <- tdata$plot[tspec$species == specNames[j]]
    xspec <- xspec[tspec$species == specNames[j],]
    fspec <- fspec$matrEst[fspec$species == specNames[j]]
    dspec <- tspec$diam[tspec$species == specNames[j]]
     
    pspec <- droplevels(pspec)
    pall  <- sort(unique(pspec))
    
    dq   <- round( quantile(dspec, lseq), 1)
    dq   <- unique(dq)
    
    ttt  <- nn2(dspec, dq, k = 1)[[1]]
    www  <- which(duplicated(ttt))
    if(length(www) > 0)ttt <- ttt[-www]
    
    ntt  <- length(ttt)                          # tree yr
    ksamp <- sample(nrow(brep),nsim,replace=T)   # MCMC row
    bcols <- colnames(brep)
    if(nspec > 1)bcols <- bcols[ grep(specNames[j], bcols) ]
    rk    <- matrix(NA, ntt, nsim)
    
    plot(NULL, xlim=c(0, max(dspec,na.rm=T)),ylim=c(0,1),xlab='',
         ylab='')

    for(k in 1:nsim){
      rk[,k] <- pnorm( xspec[ttt,bcols]%*%matrix(brep[ksamp[k],bcols],ncol=1))
    }
    qr <- apply(rk, 1, quantile, c(.5, .025, .975) )
    
    .shadeInterval( dspec[ttt], t(qr[2:3,]) , col=
                                     .getColor('black', .1) )
    lines(dspec[ttt], qr[1,], col='white', lwd=5)
    lines(dspec[ttt], qr[1,], lwd=2)
      
     wp <- match(as.character(pspec),plots)
     points(dspec,fspec, pch=16, col='white', cex=1)
     points(dspec,fspec, pch=16, col=.getColor(plotCol[wp],.4),
            cex=.7)
     jplot <- plots[plots %in% pspec]
     legend('bottomright',jplot,text.col=plotCol[jplot],bty='n')
     title(specNames[j])
  }
  mtext('Diameter (cm)', side=1, line=0, outer=T)
  mtext('Probability', side=2, line=0, outer=T)
  
  if(!SAVEPLOTS){
    readline('maturation by diameter -- return to continue ')
  } else {
    dev.off( )
  }
  
  ############ dispersal
  
  udat <- ugibbs[burnin:ng,,drop=F]
  
  if(colnames(udat)[1] %in% yeGr){
    cols <- groupCol[yeGr]
  }
  if(yeGr[1] == 'all'){
    cols <- groupCol[colnames(udat)]
  }
  if(colnames(udat)[1] %in% specNames){
    cols <- specCol[specNames]
  }
  
  upars <- parameters$upars
  
  if(ncol(udat) > 1){
    
    if( is.null(rownames(upars)) )rownames(upars) <- 'mean'
    
    if(SAVEPLOTS)pdf( file=.outFile(outFolder,'dispersalCoeffs.pdf') )
    
    ncc <- ncol(udat)
    sideMar <- min( c(3, 3 + 10/ncc) )
    
    par(mfrow=c(1,1), mar=c(3,5,5,4), bty='n', oma=c(2,sideMar,1,sideMar))
 
    ylim <- range(udat)
    ylim[1] <- .75*ylim[1]
    ylim[2] <- 1.25*ylim[2]
    
    ylab1 <- expression( paste(hat(u),' (', plain(m)^2, ')') )
    ylab2 <- expression( paste('Kernel mean ', bar(d), ' (m)') )
    
    tmp <- .boxplotQuant( udat, xaxt='n', add = F, xlim = NULL, 
                          ylab=ylab1, ylim = ylim, boxwex=.6,
                          outline=F, col=.getColor(cols,.3), 
                          border=cols, xaxt='n',lty=1, boxfill=NULL )
    abline(h=upars['mean',3:4], lty=2, lwd=2, col='grey')
    .boxCoeffsLabs( tmp, names(cols), cols, cex=.8)
    
    rr <- range( pi*sqrt(udat)/2 )
    rm <- -1
    by <- 10
    if(diff(rr) < 20){
      by <- 5
      rm <- 0
    }
    if(diff(rr) < 10){
      by <- 2
      rr[1] <- rr[1] - 1
      rr[2] <- rr[2] + 1
    }
    
    rr <- round(rr, rm)
    rr <- seq(rr[1],rr[2], by=by)
    
    axis(4, at = (2*rr/pi)^2, labels = rr )
    mtext(ylab2, 4, line=0, outer=T)
    abline(v=par('usr')[1:2])
    
    if(!SAVEPLOTS){
      readline('dispersal by group -- return to continue ')
    } else {
      dev.off( )
    }
  }
  
  graphics.off()
  
  if(SAVEPLOTS)pdf( file=.outFile(outFolder,'seedShadow.pdf') )
  
  # standardize betas
  xfec  <- data$setupData$xfec
  xfecU <- data$setupData$xfecU
  
  xfecu2s <- output$data$setupData$xfecu2s
  bchainStand <- bfec%*%xfecu2s
  
  # common diameter across groups
  
  qd     <- .75
  
  dcol   <- grep('diam',colnames(xfec))
  
  diam60 <- xfec[,dcol]
  diam60 <- diam60[diam60 != 0]
  diam60 <- quantile(diam60, qd)
  
  dtcol <- grep('diam',colnames(tdata))
  lab60 <- tdata[,dtcol]
  lab60 <- lab60[lab60 != 0]
  lab60 <- quantile(lab60, qd)
  lab60 <- signif(lab60,1)
  
  # each group at mean for other predictors
  xbar   <- numeric(0)
  rnames <- character(0)
  ucol   <- numeric(0)
  
  for(j in 1:nspec){
    
    wj <- which(tdata$species == specNames[j])
    
    if(ngroup > 1){
      
      for(g in 1:ngroup){
        wg <- wj
        if(ngroup > 1){
  #        wg   <- which(as.character(tdata$group) %in% yeGr[g])
          wg   <- which(yrIndex[,'group'] == g)
          wg   <- intersect(wj, wg)
        }
        if(length(wg) == 0)next
        xmu  <- colMeans(xfec[drop=F,wg,])
        
        w0 <- which(xmu != 0)       #insert diameter value
        w0 <- intersect(w0, dcol)
        xmu[w0] <- diam60
        
        xbar <- rbind(xbar, xmu)
        rnames <- c(rnames, yeGr[g])
        
        wu <- which(colnames(ugibbs) == yeGr[g])
        ucol <- c(ucol, wu)
      }
      
    }else{
      xmu  <- colMeans(xfec[drop=F,wj,])
      
      w0 <- which(xmu != 0)       #insert diameter value
      w0 <- intersect(w0, dcol)
      xmu[w0] <- diam60
      
      xbar <- rbind(xbar, xmu)
      rnames <- c(rnames, specNames[j])
      ucol   <- c(ucol, 1)
    }
  }
  rownames(xbar) <- rnames
  
  nsim <- 2000
  ns  <- 100
  buffer <- 15
  
  mm <- 6*max(dpars)
  if(ngroup > 1)mm <- 6*max(dpars['mean',])
  
  dseq <- 10^seq(0,log10(mm),length=ns) - 1
  dseq <- matrix( dseq, ns, nsim)

  
  ij <- sample(burnin:ng, nsim, replace=T)
  
  
  ssList <- numeric(0)
  maxy   <- 0
  
  ff <-  bchainStand[ij,colnames(xbar)]%*%t(xbar)
  ff[ff > log(maxF)] <- log(maxF)
  maxy <- numeric(0)
  
  kss <- 1:(ns-buffer)
  keepSeq <- c(rev(kss), kss) 
  dss  <- c(-rev(dseq[kss,1]),dseq[kss,1])
  
  for(k in 1:nrow(xbar)){
    
    uj <- ugibbs[ij,ucol[k]]
    sj <- sgibbs[ij,'sigma']

    kj <- uj/pi/(uj + dseq^2)^2
    kj <- kj*matrix(exp(ff[,k] + sj/2), ns, nsim, byrow=T)
    kj[!is.finite(kj)] <- NA
    qj <- t( apply(kj, 1, quantile, c(.5, .05, .95), na.rm=T ) )

    for(m in 1:3)qj[,m] <- runmed(qj[,m], k = buffer, endrule='constant')
    maxy <- c(maxy, max( qj[,1] ))
    
    ssList <- append(ssList, list( qj[keepSeq,] ) )
  }
  
  maxy[maxy > 1e+10] <- 1e+10
  maxy[maxy < 1e-10] <- 1e-10
  names(ssList) <- rownames(xbar)
  
  par(bty='n', mar=c(5,5,1,1))
  
  labSeeds <- expression( paste('Seeds (', plain(m)^-2, ')') )
  
  rmax <- diff( range(log10(maxy)) )
  SQRT <- F
  if( rmax > 3 ){
    SQRT <- T
    smax <- sqrt(max(maxy))
    plot(NULL, xlim=range(dss), ylim=c(0,smax),
         xlab='Distance (m)', ylab = labSeeds, yaxt='n')
    tt   <- sqrtSeq(1.2*smax)
    at   <- tt$at
    labs <- tt$labs
    axis(2, at = at, labels = labs)
    
  }else{
    plot(NULL, xlim=range(dss), ylim=c(0,1.2*max(maxy)),
         xlab='Distance (m)', ylab = labSeeds)
  }
  
  for(k in 1:length(ssList)){
    if(maxy[k] <= 1e-10 | maxy[k] >= 1e+10)next
    ss <- ssList[[k]][,2:3]
    
    
    if(SQRT)ss <- sqrt(ss)
    .shadeInterval( dss, ss , col=
                      .getColor(cols[names(cols)[k]], .2) )
  }
  mc <- numeric(0)
  for(k in 1:length(ssList)){
    ss <- ssList[[k]][,1]
    if( max(ss,na.rm=T) >= 1e+10)next
    if(SQRT)ss <- sqrt(ss)
    lines(dss, ss, col='white',lwd=6)
    lines(dss, ss, col=cols[names(cols)[k]], lwd=2)
    mc <- c(mc, max(ssList[[k]][,1]))
  }

  if(nspec > 1 | ngroup > 1){
    ord <- order(mc, decreasing=T)
    legend('topright',names(cols)[ord],
                                   text.col=cols[names(cols)[ord]], bty='n')
  }
  legend('topleft', paste(lab60,' cm diameter tree',sep=''), bty='n')
  
  if(!SAVEPLOTS){
    readline('seed shadow -- return to continue ')
  } else {
    dev.off( )
  }
  
  ########### random coefficients
  
  if(RANDOM){
    
    if(SAVEPLOTS)pdf( file=.outFile(outFolder,'randomCoeffs.pdf') )
    
    xrandCols  <- output$data$setupRandom$xrandCols
    formulaRan <- output$data$setupRandom$formulaRan
    
    betaFec <- parameters$betaFec[,1]
    names(xrandCols) <- rownames(parameters$betaFec)[xrandCols]
    alphaMu <- parameters$alphaMu
    intercept <- paste('species',specNames,sep='')
    if(nspec == 1)intercept <- '(Intercept)'
    slopes    <- names(xrandCols)[!names(xrandCols) %in% intercept]
    
    npp <- length(slopes)/length(intercept)
    
    npp <- 1
    
    mfrow <- c(1,1)
    if(npp > 1)mfrow <- c(2,2)
    
    par(mfrow=mfrow, bty='n', mar=c(4,4,1,1), cex=1.2)
    
    icol <- match(intercept, colnames(alphaMu))
    icol <- icol[is.finite(icol)]
    intercept <- intercept[is.finite(icol)]
    
    xlim <- range(betaFec[icol] + alphaMu[,icol])
    
    if(length(slopes) == 0){
      
      rlim <- range(alphaMu)
      breaks <- seq(rlim[1]-.1, rlim[2]+.1,length=20)
      hvals <- xvals <- numeric(0)
      for(s in 1:length(specNames)){
        ra <- alphaMu[,icol[s]]
        ra <- ra[ra != 0]
        tmp <- hist(ra, plot=F, breaks = breaks)
        hvals <- cbind(hvals,tmp$density)
        xvals <- cbind(xvals,tmp$mids)
      }
      ylim <- range(hvals, na.rm=T)
      xl <- range(betaFec,na.rm=T) + range(rlim, na.rm=T)
      xl[1] <- xl[1] - 1
      xl[2] <- xl[2] + 1
      
      plot(NA,xlim=xl,ylim = ylim, xlab='log fecundity', ylab='frequency')
      for(s in 1:length(specNames)){
        ws <- range( which(hvals[,s] > 0) )
        ss <- ws[1]:ws[2]
        xs <- betaFec[icol[s]] + xvals[ss,1]
        xs <- c(xs[1],xs,xs[length(xs)])
        ys <- c(0, hvals[ss,s], 0)
        lines( xs, ys,type='s',col=specCol[s], lwd=2)
        lines( xs, xs*0,col=specCol[s], lwd=2)
      }
    }else{

      ww <- grep(':', slopes)
      
      if(length(ww) > 0){
        slopeLab  <- matrix( unlist( strsplit( slopes, ':')), ncol=2,byrow=T)[,2]
      }else{
        slopeLab <- slopes
      }
      jcol <- match(slopes, colnames(alphaMu))
      ylim <- range(betaFec[jcol] + alphaMu[,jcol])
      
      if(length(icol) == 1)jcol <- jcol[1]
      
      plot(betaFec[icol],betaFec[jcol],col=.getColor(specCol,.7), 
           xlim=xlim, ylim=ylim, cex=.8, xlab='Intercept',ylab=slopeLab[1])
      abline(h = 0, lwd=2, lty=2, col='grey')
      abline(v = 0, lwd=2, lty=2, col='grey')
      
      for(s in 1:length(specNames)){
        points(betaFec[icol[s]] + alphaMu[,icol[s]],
               betaFec[jcol[s]] + alphaMu[,jcol[s]],
               col=.getColor(specCol[s],.3), pch=16, cex=1)
      }
    }
    legend('topright', specNames, text.col = specCol,bty='n')
    
    if(!SAVEPLOTS){
      readline('fixed plus random effects -- return to continue ')
    } else {
      dev.off( )
    }
  }
  
  ########### predicted seed counts, true fecundity
   
  if(SAVEPLOTS)pdf( file=.outFile(outFolder,'seedPrediction.pdf') )
  
  obsRowSeed <- inputs$obsRowSeed
  
  xs <- rowSums(sdata[obsRowSeed,seedNames,drop=F])
  pcols <- grep('predMean',colnames(seedPred))
  ecols <- grep('estMean',colnames(seedPred))
  
  ys <- seedPred[,pcols]
  zs <- seedPred[,ecols]
  
  ww <- which(is.finite(xs))
  xs <- xs[ww]
  ys <- ys[ww]
  zs <- zs[ww]

  ylim <- range(c(ys, zs),na.rm=T)
  xlim <- range(xs,na.rm=T) + 1

  mfrow <- c(1,2)
  title <- 'a) From posterior distribution'
  
  if(TV)mfrow <- c(1,3)
  
  par(mfrow=mfrow, mar=c(4,4,2,1), bty='l')
  
  bins <- getBins(xs, nbin=20)
  nbin <- length(bins)
  
  opt <- list(log=F, xlabel='Observed', bins = bins,
              nbin=nbin, ylabel='Predicted', col='forestgreen', 
              ylimit=ylim, xlimit = xlim, SQRT=T)
  tmp <- .plotObsPred(xs, ys, opt = opt)
  .plotLabel(title, above=T, cex=.8)
  abline(0,1,lty=2)

  opt <- list(log=F, xlabel='Observed', bins = bins, atx = tmp$atx, labx = tmp$labx,
              aty = tmp$aty, laby = tmp$laby,
              ylabel='', col='forestgreen', ylimit=ylim, xlimit = xlim,
              nbin=nbin, SQRT=T)
  tmp <- .plotObsPred(xs, zs, opt = opt)
  .plotLabel('b) From fecundity estimate', above=T, cex=.8)
  abline(0,1,lty=2)
  
  if(TV){
    xs <- inputs$trueValues$fec
    ys <- prediction$fecPred$fecEstMu
    ylim <- quantile(ys[ys > 1],c(.02,1))
    ylim[1] <- max(c(ylim[1],1))
    
    ws <- which(xs > 0 & ys > 0)
    xlim <- quantile(xs[ws],c(.02,.99))
    
    bins <- getBins(xs)
    nbin <- length(bins)
    
    opt <- list(xlimit = xlim, ylimit = ylim, bins = bins,
                nbin=nbin, log=F, xlabel='True values', 
                ylabel='Estimates', col='darkgreen', SQRT=T)
    .plotObsPred(xs, ys, opt = opt)
    .plotLabel('c) Fecundity prediction', above=T, cex=.8)
    abline(0,1,lty=2)
    abline(h=mean(xs,na.rm=T), lty=2)
  }
  
  if(!SAVEPLOTS){
    readline('prediction -- return to continue ')
  } else {
    dev.off( )
  }
  
  # true parameter values
  
  if(TV){
    
    if(SAVEPLOTS)pdf( file=.outFile(outFolder,'trueParameters.pdf') )
    
    par(mfrow=c(1,2),bty='n', mar=c(4,4,3,1))
    bfec <- output$parameters$betaFec
    brep <- output$parameters$betaRep
    plot( trueValues$betaFec, bfec[,1], ylim=range(bfec), 
          xlab='True value', ylab='Estimate', pch=3)
    abline(0,1, lwd=2, lty=2, col=.getColor('black',.4))
    suppressWarnings(
    arrows( trueValues$betaFec, bfec[,3], trueValues$betaFec, bfec[,4], lwd=2,
            angle=90, length=.05, code=3)
    )
    .plotLabel('a) Fecundity parameters', above=T, cex=1)
    
    plot( trueValues$betaRep, brep[,1], ylim = range(brep), xlab='True value', 
          ylab=' ', pch=3)
    abline(0, 1, lwd=2, lty=2, col=.getColor('black',.4))
    suppressWarnings(
    arrows( trueValues$betaRep, brep[,3], trueValues$betaRep, brep[,4], lwd=2,
            angle=90, length=.05, code=3)
    )
    .plotLabel('b) Maturation parameters', above=T, cex=1)
    
    if(!SAVEPLOTS){
      readline('parameter recovery -- return to continue ')
    } else {
      dev.off( )
    }
  }
  
  ############# predicted maps
  
  if(MAPS){
    
    treeSymbol <- fecPred$fecEstMu
    if(is.null(treeSymbol))treeSymbol <- treeData$diam
    
    graphics.off()
    
    mpp <- 1:length(plots)
    
    plotDims <- as.matrix(plotDims)
    
    for(m in mpp){
      
      xlim <- plotDims[plots[m],c('xmin','xmax')]
      ylim <- plotDims[plots[m],c('ymin','ymax')]
      
      dx <- diff(xlim)
      dy <- diff(ylim)
      ratio <- dx/dy
      
      pyr <- plotYrTable[plots[m],]
      pyr <- as.numeric( colnames(plotYrTable)[pyr > 0] )
      
      mfrow <- c(2, 2)
      if( max(c(dx,dy)) > 100){
        if(ratio > 2) mfrow <- c(2,1)
        if(ratio < .5)mfrow <- c(1,2)
      }
      
      nperPage <- prod(mfrow)
      
      yrm <- years[years %in% pyr]
      ny  <- length(yrm)
      
      k   <- 0
      add <- F
      o   <- 1:nperPage
      o   <- o[o <= nyr]
      
      seedMax <- sdata[sdata$plot == plots[m] & sdata$year %in% yrm ,seedNames]
      if(is.matrix(seedMax))seedMax <- rowSums(seedMax, na.rm=T)
      seedMax <- max(seedMax, na.rm=T) + 1
      
      while( max(o) <=  nyr & length(o) > 0 ){
        
        if(length(o) < 5)mfrow <- c(2,2)
        
        yr <- yrm[o]
        yr <- yr[is.finite(yr)]
        if(length(yr) == 0)break
        
        graphics.off()
        
        file <- paste('map_',plots[m],'_',years[o[1]],'.pdf',sep='')
        
        if(SAVEPLOTS)pdf( file=.outFile(outFolder,file) )
        
        add <- mastMap(output, treeSymbol=treeSymbol, mapPlot = plots[m], 
                       xlim = xlim, ylim = ylim, PREDICT=T, mapYears = yr, 
                       treeScale = 1, trapScale = .5, mfrow=mfrow, seedMax = seedMax,
                       COLORSCALE=F, LEGEND=T)
        
        if(add)scaleBar('m', value = 20, yadj=.07)
        
        if(!SAVEPLOTS){
          readline('predicted fecundity, seed data -- return to continue ')
        } else {
          dev.off()
        }
        
        o <- o + nperPage
        o <- o[o <= nyr]
        
        
        if(length(o) == 0){
          break
        }
        
        if(!add)next
        
  #      if(!SAVEPLOTS){
  #        readline('predicted fecundity, seed data -- return to continue ')
  #      } else {
  #        dev.off()
  #      }
      }  
    }
  }
  
  ################# year effects
  
  if( 'betaYrMu' %in% names(parameters) ){
    
    graphics.off()
    
    file <- paste('yearEffect.pdf',sep='')
    
    if(SAVEPLOTS)pdf( file=.outFile(outFolder,file) )
    
  #  YR <- T
    if('lagGroup' %in% names(inputs))AR <- T
    
    betaYrMu <- parameters$betaYrMu
    betaYrSe <- parameters$betaYrSe
    
    if(nrow(betaYrMu) == 1 & length(specNames) == 1)
      rownames(betaYrMu) <- rownames(betaYrSe) <- yeGr[1]
    
    if('betaYrRand' %in% names(parameters)){  #combine fixed and random
      
      betaYrRand   <- parameters$betaYrRand
      betaYrRandSE <- parameters$betaYrRandSE
      betaYrRandSE[is.na(betaYrRandSE)] <- 0
      bmu <- bsd <- betaYrRand*0
      
      if(!AR){
        ttab <- table( yeGr[yrIndex[,'group']], tdata$year )
      }else{
        ttab <- betaYrRand
      }
      ttab[ttab != 0] <- 1
      

      for(k in 1:nrow(betaYrRand)){
        bmu[k,] <- betaYrMu + betaYrRand[k,]
        bsd[k,] <- sqrt(betaYrSe^2 + betaYrRandSE[k,]^2)
        bmu[k,] <- bmu[k,]*ttab[rownames(betaYrRand)[k],]
        bsd[k,] <- bsd[k,]*ttab[rownames(betaYrRand)[k],]
      }
      betaYrMu <- bmu
      betaYrSe <- bsd
      betaYrSe[is.na(betaYrSe)] <- 0
    }

    
    if(AR){
      par(mfrow=c(1,2),bty='n', mar=c(5,4,2,1))
      xlab <- 'lag (yr)'
      yr <- c(1:plag)
      mr <- .5
    }
    
    
    xtt <- seq(1900,2100,by=5)                #reference xtick
    
    par(mfrow=c(1,1),bty='n', mar=c(4,4,2,5), mai=c(1,1,1,1.1))
    if(AR){
      yr <- xtick <- 1:plag
      xlab <- 'lag (yr)'
    }
    if(YR & !AR){
      yr <- xtick <- years
      xlab <- ''
      if(length(yr) > 10)xtick <- xtick[xtick %in% xtt]
    }
    mr  <- max(betaYrMu + betaYrSe, na.rm=T)
 
    
    ylim = c(-2*mr,mr)
    
    plot(NULL, xlim = range(yr), ylim = ylim, xlab = xlab, 
         ylab = 'log fecundity', xaxt='n')
    axis(1, xtick)
    abline(h = 0, lty=2, col='grey', lwd=2)
 #   abline(v = xtick, lty=2, col='grey', lwd=2)
    leg <- character(0); col <- numeric(0)
    
    if(!AR & !YR){
      loHi <- cbind( betaYrMu[1,] - 1*betaYrSe[1,],
                     betaYrMu[1,] + 1*betaYrSe[1,])
      .shadeInterval(yr,loHi,col='black',PLOT=T,add=T, trans = .3)
      lines(yr, betaYrMu[1,], col=.getColor('white',.7), lwd=5)
      lines(yr, betaYrMu[1,], col='grey', lwd=2)
    }
    col <- numeric(0)
    
    if(!is.null(betaYrRand)){
      betaYrRand   <- betaYrRand[drop=F,yeGr,]
      betaYrRandSE <- betaYrRandSE[drop=F,yeGr,]
    }
    
    for(j in 1:ngroup){
      nj <- yeGr[j]
      wj <- which(is.finite(betaYrMu[nj,]) & betaYrMu[nj,] != 0) 
   #   if(length(wj) < 3)next
      loHi <- cbind( betaYrMu[nj,wj] - 1*betaYrSe[nj,wj],
                     betaYrMu[nj,wj] + 1*betaYrSe[nj,wj])
      .shadeInterval(yr[wj],loHi,col=groupCol[nj],PLOT=T,add=T, trans = .3)
    }
    
    for(j in 1:ngroup){
      nj <- yeGr[j]
      wj <- which(is.finite(betaYrMu[nj,]) & betaYrMu[nj,] != 0) 
   #   if(length(wj) < 3)next
      lines(yr[wj], betaYrMu[nj,wj], col=.getColor('white',.7), lwd=5)
      lines(yr[wj], betaYrMu[nj,wj], col=groupCol[nj], lwd=2)
    }
    if(ngroup > 1)cornerLegend('topright', yeGr, text.col = groupCol[yeGr],
                               cex=.8, bty='n')
    if(AR){
      title('AR coefficients',adj=0, font.main=1, font.lab=1,
                cex.main=.9)
      
      if(!SAVEPLOTS){
        readline('lag effect groups -- return to continue ')
      } else {
        dev.off()
      }
      
      file <- paste('countByYr.pdf',sep='')
      if(SAVEPLOTS)pdf( file=.outFile(outFolder,file) )
    }
  
    if(YR)title('Year effects, +/- 1 se',adj=0, font.main=1, font.lab=1,
                cex.main=.9)
    
    s2p <- as.matrix( sdata[,seedNames] )
    nseed <- length(seedNames)
    
    if(!AR)par(new = T)
    ymax <- 2*max(sqrt(s2p), na.rm=T)
    plot(jitter(sdata[,'year']),sqrt( s2p[,1]), cex=.8, axes=F,
         xlab='Year', ylim=c(0,ymax), ylab=NA, pch=16, 
         col=.getColor('black',.3))
    if(nseed > 1){
      for(j in 2:nseed){
        points(jitter(sdata[,'year']),sqrt( s2p[,j]), cex=.6, 
       pch=16, col=.getColor('black',.3))
      }
    }
    ii <- match(sdata$year,years)
    ii <- rep(ii, ncol(s2p))
    jj <- rep(1, length=length(s2p))
    ssum <- .myBy(as.vector( sqrt(s2p) ), ii, jj, 
                   summat = matrix(0, length(years), 1), fun='sum')
    nsum <- .myBy(as.vector( s2p*0 + 1 ), ii, jj, 
                  summat = matrix(0, length(years), 1), fun='sum')
    
    sall <- ssum/nsum
    lines(years, sall, lwd=5, col='white')
    lines(years, sall, lwd=2, col='black')
    
    side <- 4
    if(AR)side <- 2
    
    at <- c(0:100)^2
    at <- at[at < ymax]
    
    axis( side = side, at=at, labels=at )
    mtext(side = side, line = 3, 'seed count')
    
    if(!SAVEPLOTS){
      readline('year effect groups -- return to continue ')
    } else {
      dev.off()
    }
    
    if('plotGroups' %in% names(yearEffect)){  # by region
      
      graphics.off()
      
      file <- 'yearEffectByRegion.pdf'
      if(AR)file <- 'lagEffectByRegion.pdf'
      
      if(SAVEPLOTS)pdf( file=.outFile(outFolder,file) )
      
      yeGr <- output$data$setupYear$yeGr
      if(is.null(yeGr))yeGr <- output$data$setupData$arList$yeGr
      region <- yearEffect$plotGroups
      spec   <- yearEffect$specGroups
   #   if(!'specGroups' %in% names(yearEffect))spec   <- character(0)
   #   if(!'plotGroups' %in% names(yearEffect))region <- character(0)
      
      wd <- grep('-',yeGr)
      if(length(wd) > 0){
        spec   <- columnSplit(yeGr,sep='-')[,1]
        region <- columnSplit(yeGr,sep='-')[,2]
      }else{
        if(length(spec) > 0)spec <- yeGr
        if(length(region) > 0)region <- yeGr
      }
      regs <- unique(region)
      nreg <- length(regs)
      
      if(nreg == 0){
        regs <- spec
        nreg <- length(spec)
      }
      
      par(mfrow=c(nreg,1), bty='n', mar=c(2,2,.1,2), oma=c(2,3,1,1))  
      
      for(k in 1:nreg){
        
        wk <- which(region == regs[k])
        
        plot(NULL, xlim = range(yr), ylim = c(-2,2), xlab = xlab, 
             ylab = '', xaxt='n',yaxt='n')
        if(k == nreg)axis(1, xtick)
        axis(2,c(-2,0,2))
        abline(h = 0, lty=2, col='grey', lwd=2)
     #   abline(v = xtick, lty=2, col='grey', lwd=2)
        
        wll <- numeric(0)
        
        for(j in wk){
          nj <- yeGr[j]
          wj <- which(is.finite(betaYrMu[nj,]) & betaYrMu[nj,] != 0)
     #     if(length(wj) < 3)next
          wll <- c(wll,j)
          loHi <- cbind( betaYrMu[nj,wj] - 1*betaYrSe[nj,wj],
                         betaYrMu[nj,wj] + 1*betaYrSe[nj,wj])
          .shadeInterval(yr[wj],loHi,col=groupCol[nj],PLOT=T,add=T, trans = .3)
        }
        
        for(j in wk){
          nj <- yeGr[j]
          wj <- which(is.finite(betaYrMu[nj,]) & betaYrMu[nj,] != 0)
          if(length(wj) < 3)next
          lines(yr[wj], betaYrMu[nj,wj], col=.getColor('white',.7), lwd=5)
          lines(yr[wj], betaYrMu[nj,wj], col=groupCol[nj], lwd=2)
        }
        legend('bottomleft', yeGr[wll], text.col = groupCol[yeGr[wll]],
                                   cex=1.2, bty='n',)
        .plotLabel(region[j],'topleft')
      
      
      if(AR){
        mtext('Lag',side=1,line=1,outer=T)
      }else{
        mtext('Year',side=1,line=1,outer=T)
      }
      mtext('log fecundity',side=2,line=1,outer=T)
      }
      
      if(!SAVEPLOTS){
        readline('year effect groups -- return to continue ')
      } else {
        dev.off()
      
      }
    }
  }
  
  ############# pacf
  
  nyy <- ceiling(nyr*.6)
  if(nyy > 10)nyy <- 10
  
  if(nyr > 3 & sum(pacfMat,na.rm=T) != 0){
    
    graphics.off()
    
    file <- paste('pacf.pdf',sep='')
    if(SAVEPLOTS)pdf( file=.outFile(outFolder,file) )
    
    par(mfrow=c(1,2),bty='n', mar=c(3,2,1,.5), oma=c(2,3,1,1))
    
    mr  <- .9
    ylim <- range(pacfMat[,-1],na.rm=T)
    ylim <- range( c(ylim, pacsMat[-1]), na.rm=T )
    
    plot(NULL, xlim = c(1,nyy), ylim = ylim, xaxt = 'n',
         xlab = '', ylab = '')
    axis(1, at=c(1:nyy))
    abline(h = 0, lty=2, col='grey', lwd=2)
    leg <- character(0); col <- numeric(0)
    lag <- c(0:(nyr-1))
    
    if(yeGr[1] %in% specNames & specNames[1] %in% rownames(pacfMat)){
      pacfMat <- pacfMat[drop=F,specNames,]
      cols <- specCol[specNames]
      leg  <- specNames
    }else{
      cols <- groupCol
      leg <- rownames(pacfMat)
    }
    
    pacCol <- gfun(nrow(pacfMat))
    names(pacCol) <- rownames(pacfMat)
    
    for(j in 1:nrow(pacfMat)){
      
      wj <- which( is.finite(pacfMat[j,]))
      wj <- wj[-1]                   # omit lag zero
      nj <- length(wj)
      if(nj < 2)next
      
      ac <- pacfMat[j,wj]
      
      lines(lag[wj], ac, col=pacCol[j], lwd=2)
      loHi <- cbind( ac - 1.96*pacfSe[j,wj],
                     ac + 1.96*pacfSe[j,wj])
      .shadeInterval(lag[wj],loHi,col=pacCol[j],PLOT=T,add=T, 
                     trans = .2)
      
      up <- which(loHi[,1] > 0)
      up <- up[up > 1]
      points(lag[wj[up]],ac[up],cex=1.3,pch=16,col=.getColor(pacCol[j],.5))
      
      up <- up[up < 10]
      up <- paste0( up[up > 2], collapse=', ')
      ll <- rownames(pacfMat)[j]
          leg <- c(leg,ll)
          col <- c(col,j)
    }
    if(length(leg) > 1)legend('topright', leg, text.col = pacCol, 
                              bty='n', cex=.6)
    .plotLabel('a) log Fecundity', location='topleft',above=T, cex=.9)
    
    xlim <-  c(1,nyy)
    plot(NULL, xlim = xlim, ylim = ylim, xaxt = 'n', yaxt='n',
         xlab = '', ylab = '')
    axis(1, at=c(1:nyy))
    axis(2, labels=F)
    abline(h=0, lty=2,lwd=2, col='grey')
    leg <- character(0); col <- numeric(0)
    lag <- c(0:nyy)
    
    wj <- which(pacsMat != 0)
    nj <- length(wj)
    if(nj > 2){
      
      ac <- pacsMat[wj]
      lines(lag[wj], ac, col=1, lwd=2)
      segments( lag[wj], ac*0, lag[wj], ac)
    }
    .plotLabel('b) Seed counts', location='topleft', above=T, cex=.9)
    mtext('Lag (yr)', side=1, outer=T)
    mtext('PACF', side=2, outer=T, line=1)
    
    if(!SAVEPLOTS){
      readline('partial ACF -- return to continue ')
    } else {
      dev.off()
    }
  }
  
  
  if('plotGroups' %in% names(yearEffect) & sum(pacfMat,na.rm=T) != 0){  
      # by species/plot
    
    file <- paste('yearEffectByGroup.pdf',sep='')
    
    graphics.off()
    
    if(SAVEPLOTS)pdf( file=.outFile(outFolder,file) )
    
    
    mfrow <- .getPlotLayout(nspec)
    par(mfrow=mfrow$mfrow, bty='n', mar=c(1,1,1,1), oma=c(3,3,1,3))  
    
    preg <- columnSplit(rownames(pacfMat),'-')
    
    for(k in 1:nspec){
      
      wk <- which(preg[,1] == specNames[k])
      
      ylab <- xlab <- F
      if(k %in% mfrow$left)ylab=T
      if(k %in% mfrow$bottom)xlab=T
      
      plot(NULL, xlim = c(1,nyy), ylim = ylim, xaxt = 'n', yaxt='n',
           xlab = '', ylab = '')
      axis(1, at=c(1:nyy), labels=xlab)
      axis(2,labels=ylab)
      abline(h = 0, lty=2, col='grey', lwd=2)
      leg <- character(0); col <- numeric(0)
      lag <- c(0:(nyr-1))
      
      for(j in wk){
        wj <- which( is.finite(pacfMat[j,]))
        wj <- wj[-1]                   # omit lag zero
        nj <- length(wj)
        if(nj < 2)next
        
        ac <- pacfMat[j,wj]
        
        lines(lag[wj], ac, col=plotCol[preg[j,2]], lwd=2)
        loHi <- cbind( ac - 1.96*pacfSe[j,wj],
                       ac + 1.96*pacfSe[j,wj])
        .shadeInterval(lag[wj],loHi,col=plotCol[preg[j,2]],PLOT=T,add=T, 
                       trans = .4)
        
        up <- which(loHi[,1] > 0)
        up <- up[up > 1]
        points(lag[wj[up]],ac[up],cex=1.3,pch=16,
               col=.getColor(plotCol[preg[j,2]],.5))
        
        up <- up[up < 10]
        up <- paste0( up[up > 2], collapse=', ')
        ll <- rownames(pacfMat)[j]
        #    leg <- c(leg,ll)
        #    col <- c(col,j)
      }
      .plotLabel(specNames[k],'topright') 
    }
    if( nspec == prod(mfrow$mfrow) ){
      cornerLegend('bottomright',plots,bty='n',text.col=plotCol[plots], cex=.9)
    }else{
      plot(NULL, xlim=xlim, ylim=ylim, xlab='', ylab='', axes=F)
      legend('topleft',plots,bty='n',text.col=plotCol[plots], cex=1.2)
    }
    mtext('Lag (yr)',side=1,line=1,outer=T)
    mtext('PACF',side=2,line=1,outer=T)
    
    if(!SAVEPLOTS){
      readline('partial ACF -- return to continue ')
    } else {
      dev.off()
    }
  }
  
  
  ############# eigenvalues AR
  
  if(AR){
    
    graphics.off()
    
    file <- paste('eigenAR.pdf',sep='')
    if(SAVEPLOTS)pdf( file=.outFile(outFolder,file) )
    
    par(bty='n', mai=c(1,1,1,1),mar=c(5,5,1,1))
    
    ename <-rownames(eigenMu)
    ename <- rep(ename,each=ncol(eigenMu))
    text.col <- rep(groupCol[yeGr],each=ncol(eigenMu))
    
    xlab <- expression( paste( plain(Re),' ',lambda ))
    ylab <- expression( paste( plain(Im),' ',lambda ))
    
    xseq <- seq(-1,1,length=100)
    yseq <- sqrt(1 - xseq^2)
    plot(eigenMu,xlim=c(-1.2,1.2),ylim=c(-1.1,1.1), 
         cex=.1,xlab=xlab,ylab=ylab)
    lines(xseq,yseq,lwd=2,col='grey',lt=2)
    lines(xseq,-yseq,lwd=2,col='grey',lt=2)
    lines(c(0,0),c(-1,1),col='grey',lt=2)
    lines(c(-1,1),c(0,0),col='grey',lt=2)
    text(Re(eigenMu),Im(eigenMu),ename, cex=.9, col=text.col)
    
    if(!SAVEPLOTS){
      readline('ACF eigenvalues -- return to continue ')
    } else {
      dev.off()
    }
  }
  
  ############# fecundity and seed prediction
  
  tyears <- years  <- sort(unique(sdata$year))
  tplots <- pplots <- attr(sdata$plot,'levels')
  if(AR)tyears <- sort(unique(tdata$year))
  
  if(PREDICT & length(inputs$predList$years) > 1){
    
    file <- paste('forecast.pdf',sep='')
    if(SAVEPLOTS)pdf( file=.outFile(outFolder,file) )
    
    pyears <- sort(unique(c(fecPred$year,seedPredGrid$year)))
    pplots <- attr(seedPredGrid$plot,'levels')
    tyears <- sort(unique(c(years,pyears,tyears)))
    tplots <- sort(unique(c(plots,pplots)))
    
    # prediction grid closest to traps
    
    pcol <- grep('meanM2',colnames(seedPredGrid))
    scol <- grep('seM2',colnames(seedPredGrid))
    
    kcol <- c("trapID","plot","year","trap","plotYr","plotyr",
              "drow","area","active")
    spred <- numeric(0)
    
    for(j in 1:length(pplots)){
      
      xyt <- xytrap[xytrap$plot == pplots[j],]
      
      sp <- seedPredGrid[seedPredGrid$plot == pplots[j],]
      so <- sdata[sdata$plot == pplots[j],]
      if(nrow(so) == 0)next
      sxy <- xyt[match(so$trapID,xyt$trapID),c('x','y')]
      
      spi <- as.character( columnPaste(sp$trapID,sp$year) )
      soo <- as.character( columnPaste(so$trapID,so$year) )
      spi <- match(soo,spi)
      wf  <- which(is.finite(spi))
      
      spp <- sp[spi[wf],pcol,drop=F]
      countPerM2 <- rowSums(so[wf,seedNames,drop=F])/so$area[wf]/so$active[wf]
      predPerM2  <- rowSums(spp)
      
      sall  <- cbind(so[wf,],spp,countPerM2,predPerM2)
      spred <- rbind(spred,sall)
    }
    
    #error by year
    rmse <- sqrt( (spred$countPerM2 - spred$predPerM2)^2 )
    aerr <- spred$countPerM2 - spred$predPerM2
    
    xlim <- range(spred$year,na.rm=T)
    
    maxMod <- NULL
    if(!is.null(modelYears))maxMod <- max(modelYears)
    
    pplots <- as.character(sort(unique(spred$plot)))
    
    mfrow <- .getPlotLayout(length(pplots))$mfrow
    opt <- list(log=F, xlabel='Year', POINTS=F, 
                ylabel='Residual', col='brown', add=T)
    
    par(mfrow=mfrow,bty='n', mar=c(3,3,1,1), oma=c(2,2,0,2))
    
    xseq <- c(0,2^c(0:15))[-2]
    
    ylabel <- expression( paste('Count (', plain(m)^-2,')') )
    zlabel <- expression( bar(y) %+-% plain(sd) )
    
    for(j in 1:length(pplots)){
      
      wj    <- which(spred$plot == pplots[j])
      obs   <- spred$year[wj]
      yMean <- spred$predPerM2[wj]
      yObs  <- spred$countPerM2[wj]
      
      if( max(yMean, na.rm=T) == 0 | length(yMean) < 3 )next
      
      tj <- by(yMean, obs, quantile, probs=pnorm(c(0,-1,1)), na.rm=T)
      cj <- names(tj)
      tj <- matrix( unlist(tj), ncol=3, byrow=T )
      rownames(tj) <- cj
      tj <- sqrt(tj)     
      ww <- which(is.finite(tj[,1]))
      yj <- as.numeric(rownames(tj))
      
      omu <- by(yObs, obs, quantile, probs=pnorm(c(0,-1,1)), na.rm=T)
      cj <- names(omu)
      oj <- matrix( unlist(omu), ncol=3, byrow=T )
      rownames(oj) <- cj
      oj <- sqrt(oj)     # not for residuals
      
      smax <- max( c(tj,oj,4) )
      tt   <- sqrtSeq(smax)
      at   <- tt$at
      labs <- tt$labs
    
    #  xlim <- range(obs,na.rm=T)
      ylim <- range(at)
      
      plot(NULL,xlim=xlim,ylim=ylim, ylab='', xlab='',yaxt='n')
      axis(2,at=at, labels = labs)
             
      if(!is.null(maxMod)){
        rect(maxMod+.5,ylim[1],xlim[2]+.5,ylim[2],col='wheat', border='wheat')
      }
      
      .shadeInterval(yj,loHi=tj[ww,2:3], col=.getColor('grey', .8))
      abline(h=0,lty=2,lwd=4,col='white')
      
      .shadeInterval(yj,loHi=oj[ww,2:3], col=.getColor('green', .3))
 
      lines(yj, tj[ww,1], col='white', lwd=8)
      lines(yj, tj[ww,1], lwd=3)
      lines(yj, oj[ww,1], col=.getColor('white',.5), lwd=8)
      lines(yj, oj[ww,1], col=.getColor('forestgreen',.7),lwd=3)
      
      points(jitter(obs),sqrt(yObs),pch=16,col=.getColor('forestgreen',.2))
      
      .plotLabel(tplots[j],'topleft')
    }
    mtext('Year',side=1,line=0,outer=T, cex=1.4)
    mtext(ylabel,side=2,line=0,outer=T, cex=1.4)
    mtext(zlabel,side=4,line=0,outer=T, cex=1.4)
    
    if(!SAVEPLOTS){
      readline('observed (green), predicted (black), shaded forecast (if modelYears) -- return to continue ')
    } else {
      dev.off()
    }
  }    
  
  yfun    <- colorRampPalette( c('tan', 'brown','turquoise','steelblue') )
  yearCol <- yfun(nyr)
  names(yearCol) <- tyears
  
  ########## tree correlations over years
  
  graphics.off()
  
  file <- paste('treeCor.pdf',sep='')
  
  if(SAVEPLOTS)pdf( file=.outFile(outFolder,file) )
  
  breaks <- seq(-1.1,1.1,by=.1)
  ylim <- c(0,5)
  nplot <- length(plots)
  
  ppt <- character(0)
  
  for(j in 1:nplot){
    
    wjk <- tdata$dcol[ tdata$plot == plots[j] ]
   njk <- length(unique(wjk))
      if(njk < 2)next
   ojk <- omegaE[wjk,wjk]
   ojk[is.na(ojk)] <- 0
   if(max(ojk) == 0)next
   
   ppt <- c(ppt,plots[j])
  }
  
  npp <- length(ppt)
  
  mfrow <- .getPlotLayout(npp)
  par(mfrow=mfrow$mfrow, mar=c(1,1,1,2), oma=c(3,3,1,1), bty='n')
  
  for(j in 1:npp){
    
    jk <- 0
    sk <- character(0)
    ek <- numeric(0)
    
    for(k in 1:nspec){
      
      wjk <- tdata$dcol[ tdata$species == specNames[k] &
                           tdata$plot == ppt[j] ]
      njk <- length(unique(wjk))
      if(njk < 2)next

      wjk <- sort(unique(wjk))
      ojk <- omegaE[wjk,wjk]
      
      ojk[is.na(ojk)] <- 0
      oj <- ojk
      diag(oj) <- 0
      rs <- which( rowSums(oj) == 0 )
      diag(oj) <- diag(ojk)
      if(length(rs) > 0)oj <- oj[-rs,-rs]
      
      oj[oj > .95] <- .95
      oj[oj < -.95] <- -.95
      
      diag(oj) <- 1
      
      if(length(oj) < 2)next
      jk <- jk + 1
      sk <- c(sk,specNames[k])
      
      oj <- oj[lower.tri(ojk)]
      ovec <- hist(oj, breaks = breaks, plot=F)$density
      
      tmp <- .getPoly(breaks[-1],ovec)
      if(jk == 1){
        plot(tmp[1,], tmp[2,],type='s',lwd=2,
             col=.getColor(specCol[specNames[k]],.3),
             xlab='', ylab='', ylim=ylim, xaxt='n', yaxt='n')
        axis(1, at=c(-1,0,1), labels=c(-1,0,1))
        axis(2, labels=T)
      }
      
      tmp <- .getPoly(breaks[-1],ovec)
      polygon( tmp[1,], tmp[2,], col=.getColor(specCol[specNames[k]],.3), lwd=2, 
               border=specCol[specNames[k]])
    }
    if(length(sk) == 0)next
    .plotLabel(ppt[j],'topleft')
    legend('topright',sk,text.col=specCol[sk],bty='n')
  }
  
  mtext('Correlation', side=1, outer=T, line=1)
  mtext('Density', side=2, outer=T, line=1)
  
  if(!SAVEPLOTS){
    readline('tree correlation in time -- return to continue ')
  } else {
    dev.off()
  }
  
  ################# spatio-temporal correlation
  
  if(SPACETIME){
    
    # trees/sites ordered by similarity at zero lag
    
    mvs <- suppressWarnings(
      meanVarianceScore(output, Q = pnorm(c(0, -1, 1)), nsim=1, LAGMAT = T,
                        ktree = 30, maxSite = 30, CLOSE = F)
    )
    
    treeCov <- mvs$lagCanopy
    trapCov <- mvs$lagGround
    nkk <- length(treeCov)
    plotk <- names(treeCov)
    
    graphics.off()
    
    if(nkk > 0){
      
      col2 <- colorRampPalette(c("#67001F", "#B2182B", "#D6604D", "#F4A582",
                                 "#FDDBC7", "#FFFFFF", "#D1E5F0", "#92C5DE",
                                 "#4393C3", "#2166AC", "#053061"))
      for(k in 1:nkk){
        
        wpt <- which(names(treeCov) == plotk[k])
        wpc <- which(names(trapCov) == plotk[k])
        
        if(length(wpt) == 0 | length(wpc) == 0)next
        
        tvar <- treeCov[[wpt]]
        cvar <- trapCov[[wpc]]
        
        if(length(tvar) < 2 | length(cvar) < 2) next
        
        tmp  <- columnSplit(colnames(tvar),'_')
        
        klag <- as.numeric(tmp[,ncol(tmp)])
        tvar[tvar < -1] <- 0
        tvar[tvar > 1] <- 0
        
        if(nrow(tvar) < 2)next
        
        km <- max(klag)
        if(km > 5)km <- 5
        
        
        file <- paste('spaceTime',plotk[k],'.pdf',sep='')
        if(SAVEPLOTS)pdf( file=.outFile(outFolder,file) )
        
        
        par(mfrow=c(2,km+1), mar=c(2,1,1,1), oma=c(2,1,2,1), bty='n',xpd=T)
        
        order  <- 'hclust'
        
        for(i in 0:km){
          
          wl <- which(klag == i)
          stree <- tvar[,wl]
          
          if(i > 0){
            order <- 'original'
            colnames(stree) <- rownames(stree)
            stree <- stree[rnames,rnames]
          }
          
          #    print(dim(stree))
          
          tmp <- corrplot(stree, is.corr=T, method='color', col=rev(col2(200)),
                          tl.pos='n', cl.length=3, cl.lim=c(-1,1), type='lower',
                          order=order, cl.pos='n')
          rnames <- rownames(tmp)
          tlab <- paste('lag',i)
          title(main = list(tlab, cex = 1.5,
                            font = 1), line = -2, adj=1)
          if(i == 0)title(main = list("Canopy", cex = 1.5,
                                      font = 3))
        }
        
        
        tmp  <- columnSplit(colnames(cvar),'_')
        klag <- as.numeric(tmp[,ncol(tmp)])
        cvar[cvar < -1] <- 0
        cvar[cvar > 1] <- 0
        
        order  <- 'hclust'
        
        for(i in 0:km){
          wl <- which(klag == i)
          if(length(wl) == 0)next
          stree <- cvar[,wl]
          
          if(ncol(stree) < nrow(stree)){
            sc <- columnSplit(colnames(stree),'_')[,1]
            w0 <- which(!rownames(stree) %in% sc)
            c1 <- c(1:(w0-1))
            c2 <- c((w0+1):ncol(stree))
            c1 <- c1[c1 > 0]
            c2 <- c2[c2 < nrow(stree) & c2 > w0]
            stree <- cbind( stree[,c1], 0, stree[,c2] )
            colnames(stree)[w0] <- paste(rownames(stree)[w0],i,sep='_')
          }
          
          if(i > 0){
            order <- 'original'
            colnames(stree) <- rownames(stree)
            stree <- stree[rnames,rnames]
          }
          
          tmp <- corrplot(stree, is.corr=T, method='color', col=rev(col2(200)),
                          tl.pos='n', cl.length=3, cl.lim=c(-1,1), type='lower',
                          cl.pos='n')
          rnames <- rownames(tmp)
          if(i == 0)title(main = list("Forest floor", cex = 1.5,
                                      font = 3))
        }
        mtext(plotk[k], 1, outer=T, line=0)
        
        
        if(!SAVEPLOTS){
          readline('tree-time (above), space-time (below) -- return to continue ')
        } else {
          dev.off()
        }
      }
    }
    
    ############## score by scale
    
    darea <- 100
    
    mvs <- suppressWarnings(
      meanVarianceScore(output, Q = pnorm(c(0, -1, 1)), ktree = 20, 
                        nsim=100, LAGMAT = T, CLOSE = T)
    )
    
    graphics.off()
    
    
    scoreT <- scoreS <- scoreTse <- scoreSse <- numeric(0)
    pname  <- character(0)
    
    for(k in 1:length(plots)){
      
      wt <- which(names(mvs$scoreTree) == plots[k])
      ws <- which(names(mvs$scoreSeed) == plots[k])
      
      if(length(wt) == 0 | length(ws) == 0)next
      
      dtree <- mvs$scoreTree[[wt]]
      dseed <- mvs$scoreSeed[[ws]]
      dtreeSE <- mvs$scoreTreeSe[[wt]]
      dseedSE <- mvs$scoreSeedSe[[ws]]
      
      
      if(length(dtree) > 2 & length(dseed) > 2){
        scoreT <- append(scoreT, list(dtree))
        scoreS <- append(scoreS, list(dseed))
        scoreTse <- append(scoreTse, list(dtree))
        scoreSse <- append(scoreSse, list(dseed))
        
        pname  <- c(pname,plots[k])
      }
    }
      
    names(scoreT) <- names(scoreS) <- names(scoreTse) <- names(scoreSse) <- pname
                                      
                                      
    
    pk <- names(scoreT)[k]
    dss <- scoreT
    ylab  <- 'Number of hosts'
    title <- 'Canopy'
    file  <- 'resourceScoreCanopy.pdf'
    carea <- 1
    q <- seq(0, 1, length=15)^1
    cols <- .getColor('black',q)
    
    npp <- length(scoreT)
    
    for(j in 1:2){
      
      if(j == 2){
        dss <- scoreS
        file <- 'resourceScoreGround.pdf'
        yy <- as.numeric(rownames(dk))
        ylab <- 'Distance (m)'
        title <- 'Forest floor'
      }
      
      zscale <- range(sapply(dss, range, na.rm=T))
      cseq   <- seq(zscale[1], zscale[2], length=30)
      
      xscale <- max(sapply(dss, ncol))
      yscale <- max(sapply(dss, nrow))
      if(yscale < 100)yscale <- 100
      if(j == 2 & yscale < 20)yscale <- 20
      
      xlim <- log(1 + c(0, xscale))
      ylim <- log(1 + c(0, yscale+10))
      
      if(SAVEPLOTS)pdf( file=.outFile(outFolder,file) )
      
      mff <- .getPlotLayout(length(dss))
      
      par(mfrow=mff$mfrow, bty='l',mar=c(4,3,.1,.1), oma=c(4,4,1,1))
      
      for(k in 1:npp){
        
        dk <- dss[[k]]
        xx <- columnSplit(colnames(dk),'_')[,2]
        xx <- as.numeric(xx)
        yy <- c(1:nrow(dk))
        
        if(j == 2){
          aa    <- as.numeric(rownames(dk))
      #    carea <- median(diff(aa))^2
          
          yy <- (pi*aa^2)/10000
        }
      
        levels <- quantile(dk, q )
        levels <- sort(unique(levels))
        
        ytick <- c(1, 5, c(1:5)*10 )
        lx <- log(xx + 1)
        ly <- log(yy + 1)
        ltick <- log(ytick + 1) 
        
        if(j == 1)ylim[1] <- ly[1]/2
        if(j == 2)ylim[1] <- diff(ly[1:2])/2
        
        plot(NA, axes = F, xlim=xlim, ylim=ylim, xlab='', ylab='')
        contour(lx, ly, t(dk), levels=levels, col=cols, labels='', add=T,
                axes=F)
        .filled.contour(lx, ly, t(dk), levels=levels, col=cols)
        
        if(length(yy) > 10){
     #     yy <- c(0,yy)
     #     ly <- c(0,ly)
          bb <- ceiling(length(yy)/10)
          ss <- seq(1, length(yy), by=bb)
          ytick <- ytick[ss]
          ltick <- ltick[ss]
        }
        
        xlabs <- ylabs <- F
        
        if(k %in% mff$bottom)xlabs <- xx
        
        if(j == 2){
          ytick <- c(100, 1000, 10000, 50000, 100000)/10000
          ltick <- log(ytick + 1)
   #       ly[1] <- diff(ly[1:2])/2
        }
        
        if(k %in% mff$left){
          ylabs <- ytick
          if(j == 2)ylabs <- c('100 m2', '1000 m2', '1 ha', '5 ha','10 ha')
        }
        
        axis(1, at=lx, labels=xlabs)
        axis(2, at=ltick, labels=ylabs)
        
        .plotLabel(names(dss)[k],'topright')
      }
      mtext('Years', 1, outer=T, line=1)
      mtext(ylab, 2, outer=T, line=1)
      
      endLabs <- signif(range(dk),1)
      
      clist <- list( kplot=1, ytick=NULL, text.col = 'black',
                     cols=cols, labside='right', text.col=col,
                     bg='grey', endLabels=endLabs) 
      cornerScale( clist )
      
      if(!SAVEPLOTS){
        readline('score by scale -- return to continue ')
      } else {
        dev.off()
      }
    }
    
    # plot comparison
    
    nhost <- 5
    nhy   <- 2
    tmat  <- matrix(NA, npp, 3)
    colnames(tmat) <- c('mu','lo','hi')
    rownames(tmat) <- names(scoreT)
    smat <- tmat
    
    for(k in 1:npp){
      
      nn  <- nhost
      
      dkk <- scoreT[[k]]
      if(nn > nrow(dkk))nn <- nrow(dkk)
      skk <- scoreTse[[k]]
      mu  <- dkk[nn,nhy-1]
      ss  <- skk[nn,nhy-1]
      tmat[k,] <- c(mu, mu + ss*c(-1,1))
      
      nn  <- nhost
      dkk <- scoreS[[k]]
      if(nn > nrow(dkk))nn <- nrow(dkk)
      skk <- scoreSse[[k]]
      mu  <- dkk[nn,nhy-1]
      ss  <- skk[nn,nhy-1]
      smat[k,] <- c(mu, mu + ss*c(-1,1))
    }
    
    file <- 'scoreByPlot.pdf'
    if(SAVEPLOTS)pdf( file=.outFile(outFolder,file) )
    
    par(mfrow=c(1,1), mar=c(4,4,2,2), bty='n',xpd=F)
    xlim <- range(tmat) + c(-.1,.2)
    ylim <- range(smat) + c(-.1,.2)
    
    ylab <- paste('Forest floor, ', nhost*darea,'m2')
    xlab <- paste('Canopy score,',nhost,'host trees')
 
    plot(NA, xlim=xlim, ylim=ylim, xlab=xlab,
         ylab = ylab)
    points(tmat[,1],smat[,1])
    segments(tmat[,1], smat[,2], tmat[,1], smat[,3])
    segments(tmat[,2], smat[,1], tmat[,3], smat[,1])
    text(tmat[,1]+.1, smat[,1]+.1,names(mvs$scoreTree), pos=4)
    abline(0,1,col='grey',lwd=2,lty=2)
    .plotLabel(paste(nhy,' yr'), 'bottomright')
    
    if(!SAVEPLOTS){
      readline('score by plot -- return to continue ')
    } else {
      dev.off()
    }
    
    ################ entropy
    
    entropy <- mvs$entropy
    
    if(!is.null(entropy)){
      if(nrow(entropy) > 4){
        
        graphics.off()
        
        file <- 'entropy.pdf'
        if(SAVEPLOTS)pdf( file=.outFile(outFolder,file) )
        
        par(mfrow=c(1,2), mar=c(4,4,1,.5), oma=c(1,1,1,1), bty='n',xpd=F)
        
        entropy[!is.finite(entropy)] <- NA
        
        xl <- range(entropy[,1], na.rm=T )
        dx <- .3*diff(xl)
        xl[2] <- xl[2] + dx
        yl <- range(entropy[,1])
        
        we <- grep('tree-tree',rownames(entropy))
        wr <- grep('site-site',rownames(entropy))
        rnames <- unlist( strsplit(rownames(entropy)[we],'_tree-tree') )
        
        xl <- range(entropy[we,1], na.rm=T ) + c(-1,1)
        dx <- .3*diff(xl)
        xl[2] <- xl[2] + dx
        yl <- range(entropy[wr,1], na.rm=T) + c(-1,1)
        
        
        plot(entropy[we,1],entropy[wr,1], xlim=xl, ylim=yl, xlab='', 
             ylab='Forest floor entropy', cex=.01, yaxt='n')
        axis(2, line=1)
        abline(0,1,lty=2)
        
        par(new=F,xpd=T)
        text(entropy[we,1],entropy[wr,1],rnames)
        mtext('Canopy entropy',1,outer=T,line=-1)
        .plotLabel('a) Spatial', above=T)
        
        we <- grep('tree-lag',rownames(entropy))
        wr <- grep('site-lag',rownames(entropy))
        rnames <- unlist( strsplit(rownames(entropy)[we],'_tree-lag') )
        
        xl <- range(entropy[we,1], na.rm=T ) + c(-1,1)
        dx <- .3*diff(xl)
        xl[2] <- xl[2] + dx
        yl <- range(entropy[wr,1], na.rm=T) + c(-1,1)
        
        plot(entropy[we,1],entropy[wr,1], xlim=xl, ylim=yl, xlab='', ylab='', 
             cex=.01, yaxt='n')
        axis(2, line=1)
        # abline(0,1,lty=2)
        
        par(new=F,xpd=T)
        text(entropy[we,1],entropy[wr,1],rnames)
        .plotLabel('b) Temporal', above=T)
        
        par(new=T,xpd=F)
        
        if(!SAVEPLOTS){
          readline('mean-entropy -- return to continue ')
        } else {
          dev.off()
        }
      }
    }
  }
}

corrLabs <- function(xvar, nsite=5){
  
  tname <- columnSplit(rownames(xvar),'-')[,2]
  tall  <- unique(tname)
  if(length(tall) > nsite){
    ttab <- table(tname)
    most <- names(ttab)[order(ttab,decreasing=T)]
    tall <- most[1:nsite]
  }
  ws   <- which(tname %in% tall)
  xvar <- xvar[ws,ws]
  
  ttmp  <- columnSplit(rownames(xvar),'-')
  tline <- which( diff(rev(as.numeric(ttmp[,3]))) != -1 )
  tname <- unique(ttmp[,2])
  
  vx <- rbind(ncol(xvar) - tline+.5, ncol(xvar) - tline+.5)
  vy <- rbind(.5, tline-.5)
  hx <- rbind(.5, ncol(xvar)-tline+.5)
  hy <- rbind(tline-.5, tline-.5)
  
  tx <- c(ncol(xvar),hx[2,])
  dx <- -diff(tx)/2
  dx <- c(dx,dx[length(dx)])
  tx <- tx - dx*.7
  ty <- c(vy[2,],ncol(xvar)) - dx*.4
  
  list(vx = vx, vy = vy, hx = hx, hy = hy, tx = tx, ty = ty,
       xvar = xvar, tname = tname, tline = tline)
}


t2tr <- function(TV, R){
  
  # TV = nT x nT tree-time covariance
  # R  = nT X R seed type matrix
  
  rvec <- matrix(t(R),ncol=1)
  
  RR <- crossprod(t(rvec))    # R x R blocks in RnT x RnT matrix
  
  ti <- rep(1:ncol(TV), each=ncol(R))
  TT <- TV[ti,]
  TT <- TT[,ti]
  TR <- TT*RR
  
}

getCvol <- function(Fmat, Smat, tree, seed, rCons){
  
  # rCons - species by seed preference or seed mass
  # nrep = 1, because x,y currently fixed
  # on a per-tree and per-location basis
  
  yr <- sort(unique(tree$year))
  
  Rvec <- rCons[tree$species,1]
  Fmat <- t( t(Fmat)*Rvec )
  
  # total
  
  Tv   <- cov(Fmat)
  w0    <- which(diag(Tv) == 0)
  if(length(w0) > 0){
    Fmat <- Fmat[,-w0]
    Tv <- cov(jitter(Fmat))
    tree <- tree[-w0,]
  }
  FF     <- mean(Fmat)
  entT   <- nrow(Tv)/2*(1 + log(2*pi)) + determinant(Tv)$modulus/2
  
  C   <- cov(Smat)
  w0  <- which(diag(C) == 0)
  if(length(w0) > 0){
    Smat <- Smat[,-w0]
    C <- cov(jitter(Smat))
    seed <- seed[-w0,]
  }
  M    <- mean(Smat)
  entC <- nrow(C)/2*(1 + log(2*pi)) + determinant(C)$modulus/2
  rownames(C) <- colnames(C) <- columnPaste(seed$trapID,seed$year)
  
  entT <- entT/nrow(Tv)
  entC <- entC/nrow(C)
  
  # per year
  
  ec <- et <- yri <- numeric(0)
  
  for(k in 1:length(yr)){
    
    wt <- which(tree$year == yr[k])
    ws <- which(seed$year == yr[k])
    
    if(length(wt) < 3 | length(ws) < 3)next
    
    Tk    <- cov(Fmat[,wt])
    w0    <- which(diag(Tk) == 0)
    if(length(w0) > 0)Tk <- Tk[-w0,-w0]
    entTk <- nrow(Tk)/2*(1 + log(2*pi)) + determinant(Tk)$modulus/2
    
    Ck    <- cov(Smat[,ws])
    w0    <- which(diag(Ck) == 0)
    if(length(w0) > 0)Ck <- Ck[-w0,-w0]
    entCk <- nrow(Ck)/2*(1 + log(2*pi)) + determinant(Ck)$modulus/2
    
    entCk <- entCk/nrow(Ck)    #per dimension
    entTk <- entTk/nrow(Tk)
    
    ec  <- c(ec,entCk)
    et  <- c(et,entTk)
    yri <- c(yri,yr[k])
  }
  
  eyr <- rbind(et,ec)
  colnames(eyr) <- yri

  ent <- cbind( rbind(sum(FF),sum(M)),rbind(entT,entC) , eyr)
  colnames(ent)[1:2] <- c('mean_seeds','E_per_dim')
  rownames(ent) <- c('T','C')
  
  list(Tvar = Tv, Cvar = C, ent = ent)
}

.getPoly <- function(x,y){
  dx <- diff(x)
  xx <- c(x[1] - dx[1]/2, x[-1] - dx/2)
  xx <- rep(xx,each=2)
  yy <- rep(y,each=2)
  yy <- c(0,yy,0)
  xx <- c(xx,xx[length(xx)],xx[length(xx)])
  rbind(xx,yy)
}

.chainPlot <- function(mat, burnin, label, refVals = NULL,
                       SAVEPLOTS=F, outFolder='', ylim = NULL){
  
  cseq <- 1:nrow(mat)
  if(length(cseq) > 2000)cseq <- round(seq(1,length(cseq),length=1000))
  
  if(SAVEPLOTS){
    fileName <- .replaceString(label,',','')
    fileName <- .replaceString(label,' ','')
    fileName <- paste(fileName,'.pdf',sep='')
    pdf( file=.outFile(outFolder,fileName) )
  }
  
  cnames <- .coeffNames( colnames(mat) )
  colnames(mat) <- cnames
  
  mfrow <- .getPlotLayout(length(cnames))
  par(mfrow=mfrow$mfrow, bty='n', mar=c(2,2,2,2), oma=c(2,3,1,1)) 
  
  cex <- 1/(1 + mfrow$mfrow[2])^.1
  
  ng <- nrow(mat)
  
  cseq <- 1:ng
  burnline <- burnin
  ss   <- burnin:ng
  if(nrow(mat) > 1000){
    cseq <- seq(1,ng,length=1000)
    burnline <- burnin/ng*1000
    ss <- cseq[ cseq > burnin ]
  }
  
  NEWY <- F
  if(is.null(ylim))NEWY <- T
  
  for(j in 1:ncol(mat)){
    
    xlabels <- F
    
    if(j %in% mfrow$bottom)xlabels <- T
    
    mj <- mat[,j]
    
    if(NEWY){
      ylim <- range(mj)
      if(!is.null(refVals) & NEWY)ylim <- range(c(refVals[j], mj))
    }
    
    plot(mj[cseq], type='l', ylim = ylim, xaxt='n', xlab='', ylab='')
    if(xlabels){
      axis(1, at = c(0, burnline, 1000), labels = c(0, burnin, ng))
    }else{
      axis(1, at = c(0, burnline, 1000), labels = F)
    }
    
    q <- quantile( mj[ss], c(.025, .5, .975) )
    for(k in 1:3){
      segments(burnline, q[k], 1000, q[k], col='white',lwd=1)
      segments(burnline, q[k], 1000, q[k], lty=2)
    }
    segments(burnline, q[1], burnline, q[3], col='white',lwd=1)
    segments(burnline, q[1], burnline, q[3], lty=2)
    if(!is.null(refVals))abline(h = refVals[j], col='blue')
    .plotLabel(cnames[j], above=T, cex=cex)
  }
  mtext('Iteration', outer=T, side=1, line=1)
  mtext('Parameter value', outer=T, side=2, line=1)
    
  if(!SAVEPLOTS){
    lab <- paste(label, ' -- return to continue')
    readline(lab)
  } else {
    dev.off( )
  }
}
 
.outFile <- function(outFolder=character(0),file){
  paste(outFolder,file,sep='/')
}

.plotLabel <- function(label,location='topleft',cex=1.3,font=1,
                       above=F,below=F,bg=NULL){
  
  if(above){
    adj <- 0
    if(location == 'topright')adj=1
    title(label,adj=adj, font.main = font, font.lab =font, cex.main=cex)
    return()
  }
  if(below){
    adj <- 0
    if(location == 'bottomright')adj=1
    mtext(label,side=1,adj=adj, outer=F,font.main = font, font.lab =font,cex=cex)
    return()
  }
  
  if(is.null(bg)){
    tmp <- legend(location,legend=' ',bty='n')
  } else {
    tmp <- legend(location,legend=label,bg=bg,border=bg,text.col=bg,bty='o')
  }
  
  xt <- tmp$rect$left # + tmp$rect$w
  yt <- tmp$text$y
  
  pos <- 4
  tmp <- grep('right',location)
  if(length(tmp) > 0)pos <- 2
  
  XX <- par()$xlog
  YY <- par()$ylog
  
  if(XX)xt <- 10^xt
  if(YY)yt <- 10^yt
  
  text(xt,yt,label,cex=cex,font=font,pos=pos)
}

.boxCoeffs <- function(chain, snames, xlab = "", ylab='Coefficient',
                       addSpec = 'species', ylim=NULL, cols = NULL,
                       xaxt = 's', yaxt = 's'){
  
  nspec  <- length(snames)
  cnames <- colnames(chain)
  xn     <- character(0)
  vnames <- numeric(0)
  iname  <- character(0)
  gnames <- paste(addSpec,snames,sep='')

  for(j in 1:nspec){
    ij     <- which(cnames == gnames[j])
    if(length(ij) > 0)iname <- 'intercept'
    wj     <- grep(gnames[j],cnames)
    if(length(wj) > 0)vnames <- rbind(vnames, wj)
    
    wk <- grep(':',cnames[wj])
    if(length(wk) > 0){
      xn <- matrix( unlist( strsplit(cnames[wj[wk]],':')),
                    ncol=2, byrow=T)[,2]
    }
  }
    
  rownames(vnames) <- snames[vnames[,1]]
  colnames(vnames) <- c(iname,xn)
  nv <- ncol(vnames)
  
  nss <- nrow(vnames)
  
  atvals <- c(1:nss)/(nss + 1)
  atvals <- atvals - mean(atvals)
  sseq   <- c(1:nv)
  xlim   <- c(.5, nv +.5)
  
  if(is.null(ylim)){
    ylim <- range(chain)
    ylim[1] <- ylim[1] - diff(ylim)*.25
  }
  
  add <- F
  if(is.null(cols))cols <- seq(1:nss)
  
  xlabel <- ''
  
  stats <- numeric(0)
  
  for(j in 1:nv){
    
    jcol <- vnames[,j]
    if(j > 1)add <- T
    chainj <- chain[,jcol, drop=F]
    mj     <- mean(chainj)
    sj     <- sd(chainj)
    chainj <- chainj #/sj
    
    if(j == nv)xlabel <- xlab
    
    colj <- cols[j]
    
    if(any(colnames(chainj) %in% snames))colj <- cols[colnames(chainj)]
    
    wi <- grep(':',colnames(chainj))
    if(length(wi) > 0){
      tt <- columnSplit(colnames(chainj)[wi],':')[,1]
      colj <- cols[tt]
    }
    
    boxPars <- .boxplotQuant( chainj, xaxt=xaxt, add = add,
                          at = atvals + j, xlim = xlim,
                          outline=F, ylim=ylim,
                          col=.getColor(colj, .5), 
                          border=colj, lty=1,
                          ylab = ylab, yaxt = yaxt)
    stats <- cbind(stats,boxPars$stats)
 #   text(j,ylim[1],colnames(vnames)[j],pos=3)
  }
  .plotLabel(xlab,'topleft',above=T)
 # legend('topright',snames, text.col=1:nspec, bty='n')
  abline(h = c(0), lwd=1, col=.getColor('grey',.6))
  boxPars$stats <- stats
  invisible(boxPars)
}

.coeffNames <- function(cvec){
  
  # clean names  for coefficients
  
  fnames  <- .replaceString(cvec, '(Intercept)', 'intercept' )
  fnames  <- .replaceString(fnames, 'I(', '' )
  fnames  <- .replaceString(fnames, '))', ')' )
  fnames  <- .replaceString(fnames, 'species','')
  fnames
}


.fixNamesVector <- function(vnames, data, MODE='keep'){
  
  wgg <- match(vnames,names(data))
  
  for(k in wgg){
    data[[k]] <- .fixNames(data[[k]], all=T, MODE)$fixed
  }
  data
}

.fixNames <- function(cvec, all=F, MODE='keep'){
  
  # MODE == 'character', 'factor', or 'keep' (return same mode)
  
  cdup <- numeric(0)
  
  FACT <- F
  if(is.factor(cvec)){
    FACT <- T
    cvec <- as.character(cvec)
  }
  if(is.numeric(cvec)){
    cvec <- as.character(cvec)
    wdot <- grep('.',cvec)
    if(length(wdot) > 0)cvec <- .replaceString(cvec,'.','dot')
  }
  if(all) cvec <- .replaceString(cvec, '_','')
  cvec <- .replaceString(cvec, '-','')
  cvec <- .replaceString(cvec, ' ','')
  cvec <- .replaceString(cvec, "'","")
  cvec <- .replaceString(cvec, ".","dot")
  if(FACT | MODE == 'factor'){
    cvec <- as.factor(cvec)
    droplevels(cvec)
  }
  
  cvec <- .replaceString(cvec, 'acerPenn', 'acerPens')
  
  wd <- which(duplicated(cvec))
  if(length(wd) > 0)cdup <- wd
    
  list(fixed = cvec, dups = cdup)
}

.setupR <- function(sdata, tdata, seedNames, specNames, minDiam, unknown='UNKN'){
  
  SAMPR <- T
  
  if(!'specPlot' %in% colnames(tdata))
    tdata$specPlot <- columnPaste(tdata$species,tdata$plot)
  
  if( length(seedNames) == 1 )SAMPR <- F
  
  plots  <- sort(unique(as.character(tdata$plot)))
  priorR <- numeric(0)
  
  for(j in 1:length(plots)){
    
    wj   <- which(tdata$plot == plots[j])
    if(length(wj) == 0)wj <- which(tdata$plot == plots[j])
    jtab <- table(tdata$species[wj])
    jtab <- jtab[jtab > 0]
    
    ws   <- which(sdata$plot == plots[j])
    stab <- colSums(sdata[drop=F,ws,seedNames], na.rm=T)
    sname <- names(stab)
    
    unkn <- grep(unknown, sname)
    if(length(unkn) > 1)stop(paste('more than one seedNames with string',unknown))
    intr <- names(jtab)
    stab <- stab[ c(which(names(stab) %in% names(jtab)), unkn) ] 
    
    
    jmat <- matrix(0, length(jtab), length(seedNames))
    rownames(jmat) <- names(jtab)
    colnames(jmat) <- seedNames
    
    
    if(ncol(jmat) == 1){
      jmat[,1] <- 1
      rownames(jmat) <- paste(names(jtab),plots[j],sep='-')
      priorR <- rbind(priorR,jmat)
      next
    }
    
    knownSeed <- seedNames
    unkn <- grep(unknown,names(stab))
    
    if(length(unkn) == 0){
      for(k in 1:length(seedNames)){
        wk <- grep(seedNames[k],specNames)
        if(length(wk) > 0)jmat[wk,seedNames[k]] <- 1
      }
    }else{
      treeFrac <- jtab/sum(jtab)                           #tree fraction
      sname     <- sname[sname %in% names(treeFrac)]   #known types
      if(length(sname) > 0){
        toKnown <- stab[sname]
      }else{
        toKnown <- rep(0, length(treeFrac))
        names(toKnown) <- names(treeFrac)
      }
   #   if(length(unkn) > 0)toKnown <- toKnown[-unkn]
    #  toKnown   <- toKnown[toKnown > 0]
      toUnkn    <- treeFrac*stab[unkn] #each spec to unkn
  #    tname     <- intersect(names(toKnown),names(toUnkn))
  #    if(length(tname) > 0){
        fraction  <- toKnown/(toKnown + toUnkn)   # fraction to known types
        
        wna <- which(!is.finite(fraction))
        if(length(wna) > 0)fraction[wna] <- 0
        
        if( length(fraction) == 0 )fraction <- 0
        unn       <- 1 - fraction
        
        fraction <- fraction[fraction > 0]
        if(length(fraction) > 0)
          jmat[ cbind(names(fraction),names(fraction)) ] <- fraction
        jmat[,names(stab)[unkn]] <- unn
    #  }else{
    #    jmat[,unkn] <- 1 
    #  }
    }
    
    rownames(jmat) <- paste(names(jtab),plots[j],sep='-')
    priorR <- rbind(priorR,jmat)
  }
  priorRwt <- priorR*10
  
  seedCount <- as.matrix(sdata[,seedNames,drop=F])
  posR <- which(priorR > 0)
  
  tt <- columnSplit(rownames(priorR),'-')
  
  attr(priorR,'species') <- tt[,1]
  attr(priorR,'plot')    <- tt[,2]
  
  return( list(SAMPR = SAMPR, R = priorR, priorR = priorR, priorRwt = priorRwt,
               seedCount = seedCount, posR = posR, tdata = tdata) )
}

setupZ <- function(tdata, ntree, years, minDiam, maxDiam, maxF){
  
  if(is.na(minDiam) | minDiam == 0)minDiam <- 1
  nyr <- length(years)
  
  # initialize repro
  tdata$repr[tdata$repr < .5] <- 0
  tdata$repr[tdata$repr >= .5] <- 1
  iy   <- match(tdata$year, years)
  tdata$repr[tdata$diam < minDiam & is.na(tdata$repr)] <- 0
  tdata$repr[tdata$diam > maxDiam & is.na(tdata$repr)] <- 1
  zknown <- dmat <- matrix(NA, max(tdata$dcol), length(years))
  zknown[ cbind(tdata$dcol, iy) ] <- tdata$repr
  dmat[ cbind(tdata$dcol, iy) ] <- tdata$diam
  
  tmp <- .getRepro(zknown, dmat, minDiam, maxDiam)
  last0first1 <- tmp$last0first1
  zmat  <- tmp$rmat
  matYr <- tmp$matYr
  
  tmp <- .propZ(zmat, last0first1, matYr)
  zmat <- tmp$zmat
  matYr <- tmp$matYr
  
  # attributes for zmat
  tt <- matrix(NA, ntree, nyr)
  tt[ cbind(tdata$dcol, iy) ] <- tdata$treeID
  tt <- rowMeans(tt,na.rm=T)
  tt <- attr(tdata$treeID,'levels')[tt]
  pt <- matrix( unlist( strsplit(tt, '-') ), ncol=2, byrow=T)
  
  attr(zmat,'plot') <- pt[,1]
  attr(zmat,'treeID') <- tt
  
  tyindex <- cbind(tdata$dcol, iy)  #tree-yr index
  z    <- zmat[ tyindex ]
  
  
  if( is.null(tdata$fecMax) ){
    tdata$fecMax <- 1 + (maxF - 1)*z
    tdata$fecMin <- z
  }
  fmin <- z + 1e-4
  fmax <- 1 + z*maxF
  
  fmin[is.finite(tdata$fecMin)] <- ((z+1e-4)*tdata$fecMin)[is.finite(tdata$fecMin)]
  fmax[is.finite(tdata$fecMax)] <- ((z+1e-4)*tdata$fecMax)[is.finite(tdata$fecMax)]
  
  fmin[fmin < 1e-4] <- 1e-4
  fmax[fmax < 1] <- 1
  tdata$fecMax <- fmax
  tdata$fecMin <- fmin
  
  list(z = z, zmat = zmat, matYr = matYr, 
       last0first1 = last0first1, tdata = tdata)
}

getPredGrid <- function(predList, tdata, sdata, xytree, xytrap, group){
  
  cat('\npredList:\n')
  print(predList)
  
  mapMeters  <- predList$mapMeters
  mapPlot    <- predList$plots
  mapYear    <- predList$years
  
  if(is.null(mapMeters)){
    mapMeters <- 5
    cat('\nMissing mapMeters for prediction grid set to 5 m\n')
  }
  
  plotYrComb <- table( tdata$plot, tdata$year )
  plotYrComb <- plotYrComb[drop=F,rownames(plotYrComb) %in% mapPlot,]
  plotYrComb <- plotYrComb[,colnames(plotYrComb) %in% mapYear, drop=F ]
  npred      <- nrow(plotYrComb)
  
  if(sum(plotYrComb) == 0){
    cat('\nPlot-years in predList missing from data\n')
    return( list(seedPred = NULL, distPred = NULL) )
  }
  
  if(length(mapMeters) == 1 & npred > 1)mapMeters <- rep( mapMeters, npred )
  
  print(mapMeters)
  
  seedPred <- numeric(0)
  drowj <- drowTot <- 0
  
  distPred <- grp <- numeric(0)
  treeid   <- trapid <- numeric(0)
  
  gridSize <- rep(0, npred)
  names(gridSize) <- rownames(plotYrComb)
  
  for(j in 1:npred){
    
    wj <- which(plotYrComb[j,] > 0)
    pj <- rownames(plotYrComb)[j]
    wy <- as.numeric(colnames(plotYrComb)[wj])
  
    jx <- range( c(xytree$x[xytree$plot == pj], xytrap$x[xytrap$plot == pj]) ) 
    jy <- range( c(xytree$y[xytree$plot == pj], xytrap$y[xytrap$plot == pj]) )
    jx[1] <- jx[1] - mapMeters[j]/2
    jx[2] <- jx[2] + mapMeters[j]/2
    jy[1] <- jy[1] - mapMeters[j]/2
    jy[2] <- jy[2] + mapMeters[j]/2
    
    
    
    sx <- seq(jx[1],jx[2],by=mapMeters[j])
    sy <- seq(jy[1],jy[2],by=mapMeters[j])
    
    jgrid <- expand.grid(x = sx, y = sy)
    gridSize[j] <- nrow(jgrid)
    
    yrj   <- rep( wy, nrow(jgrid) )
    
    jseq  <- rep(1:nrow(jgrid), each=length(wy) )
    jgrid <- jgrid[ jseq ,]
    
    dgrid <- drowTot + jseq
    drowTot <- max(dgrid)
    
    trapID <- paste(pj,'-g',dgrid,sep='')
    
    dj  <- data.frame(trapID = trapID, year = yrj, jgrid, drow=0, dgrid = dgrid)
    dj$plot  <- pj
    
    #includes trap years from data
    wmatch <- which(sdata$plot %in% pj)
 
    id  <- as.character(unique( sdata$trapID[wmatch] ) )
    
    sdd <- sdata[match(id,sdata$trapID),]
    
    xy  <- xytrap[match(sdd$trapID,xytrap$trapID),c('x','y')]
  #  ii  <- rep( c(1:length(wy)), each=nrow(xy))
    jj  <- rep( c(1:nrow(sdd)), each=length(wy))
    yrj <- rep( wy, nrow(sdd) )
    drAll  <- sort(unique(sdd$drow))
    dgrid  <- drowTot + c(1:length(drAll))
    dgrid  <- dgrid[ match(sdd$drow[jj],drAll) ]
    drowTot <- max(dgrid)
    
    tj <- data.frame(trapID = sdd$trapID[jj], year = yrj, 
                     x = xy[jj,'x'], y = xy[jj,'y'],
                     drow = sdd$drow[jj], dgrid=dgrid, plot = sdd$plot[jj])
    
    dj <- rbind(dj,tj)
    
    seedPred <- rbind(seedPred, dj)
    
    tj      <- which(xytree$plot == pj)
    xy1     <- xytree[tj,]
    treeid  <- c(treeid,xytree$treeID[tj])
    
    grr <- match( xytree$treeID[tj], tdata$treeID )
    grp <- c(grp, group[grr])
    
    wk <- which(!duplicated(dj$dgrid))
    tid <- dj$trapID[wk]
    
    kgrid <- dj[!duplicated(dj$dgrid),c('x','y')]
    
   # tid     <- paste(pj,1:nrow(kgrid),sep='_')   # index for grid
    trapid  <- c(trapid, tid )
    da      <- .distmat(xy1[,'x'],xy1[,'y'],kgrid[,'x'],kgrid[,'y']) 
    colnames(da) <- xy1$treeID
    rownames(da) <- tid
    distPred     <- .blockDiag(distPred,da)
  }
  
  plotYrComb <- cbind(plotYrComb, mapMeters, gridSize)
  seedPred$active  <- seedPred$area <- 1 # note for 1 m2
  attr(distPred, 'group') <- grp
  
  distPred[distPred == 0] <- 100000
  distPred <- round(distPred,1)
  
  rownames(seedPred) <- NULL
  
  cat("\nPrediction grid size by plot-year\n")
  cat("\nif too large, increase init$mapMeters or init$PREDICTSEED  = F:\n")
  print( plotYrComb )
  
  list(seedPred = seedPred, distPred = distPred)
}

cleanFactors <- function(x){
  
  #fix factor levels
  
  scode <- names(x[ which(sapply( x, is.factor )) ])
  if(length(scode) > 0){
    for(j in 1:length(scode)) {
  #    x[,scode[j]] <- .fixNames(x[,scode[j]])$fixed
      x[,scode[j]] <- droplevels(x[,scode[j]])
    }
  }
  x
}
  
firstClean <- function(tdata, sdata, xytree, xytrap,
                       specNames, seedNames, minDiam){
  
  facLevels <- character(0)
  
  if('repr' %in% names(tdata)){
    wf <- which(is.finite(tdata$repr))
    tdata$repMu <- rep(.5, nrow(tdata))
    tdata$repSd <- rep(1, nrow(tdata))
    if(length(wf) > 0){
      tdata$repMu[wf] <- tdata$repr[wf]
      tdata$repSd[wf] <- .1
    }
    rr <- range(tdata$repr,na.rm=T)
    
    if(rr[1] == 1)warning('every tree declared to be reproductive')
    if(rr[2] == 0){
      warning('every tree declared to be immature')
      tdata$repr[tdata$diam < minDiam] <- NA
    }
  }
  
  rdiam <- range(tdata$diam,na.rm=T)
  if(rdiam[1] <= 0)stop('some diameters <= 0')
  if(rdiam[2] > 500)warning('some diameters > 5 m')
  
  if(!'active' %in% names(sdata))sdata$active <- 1
  
  #coerce characters
  if( is.character(sdata$area) ){
    warning('seedData$area coerced to numeric')
    sdata$area <- as.numeric(sdata$area)
  }
  if( is.character(sdata$active) ){
    warning('seedData$active coerced to numeric')
    sdata$active <- as.numeric(sdata$active)
  }
  if(max(sdata$active, na.rm=T) > 1 | min(sdata$active, na.rm=T) < 0){
    warning('seedData$active outside (0, 1)')
  }
  if( is.character(tdata$diam) ){
    warning('tdata$diam coerced to numeric')
    tdata$diam <- as.numeric(tdata$diam)
  }
  
  tdata  <- .fixNamesVector(vnames=c('plot','tree','species'), 
                           data = tdata, MODE='factor')
  xytree <- .fixNamesVector(vnames=c('plot','tree'), 
                           data = xytree, MODE='factor')
  sdata <- .fixNamesVector(vnames=c('plot','trap'), 
                            data = sdata, MODE='factor')
  xytrap <- .fixNamesVector(vnames=c('plot','trap'), 
                            data = xytrap, MODE='factor')
  
  
  #colnames(tdata) <- .fixNames(colnames(tdata), all=T)$fixed
  colnames(sdata) <- .fixNames(colnames(sdata), all=T)$fixed
  tdata$species   <- .fixNames(tdata$species, all=T)$fixed
  specNames <- .fixNames(specNames, all=T)$fixed
  seedNames <- .fixNames(seedNames, all=T)$fixed
  
  wf <- which( sapply(tdata, is.factor) )
  wf <- wf[!names(wf) %in% c( 'plot', 'tree', 'species')]
  if(length(wf) > 0){
    for(k in wf){
      names(tdata)[k] <- .fixNames( names(tdata)[k], all=T )$fixed
      kl <- .fixNames( levels(tdata[[k]]), all=T )$fixed
      facLevels <- append(facLevels, list( kl ) )
      tdata[[k]] <- .fixNames(tdata[[k]], all=T)$fixed
    }
    names(facLevels) <- names(wf)
  }
  wf <- which( sapply(sdata, is.factor) )
  if(length(wf) > 0)for(k in wf)sdata[[k]] <- .fixNames(sdata[[k]], all=T)$fixed
  wf <- which( sapply(xytree, is.factor) )
  if(length(wf) > 0)for(k in wf)xytree[[k]] <- .fixNames(xytree[[k]], all=T)$fixed
  wf <- which( sapply(xytrap, is.factor) )
  if(length(wf) > 0)for(k in wf)xytrap[[k]] <- .fixNames(xytrap[[k]], all=T)$fixed
  
  plots <- sort(unique(as.character(tdata$plot)))
  years <- sort(unique(tdata$year))
  
  xytrap$trapID <- columnPaste(xytrap$plot,xytrap$trap)
  sdata$trapID  <- columnPaste(sdata$plot,sdata$trap)
  xytree$treeID <- columnPaste(xytree$plot,xytree$tree)
  tdata$treeID  <- columnPaste(tdata$plot,tdata$tree)
  
  xytree <- .cleanRows(xytree, 'treeID')
  xytrap <- .cleanRows(xytrap, 'trapID')
  
  ww <- which(!tdata$treeID %in% xytree$treeID)
  if(length(ww) > 0)tdata <- tdata[-ww,]
  
  ww <- which(!sdata$trapID %in% xytrap$trapID)
  if(length(ww) > 0)sdata <- sdata[-ww,]
  
  tdata  <- cleanFactors(tdata)
  sdata  <- cleanFactors(sdata)
  xytree <- cleanFactors(xytree)
  xytrap <- cleanFactors(xytrap)
  
  ty <- with(tdata, table(treeID, year) )
  
  if(max(ty) > 1){
    wm <- which(ty > 1,arr.ind=T)
    cw <- unique(rownames(wm))
    cy <- paste0( cw , collapse=', ')
    print( paste( 'trees with duplicated years:', cy ) )
    cat( '\nduplicated tree-years removed\n' )
    tdata <- tdata[!as.character(tdata$treeID) %in% cw,]
  }
  
  tdata    <- tdata[as.character(tdata$plot) %in% plots & 
                      tdata$year %in% years,]
  sdata <- sdata[as.character(sdata$plot) %in% plots & 
                         sdata$year %in% years,]
  xytree   <- xytree[as.character(xytree$plot) %in% plots,]
  xytrap   <- xytrap[as.character(xytrap$plot) %in% plots,]
  tdata$plot  <- droplevels(tdata$plot)
  sdata$plot  <- droplevels(sdata$plot)
  xytree$plot <- droplevels(xytree$plot)
  xytrap$plot <- droplevels(xytrap$plot)
  
  # fix names
  tdata$species <- as.factor( .fixNames(as.character(tdata$species), all=T)$fixed )
  tdata$species <- droplevels( tdata$species )
  
  specNew <- .fixNames(specNames, all=T)$fixed
  seedNew <- .fixNames(seedNames, all=T)$fixed
  
  ww <- which(!specNames %in% specNew)
  if(length(ww) > 0){
    ts <- as.character(tdata$species)
    for(k in 1:length(ww)){
       ts[ ts == specNames[ww[k]] ] <- specNew[ww[k]]
    }
    tdata$species <- as.factor(ts)
    tdata$species <- droplevels( tdata$species )
    specNames <- sort( unique( specNew ) )
  }
  
  ww <- which(!seedNames %in% seedNew)
  if(length(ww) > 0){
    
    mm <- match(seedNames[ww], colnames(sdata))
    sdata <- sdata[,-mm]
    seedNames <- sort( unique( seedNew ) )
  }
           
  colnames(sdata) <- .fixNames(colnames(sdata), all=T)$fixed
  if('species' %in% colnames(xytree))
    xytree$species <- as.factor( .fixNames(as.character(xytree$species), all=T)$fixed )
  
  
  if(!'diam' %in% colnames(tdata))stop("include 'diam' column in treeData")
  
  tdata <- tdata[tdata$year %in% sdata$year,]
  
  tdata <- tdata[,!names(tdata) %in% c('treeID','plotYr')]
  sdata <- sdata[,!names(sdata) %in% c('trapID','plotYr')]
  
  tdata$year <- as.numeric(tdata$year)
  sdata$year <- as.numeric(sdata$year)
  

  
  list(tdata = tdata, sdata = sdata, xytree = xytree, xytrap = xytrap,
       specNames = specNames, seedNames = seedNames, plots = plots, 
       years = years, facLevels = facLevels)
}

.setupRandom <- function(randomEffect, tdata, xfec, xFecNames, specNames){
  
  nspec      <- length(specNames)
  formulaRan <- randomEffect$formulaRan
  if(nspec > 1)formulaRan <- .specFormula(randomEffect$formulaRan)
  xx        <- .getDesign(formulaRan, tdata)$x
  if(nspec > 1)xx <- xx[,grep('species',colnames(xx)),drop=F]  # CHECK for 1 spp
  xrandCols  <- match(colnames(xx),colnames(xfec))
  Qrand      <- length(xrandCols)
  reI        <- as.character(tdata[,randomEffect$randGroups])
  rnGroups   <- sort(unique(reI))
  reIndex    <- match(reI, rnGroups)
  reGroups   <- sort(unique(reIndex))
  nRand      <- length(reGroups)
  Arand      <- priorVA <- diag(1, Qrand)
  dfA        <- ceiling( Qrand + sqrt(nRand) )
  alphaRand  <- matrix(0, nRand, Qrand)
  colnames(alphaRand) <- xFecNames[xrandCols]
  rownames(alphaRand) <- rnGroups
  
 # formulaRan <- .specFormula(randomEffect$formulaRan, NOINTERCEPT=T)
 # xunstand   <- model.matrix(formulaRan,tdata)
  
  XX <- crossprod(xx)
  diag(XX) <- diag(XX) + .00000001
  xrands2u <- solve(XX)%*%crossprod(xx,xfec[,xrandCols]) 
  xrands2u[abs(xrands2u) < 1e-8] <- 0
  
  list(formulaRan = formulaRan, xrandCols = xrandCols, Qrand = Qrand,
       rnGroups = rnGroups, reIndex = reIndex, reGroups = reGroups, 
       Arand = Arand, dfA = dfA, alphaRand = alphaRand, priorVA = priorVA,
       xrands2u = xrands2u)
}

getPlotDims <- function(xytree, xytrap){
  
  plots <- sort( unique(as.character(xytree$plot) ) )
  npp   <- length(plots)
  
  pdims <- numeric(0)
  
  for(j in 1:npp){
    
    jx <- range( c(xytree$x[xytree$plot == plots[j]],
                   xytrap$x[xytrap$plot == plots[j]]) )
    jy <- range( c(xytree$y[xytree$plot == plots[j]],
                   xytrap$y[xytrap$plot == plots[j]]) )
    jx[1] <- floor( jx[1] - 1 )
    jx[2] <- ceiling( jx[2] + 1 )
    jy[1] <- floor( jy[1] - 1 )
    jy[2] <- ceiling( jy[2] + 1 )
    
    area <- diff(jx)*diff(jy)/10000
    
    pdims <- rbind( pdims, c(jx, jy, area) )
  }
  colnames(pdims) <- c('xmin','xmax','ymin','ymax','area')
  rownames(pdims) <- plots
  pdims
}
  
.orderChain <- function(xchain, snames){
  
  if(!snames[1] %in% colnames(xchain))return(xchain)
  
  ns <- length(snames)
  
  mnames <- .coeffNames(colnames(xchain))
  first  <- mnames[1:ns]
  tmp    <- grep('_',first)
  
  if(length(tmp) > 0){
    first <- matrix( unlist(strsplit(first,'_')),ncol=2,byrow=T)[,2]
  }
  
  orr    <- match(snames,first)
  if(is.na(orr[1]))return(xchain)
  
  newChain <- xchain*0
  
  k <- orr
  m <- 1:ns
  while(max(k) <= ncol(xchain)){
    
    newChain[,m] <- xchain[,k]
    colnames(newChain)[m] <- colnames(xchain)[k]
    
    m <- m + ns
    k <- k + ns
    
  }
  newChain
}

mastif <- function(inputs, formulaFec=NULL, formulaRep=NULL, 
                 predList=NULL, yearEffect = NULL, randomEffect = NULL, 
                 modelYears = NULL, ng = NULL, burnin = NULL){   
  data  <-  NULL
  
  if( class(inputs) == 'mastif' ){
    data   <- inputs$data
    if( !is.null(modelYears) ){
      inputs$inputs$tdataOut <- inputs$prediction$tdataOut
      inputs$inputs$sdataOut <- inputs$prediction$sdataOut
    }
    inputs <- inputs$inputs
    class(inputs) <- 'mastif'
  }
  if(!'specNames' %in% names(inputs))stop('missing specNames')
  if(!'seedNames' %in% names(inputs))stop('missing seedNames')

  if(!'treeData' %in% names(inputs))stop('missing treeData')
  if(!'seedData' %in% names(inputs))stop('missing seedData')
  if(!'xytree' %in% names(inputs))stop('missing xytree')
  if(!'xytrap' %in% names(inputs))stop('missing xytrap')
  if(is.null(ng))stop("supply no. MCMC steps, 'ng'")
  if(is.null(burnin))stop("supply 'burnin'")
  
  .mast(inputs, data, formulaFec, formulaRep, predList, yearEffect, 
        randomEffect, modelYears, ng, burnin) 
}
  
.mast <- function(inputs, data, formulaFec, formulaRep, predList, yearEffect,
                  randomEffect, modelYears, ng, burnin){
  
  upar <- xytree <- xytrap <- specNames <- treeData <- seedData <-
  seedNames <- arList <- times <- xmean <- xfecCols <- xrepCols <-
  groupByInd <- dfA <- xrands2u <- lagGroup <- lagMatrix <- xfecs2u <-
  xreps2u <- Qrand <- xfecU <- xrepU <- seedMass <- plotDims <- plotArea <-
  tdataOut <- sdataOut <- specPlots <- NULL

  priorDist <- 20; priorVDist <- 10; maxDist <- 40; minDist <- 2
  minDiam <- 2; maxDiam <- 80; sigmaMu <- 1
  sigmaWt <- nrow(inputs$treeData)/10
  plag    <- 1; maxF   <- 1e+8; priorVtau <- 6
  ug <- plag <- priorTauWt <- NULL
  alphaRand <- Arand <- priorB <- priorIVB <- betaPrior <- NULL
  RANDOM <- YR <- AR <- ARSETUP <- F
  PREDSEED <- T
  if(is.null(predList))PREDSEED <- F
  betaYr <- betaLag <- yeGr <- plots <- years <- NULL
  facLevels <- character(0)
  ngroup <- 1
  
  nng <- ng
  
  if( class(inputs) == 'mastif' ){
    
    ARSETUP <- T
    
    ww <- c(1:length(inputs))
    if(!is.null(formulaFec))ww <- ww[ names(inputs)  != 'formulaFec' ] 
    if(!is.null(formulaRep))ww <- ww[ names(inputs)  != 'formulaRep' ] 
    ww <- ww[ !names(inputs) %in% c('burnin', 'ng', 'predList') ]
    
    for(k in ww)assign( names(inputs)[k], inputs[[k]] )
    
    for(k in 1:length(data$setupData)){
      assign( names(data$setupData)[k], data$setupData[[k]] )
    }
    if('arList' %in% names(data)){
      for(k in 1:length(data$arList))
        assign( names(data$arList)[k], data$arList[[k]] )
      AR <- T
    }
    if('setupRandom' %in% names(data)){
      for(k in 1:length(data$setupRandom))
        assign( names(data$setupRandom)[k], data$setupRandom[[k]] )
      RANDOM <- T
    }
    if('setupYear' %in% names(data)){
      for(k in 1:length(data$setupYear)){
        if(names(data$setupYear)[k] == 'yrIndex')next
        assign( names(data$setupYear)[k], data$setupYear[[k]] )
      }
      YR <- T
    }
    ug <- upar <- inputs$parameters$upars[,1]
    years <- sort(unique(tdata$year))
    nyr <- years
    if(!is.null(predList)){
      predList$years <- predList$years[predList$years %in% years]
      if('plots' %in% names(predList))
        predList$plots <- .fixNames(predList$plots, all=T)$fixed
    }
    yrIndex <- yrIndex[,!duplicated(colnames(yrIndex))]
    
  }else{ 
    data <- numeric(0)
    for(k in 1:length(inputs))assign( names(inputs)[k], inputs[[k]] )
    seedNames <- sort(seedNames)
    specNames <- sort(specNames)
    years <- sort(unique(c(treeData$year,seedData$year)))
  }
  
  if(length(which(duplicated(c(seedNames)))) > 0)stop('duplicate seedNames')
  if(length(which(duplicated(c(specNames)))) > 0)stop('duplicate specNames')
  
  if(!is.null(randomEffect)){
    if('randGroups' %in% names(randomEffect)){
      randomEffect$randGroups <- 'treeID'
      randGroups <- 'treeID'
    }
  }
  
  if(!is.null(predList)){
    PREDSEED <- T
    if(!'plots' %in% names(predList))stop('predlist must include plots')
    predList$plots <- .fixNames(predList$plots, all=T)$fixed
    if(!'mapGrid' %in% names(predList))predList$mapGrid <- 5
  }
  if(!is.null(yearEffect)){
    plotGroups <- yearEffect$plotGroups
    specGroups <- yearEffect$specGroups
    YR <- T
    if('p' %in% names(yearEffect)){
      YR <- F
      AR <- T
      plag <- yearEffect$p
    }
  }
  
  ng <- nng
  
  treeData <- .fixNamesVector(vnames = c('plot','tree','species'), data=treeData, 
                              MODE='character')
  xytree   <- .fixNamesVector(vnames = c('plot','tree'), data=xytree, 
                              MODE='character')
  seedData <- .fixNamesVector(vnames = c('plot','trap'), data=seedData, 
                              MODE='character')
  xytrap   <- .fixNamesVector(vnames = c('plot','trap'), data=xytrap, 
                              MODE='character')
  
  colnames(seedData) <- .fixNames(colnames(seedData), all=T)$fixed
  
  ws <- which(colnames(seedData) %in% seedNames)
  colnames(seedData)[ws] <- .fixNames(colnames(seedData)[ws], all=T)$fixed
    
  specNames   <- .fixNames(specNames, all=T)$fixed
  treeData$species <- .fixNames(treeData$species, all=T)$fixed
  
  
  
  plots <- sort(unique(as.character(treeData$plot)))
  
  tmp <- checkPlotDims(plots, years, xytree, xytrap, plotDims, plotArea)
  plotDims <- tmp$plotDims
  plotArea <- tmp$plotArea

  
  priorU  <- (2*priorDist/pi)^2
  priorVU <- (2/pi)^2*priorVDist^2
  priorVU <- priorVU            
  maxU    <- (2*maxDist/pi)^2
  minU    <- (2*minDist/pi)^2

  umean <- priorU
  propU <- priorU/100
  uvar  <- priorVU
  if(is.null(ug))ug <- priorU

  if(is.null(priorTauWt))priorTauWt <- nrow(treeData)
  tau1 <- priorTauWt
  tau2 <- priorVU*(tau1 - 1)
  
  
  nspec <- length(specNames)
  if(is.null(yeGr))yeGr <- specNames[1]
  
  if(!is.null(betaPrior)){
    if(!is.null(yearEffect) | !is.null(randomEffect))
      stop('if betaPrior included, do not include yearEffect or randomEffect')
    plag <- yearEffect$p
  }
  
  
  if(nspec > 1){
    fc <- as.formula( .replaceString( as.character(formulaFec), 'species *','') )
    fr <- as.formula( .replaceString( as.character(formulaRep), 'species *','') )
    formulaFec <- .specFormula(fc)
    formulaRep <- .specFormula(fr)
  }
  
  tdata <- treeData
  sdata <- seedData
  
  if( !is.null(modelYears) ){
    
    inputs$modelYears <- modelYears
    
    wtree <- which(tdata$year %in% modelYears)
    wtrap <- which(sdata$year %in% modelYears)
    
    tdataOut <- tdata[-wtree,]
    sdataOut <- sdata[-wtrap,]
    sdataOut$seedM2  <- round(rowSums( as.matrix(sdataOut[,seedNames]) )/
                                sdataOut$active/sdataOut$area,1)
    
    xy <- xytrap[ match(sdataOut$trapID, xytrap$trapID),c('x','y')]
    sdataOut <- cbind(xy, sdataOut)
  }
  
  ##################
  rm(treeData)
  rm(seedData)
  ##################
 
  if( !ARSETUP ){ 
    
    tmp <- firstClean(tdata, sdata, xytree, xytrap, specNames, seedNames, minDiam)
    for(k in 1:length(tmp))assign( names(tmp)[k], tmp[[k]] )
    
    tmp <- .setupData(formulaFec, formulaRep, tdata, sdata, facLevels,
                      xytree, xytrap, specNames, seedNames, AR, YR, yearEffect, 
                      minDiam, maxDiam)
    for(k in 1:length(tmp))assign( names(tmp)[k], tmp[[k]] ) 
    data      <- append(data, list(setupData = tmp))
    
    # enough seed?
    ntype <- length(seedNames)
    traps <- sort(unique(sdata$trapID))
    ii <- match(as.character(sdata$plot),plots)
    ii <- rep(ii,ntype)
    jj <- rep(c(1:ntype),each=nrow(sdata))
    tt <- .myBy(as.vector(unlist(sdata[,seedNames])),ii,jj,fun='sum')
    rownames(tt) <- plots
    colnames(tt) <- seedNames
    cat('\nSeed count by plot:\n')
    print(tt)
    if(sum(tt) < 10)stop('not enough seeds')

    if(AR){
      data      <- append(data, list(arList = arList))
      for(k in 1:length(arList))assign( names(arList)[k], arList[[k]] ) 
      years    <- c(min(years) - plag):(max(years) + plag)
      nyrAR    <- length(times)
      groupYr  <- yrIndex[,'group'] # for AR, group and groupYr identical
      yrIndex  <- cbind(yrIndex,groupYr)
      preYr    <- years[-c(1:plag)]

      # check order:
      # plot(tdata$dcol,yrIndex$group)      
      # points(1:ntree, groupByInd, col=2, cex=.6)
      # points(tdata$dcol[lagMatrix[,1]],lagGroup,col=3,cex=.3) # must be > p yr
      # plot(tdata$dcol[match(tdata$treeID,colnames(distall))], cex=.2)
    }

    if(YR){
      yrTmp <- yrIndex
      tmp <- .setupYear(yearEffect, tdata, years)
      for(k in 1:length(tmp))assign( names(tmp)[k], tmp[[k]] ) 
      data     <- append(data, list(setupYear = tmp))
      
      wn <- which(!colnames(yrTmp) %in% colnames(tmp$yrIndex))
      if(length(wn) > 0)yrIndex <- cbind(yrTmp[,wn],tmp$yrIndex)
      
      betaYrR  <- betaYr
      betaYrF  <- betaYrR[drop=F,1,]
    }
    if(YR | AR){
      dcol     <- tdata$dcol
      yrIndex <- as.matrix( cbind(yrIndex, dcol ) )
    }
  } ################
  
  if(!'dcol' %in% colnames(yrIndex)){
    dcol <- tdata$dcol
    yrIndex <- cbind(yrIndex,dcol)
  }
  if(is.null(upar))upar <- ug
  
  if( is.list(yrIndex) )yrIndex <- as.matrix(yrIndex)
  
  nseed     <- nrow(sdata)
  nplot     <- length(plots)
  n         <- nrow(xfec)
  ntobs     <- table(tdata$plot); nsobs <- table(sdata$plot)
  ttab      <- table(tdata$plot, tdata$year)
  wtab      <- which(ttab > 0, arr.ind=T) 
  ntree     <- nrow(xytree)
  ntrap     <- nrow(xytrap)
  obsYr     <- sort(unique(sdata$year))
  
  if(!is.null(modelYears))      
    obsYr <- intersect(obsYr,modelYears)  #also, plotYears? yrIndex?
  
  
  obsTimes <- match(obsYr,years)
  
  RANDYR <- F
  species <- match( as.character(tdata$species), specNames )
  yrIndex <- cbind(yrIndex,species)

  tdata  <- cleanFactors(tdata)
  plotYears <- attr(tdata$plotYr, 'levels' )
  nyr    <- nacc <- length(years)
  ngroup <- length(yeGr)
  nspec  <- length(specNames)
  
  yrIndex <- yrIndex[,!duplicated(colnames(yrIndex))]
  
  UNSTAND <- T
  if(length(xmean) == 0)UNSTAND <- F
  
  Qfec <- ncol(xfec)
  Qrep <- ncol(xrep)
  xFecNames <- colnames(xfec)
  xRepNames <- colnames(xrep)
  
  nSpecPlot <- max(yrIndex[,'specPlot'])
  
  pacfMat <- matrix(0, nSpecPlot, nacc)
  rownames(pacfMat) <- as.character(specPlots)
  colnames(pacfMat) <- paste('lag',0:(nacc-1),sep='-')

  pacf2 <- pacfMat
  acfMat <- pacfMat
  pacN <- pacfMat
  
  if(YR | AR){
    muyr <- rep(0,nrow(tdata))
    if(ngroup > 1)RANDYR <- T
  }
  
  if( !ARSETUP ){        
    # random effects
    rnGroups <- reIndex <- reGroups <- priorVA <- NULL
    alphaRand <- Arand <- NULL
    if( !is.null(randomEffect) ){
      RANDOM     <- T
      tmp <- .setupRandom(randomEffect, tdata, xfec, xFecNames, specNames) 
      for(k in 1:length(tmp))assign( names(tmp)[k], tmp[[k]] ) 
      data <- append(data, list(setupRandom = tmp))
      if(length(Arand) == 1)ONEA <- T
    }
  }           
 
  obsRows <- which(tdata$obs == 1)
  
  ONEF <- ONER <- ONEA <- F  # intercept only
  if(ncol(xfec) == 1)ONEF <- T  # intercept only
  if(ncol(xrep) == 1)ONER <- T  # intercept only
  
  seedPredGrid <- distPred <- NULL
  nseedPred <- 0
  
  rownames(yrIndex) <- NULL
  
  
  
  # species to seed type
  
  tmp <- .setupR(sdata, tdata, seedNames, specNames, minDiam)
  R         <- tmp$R
  priorR    <- tmp$priorR
  priorRwt  <- tmp$priorRwt
  SAMPR     <- tmp$SAMPR
  seedCount <- tmp$seedCount
  posR      <- tmp$posR
  tdata     <- tmp$tdata
  
  if(PREDSEED){
    
    if(AR)predList$years <- years         # may not be desireable 
    
    if(is.character(predList$years)) predList$years <- 
                            as.numeric(predList$years)
    
    tmp <- getPredGrid(predList, tdata, sdata, xytree, xytrap, 
                       group = yrIndex[,'group'])
    
    seedPredGrid <- tmp$seedPred
    distPred     <- tmp$distPred
    nseedPred    <- nrow(seedPredGrid)
    if( !is.null(nseedPred) ){ 
      treeRows <- which( tdata$plot %in% predList$plots &
                           tdata$year %in% predList$years )
      predDcol <- match( tdata$treeID[treeRows], colnames(distPred) )
      predTree <- data.frame( treeID = tdata$treeID[treeRows], row = treeRows, 
                              dcol = predDcol, 
                              species = tdata$species[treeRows],
                              year = tdata$year[treeRows])
    }else{
      PREDSEED <- F
    }
    tdatPred <- data.frame(tdata[predTree[,'row'],c('species','specPlot',
                                                    'year')])
    tdatPred$dcol <- predTree[,'dcol']
  #  colnames(tdatPred)[3] <- 'dcol'
    sdatPred <- data.frame(year = seedPredGrid[,'year'],
                           drow= seedPredGrid[,'dgrid'])
  }
  
  ti <- match(colnames(distall),as.character(tdata$treeID))
  attr(distall, 'group') <- yrIndex[ti,'group'] 
  
  if(YR){
    cnames  <- paste('yr',1:nyr,sep='-')
    sygibbs <- matrix(0,ng, nyr)
    colnames(sygibbs) <- cnames        # variance, random years
  }
  
  if(AR)Alag   <- diag(.1, plag, plag)
  
  if( (AR | YR) & ngroup > 1){
    if(length(ug) == 1)ug <- rep(ug, ngroup)
    names(ug) <- yeGr
    propU <- rep(propU, ngroup)
  }
  

  tmp <- setupZ(tdata, ntree, years, minDiam, maxDiam, maxF)
  z           <- tmp$z
  zmat        <- tmp$zmat
  matYr       <- tmp$matYr 
  last0first1 <- tmp$last0first1
  tdata       <- tmp$tdata       # fecMin, fecMax included
  fecMin <- tdata$fecMin
  fecMax <- tdata$fecMax
  
  
  MISS <- F
  missSeed <- which(is.na(seedCount),arr.ind=T)
  if(length(missSeed) > 0){
    MISS <- T
    sdata[,seedNames][missSeed] <- 0
  }
  
  
  
  # initial maturation fecundity
  
  if('lastRepr' %in% names(tdata)){
    z <- tdata$lastRepr
    zmat[ yrIndex[,c('dcol','year')] ] <- z
  }
  if('lastFec' %in% names(tdata)){
    fg <- tdata$lastFec
  }else{
    tmp <- .initEM(last0first1, yeGr, distall, ug, priorU, minDiam, 
                   minU, maxU, tdata, sdata, specNames, seedNames, R, 
                   SAMPR, years, plotYears, z, xfec, zobs=tdata$repr)
    fg <- tmp
  }
  fg[fg > (maxF - 1)] <- maxF - 1
  fecMax[fecMax > maxF] <- maxF
  
  propF <- fg/10
  propF[propF < .001] <- .001
  propF[propF > .1*maxF] <- .1*maxF
  
  
  # reg variance
  sg <- sigmaMu
  s1 <- sigmaWt
  s2 <- sigmaMu*(s1 - 1)
  
  
  bgFec <- matrix(0, Qfec, 1)
  bgRep <- matrix(0, Qrep, 1)
  rownames(bgFec) <- xFecNames
  rownames(bgRep) <- xRepNames
  
  if(!is.null(betaPrior)){
    betaPrior <- .getBetaPrior(betaPrior, bgFec, bgRep, specNames)
  }
  
  
  obsRowSeed <- which(sdata$year %in% obsYr)
  
  .updateBeta <- .wrapperBeta(priorB, priorIVB, minU, maxU,
                              priorU, priorVU, SAMPR, obsRows, obsYr, obsRowSeed,
                              sdata[,c('year','drow','active','area', seedNames)], 
                              tdata[,c('specPlot','year','dcol')],
                              seedNames, xfecCols,
                              xrepCols, last0first1, ntree, nyr,
                              betaPrior, years, distall, YR, AR, yrIndex,
                              RANDOM, reIndex, xrandCols, RANDYR, tau1, tau2)
  
  predYr  <- sort( unique(tdata$year) )
  
  tcols <- c('specPlot','species','dcol','year','plotYr','plotyr')
  if(AR)tcols <- c(tcols,'times')

  
  updateProp <- c( 1:1000, seq(1001, 10000, by=100) )
  updateProp <- updateProp[updateProp < .9*ng]
  
  .updateFecundity <- .wrapperStates( maxF, SAMPR, RANDOM, obsTimes, plotYears, 
                                      sdata, tdat = tdata[,tcols], seedNames,
                                      last0first1, distall, YR, AR,
                                      obsRows, obsYr, predYr, obsRowSeed,
                                      ntree, years, nyr, xrandCols, reIndex, 
                                      yrIndex, plag, groupByInd, RANDYR, updateProp)
  
  bfgibbs  <- matrix(0, ng, Qfec); colnames(bfgibbs) <- xFecNames
  brgibbs  <- matrix(0, ng, Qrep); colnames(brgibbs) <- xRepNames
  bygibbsF <- bygibbsR <- NULL

  ugibbs <- matrix(0, ng, ngroup)
  sgibbs <- matrix(0, ng, 3)
  colnames(ugibbs) <- yeGr
  colnames(sgibbs) <- c('sigma','rmspe','deviance')
  
  ncols <- nyr
  if(AR){
    ncols <- plag
    cnames <- paste('lag',c(1:plag),sep='-')
  }
  
  betaYrF  <- betaYrR <- matrix(0, 1, ncols)
  if( YR ) sgYr <- rep(1, ncols)
  if(YR | AR){
    betaYrF  <- matrix(0, 1, ncols)
    betaYrR  <- matrix(0, ngroup, ncols)
    rownames(betaYrR) <- yeGr
    colnames(betaYrF) <- colnames(betaYrR) <- cnames
    bygibbsF <- matrix(NA, ng, length(betaYrF))
    bygibbsR <- matrix(0, ng, length(betaYrR))
    colnames(bygibbsF) <- colnames(betaYrF)
    colnames(bygibbsR) <- .multivarChainNames(yeGr, colnames(betaYrR))
    bygibbsN <- bygibbsR
  }
  if(AR){
    Gmat  <- rbind(0, cbind( diag(plag-1), 0 ) )
    eigenMat <- eigen1 <- eigen2 <- betaYrR*0
  }
  if(ngroup > 1){
    priorUgibbs <- matrix(0,ng,2)
    colnames(priorUgibbs) <- c('mean','var')
  }
  
  if(RANDOM){
    agibbs <- matrix(NA, ng, length(Arand))
    colnames(agibbs) <- .multivarChainNames(xFecNames[xrandCols],
                                            xFecNames[xrandCols])
    asum <- asum2 <- alphaRand*0
  }
  
  colnames(brgibbs) <- xRepNames
  
  if(SAMPR){
    rgibbs <- matrix(0, ng, length(R))
    colnames(rgibbs) <- .multivarChainNames(rownames(R), colnames(R) )
    rgibbs <- rgibbs[,posR]
  }
  
  accept <- rep(0, length(plotYears))
  
  cat('\nMCMC\n')
  
  ug <- .tnorm( length(ug), minU, maxU, ug, 5)
  
  pars  <- list(fg = fg, fecMin = tdata$fecMin, fecMax = tdata$fecMax, 
                ug = ug, umean = umean, uvar = uvar,
                sg = sg, bgFec = bgFec, bgRep = bgRep,
                betaYrR = betaYrR*0, betaYrF = betaYrF, alphaRand = alphaRand, 
                Arand = Arand, R = R)
  
  mufec <- xfec%*%bgFec
  muyr  <- muran <- mufec*0
  
  # draw from probit
  wlo <- rep(-Inf, length(z))
  whi <- rep(Inf, length(z))
  whi[z == 0] <- 0
  wlo[z == 1] <- 0
  w <- .tnorm(length(z), wlo, whi, tdata$repMu, 1)
  
  tmp <- .updateBeta(pars, xfec, xrep, R, propU, propF, w, 
                     z, zmat, matYr, muyr)
  ug    <- pars$ug    <- tmp$ug
  umean <- pars$umean <- tmp$umean
  uvar  <- pars$uvar  <- tmp$uvar
  bgFec <- pars$bgFec <- tmp$bgFec
  bgRep <- pars$bgRep <- tmp$bgRep
  propU <- tmp$propU
  
  # tree correlation
  nSpecPlot <- max(yrIndex[,'specPlot'])
  fmat <- matrix(0,ntree,nyr)
  
  wwi <-  match( unique(tdata$dcol),tdata$dcol )
  rownames(fmat) <- tdata$treeID[ wwi ]
  
  if(is.null(seedMass)){
    seedMass <- matrix(1,length(specNames),1)
    rownames(seedMass) <- specNames
  }
  
  if(MISS){
    lm <- .getLambda(tdata[obsRows,c('specPlot','year','dcol')],
                     sdata[obsRowSeed,c('year','drow')],
                     AA=1, ug, fg[obsRows], z[obsRows], R, 
                     SAMPR, distall, obsYr, PERAREA=F, SPECPRED=F)
    sdata[obsRowSeed,seedNames][missSeed] <- round(lm[missSeed])
  }
  
  
  omegaE <- omegaN <- matrix(0,ntree,ntree)
  
  if(PREDSEED){  #predict species, not seedNames
    specPredSum <- specPredSum2 <- matrix(0, nseedPred, nspec)
    colnames(specPredSum) <- specNames
  }
  
  gupdate <- c(5, 10, 15, 20, 25, 50, 100, 200, 500, 800, 1200,2000)  
  g1      <- 1
  yupdate <- sort( sample(burnin:ng, 50, replace=T) )
  yupdate <- unique(yupdate)
  pupdate <- burnin:ng
  if(length(pupdate) > 1000)pupdate <- 
    unique( round( seq(burnin, ng, length=100) ))
  
 # TvarEst <- matrix(0,length(pupdate),length(obsRows))
 # colnames(TvarEst) <- columnPaste(tdata$treeID[obsRows], tdata$year[obsRows])
 # TvarPred <- TvarEst
  SvarEst  <- SvarPred <- matrix(0,length(pupdate),length(obsRowSeed))
  colnames(SvarEst)   <- rownames(sdata)[obsRowSeed]
  specSum <- specSum2 <- matrix(0,length(obsRowSeed),nspec)
  colnames(specSum) <- colnames(specSum2) <- specNames
  rownames(specSum) <- rownames(specSum2) <- rownames(sdata)[obsRowSeed]
  
  ntoty <- rmspe <- deviance <- 0
  ntot  <- 0
  zest <- zpred <- fest <- fest2 <- fpred <- fpred2 <- fg*0 # fecundity 
  sumDev <- ndev <- 0   #for DIC
 # spred <- spred2 <- sest <- sest2 <-
 #        matrix(0, nseed, length(seedNames)) # seed predictions
  
  maxFLog <- log(maxF)
  
  activeArea <- sdata$active*sdata$area
  
  nPlotYr <- max(tdata$plotyr)
  acceptRate <- nPlotYr/5
  
  zdmat <- matrix(1, ntree, nyr)                    # this needs ztrue!
  zdmat[ yrIndex[,c('dcol','year')] ] <- tdata$diam
  zdmat[zdmat < minDiam] <- 0                       # 1 is unknown
  rpred <- zdmat*0
  
  pbar <- txtProgressBar(min=1,max=ng,style=1)
  
  logScoreStates <- logScoreFull <- obsRowSeed*0
  
  nprob <- 0
  
  epsilon <- (log(fg) + log(maxF+1 - fg))/100000
  
  for(g in 1:ng){ 
    
    pars$fg <- fg
    yg      <- log(fg)      
    mufec   <- xfec%*%bgFec
    
    if(RANDOM){
      yg <- yg - mufec
      if(YR){
        yg <- yg - betaYrF[yrIndex[,'year']] 
        if(RANDYR)yg <- yg - betaYrR[yrIndex[,c('group','year')]] 
      }
      if(AR)yg <- yg - muyr
      tmp <- .updateAlphaRand(yg[obsRows], xfec[drop=F,obsRows,], z[obsRows], sg, 
                              reIndex[obsRows], reGroups,
                              xrandCols, Arand, priorVA, dfA)
      alphaRand  <- pars$alphaRand <- tmp$alphaRand%*%xrands2u
      
      Arand      <- pars$Arand <- xrands2u%*%tmp$Arand%*%xrands2u
      agibbs[g,] <- as.vector(Arand)                  
      muran <- xfec[,xrandCols]*alphaRand[reIndex,]
      if(length(Arand) > 1)muran <- rowSums( muran )
    }
    
    if(YR){
      yg <- log(fg) - mufec
      if(RANDOM)yg <- yg - muran
      
      tmp      <- .updateBetaYr(yg, z, sg, sgYr, betaYrF, betaYrR, yrIndex, yeGr,
                                RANDYR, tdata$obs)
      betaYrF  <- pars$betaYrF <- tmp$betaYrF
      betaYrR  <- pars$betaYrR <- tmp$betaYrR
      sgYr     <- pars$sgYr    <- tmp$sgYr
      wfinite  <- tmp$wfinite
      
      bygibbsF[g,] <- betaYrF
      bygibbsR[g,] <- betaYrR
      bygibbsN[g,wfinite] <-  1
      sygibbs[g,]  <- sgYr
      muyr <- betaYrF[yrIndex[,'year']] 
      if(RANDYR) muyr <- muyr + betaYrR[ yrIndex[,c('group','year')] ]
    }
    
    if(AR){
      yg <- log(fg) 
      mu <- mufec
      if(RANDOM)mu <- mu + muran
      
      tmp <- .updateBetaAR_RE(betaYrF, betaYrR, Alag, yg, mu, z, 
                              lagGroup, lagMatrix, plag, ngroup, sg)
      betaYrF <- pars$betaYrF <- tmp$betaYrF
      betaYrR <- pars$betaYrR <- tmp$betaYrR
      wfinite <- which(betaYrR != 0)
      Alag    <- tmp$Alag
      muyr <- tmp$ylag
      bygibbsF[g,] <- betaYrF
      bygibbsR[g,] <- betaYrR
      bygibbsN[g,wfinite] <-  1
    }
    
    wlo <- 10*(z - 1)
    whi <- 10*z
    w   <- .tnorm(length(z), wlo, whi, xrep%*%bgRep, 1)
      
    tmp <- .updateBeta(pars, xfec, xrep, R, propU, propF, w, 
                       z, zmat, matYr, muyr)
    ug    <- pars$ug    <- tmp$ug
    umean <- pars$umean <- tmp$umean
    uvar  <- pars$uvar  <- tmp$uvar
    bgFec <- pars$bgFec <- tmp$bgFec
    bgRep <- pars$bgRep <- tmp$bgRep
    propU <- tmp$propU
    

    tmp <- .updateFecundity(g, pars, xfec, xrep, propF, z, zmat, matYr, muyr,
                            epsilon = epsilon, pHMC = .03)

   # if(is.na(max(tmp$fg)))stop('fecundity error')
    
    fg    <- pars$fg <- tmp$fg
    pars$fecMin <- tmp$fecMin
    pars$fecMax <- tmp$fecMax
    
    z     <- tmp$z
    zmat  <- tmp$zmat
    matYr <- tmp$matYr
    propF <- tmp$propF
    epsilon <- tmp$epsilon

    muf <- xfec%*%bgFec
    if(YR | AR)muf <- muf + muyr
    if(RANDOM)muf  <- muf + muran
    
    sg <- pars$sg <- .updateVariance(log(fg[obsRows][z[obsRows] == 1]), 
                                     muf[obsRows][z[obsRows] == 1], s1, s2)
    if(sg > 50)sg <- pars$sg <- 50
      
    if(SAMPR){
      tmp  <- .updateR(ug, fg[obsRows], z[obsRows], SAMPR, distall, 
                     sdata, seedNames, tdata[obsRows,], R, priorR, 
                     priorRwt, obsYr, posR, plots)
      pars$R <- R <- tmp
      rgibbs[g,]  <- R[posR]
    }
    
    # prediction
 
    rpred[ yrIndex[,c('dcol','year')] ] <- rnorm(nrow(tdata),xrep%*%bgRep)
    
    Fpred <- pnorm(rpred)
    Fpred[rpred == 0] <- 0
    dvals <- matrix( runif(ntree, 0, 1), ntree, nyr)
    Fpred[Fpred > dvals] <- 1
    Fpred[zdmat == 0] <- 0
    Fpred[Fpred < 1] <- 0
    
    zw <- Fpred[yrIndex[,c('dcol','year')] ]
    
    zw <- rnorm(nrow(tdata),xrep%*%bgRep)
    zw[zw <= 0] <- 0
    zw[zw > 0] <- 1
    
    flo <- fhi <- zw*0
    flo[zw == 0] <- -maxFLog
    fhi[zw == 1] <- maxFLog
    
    ymu <- .tnorm( length(muf), flo, fhi, muf, sqrt(sg) )
    fmu <- exp( ymu )
 
    # save unstandardized
    bfSave <- bgFec
    brSave <- bgRep
    
    if(UNSTAND){
      if(length(xfecs2u) > 0)bfSave <- xfecs2u%*%bgFec
      if(length(xreps2u) > 0)brSave <- xreps2u%*%bgRep
    }
    
    bfgibbs[g,] <- bfSave
    brgibbs[g,] <- brSave
    ugibbs[g,]  <- ug
    
    
    if( ngroup > 1 )priorUgibbs[g,] <- c(umean,uvar)
    
    
    if(g %in% gupdate){
      gi <- g1:g
      g1 <- g
      if(length(ug) > 1){
        propU  <- apply(ugibbs[gi,],2,sd) + .00001
      }else{
        propU <- sd(ugibbs[gi,]) + .00001
      }
      ng1 <- length(gi)
    }
    
    
    if(g %in% pupdate){
      
      nprob <- nprob + 1
      
      # estimated fecundity per m2
      
      lm <- .getLambda(tdata[obsRows,c('specPlot','year','dcol')],
                       sdata[obsRowSeed,c('year','drow')],
                       AA=1, ug, fg[obsRows], z[obsRows], R, 
                       SAMPR, distall, obsYr, PERAREA=T, SPECPRED=T)  # per m^2
      lm[lm < 1e-9] <- 1e-9
      pm    <- matrix( rpois(length(lm), lm), nrow(lm), ncol(lm) )
      specSum  <- specSum + pm
      specSum2 <- specSum2 + pm^2
      
      #estimated fecundity per trap
      lf <- .getLambda(tdata[obsRows,c('specPlot','year','dcol')],
                       sdata[obsRowSeed,c('year','drow')],
                       activeArea[obsRowSeed], ug, fg[obsRows], z[obsRows], R, 
                       SAMPR, distall, obsYr, PERAREA=F)  # per trap
                       
      lf[lf < 1e-9] <- 1e-9
      pf    <- matrix( rpois(length(lf), lf), nrow(lf), ncol(lf) )
      
      if(MISS){
        sdata[obsRowSeed,seedNames][missSeed] <- pf[missSeed]
      }
      
      SvarEst[nprob,] <- rowSums(pf)
      
      # predicted fecundity per trap
      
      la    <- .getLambda(tdata[obsRows,c('specPlot','year','dcol')], 
                          sdata[obsRowSeed,c('year','drow')],
                          AA=activeArea[obsRowSeed], ug, fmu[obsRows], zw[obsRows], 
                          R, SAMPR, distall, obsYr, PERAREA=F)     # per trap
      la[la < 1e-9] <- 1e-9
      pg    <- matrix( rpois(length(la), la), nrow(lf), ncol(la) )
      rmspe <- sqrt( mean( (seedCount[obsRowSeed,] - pg)^2, na.rm=T ) )
      
      SvarPred[nprob,] <- rowSums(pg)
      
      # deviance from predicted fecundity
      
      dev   <- dpois(seedCount[obsRowSeed,], la, log=T)
      deviance <- -2*sum( dev )
      
      sumDev <- sumDev - deviance
      ndev <- ndev + 1
      
      # score from predicted fecundity
      
      logScoreStates <- logScoreStates - dpois(seedCount[obsRowSeed,], lf, log=T)
      logScoreFull   <- logScoreFull - dev
      
      if(PREDSEED){
      
        ls <- .getLambda(tdat = tdatPred, sdat = sdatPred, AA = 1,
                         ug, fmu[predTree[,'row']], zw[predTree[,'row']], R, 
                         SAMPR, distPred, predList$years, PERAREA=T,
                         SPECPRED=T)   # per m^2
        ls <- ls + 1e-9
        ps <- matrix( rpois(nseedPred*ncol(ls), ls), nseedPred, nspec )
        
        specPredSum  <- specPredSum + ps
        specPredSum2 <- specPredSum2 + ps^2
      }
    }
    
    sgibbs[g,]  <- c(sg, rmspe, deviance)
    
    if(g %in% yupdate){
      
      ntoty  <- ntoty + 1
      
      # by group
      if(AR){
        for(j in 1:ngroup){
          Gmat[1,] <- betaYrF + betaYrR[j,]
          eigenMat[j,] <- eigen(Gmat)$values
        }
        eigen1 <- eigen1 + eigenMat
        eigen2 <- eigen2 + eigenMat^2
      }
      
      fecRes <- fg
      yRes <- yg
      yRes <- yRes - mean(yRes)           #tree autocorrelation
      
        
      for(m in 1:nSpecPlot){
        
        fmat <- fmat*0
        
        wm   <- which(yrIndex[,'specPlot'] == m & tdata$obs == 1 & fg > 1)
        dm   <- tdata$dcol[wm]
        ym   <- yrIndex[wm,'year']
        dr   <- unique(dm)
        if(length(dr) < 2)next
        
        fmat[ cbind(dm,ym) ] <- fecRes[wm]
        fm  <- fmat[unique(dm),unique(ym)]
        
        # between trees
        ff <- suppressWarnings( cor(t(fm), use="pairwise.complete.obs") ) 
        wf <- which(is.finite(ff))
        
        omegaE[dr,dr][wf] <- omegaE[dr,dr][wf] + ff[wf]
        omegaN[dr,dr][wf] <- omegaN[dr,dr][wf] + 1
        
        acm <- acfEmp(yRes[wm], tdata$dcol[wm], yrIndex[wm,'year'])
        wfin <- which(is.finite(acm))
        if(length(wfin) == 0)next
        
        pa <- try( pacf(acm[wfin]), T )
        
        if( !inherits(pa,'try-error') ){
          pacfMat[m,wfin] <- pacfMat[m,wfin] + pa
          pacf2[m,wfin]   <- pacf2[m,wfin] + pa^2
          acfMat[m,wfin] <- acm[wfin]
          pacN[m,wfin] <- pacN[m,wfin] + 1
        }
      }
    }
     
    if(g >= burnin){
      
      ntot <- ntot + 1
      
      zest  <- zest + z
      fz    <- fg*z
      fest  <- fest + fz  # only add to mature 
      fest2 <- fest2 + fz^2
      
      zpred  <- zpred + zw
      fz     <- fmu*zw
      fpred  <- fpred + fz
      fpred2 <- fpred2 + fz^2
      
      if(RANDOM){
        if(length(xrands2u) > 0)bfSave <- alphaRand%*%xrands2u
        asum <- asum + bfSave
        asum2 <- asum2 + bfSave^2
      }
    }
    
    setTxtProgressBar(pbar,g)
  } ###########################################################
   
  
  # to re-initialize
  tdata$lastFec  <- fg
  tdata$lastRepr <- z
    
  # fecundity
  matrEst  <- zest/ntot
  matrPred <- zpred/ntot
  
  fecEstMu <- fest/zest     #  E[fecundity]|z = 1
  fecEstSe <- fest2/zest - fecEstMu^2
  fecEstSe[fecEstSe < 0] <- 0
  fecEstSe <- sqrt(fecEstSe)
  
  fecPredMu <- fpred/zpred     #  pred|z = 1
  fecPredSe <- fpred2/zpred - fecPredMu^2
  fecPredSe[fecPredSe < 0] <- 0
  fecPredSe <- sqrt(fecPredSe)
  
  fecPred <- tdata[,c('plot','treeID','species','year','diam','dcol')]
  if(YR | AR){
    ygr <- yrIndex
    colnames(ygr) <- paste('ind_',colnames(ygr),sep='')
    fecPred <- cbind(fecPred, ygr)
  }
  mpm <- round( cbind(matrEst, matrPred), 3)
  
  fpm <- cbind(fecEstMu, fecEstSe, fecPredMu, fecPredSe)
  fpm[is.na(fpm)] <- 0
  fpm <- round(fpm,1)
  
  fecPred <- cbind(fecPred, mpm, fpm)
  
  scols <- c('plot','year','trapID','drow','area','active')
  countPerTrap <- rowSums(sdata[obsRowSeed,seedNames,drop=F])

  seedEst <- cbind( colMeans(SvarEst), apply(SvarEst, 2, sd) )
  colnames(seedEst) <- c('estMeanTrap', 'estSeTrap')
  seedPred <- cbind( colMeans(SvarPred), apply(SvarPred, 2, sd) )
  colnames(seedPred) <- c('predMeanTrap', 'predSeTrap')
  
  seedSpecMu <- specSum/nprob
  seedSpecSe <- specSum2/nprob - seedSpecMu^2
  seedSpecSe[seedSpecSe < 0] <- 0
  seedSpecSe <- sqrt(seedSpecSe)
  
  
  seedPred <- data.frame( cbind( sdata[obsRowSeed,scols], countPerTrap,
                                 signif(seedEst, 3),signif(seedPred, 3)) )
  
  nss <- length(obsRowSeed)
  mss <- matrix(seedMass[specNames,],nss,nspec,byrow=T)
  massMu <- seedSpecMu[,specNames,drop=F]*mss
  massSe <- sqrt( seedSpecSe[,specNames,drop=F]^2*mss^2 ) 
  
  m1 <- paste(colnames(massMu),'meanM2',sep='_')
  m2 <- paste(colnames(massSe),'sdM2',sep='_')
  
  MM <- F
  if(!all(seedMass == 1))MM <- T
  
  if( MM ){
    m1 <- paste(colnames(massMu),'meanGmM2',sep='_')
    m2 <- paste(colnames(massSe),'sdGmM2',sep='_')
  }
  colnames(massMu) <- m1
  colnames(massSe) <- m2
  seedPred <- cbind(seedPred, signif(massMu, 3), signif(massSe,3))
  
  
  if(PREDSEED){
    scols <- c('plot','trapID','year','x','y','drow','dgrid','area','active')
    specMu <- specPredSum/nprob
    specSe <- specPredSum2/nprob - specMu^2
    colnames(specMu) <- paste(colnames(specMu),'_meanM2',sep='')
    colnames(specSe) <- paste(colnames(specSe),'_sdM2',sep='')
    
    
    preds <- signif(cbind(specMu, specSe), 3)
    
    seedPredGrid <- data.frame( cbind( seedPredGrid[,scols], preds ) )
    treePredGrid <- cbind(predTree, fecPred[predTree$row,])
    seedPredGrid$plot <- as.factor(seedPredGrid$plot)
    treePredGrid$plot <- as.factor(treePredGrid$plot)
    
    # out-of-sample
    if(!is.null(modelYears)){
      
      sdataOut$plotTrapYr <- columnPaste(sdataOut$trapID,sdataOut$year)
      sdataOut$fore <- sdataOut$year - max(modelYears)

      tdataOut$plotTreeYr <- columnPaste(tdataOut$treeID,tdataOut$year)
      treePredGrid$plotTreeYr <- columnPaste(treePredGrid$treeID,treePredGrid$year)
      
      fec <- matrix(NA,nrow(tdataOut),4)
      colnames(fec) <- colnames(fpm)
      ww <- which(tdataOut$plotTreeYr %in% treePredGrid$plotTreeYr)
      qq <- match(tdataOut$plotTreeYr[ww], treePredGrid$plotTreeYr)
      fec[ww,] <- fpm[qq,]
      tdataOut <- cbind(tdataOut,fec)
    }
  }
  
  kg <- burnin:ng
  
  betaFecMu <- matrix( colMeans(bfgibbs[drop=F,kg,]), Qfec )
  betaFecSe <- matrix( apply(bfgibbs[drop=F,kg,],2,sd), Qfec )
  betaFecCI <- t( apply( bfgibbs[drop=F,kg,],2,quantile,c(.025,.975) ) )
  rownames(betaFecMu) <- rownames(betaFecSe) <- rownames(betaFecCI) <- xFecNames
  betaFec <- cbind(betaFecMu, betaFecSe, betaFecCI)

  betaRepMu <- matrix( colMeans(brgibbs[drop=F,kg,]), Qrep )
  betaRepSe <- matrix( apply(brgibbs[drop=F,kg,],2,sd), Qrep )
  betaRepCI <- t( apply( brgibbs[drop=F,kg,],2,quantile,c(.025,.975) ) )
  rownames(betaRepMu) <- rownames(betaRepSe) <- rownames(betaRepCI) <- xRepNames
  betaRep <- cbind(betaRepMu, betaRepSe, betaRepCI)
  colnames(betaFec)[1:2] <- colnames(betaRep)[1:2] <- c('estimate','se')
  
  # seed, fecundity acf
  seedRes <- t( t(seedCount) - colMeans(seedCount, na.rm=T) )
  ww <- colSums(seedCount, na.rm=T)
  seedRes <- seedRes[,ww > 0,drop=F]
  
  ii <- rep(sdata$drow,ncol(seedRes))
  yy <- match(sdata$year,years)
  jj <- rep(yy, ncol(seedRes))
  
  acm  <- acfEmp(as.vector(seedRes), ii, jj)
  wfin <- which(is.finite(acm))
  pacsMat <- acm*0
  pacsMat[wfin] <- pacf(acm[wfin])
  
  acfMat  <- acfMat/ntoty
  pacfMat <- pacfMat/pacN
  pacfSe  <- pacf2/pacN - pacfMat^2
  
  omegaE <- omegaE/omegaN
  
  ################ years/lags
  
  if(YR | AR){
    ncol <- nyr
    ccol <- years
    if(AR){
      ncol <- plag
      ccol <- colnames(bygibbsF)
    }
    betaYrMu <- matrix( colMeans(bygibbsF[kg,], na.rm=T), nrow(betaYrF), ncol )
    betaYrSe <- matrix( apply(bygibbsF[kg,], 2, sd, na.rm=T), nrow(betaYrF), ncol )
    colnames(betaYrMu) <- colnames(betaYrSe) <- ccol
    betaYrRand <- betaYrRandSE <- betaYrMu*0
    if(RANDYR){
      brsum <- matrix( colSums(bygibbsR[kg,], na.rm=T), 
                       nrow(betaYrR), ncol )
      brn <- matrix( colSums(bygibbsN[kg,], na.rm=T), 
                     nrow(betaYrR), ncol )
      betaYrRand <- brsum/brn
      brn2 <- matrix( colSums(bygibbsR[kg,]^2, na.rm=T), 
                      nrow(betaYrR), ncol )
      ser <- sqrt( brn2/brn - betaYrRand^2 )
      betaYrRandSE <- ser
      betaYrRand[!is.finite(betaYrRand)] <- 0
      rownames(betaYrRand) <- rownames(betaYrRandSE) <- yeGr
      colnames(betaYrRand) <- colnames(betaYrRandSE) <- ccol
    }
  }
  
  ################ REs
  
  if(RANDOM){
    alphaMu <- asum/ntot
    av    <- asum2/ntot - alphaMu^2
    av[av < 0] <- 0
    alphaSe <- sqrt(av)
    if(ONEA){
      aMu <- mean(agibbs)
      aSe <- sd(agibbs)
    }else{
      colnames(alphaMu) <- colnames(alphaSe) <- xFecNames[xrandCols]
      rownames(alphaMu) <- rownames(alphaSe) <- as.character()
      aMu <- matrix( apply(agibbs, 2, mean), Qrand, Qrand )
      aSe <- matrix( apply(agibbs, 2, sd), Qrand, Qrand )
      colnames(aMu) <- rownames(aMu) <- colnames(aSe) <- rownames(aSe) <-
        colnames(alphaMu)
      names(xrandCols) <- xFecNames[xrandCols]
    }
  }
  
  mmu <- 1
  if(SAMPR){
    
    tmu <- apply(rgibbs[kg,], 2, mean)
    tse <- apply(rgibbs[kg,], 2, sd)
    
    mmu <- mse <- priorR*0
    mmu[posR] <- tmu
    mse[posR] <- tse
    
    attr(mmu,'posR') <- attr(mse,'posR') <- posR
    
    colnames(mmu) <- colnames(mse) <- seedNames
    rownames(mmu) <- rownames(mse) <- rownames(R)
  }
  
  #################### dispersal
  
  dgibbs <- pi*sqrt(ugibbs[kg,])/2
  
  if( ngroup > 1){
    uByGroup <- colMeans(ugibbs[kg,])
    dByGroup <- colMeans(dgibbs)
    priorDgibbs <- pi*sqrt(priorUgibbs[kg,])/2
  
    upars <- t( rbind( uByGroup, apply(ugibbs[kg,], 2, sd) ) )
    upars <- cbind(upars, t(apply(ugibbs[kg,],2,quantile,c(.025,.975))) )
  
    dpars <- t(rbind( dByGroup, apply(dgibbs, 2, sd) ) )
    dpars <- cbind(dpars, t(apply(dgibbs,2,quantile,c(.025,.975))) )
  
    uall <- cbind( colMeans(priorUgibbs[kg,]), apply(priorUgibbs[kg,], 2, sd),
                  t( apply(priorUgibbs[kg,],2,quantile,c(.025,.975))) )
    dall <- cbind( colMeans(priorDgibbs), apply(priorDgibbs, 2, sd),
                  t( apply(priorDgibbs,2,quantile,c(.025,.975))) )
  
  colnames(uall)[1:2] <- colnames(dall)[1:2] <- c('global mean','global var')
  upars <- rbind(upars,uall)
  dpars <- rbind(dpars,dall)
 
  }else{
    uu <- quantile(ugibbs[kg,],c(.025,.975))
    dd <- quantile(dgibbs,c(.025,.975))
    upars <- matrix( c( mean(ugibbs[kg,]), sd(ugibbs[kg,]), uu), 1)
    dpars <- matrix( c( mean(dgibbs), sd(dgibbs), dd), 1)
    uByGroup <- upars[1]
  }
  colnames(upars) <- colnames(dpars) <- c('estimate','std err','2.5%','97.5%')
  rownames(upars)[rownames(upars) == 'umean'] <- 'mean'
  rownames(upars)[rownames(upars) == 'var'] <- 'variance'
  rownames(dpars) <- rownames(upars)
  
  su <- cbind( colMeans(sgibbs,na.rm=T), apply(sgibbs,2,sd,na.rm=T), 
               t(apply(sgibbs,2,quantile,c(.025,.975),na.rm=T)) )
  colnames(su)[1:2] <- c('estimate','std err')
  
  # coefficients are saved unstandardized 
  
  if(UNSTAND){
    xfec <- xfecU
    xrep <- xrepU
  }
  
  zp <- as.vector( pnorm(xrep%*%betaRepMu) )
  fp <- xfec%*%betaFecMu
  if(YR){
    byr <- betaYrMu[ yrIndex[,'year'] ] 
    if(RANDYR)byr <- byr + betaYrRand[ yrIndex[,c('group','year')] ]
    byr[is.na(byr)] <- 0
    fp <- fp + byr
  }
  
  ############### AR
  
  if(AR){

    eigenMu <- eigen1/ntoty
    eigenSe <- sqrt( eigen2/ntoty - eigenMu^2 )
    
    yg   <- log(fecPred$fecPredMu)      
    yg[!is.finite(yg)] <- 0
    ylag <- yg*0
    mu   <- fp
    
    nl <- nrow(lagMatrix)
    
    zlag <- matrix(zp[lagMatrix[,-1]],nl,plag)
    zlag[zlag < .5] <- 0
    zlag[zlag > 0] <- 1
    xm <- matrix( yg[ lagMatrix[,-1]],nl,plag)*zlag
    ylag[lagMatrix[,1]] <- xm%*%t(betaYrMu)
    
    # random effects
    if(RANDYR){
      for(m in 1:ngroup){
        tg <- which(lagGroup == m)
        if(length(tg) == 0)next
        ylag[lagMatrix[tg,1]] <- xm[tg,]%*%betaYrRand[m,]
      }
    }
    ylag[!is.finite(ylag)] <- 0
    fp <- mu + ylag
  }
  
  ################ fit
  
  fp <- exp(fp + su['sigma','estimate']/2)
  
  meanDev <- sumDev/ndev
  zp[zp >=.5] <- 1
  zp[zp < 1] <- 0
  
  la <- .getLambda(tdata[obsRows,c('specPlot','year','dcol')], 
                   sdata[obsRowSeed,c('year','drow')], activeArea[obsRowSeed],
                   uByGroup, as.vector(fp[obsRows]), as.vector(zp[obsRows]), 
                   mmu, SAMPR, distall, obsYr, PERAREA=F)
  la <- la + 1e-9
  pd  <- meanDev - 2*sum(dpois(seedCount, la, log=T))
  DIC <- pd + meanDev
  
  RMSPE <- mean(sgibbs[(1+burnin):ng,'rmspe'])
  
  logScoreStates <- logScoreStates/nprob
  logScoreFull   <- logScoreFull/nprob
  
  
  fit <- list( DIC = round(DIC), scoreStates = mean(logScoreStates),
               predictScore = signif( mean(logScoreFull), 3), 
               RMSPE = signif(RMSPE,3) )
  
  inputs$treeData   <- tdata
  inputs$seedData   <- sdata
  inputs$formulaFec <- formulaFec
  inputs$formulaRep <- formulaRep
  inputs$upar       <- ug
  inputs$ng         <- ng
  inputs$burnin     <- burnin
  inputs$plotDims   <- plotDims
  inputs$plotArea   <- plotArea
  inputs$specNames  <- specNames
  inputs$seedNames  <- seedNames
  inputs$priorR     <- priorR
  inputs$priorRwt   <- inputs$priorRwt
  inputs$xytree     <- xytree
  inputs$xytrap     <- xytrap
  inputs$yrIndex    <- yrIndex
  inputs$obsRowSeed <- obsRowSeed
  inputs$obsRows    <- obsRows
  inputs$maxF       <- maxF
  if(!is.null(plotDims))inputs$plotDims <- plotDims
  if(!is.null(seedMass))inputs$seedMass <- seedMass
  if(!is.null(yearEffect) & !'yearEffect' %in% names(inputs))
               inputs$yearEffect <- yearEffect

  chains <- list(bfec = .orderChain(bfgibbs, specNames), 
                 brep = .orderChain(brgibbs, specNames), 
                 ugibbs = ugibbs, sgibbs = sgibbs)
  parameters <- list( betaFec = signif(betaFec, 3),  
                      betaRep = signif(betaRep, 3),
                      sigma = signif(su), acfMat = acfMat, 
                      upars = signif(upars, 4), dpars = signif(dpars,4),
                      pacfMat = pacfMat,
                      pacfSe = pacfSe, pacsMat = pacsMat, omegaE = omegaE,
                      omegaN = omegaN)
  
  if(SAMPR){
    chains <- append(chains, list(rgibbs = rgibbs))
    parameters <- append(parameters, list(rMu = signif(mmu,3), 
                                          rSe = signif(mse,3)))
  }
  prediction <-  list( fecPred = fecPred, seedPred = seedPred, 
                      # entropy = entropy, sdResource = sdResource,
  #                     logScoreStates = logScoreStates,
                       predictScore   = signif( sum(logScoreFull), 3))
  if(ngroup > 1)chains <- append(chains, list(priorUgibbs = priorUgibbs))
  if(YR){
    chains <- append(chains, list(sygibbs = sygibbs))
  }
  if(AR){
    parameters <- append(parameters, 
                         list(eigenMu = eigenMu, eigenSe = eigenSe))
    prediction <- append(prediction, 
                         list(tdataOut = tdataOut, sdataOut = sdataOut))
  }
  if(AR | YR){
    if(yeGr[1] %in% specNames){
      tmp <- .orderChain(bygibbsR, specNames)
    }
    
    chains     <- append(chains, list(bygibbsF = bygibbsF, bygibbsR = bygibbsR) )
    parameters <- append(parameters, list(betaYrMu = signif(betaYrMu, 3),
                                          betaYrSe = signif(betaYrSe, 3)))
    
    if(RANDYR)parameters <- append(parameters, 
                                   list(betaYrRand = signif(betaYrRand,3),
                                        betaYrRandSE = signif(betaYrRandSE,3)))
  }
  if(PREDSEED) {
    prediction <- append(prediction, 
                         list(seedPredGrid = seedPredGrid,
                              treePredGrid = treePredGrid) )
    if(!'predList' %in% names(inputs))inputs <- append(inputs, list(predList = predList))
  }
  if(RANDOM){
    inputs$randomEffect <- randomEffect
    chains     <- append(chains, list(agibbs = .orderChain(agibbs, specNames)) )
    parameters <- append(parameters, 
                         list(alphaMu = alphaMu, alphaSe = alphaSe,
                              aMu = aMu, aSe = aSe) )
  }
  data$setupData$distall <- distall
  
  chains     <- chains[ sort( names(chains) )]
  inputs     <- inputs[ sort( names(inputs) )]
  data       <- data[ sort(names(data)) ]
  inputs$ng  <- ng
  inputs$burnin <- burnin
  parameters <- parameters[ sort( names(parameters) )]
  
  out <- list(inputs = inputs, chains = chains, data = data, fit = fit, 
              parameters = parameters, prediction = prediction )
  
  class(out) <- 'mastif'
  out
} 

meanVarianceScore <- function(output, ktree = 30, maxSite = 100,
                              LAGMAT=F, Q = c(.5, .025, .975), nsim = 1,
                              CLOSE=F){
  
  # ntree  - no. trees visited
  # cspace - no. m2 visited
  # cyr    - no. years
  # CLOSE = T: sample small neighborhoods (selected to be close)
  # CLOSE = F: random neighborhoods
  
  lagCanopy <- lagGround <- NULL
  
  sdata   <- output$inputs$seedData
  tdata   <- output$inputs$treeData
  plots   <- output$data$setupData$plots
  years   <- output$data$setupData$years
  xytrap  <- output$inputs$xytrap
  xytree  <- output$inputs$xytree
  specNames <- output$data$setupData$specNames
  seedNames <- output$data$setupData$seedNames
  seedMass  <- output$inputs$seedMass
  yrIndex   <- output$data$setupData$yrIndex
  nplot     <- length(plots)
  nyr       <- length(years)
  ntree     <- nrow(xytree)
  ntrap     <- nrow(xytrap)
  obsYr     <- sort(unique(sdata$year))
  plotDims  <- output$inputs$plotDims
  predList  <- output$inputs$predList
  
  fecPred      <- output$prediction$fecPred
  seedPred     <- output$prediction$seedPred
  seedPredGrid <- output$prediction$seedPredGrid
  
  fmat <- matrix(0, ntree, nyr)
  colnames(fmat) <- obsYr
  rownames(fmat) <- tdata$treeID[ match( unique(tdata$dcol),tdata$dcol ) ]
  
  if(is.null(seedMass)){
    seedMass <- matrix(1,length(specNames),1)
    rownames(seedMass) <- specNames
  }
  
  scoreT <- scoreS <- deltaT <- deltaS <- numeric(0)
  scoreTse <- scoreSse <- deltaTse <- deltaSse <- numeric(0)
  treeCor <- trapCor <- numeric(0)
  entropy <- domain <- numeric(0)
  resourceScore <- resourceMean <- totalVar <- numeric(0)
  win <- floor(nyr/2)
  if(win > 10)win <- 10
  
  GRID <- SMASS <- F
  
  trapID <- as.character(seedPred$trapID)
  allTraps   <- unique(trapID)
  seedPred$x <- xytrap$x[match(trapID,xytrap$trapID)]
  seedPred$y <- xytrap$y[match(trapID,xytrap$trapID)]
  meanCols <- grep('_meanGmM2', colnames(seedPred))
  if(length(meanCols) == 0){
    meanCols <- grep('_meanM2', colnames(seedPred))
    rs       <- columnSplit(colnames(seedPred)[meanCols],'_meanM2')[,1]
    rSeedMat <- seedPred[,meanCols,drop=F]*
                matrix(seedMass[rs,1],nrow(seedPred),length(rs),byrow=T)
    colnames(rSeedMat) <- paste(specNames,'_meanGmM2',sep='')
    seedPred <- cbind(seedPred, rSeedMat)
  }
    
  sgrid <- seedPred
  GRID  <- F
  
  #replace with prediction grid if all years predicted
  
  if(!is.null(seedPredGrid)){
    
    predPlots <- predList$plots
    predYears <- predList$years
    
    oldGrid <- seedPred
    newGrid <- numeric(0)
    wcol <- which(colnames(seedPred) %in% colnames(seedPredGrid))
    massCols  <- grep('_meanGmM2', colnames(seedPred))
    areaCols <- grep('_meanM2', colnames(seedPredGrid))
    
    pcols <- c("plot","year","trapID", "drow","area","active", "x","y")           
    
    for(m in 1:nplot){
      
      wn <- which(seedPred$plot %in% plots[m])
      wm <- which(seedPredGrid$plot %in% plots[m])
      myr <- sort(unique(seedPred$year[wn]))           #all obs years in pred grid?
      pyr <- sort(unique(seedPredGrid$year[wm]))
      AYR <- all(myr %in% pyr)
      
      if( AYR & length(wm) > 0 ){
        mgrid <- seedPredGrid[wm,colnames(seedPred)[wcol]]
        acc   <- seedPredGrid[wm,areaCols, drop=F]
        
        mgrid <- cbind(mgrid, acc)
        mcols <- grep('_meanM2', colnames(mgrid))
        rs    <- columnSplit(colnames(mgrid)[mcols],'_meanM2')[,1]
        rSeedMat <- mgrid[,mcols,drop=F]*
          matrix(seedMass[rs,1],nrow(mgrid),length(rs),byrow=T)
        colnames(rSeedMat) <- paste(specNames,'_meanGmM2',sep='')
        mgrid <- cbind(mgrid[,-mcols], rSeedMat)
        
        newGrid <- rbind(newGrid, mgrid)
        if(length(oldGrid) > 0){
          wo <- which(oldGrid$plot %in% plots[m])
          if(length(wo) > 0)oldGrid <- oldGrid[-wo,]
        }
      }
    }
    if(length(newGrid) > 0){
      rownames(oldGrid) <- NULL
      pnn   <- match(pcols, colnames(oldGrid))
      mc  <- grep('_meanGmM2', colnames(oldGrid))
      sgrid <- rbind(newGrid,oldGrid[,c(pnn, mc)])
    }
    GRID  <- T
  }
  
  meanCols <- grep('_meanGmM2', colnames(sgrid))
 
  trapID  <- as.character(sgrid$trapID)
  allTrap <- sort(unique(trapID))
  ntrap   <- length(allTrap)
  smat <- matrix(0, ntrap, nyr)
  rownames(smat) <- allTrap
    
  colnames(smat) <- obsYr
  meanNames <- c('trees_PerTree','sites_PerSite','trees_PerYr', 'sites_PerYr')
  eNames <- c('tree-tree','site-site','tree-lag','site-lag')
  scoreNames <- c('gmTree','gmM2',eNames)
  
  kk   <- 0
  if(nsim > 1){
    pbar <- txtProgressBar(min=1,max=nplot*nsim,style=1)
    cat('\nScore\n')
  }
  
  if(win > 1){
    
    rjtree <- rjtrap <- rjall <- character(0)
    
    darea <- 100
    ddist <- round( sqrt(darea), 0)
    rseq  <- c(0, seq(10, 2000, by = ddist))
    
    for(m in 1:nplot){
      
      wp <- fecPred$plot == plots[m]
      wo <- tdata$obs == 1
      wm <- fecPred$matrEst > .3
      wy <- fecPred$year %in% obsYr
      
      wk <- wp&wo&wm&wy
      wc <- which(wk)
      
      if(length(wc) == 0)next
      
      yrm <- sort(unique(c(fecPred$year[wp],sgrid$year[sgrid$plot == plots[m]])))
      yrm <- yrm[yrm %in% obsYr]
      
      dm   <- tdata$dcol[wc]
      dr   <- unique(dm)             # unique trees
      
      if(length(dr) <= 1)next
      
      
      ym   <- fecPred$year[wc]
      ym   <- match(ym,obsYr)
      
      ytab  <- table(fecPred$year[wc])
      ykeep <- as.numeric(names(ytab))[ytab > 0]
      ckeep <- match(ykeep, obsYr)
      
      emat <- rmat <- vmat <- matrix(0, nsim, 4)
      cmat <- matrix(0, nsim, 6)
      size <- matrix(0, nsim, 5)
      
      #   ctree <- length(dr)
      tseq  <- c(1:length(dr))
      
      scols <- paste('0_',c(0:win),sep='')
      
      ssmat <- matrix(0, length(rseq), length(scols))
      colnames(ssmat) <- scols
      rownames(ssmat) <- rseq
      
      tdmat <- matrix(0, length(tseq),length(scols))
      colnames(tdmat) <- scols
      rownames(tdmat) <- tseq
      
      sdmat <- nsdmat <- ssmat <- nstmat <- ssmat*0
      tdmat <- ntdmat <- ttmat <- nttmat <- tdmat*0
      
      sdmat2 <- ssmat2 <- ssmat*0
      tdmat2 <- ttmat2 <- tdmat*0
      
      rjtree <- c(rjtree,plots[m])
      
      for(k in 1:nsim){
        
        kk <- kk + 1
        if(nsim > 1)setTxtProgressBar(pbar, kk)
        
        
        entTtree <- entTseed <- entYtree <-  entYseed <- rtot <- rTtree <- 
          rTseed <- rYtree <- rYseed <- Tk <- Ts <- Yk <- Ys <- carea <- 
          varTk <- varTs <- varYk <- varYs <- NA
        
        fcor <- NULL

        if(length(dr) > 1){   # canopy 
          
          fmat <- fmat*0

          yall <- sort(unique(fecPred$year[wc]))
          nyrk <- length(yall)
          
          tkeep <- ckeep
          
          fmat[ cbind(dm,ym) ] <- fecPred$fecEstMu[wc]
          rkeep <- which(rowSums(fmat) > .99)
          if(length(rkeep) > ktree)rkeep <- sample(rkeep, ktree)
          fm    <- fmat[drop=F,rkeep,tkeep]
          
          if(length(fm) > 2 & length(tkeep) > 2){
            
            nyrk <- ncol(fm)
            ntrk <- nrow(fm)
            
            sm   <- tdata$species[ match(rownames(fm),tdata$treeID) ]
            sm   <- as.character( sm )
            rvec <- seedMass[sm, 1]
            
            if(k == 1){
              fcor <- makeCrossCov( tmat = fm, win = win, MATRIX=T, COR=T )[[1]]
              fcor[!is.finite(fcor)] <- 0
            }
            
            if(LAGMAT){
              
              ff <- sample(nrow(fm))
              
              tt <- makeCrossCov( tmat = fm[drop=F,ff,], win = win, MATRIX=T, COR=F )
              mcov <- tt$lagCov
              
              if(length(mcov) > 1){
                
                tmu  <- tt$lagMean
                tmean <- mean(tmu[,1])
                
                gg <- grep( paste(rownames(mcov)[1],'_-',sep=''), colnames(mcov) )
                gg <- c(grep(paste(rownames(mcov)[1],'_0', sep=''), 
                             colnames(mcov) ), gg)
                tmu <- tmu[,gg,drop=F]
                
                scov <- mcov[,gg,drop=F]
                vars <- mcov[ cbind(rownames(mcov), paste(rownames(mcov),'_0',sep='')) ]
                vars <- matrix(vars, nrow(scov), ncol(scov))
                
                scov[-1] <- scov[-1]*2   #covariances count twice
                
                wtrow <- matrix(1:nrow(scov),nrow(scov),ncol(scov)) 
                wtcol <- matrix(1:ncol(scov),nrow(scov),ncol(scov), byrow=T) 
                
                v1 <- wtcol*vars[1]   # count diagonal elements
                v2 <- wtrow*vars
                
                scov[1] <- 0
                scov <- scov + v1 + v2
                
                scum0 <- t(apply(scov, 1, cumsum))
                scum  <- apply(t(scum0), 1, cumsum)
                if(!is.matrix(scum))scum <- scum0*0 + scum
                
                
                tcum0 <- t(apply(tmu, 1, cumsum))
                tcum <- apply(t(tcum0), 1, cumsum)
                if(!is.matrix(scum))tcum <- tcum0*0 + tcum
                
                
                tscore <- log(tcum) - 1/2*log(scum)
                rscore <- log( sum(tmu) ) - 1/2*log( length(tmu)*var(as.vector(fm)) )
                delta  <- tscore - rscore
                
                ir <- 1:nrow(delta)
                ic <- 1:ncol(delta)
                
                tdmat[ir,ic]  <- tdmat[ir,ic] + delta
                tdmat2[ir,ic] <- tdmat2[ir,ic] + delta^2
                ntdmat[ir,ic] <- ntdmat[ir,ic] + 1
                ttmat[ir,ic]  <- ttmat[ir,ic] + tscore
                ttmat2[ir,ic] <- ttmat2[ir,ic] + tscore^2
                nttmat[ir,ic] <- nttmat[ir,ic] + 1
              }
            }
            
            names(rvec) <- sm
            rmm   <- fm*rvec          
            wm    <- which(rmm > 0)
            rtree <- mean( rmm[wm] ) # mean per reproductive tree
            
            Tk  <- cov(t(rmm))       # tree cov
            Yk  <- cov(rmm)          # year cov
            
            varTk <- sum(Tk)
            varYk <- sum(Yk)
            
            rTtree <- rYtree <- sum(rmm)
            
            if(!is.na(max(Tk))){
              tmp        <- var2score(rTtree, varTk, Tk) # var between trees 
              entTtree   <- tmp$entropy
            }
            if(!is.na(max(Yk))){
              tmp        <- var2score(rYtree, varYk, Yk) # var between years
              entYtree   <- tmp$entropy
            }
          } #end canopy
          
          lagCanopy <- append(lagCanopy, list(fcor))
          
        }
        
        ############# seed traps or seed prediction grid
        
        smat <- smat*0
        
        wy <- unique(sgrid$year[sgrid$plot == plots[m]])
        wy <- wy[wy %in% obsYr]
        wm <- which(sgrid$plot == plots[m] & 
                      sgrid$year == wy[1])
        ix <- sample(wm, 1)
        wo <- wm[wm != ix]
        
        dist <- NULL
        scor <- NULL
        
        if(CLOSE){
          distSite <- .distmat( sgrid$x[ix], sgrid$y[ix], 
                                sgrid$x[wo], sgrid$y[wo] )
          oo <- order(distSite)
          #     if(length(oo) > cspace)oo <- oo[1:(cspace-1)]
          wm <- c(ix,wo[oo])
          dist <- c(0,distSite[oo])
        }else{
          #     if(length(wm) > cspace)wm <- sample(wm,cspace)
        }
        
        trapIDS <- as.character( sgrid$trapID[wm] )
        if(!is.null(dist))names(dist) <- trapIDS
        wm <- which( as.character(sgrid$trapID) %in% trapIDS &
                       sgrid$year %in% yrm &
                       sgrid$plot == plots[m])
        j <- match(sgrid$year[wm], obsYr)
        i <- match(as.character(sgrid$trapID[wm]),rownames(smat))
        
        smat[ cbind(i,j) ] <- rowSums(sgrid[wm,meanCols,drop=F])
        
        yall <- sort(unique(j))
        nyrj <- length(yall)
        
        wsp <- columnSplit(rownames(smat),'-')[,1]
        rkeep <- which(wsp == plots[m])
        
        #    rkeep <- which(rowSums(smat) > 0)
        rmm <- smat[drop=F,rkeep, as.character(yrm)]
        
        crr <- which(colSums(rmm) > 0)
        rmm <- rmm[,crr]
        
        
        gind  <- match(rownames(rmm),sgrid$trapID)
        carea <- diff(range(sgrid$x[gind]))*
          diff(range(sgrid$y[gind]))
        if(carea == 0 | !is.finite(carea))carea <- 100
        
        
        if(LAGMAT & length(crr) > 2){
          
          if(!is.null(dist)){
            attr(rmm,'distance') <- dist
            i  <- findInterval(dist,rseq)
            ii <- rep( i, ncol(rmm) )
            jj <- rep( 1:ncol(rmm), each=nrow(rmm) )
            tmp <- .myBy(as.vector(rmm), ii, jj, fun='mean')*darea  # per m2
            rownames(tmp) <- rseq[ min(i):max(i) ]
            colnames(tmp) <- colnames(rmm)
            tmp <- tmp[,colSums(tmp) != 0, drop=F]
            tmp <- tmp[drop=F, rowSums(tmp) != 0,]
          }else{
            tmp <- rmm
            if(nrow(tmp) > maxSite)tmp <- tmp[sample(nrow(rmm),maxSite),]
          }
          
          tt   <- makeCrossCov( tmat = tmp, win = win, MATRIX=T, COR=F )
          mcov <- tt$lagCov
          tmu  <- tt$lagMean
          
          zrows <- rownames(mcov)
          zcols <- paste(rownames(mcov),'_0',sep='')
          
          if(!is.null(dist)){
            
            p1 <- paste('^',rownames(mcov)[1],'_',sep='')
            
            gg <- grep(p1,colnames(mcov))
            scov <- mcov
            #      if(length(gg) > 0){
            scov <- mcov[,gg, drop=F]
            tmu <- tmu[,gg, drop=F]
            #      }
            gg <- grep('-',colnames(scov))
            hh <- grep('_0',colnames(scov))
            gg <- sort(c(gg,hh))
            scov <- scov[,gg, drop=F]
            tmu <- tmu[,gg, drop=F]
            
          }else{
            gg <- grep( paste(rownames(mcov)[1],'_-',sep=''), 
                        colnames(mcov) )
            gg <- c(grep(paste(rownames(mcov)[1],'_0', sep=''), 
                         colnames(mcov) ), gg)
            scov <- mcov[,gg]
            tmu <- tmu[,gg]
            
            wz    <- which(!zcols %in% colnames(mcov))
            if(is.matrix(scov) & length(wz) > 0){
              zcols <- zcols[-wz]
              zrows <- zrows[-wz]
              scov  <- scov[-wz,]
              tmu   <- tmu[-wz,]
            }
            #   vars <- mcov[ cbind(zrows, zcols) ]
          }
          
          
          if(length(mcov) > 2)vars <- mcov[ cbind(zrows, zcols) ]
          
          if(is.matrix(scov)){
            
            vars <- matrix(vars, nrow(scov), ncol(scov))
            
            scov[-1] <- scov[-1]*2   #covariances count twice
            
            wtrow <- matrix(1:nrow(scov),nrow(scov),ncol(scov)) 
            wtcol <- matrix(1:ncol(scov),nrow(scov),ncol(scov), byrow=T) 
            
            v1 <- wtcol*vars[1]   # count diagonal elements
            v2 <- wtrow*vars
            
            scov[1] <- 0
            scov <- scov + v1 + v2
            
            scum <- t(apply(scov, 1, cumsum))
            scum <- apply(t(scum), 1, cumsum)
            
            tcum <- t(apply(tmu, 1, cumsum))
            tcum <- apply(t(tcum), 1, cumsum)
            
            tscore <- log(tcum) - 1/2*log(scum)
            rscore <- log( sum(tmu) ) - 1/2*log( length(tmu)*var(as.vector(tmp)) )
            delta  <- tscore - rscore
            
            delta  <- vec2mat(delta)
            tscore <- vec2mat(tscore)
            
            if(!is.null(dist)){
              if(is.null(rownames(delta)))
                rownames(delta) <- rownames(tscore) <- rownames(scov)
              rownames(delta) <- rownames(tscore) <- 
                columnSplit(rownames(delta),'_')[,1]
              colnames(delta) <- colnames(tscore) <- 
                colnames(ssmat)[1:ncol(delta)]
            }else{
              rownames(delta) <- rownames(tscore) <- 
                rownames(sdmat)[1:nrow(delta)]
              ccc <- columnSplit(colnames(delta),'_')[,2]
              ccc <- .replaceString(ccc,'-','')
              ccc <- paste('0_',ccc,sep='')
              colnames(delta) <- colnames(tscore) <- ccc
            }
            
            #      ir <- 1:nrow(delta)
            ic <- 1:ncol(delta)
            rd <- rownames(delta)
            
            sdmat[rd, ic]  <- sdmat[rd, ic] + delta
            sdmat2[rd, ic] <- sdmat2[rd, ic] + delta^2
            nsdmat[rd, ic] <- nsdmat[rd, ic] + 1
            ssmat[rd, ic]  <- ssmat[rd, ic] + tscore
            ssmat2[rd, ic] <- ssmat2[rd, ic] + tscore^2
            nstmat[rd, ic] <- nstmat[rd, ic] + 1
          }
          
          
          if(k == 1){
            scor <- makeCrossCov( tmat = tmp, win = win, MATRIX=T, COR=T )[[1]]
            scor[!is.finite(scor)] <- 0
            
            rjtrap   <- c(rjtrap, plots[m])
          }
        }
        
        if(k == 1)lagGround <- append(lagGround, list(scor))
        
        if(is.matrix(rmm) & length(rmm) > 1){
          
          rtot <- sum(rmm)
          
          Ts <- cov(t(rmm))  # site cov
          Ys <- cov(rmm)     # year cov
          
          varTs <- sum(Ts)
          varYs <- sum(Ys)
          
          rTseed <- rYseed <- sum(rmm)
          
          tmp        <- var2score(rTseed, varTs, Ts)
          entTseed   <- tmp$entropy
          tmp        <- var2score(rYseed, varYs, Ys)
          entYseed   <- tmp$entropy
          
          n1 <- n2 <- n3 <- n4 <- NA
          if(is.matrix(Tk))n1 <- nrow(Tk)
          if(is.matrix(Ts))n2 <- nrow(Ts)
          if(is.matrix(Yk))n3 <- nrow(Yk)
          if(is.matrix(Ys))n4 <- nrow(Ys)
          
          emat[k,] <- c(entTtree, entTseed, entYtree, entYseed)
          cmat[k,] <- c(rtree, rtot)
          rmat[k,] <- c(rTtree, rTseed, rYtree, rYseed)
          size[k,] <- c(n1,n2,n3, n4, round(carea)) 
          vmat[k,] <- c(varTk, varTs, varYk, varYs)
        }
      }
     
    if( LAGMAT){
      deltaTree <- tdmat/ntdmat
      scoreTree <- ttmat/nttmat
      deltaSeed <- sdmat/nsdmat
      scoreSeed <- ssmat/nstmat
      
      deltaTrSe <- tdmat2/ntdmat - deltaTree^2
      scoreTrSe <- ttmat2/nttmat - scoreTree^2
      deltaSdSe <- sdmat2/nsdmat - deltaSeed^2
      scoreSdSe <- ssmat2/nstmat - scoreSeed^2
      
      wr <- which(rowSums(deltaTree,na.rm=T) != 0)
      wc <- which(colSums(deltaTree,na.rm=T) != 0)
      
      deltaTree <- deltaTree[wr,wc]
      deltaTrSe <- sqrt(deltaTrSe[wr,wc])
      scoreTree <- scoreTree[wr,wc]
      scoreTrSe <- sqrt(scoreTrSe[wr,wc])
      
      wr <- which(rowSums(deltaSeed,na.rm=T) != 0)
      wc <- which(colSums(deltaSeed,na.rm=T) != 0)
      
      deltaSeed <- deltaSeed[wr,wc]
      scoreSeed <- scoreSeed[wr,wc]
      deltaSdSe <- sqrt(deltaSdSe[wr,wc])
      scoreSdSe <- sqrt(scoreSdSe[wr,wc])
      
      
      deltaT <- append(deltaT, list(deltaTree))
      deltaS <- append(deltaS, list(deltaSeed))
      scoreT <- append(scoreT, list(scoreTree))
      scoreS <- append(scoreS, list(scoreSeed))
      deltaTse <- append(deltaTse, list(deltaTrSe))
      deltaSse <- append(deltaSse, list(deltaSdSe))
      scoreTse <- append(scoreTse, list(scoreTrSe))
      scoreSse <- append(scoreSse, list(scoreSdSe))
      
    }
    
    emat[!is.finite(emat)] <- NA
    
    ee <- signif(t( apply(emat, 2, quantile, Q,na.rm=T )),3)
    rownames(ee) <- paste(plots[m],eNames,sep='_')
    vv <- t( apply(vmat, 2, quantile, Q,na.rm=T ))
    rownames(vv) <- paste(plots[m],eNames,sep='_')
    
    rjall <- c(rjall,plots[m])
    
    ii     <- apply(size,2,mean)
    domain <- rbind(domain, ii)
    
    entropy  <- rbind( entropy, ee )
    totalVar <- rbind( totalVar, vv)
    
    #   if(LAGMAT){
    #     treeCor <- append(treeCor, list(signif(fcor,3)))
    #     trapCor <- append(trapCor, list(signif(scor,3)))
    #   }
    
  } #########end plot loop
  
    if(LAGMAT){
      names(lagCanopy) <- rjtree
        names(lagGround) <- rjtrap
      names(scoreT) <- rjtree
        names(scoreS) <- rjtrap
      names(deltaT) <- rjtree
        names(deltaS) <-  rjtrap
    }
    
    rownames(domain) <- rjall
    colnames(domain) <- c('trees','sites','treeYr','siteYr','areaM2')
    
    tpr <- paste(rjtree,'tree-tree',sep='_')
    VtreeXtree <- totalVar[drop=F,tpr,]
    tpr <- paste(rjtree,'site-site',sep='_')
    VsiteXsite <- totalVar[drop=F,tpr,]
    tpr <- paste(rjtree,'tree-lag',sep='_')
    VtreeXlag <- totalVar[drop=F,tpr,]
    tpr <- paste(rjtree,'site-lag',sep='_')
    VsiteXlag <- totalVar[drop=F,tpr,]
    
    rownames(VtreeXtree) <- columnSplit(rownames(VtreeXtree),'_tree-tree')[,1]
    rownames(VsiteXsite) <- columnSplit(rownames(VsiteXsite),'_site-site')[,1]
    rownames(VtreeXlag)  <- columnSplit(rownames(VtreeXlag),'_tree-lag')[,1]
    rownames(VsiteXlag)  <- columnSplit(rownames(VsiteXlag),'_site-lag')[,1]
    
  }else{
    domain <- entropy <- 
      VtreeXtree <- VsiteXsite <- VtreeXlag <- VsiteXlag <- NULL
  }
  
  tmp <- list(domain = domain, entropy = entropy, 
       treeXtreeVar = VtreeXtree, siteXsiteVar = VsiteXsite, 
       treeXlagVar = VtreeXlag, siteXlagVar = VsiteXlag)
  if(LAGMAT){
    tmp$lagCanopy <- lagCanopy
    tmp$lagGround <- lagGround
    tmp$scoreSeed <- scoreS
    tmp$scoreTree <- scoreT
    tmp$deltaSeed <- deltaS
    tmp$deltaTree <- deltaT
    tmp$scoreSeedSe <- scoreSse
    tmp$scoreTreeSe <- scoreTse
    tmp$deltaSeedSe <- deltaSse
    tmp$deltaTreeSe <- deltaTse
  }
  for(k in 1:length(tmp)){
    if(is.null(tmp[[k]]))next
    kcol <- which(sapply(tmp[[k]],is.numeric))
    jcol <- which(sapply(tmp[[k]],is.factor))
    kcol <- intersect(kcol,!jcol)
    for(j in kcol){
      tmp[[k]][,j][ is.finite( tmp[[k]][,j] ) ] <- NA
    }
  }
    
  tmp
    
}





 
plotMeanVarScoreByDim <- function(x, mu, lo=NULL, hi=NULL, xlab='',
                                  ylab='score', ylim=NULL, YLOG=F,
                                  LINES = T, UNLOGY=F){
  BARS <- F
  if(!is.null(lo))BARS <- T
  if(UNLOGY){
    mu <- exp(mu)
    if(BARS){
      lo <- exp(lo)
      hi <- exp(hi)
    }
  }
  
  if(is.null(ylim)){
    if(!is.null(lo)){
      ylim <- range(c(lo, hi), na.rm=T)
    }else{
      ylim <- range(mu, na.rm=T)
    }
  }
  xlim <- range(x, na.rm=T)
  xlim[2] <- xlim[2] + .2
  xlim[1] <- xlim[1] - .2
  
  gfun <- colorRampPalette( c("forestgreen","#8DD3C7", "#BEBADA", "#FB8072",
                              "#80B1D3", "#FDB462", "brown") )
  plotCol <- gfun(ncol(mu))
  
  #  xs <- 1 + seq(-.2,.2,length=ncol(mu))
  
  if(!YLOG){
 #   plot(NULL, xlim = xlim, ylim = ylim, xlab=xlab, ylab=ylab, log='x')
    plot(NULL, xlim = xlim, ylim = ylim, xlab=xlab, ylab=ylab)
  }else{
    if(ylim[1] == 0)ylim[1] <- .01
    plot(NULL, xlim = xlim, ylim = ylim, xlab=xlab, ylab=ylab, log='xy')
  }
  for(j in 1:ncol(mu)){
    xs <- x[,j]
    xs[xs < xlim[1]] <- xlim[1]
  #  points(xs, mu[,j], col=plotCol[j], lwd=2)
    lohi <- cbind(lo[,j], hi[,j])
    .shadeInterval(xs, loHi=lohi, add=T, col=.getColor(plotCol[j],.1))
  #  if(LINES)lines(xs, mu[,j], col=plotCol[j], lwd=2)
  #  if(BARS){
  #    suppressWarnings(arrows(xs, lo[,j], xs, hi[,j], col=plotCol[j],
  #                 angle=90, length=.04, code=3))
  #  }
  }
  for(j in 1:ncol(mu)){
    xs <- x[,j]
    xs[xs < xlim[1]] <- xlim[1]
    lines(xs, mu[,j], col='white', lwd=4)
#    lines(xs, mu[,j], col=plotCol[j], lwd=2)
#  }
  }
  for(j in 1:ncol(mu)){
    xs <- x[,j]
    xs[xs < xlim[1]] <- xlim[1]
#    lines(xs, mu[,j], col='white', lwd=4)
    lines(xs, mu[,j], col=plotCol[j], lwd=2)
    #  }
  }
}

var2score <- function(rmean, totVr, rvar){
  
  # rmean - mean over sites or years
  # rvar  - covariance matrix
  # totVr - total variance
  # ndim  - dimension of covariance matrix
  
  if(length(rvar) < 4)return( list(score = NA, entropy = NA) )
  
  ndim <- nrow(rvar)
  
  score <- log(rmean) - 1/2*suppressWarnings(log(totVr))
  dt    <- determinant(rvar)$modulus
  if(!is.finite(dt)){
    ev <- eigen(rvar)$values
    dt <- sum( log(ev[ev > 0]) )
  }
  entr  <- ndim/2*(1 + log(2*pi)) + dt/2
  entr  <- entr/ndim
  list(score = score, entropy = entr)
}

plotMeanVarScore <- function(MVS, x1 = 'treeXtreeMu', xlim = NULL,
                             ylim = NULL, xlab=NULL, ylab=NULL, BARS=T, TEXT=F){
  
 # rmean <- MVS$resourceMean
 # score <- MVS$resourceScore
  pl    <- MVS$plots
  
  if(is.null(ylab))ylab <- 'Score'
  
  if(x1 == 'treeXtreeMu'){
    x <- MVS$treeXtreeMu
    y <- MVS$treeXtreeScore
    if(is.null(xlab))xlab <- 'g/tree'
  }
  if(x1 == 'treeXlagMu'){
    x <- MVS$treeXlagMu
    y <- MVS$treeXlagScore
    if(is.null(xlab))xlab <- 'g/tree'
  }
  if(x1 == 'siteXsiteMu'){
    x <- MVS$siteXsiteMu
    y <- MVS$siteXsiteScore
    if(is.null(xlab))xlab <- 'g/m2'
  }
  if(x1 == 'siteXlagMu'){
    x <- MVS$siteXlagMu
    y <- MVS$siteXlagScore
    if(is.null(xlab))xlab <- 'g/m2'
  }
  
  cex <- 1
  if(TEXT)cex = .01
  
  plots <- rownames(x)
  nplot <- length(plots)
  gfun <- colorRampPalette( c("forestgreen","#8DD3C7", "#BEBADA", "#FB8072",
                              "#80B1D3", "#FDB462", "brown") )
  plotCol <- gfun(nplot)
  names(plotCol) <- plots
  
  if(is.null(ylim))ylim <- range(y,na.rm=T)
  if(is.null(xlim))xlim <- range(x,na.rm=T)
  plot(x[,1], y[,1], ylim = ylim, xlab='g/tree',
       ylab='', log='x', xlim = xlim, col=plotCol, cex=cex)
  if(TEXT)text(x[,1], y[,1], plots, col=plotCol)
  if(BARS){
    segments(x[,1], y[,2],x[,1], y[,3], col=plotCol)
    segments(x[,2], y[,1],x[,3], y[,1], col=plotCol)
  }
  invisible(plotCol)
}

crossCovSetup <- function(tmat, win){
  
  nt <- ncol(tmat)
  if(win > nt/2)win <- floor(nt/2)
  
  ni    <- nrow(tmat)
  lead  <- -c(-win:win)
  ntt   <- length(lead)
  mgrid <- as.matrix( expand.grid(1:ntt,1:ni,1:ni ) )
  colnames(mgrid) <- c('t','i1','i2')
  
  ld    <- lead[mgrid[,'t']]
  mgrid <- cbind(ld,mgrid)
  colnames(mgrid)[1] <- 'lead'
  
  cc    <- columnPaste(mgrid[,'i1'],mgrid[,'i2'])
  mdex  <- columnPaste(mgrid[,'t'],cc)
  keep  <- which(mgrid[,'i2'] >= mgrid[,'i1'])
  mgrid <- mgrid[keep,]
  mdex  <- mdex[keep]
  
  mgrid <- as.data.frame(mgrid)
  
  rn <- rownames(tmat)
  if(!is.null(rn)){
    mgrid$ID1 <- rn[mgrid[,'i1']]
    mgrid$ID2 <- rn[mgrid[,'i2']]
  }
    
  list(win = win, ntt = ntt, ni = ni, mgrid = mgrid, 
       mdex = mdex, lead = lead)
}
  
makeCrossCov <- function(tmat, win = 5, MATRIX=F, COR=F ){
  
  # tmat - responses by time matrix
  # cross covariance of each row against population
  # MATRIX  - n by n*lag matrix
  # !MATRIX - n*n*lag vector
  # COR - correlation matrix
  
  tiny <- 1e-8
  
  nt  <- ncol(tmat)
  mid <- round(nt/2)
  rt  <- mid + c(-win,win)
  
  if(rt[2] <= rt[1])return(NULL)
  
  imat <- sweep(tmat, 1, rowMeans(tmat), '-')
  imat[imat == 0] <- tiny                        # no variation
  
  if(win > (ncol(tmat)-1))win <- ncol(tmat) - 1
  
  tmp <- crossCovSetup(tmat, win)
  ntt <- tmp$ntt
  ni  <- tmp$ni
  mgrid <- tmp$mgrid
  mdex  <- tmp$mdex
  lead  <- tmp$lead
  win   <- tmp$win
  
  crossCov <- rep(0,length(mdex))
  names(crossCov) <- mdex
  totSeed <- crossCov
  
  if(COR)ivar <- apply(imat,1,var)*(nt-1)/nt
  
  for(i in 1:(win+1)){
    
    ii <- 1:(nt - win + i - 1)
    pp <- (win - i + 2):nt
    
    ii  <- ii[ii > 0]
    pp  <- pp[pp <= nt]
    ldd <- pp[1] - ii[1]
    
    for(m in 1:ni){
      wm   <- which(mgrid[,'i1'] == m & mgrid[,'lead'] == ldd )
      mdx  <- mgrid[drop=F,wm,]
      
      tres <- rowMeans( imat[mdx[,'i1'],ii,drop=F]*imat[mdx[,'i2'],pp,drop=F] )
      
      tm2 <- tmat[mdx[,'i2'],pp,drop=F]
      if(ldd == 0)tm2[ rownames(tm2) == rownames(tmat[mdx[,'i1'],]) ] <- 0
      tsum <- rowMeans( tmat[mdx[,'i1'],ii,drop=F] + tm2 )
      
      if(COR)tres <- tres/sqrt(ivar[mdx[1,'i1']]*ivar[mdx[,'i2']])
      crossCov[ wm ] <- tres
      totSeed[ wm ]  <- tsum
        
      if(ldd != 0){
        wm   <- which(mgrid[,'i1'] == m & mgrid[,'lead'] == -ldd )
        mdx  <- mgrid[drop=F,wm,]
        tres <- rowMeans( imat[mdx[,'i1'], pp,drop=F]*imat[mdx[,'i2'], ii,drop=F] )
        if( identical(mdx[,'i1'], mdx[,'i2']) ){
          tres <- rowMeans( tmat[mdx[,'i1'],ii,drop=F] )
        }else{
          tsum <- rowMeans( tmat[mdx[,'i1'], pp,drop=F] + tmat[mdx[,'i2'], ii,drop=F] )
        }
        
        if(COR)tres <- tres/sqrt(ivar[mdx[1,'i1']]*ivar[mdx[,'i2']])
        crossCov[ wm ] <- tres
        totSeed[ wm ] <- tsum
      }
    }
  }
  
  totSeed <- round(totSeed, 2)
  
  covMu <- cbind(mgrid, crossCov, totSeed)
  
  
  
#  wt  <- which(covMu[,'crossCov'] != 0)
#  covMu <- covMu[wt,]
  
  if(!MATRIX)return( lagCov = covMu, lagMean = NULL )
  
  t2 <- covMu$t[covMu$lead == 0][1]
  
  rn <- columnPaste(covMu[,'ID1'],covMu[,'lead'],'_')
  cn <- columnPaste(covMu[,'ID2'],t2,'_')
  rt <- rn[!duplicated(rn)]
  ct <- cn[!duplicated(cn)]
  
  stmat <- ttmat <- matrix(0,length(rt),length(ct))
  rownames(stmat) <- rownames(ttmat) <- rt
  colnames(stmat) <- colnames(ttmat) <- ct
  cii <- cbind(rn,cn)
  stmat[ cii ] <- covMu[,'crossCov']
  stmat <- t(stmat)
  
  ttmat <- stmat*0
  tmean <- rowMeans(tmat)
  rm    <- columnSplit(ct,'_')[,1]
  ttmat[1:nrow(ttmat),] <- tmean[rm]
  
 # css <- columnSplit(rownames(stmat),'_')
 # rnn <- css[,1]
 # if( ncol(css) > 2 ){
 #   for(k in 2:(ncol(css)-1)){
 #     rnn <- columnPaste(rnn,css[,k])
 #   }
 # }
  rnn <- columnSplit(rownames(stmat),'_', LASTONLY=T)[,1]
  rownames(stmat) <- rnn
  list(lagCov = stmat, lagMean = ttmat)
}
  
.updateBetaAR <- function(betaYr, yg, mu, z,lagGroup, plag, ngroup, sg)  {
  
  # AR(p), fixed groups
  
  ylag <- yg*0
  
  for(m in 1:ngroup){
    
    lmm <- lagGroup[[m]]
    lmm <- lmm[drop=F,z[lmm[,plag+1]] == 1,]
    if(nrow(lmm) <= (ncol(lmm)+5))next
    
    ym  <- yg[lmm[,1]] - mu[lmm[,1]]
    xm  <- matrix( yg[ lmm[,-1] ], nrow(lmm) )
    V   <- solve( crossprod(xm) )*sg
    v   <- crossprod(xm,ym)/sg
    tmp <- rmvnormRcpp(1,V%*%v,V) 
    whi <- which( abs(tmp) > 1)
    if(length(whi) > 0){
      tmp <- .tnormMVNmatrix(tmp,tmp,V,tmp*0 - 1,tmp*0 + 1, 
                             whichSample=c(1:length(tmp)) )
    }
    betaYr[m,] <- tmp
    
    lmm <- lagGroup[[m]]
    xm  <- matrix( yg[ lmm[,-1] ], nrow(lmm) )
    ylag[lmm[,1]] <- xm%*%t(tmp)
  }
  
  wfinite <- which(!betaYr == 0)
  list( betaYr = betaYr, ylag = ylag, wfinite = wfinite)
}

.updateBetaAR_RE <- function(betaYrF, betaYrR, Alag,
                             yg, mu, z, lagGroup, lagMatrix, plag, ngroup, sg){
  
  # AR(p), random groups
  # betaYrF - 1 by plag fixed effects
  # betaYrR - ngroup by plag random effects
  # Alag    - random effects covariance
  
  ylag <- yg*0
  AIlag <- solve(Alag)
  
  # fixed effects
  cg   <- which(z[lagMatrix[,plag+1]] == 1)  # mature plag yr ago
  
  if(length(cg) <= plag)
    stop('too few mature individuals in AR groups--fewer groups or smaller p')
  lmm  <- lagMatrix[drop=F,cg,]
  xm   <- matrix( yg[ lmm[,-1] ], nrow(lmm))
  mvec <- yg[ lmm[,1] ] - mu[ lmm[,1] ] - 
          rowSums( betaYrR[ lagGroup[cg],]*xm )  # remove random AR effects
  v <- crossprod(xm,mvec)/sg
  V <- solve(crossprod(xm)/sg + diag(1,plag))
 
  tmp <- rmvnormRcpp(1,V%*%v,V) 
  whi <- which( abs(tmp) > 1)                               # for stability (may not be desirable)
  if(length(whi) > 0){
    tmp <- .tnormMVNmatrix(tmp,tmp,V,tmp*0 - 1,tmp*0 + 1, 
                           whichSample=c(1:length(tmp)) )
  }
  
  betaYrF <- tmp
  ylag[lmm[,1]] <- xm%*%t(betaYrF)
  
  if(ngroup == 1) return( list( betaYrF = betaYrF, betaYrR = betaYrR, ylag = ylag, Alag = Alag) )
  
  # random effects
  
  for(m in 1:ngroup){
    tg  <- lagGroup == m
    cg  <- z[lagMatrix[,plag+1]] == 1  # mature plag yr ago
    wm  <- which(tg & cg)
    
    if(length(wm) < 2){      #
      betaYrR[m,] <- 0
      next
    }
      
    if(length(wm) < plag){
      V <- AIlag
      v <- t(betaYrF)*0
    }else{
      lmm <- lagMatrix[drop=F, wm, ]
      xm  <- matrix( yg[ lmm[,-1] ], nrow(lmm) )
      mvec <- yg[lmm[,1]] - mu[lmm[,1]] - xm%*%t(betaYrF)
      v    <- crossprod(xm, mvec)/sg
      V    <- solve(crossprod(xm)/sg + Alag)
    }
    tmp <- rmvnormRcpp(1,V%*%v,V) 
    whi <- which( abs(tmp) > 1)
    if(length(whi) > 0){
      tmp <- .tnormMVNmatrix(tmp,tmp,V,tmp*0 - 1,tmp*0 + 1, 
                             whichSample=c(1:length(tmp)) )
    }
    betaYrR[m,] <- tmp 
    if(length(wm) > plag)ylag[lmm[,1]] <- ylag[lmm[,1]] + xm%*%betaYrR[m,]
  }
  
  # random effects covariance
  LL   <- crossprod(betaYrR)
  Alag <- .updateCovariance(LL, diag(1,plag), ngroup, plag+1)
  
  list( betaYrF = betaYrF, betaYrR = betaYrR, ylag = ylag, Alag = Alag)
}

acfEmp <- function(res, irow, time){
  
  # empirical (for model-based use .updateBetaAc)
  # assumes res values have unique irow and time
  # detrend
  mt <- max(time)
  st <- c(0:(mt-1))
  ni <- max(irow)
  
  resMat <- matrix(0,ni,mt)
  resMat[ cbind(irow,time) ] <- res
  xy <- t( matrix(st*t(resMat),mt) ) 
  xx <-  matrix(st^2,ni,mt,byrow=T) 
  bs <- rowSums(xy,na.rm=T)/rowSums(xx,na.rm=T)
  bi <- rowMeans(resMat,na.rm=T) - bs*mean(st)
  yy <- resMat - bi + matrix(bs,ncol=1)%*%st 
  res <- yy[ cbind(irow, time) ]
  
 # res <- res - predict.lm( lm(res ~ time) )
  
  FF  <- .myBy(res, irow, time, matrix(0,max(irow),max(time)),fun='sum')
  FF  <- crossprod(FF)
  z   <- FF[lower.tri(FF, diag=T)]
  ii  <- row(FF)-col(FF)
  ii  <- ii[lower.tri(ii, diag=T)] + 1
  tmp <- as.vector(.myBy(z,ii,ii*0+1,fun='sum'))
  tmp  <- tmp/tmp[1]

  names(tmp) <- paste('lag',st,sep='-')
  
  tmp
}

pacf <- function(xacf){
  
  nj <- length(xacf)
  st <- c(0:(nj-1))
  rmat <- diag(1,nj)
  rmat <- matrix( xacf[1 + abs( row(rmat)-col(rmat) )], nj, nj)
  tmp <- c(1, solve( toeplitz(xacf[1:(nj-1)]), xacf[2:nj] ) )
  names(tmp) <- paste('lag',st,sep='-')
  tmp
}

.mapSetup <- function(xlim,ylim,scale=NULL,widex=10.5,widey=6.5){  
  
  #scale is x per inch
  #new means not a new plot
  
  if(is.null(scale))scale <- 1
  
  px   <- diff(xlim)/scale
  py   <- diff(ylim)/scale
  
  if(px > widex){
    dx <- widex/px
    px <- widex
    py <- py*dx
  }
  if(py > widey){
    dx <- widey/py
    py <- widey
    px <- px*dx
  }
  
  par(pin=c(px,py))
  invisible( c(px,py) )
}


checkPlotDims <- function(plots, years, xytree, xytrap, plotDims, plotArea){
  
  if(is.null(plotDims)){
    plotDims <- getPlotDims(xytree, xytrap)
  }else{
    rownames(plotDims) <- .fixNames(rownames(plotDims), all=T)$fixed
    if(ncol(plotDims) != 5)stop('plotDims has 5 columns: xmin, xmax, ymin, ymax, area')
    wc <- which(!plots %in% rownames(plotDims))
    if(length(wc) > 0){
      xx <- paste(plots[wc],' missing from plotDims ')
      stop(xx)
    }
  }
  plotDims <- plotDims[plots,]
  
  if(is.null(plotArea)){
    plotArea <- matrix(plotDims[,'area'],nrow(plotDims),length(years))
    rownames(plotArea) <- plots
    colnames(plotArea) <- years
  }else{
    rownames(plotArea) <- .fixNames(rownames(plotArea), all=T)$fixed
  }
    
  list(plotDims = plotDims, plotArea = plotArea)
}



mastMap <- function(mapList, mapPlot, mapYears, treeSymbol = NULL, PREDICT = F, 
                    treeScale = 1, trapScale = 1, xlim = NULL, ylim = NULL,
                    MAPTRAPS = T, seedMax = NULL, fecMax = NULL,
                    mfrow = NULL, LEGEND = F, plotScale = 1,
                    SCALEBAR = F, scaleValue=NULL, COLORSCALE = T){
  
  # if PREDICT, needs seedPredGrid, treePredGrid
  
  if(!is.null(scaleValue))SCALEBAR <- T
  if(SCALEBAR & is.null(scaleValue))scaleValue <- 20
  
  if( class(mapList) == 'mastif' ){
    mapList <- append(mapList$inputs, mapList$prediction)
  }
  
  if('treeData' %in% names(mapList))tdata <- mapList$treeData
  if('tdata' %in% names(mapList))   tdata <- mapList$tdata
  if('seedData' %in% names(mapList))sdata <- mapList$seedData
  if('sdata' %in% names(mapList))   sdata <- mapList$sdata
  specNames <- mapList$specNames
  seedNames <- mapList$seedNames
  xytree    <- mapList$xytree
  xytrap    <- mapList$xytrap
  fecPred   <- mapList$fecPred
  seedPred  <- mapList$seedPred
  
  if(length(specNames) == 1)LEGEND <- F
  
  if(!'trapID' %in% names(xytrap))
    xytrap$trapID <- columnPaste(xytrap$plot,xytrap$trap)
  if(!'trapID' %in% names(sdata))
    sdata$trapID <- columnPaste(sdata$plot,sdata$trap)
  if(!'treeID' %in% names(xytree))
    xytree$treeID <- columnPaste(xytree$plot,xytree$tree)
  if(!'treeID' %in% names(tdata))
    tdata$treeID <- columnPaste(tdata$plot,tdata$tree)
  
  ww <- which(tdata$year %in% mapYears & 
                as.character(tdata$plot) %in% mapPlot)
  
  if(PREDICT){
    if(!'seedPredGrid' %in% names(mapList)){
      PREDICT <- F
    }else{
      seedPredGrid <- mapList$seedPredGrid
      treePredGrid <- mapList$treePredGrid
      wz <- which(seedPredGrid$year %in% mapYears & 
                    seedPredGrid$plot %in% mapPlot)
      ww <- c(ww,wz)
    }
  }
  
  if(length(ww) == 0){
    cat('\nNo obs or preds for mapPlot, mapYears\n')
    return(add=F)
  }
  
  snames <- paste(seedNames,'_meanM2',sep='')
  sdata  <- sdata[sdata$plot %in% mapPlot &
                    sdata$year %in% mapYears,]
  seedCount <- as.matrix(sdata[,seedNames])
  
  wtree <- which(tdata$plot %in% mapPlot &
                   tdata$year %in% mapYears)
  tree <- tdata[wtree,]
  if(!is.null(treeSymbol))treeSymbol <- treeSymbol[wtree]
  
  if(is.null(seedMax) & length(seedCount) > 0)
    seedMax <- max( seedCount, na.rm=T )
  
  if(!is.null(fecPred)){
    fmat <- fecPred[fecPred$plot %in% mapPlot &
                      fecPred$year %in% mapYears,]
    treeSymbol <- fmat$fecEstMu
  }else{
    if(is.null(treeSymbol))treeSymbol <- tree$diam
    fmat <- tree
    fmat$fecEstMu <- treeSymbol
  }
  
  if(is.null(fecMax))fecMax <- max( treeSymbol, na.rm=T )
  
  nspec <- length(specNames)
  
  cfun <- colorRampPalette( c('#66c2a5','#fc8d62','#8da0cb') )
  specCol <- cfun(nspec) 
  names(specCol) <- specNames
  
  xlim <- ylim <- dx <- dy <- numeric(0)
  npp  <- length(mapPlot)
  
  checkXY <- apply(xytree[,c('x','y')], 2, min)
  if(!all(is.finite(checkXY)))stop('missing x, y in xytree')
  
  checkXY <- apply(xytrap[,c('x','y')], 2, min)
  if(!all(is.finite(checkXY)))stop('missing x, y in xytrap')
  
  for(j in 1:npp){
    
    wxy1 <- which(xytree$plot == mapPlot[j])
    wxy2 <- which(xytrap$plot == mapPlot[j])
    
    xlimj <- range( c(xytree[wxy1,'x'], xytrap[wxy2,'x']) )
    ylimj <- range( c(xytree[wxy1,'y'], xytrap[wxy2,'y']) )
    if(PREDICT){
      wj <- which(seedPredGrid$plot == mapPlot[j])
      xlimj <- range( c(xlimj, seedPredGrid$x[wj]) )
      ylimj <- range( c(ylimj, seedPredGrid$y[wj]) )
    }
    dxj <- diff(xlimj)
    dyj <- diff(ylimj)
    
    dx <- c(dx, dxj)
    dy <- c(dy, dyj)
    xlim <- rbind(xlim, xlimj)
    ylim <- rbind(ylim, ylimj)
  }
  xlimit <- matrix(xlim,npp,2,byrow=F)
  ylimit <- matrix(ylim,npp,2,byrow=F)
  rownames(xlimit) <- rownames(ylimit) <- mapPlot
  
  rr  <- apply( rbind(xlimit, ylimit), 2, range)
  sc  <- max( apply(rr,1,diff) )/20
  
  opin <- par()$pin
  
  stab <- with( sdata, table(plot, year) )
  stab <- stab[drop=F, mapPlot, ]
  obs  <- stab[, colnames(stab) %in% mapYears, drop=F]
  oyr  <- colnames(obs)[ colSums(obs) > 0 ]
  pyr  <- character(0)
  pred <- NULL
  
  if(PREDICT){
    ptab <- with( seedPredGrid, table(plot, year) )
    if(!mapPlot %in% rownames(ptab)){
      cat('\nNo prediction for this plot\n')
      PREDICT <- F
    }else{
      ptab <- ptab[drop=F, mapPlot, ]
      wss  <- which(colnames(ptab) %in% mapYears)
      if(length(wss) == 0){
        cat('\nNo prediction for this plot-year\n')
        PREDICT <- F
      }else{
        pred <- ptab[, colnames(ptab) %in% mapYears, drop=F]
        pyr  <- colnames(pred)[ colSums(pred) > 0 ]
      }
    }
  }
  yr  <- sort( unique( c(oyr, pyr) ) )
  if( length(oyr) > 0 ){
    oyr <- oyr[oyr %in% yr]
    obs <- obs[, oyr, drop=F]
  }else{
    obs <- NULL
  }

  if(PREDICT){
    if( length(pyr) > 0 ){
      pyr <- pyr[pyr %in% yr]
      pred <- pred[, pyr, drop=F]
    }
  }
  
  specAll  <- table(tree$species)
  specAll  <- specAll[specAll > 0]
  specPred <- table(tree$species[tree$plot %in% rownames(pred) &
                                   tree$year %in% pyr])
  specPred  <- names(specPred)[specPred > 0]
  colList <- numeric(0)
  
  for(j in 1:npp){
    
    WOJ <- WPJ <- F
    
    jobs <- jpred <- jyr <- numeric(0)
    
    if(!is.null(obs)){
      jobs  <- obs[drop=F, mapPlot[j],]
      WOJ <- T
    }
    if(!is.null(pred)){
      jpred <- pred[drop=F, mapPlot[j],]
      WPJ <- T
    }
    jyr   <- sort(unique(as.numeric(c(colnames(jobs),colnames(jpred)))))
    njyr  <- length(jyr)
    
    if(is.null(mfrow))mfrow <- c(1,1)
    
    suppressWarnings(
      par(bty='o',mar=c(1,.4,2,.4), oma=c(3,3,1,1))
    )
    if(LEGEND)par(oma=c(4,4,1,4))
    
    if(njyr == 1){
      .mapSetup(xlimit[j,], ylimit[j,], scale = max(c(dx[j],dy[j]))/3)
    }else{
      if(prod(mfrow) == 1)mfrow <- .getPlotLayout(njyr)$mfrow
      par( mfrow=mfrow, mar=c(2,2,2,1), bty='o' )
      .mapSetup(xlimit[j,], ylimit[j,], 
                scale = max(c(dx[j],dy[j]))/plotScale/3*max(mfrow) )
    }
    
    for(k in 1:njyr){
      
      add <- WO <- WP <- F
      
      if(WOJ){
        if( yr[k] %in% colnames(obs) &
            MAPTRAPS )WO <- obs[1,colnames(obs) == yr[k]] > 0
      }
      if(WPJ){
        if( yr[k] %in% colnames(pred) )WP <- pred[,colnames(pred) == yr[k]] > 0
      }
      
      if(!WO & !WP)next
      
      if(WP){  #predicted surface, fecundity
        
        tmp <- .pmap(specNames=specNames, xytree=xytree, 
                     plot=mapPlot[j], 
                     year=jyr[k], seedPredGrid=seedPredGrid, 
                     treePredGrid = treePredGrid,
                     xlim=xlimit[j,], ylim=ylimit[j,], treeScale, trapScale,
                     sCol = specCol[specNames], add=add)
        add <- tmp$add
        tmp <- tmp[names(tmp) != 'add']
        
        if(add)colList <- append(colList, list(tmp))
        names(colList)[length(colList)] <- mapPlot[j]
      }
      
      if(WO){  #observed seed
        
        seed <- sdata[sdata$year == jyr[k],]
        sx <- xytrap$x[ match(seed$trapID, xytrap$trapID) ]
        sy <- xytrap$y[ match(seed$trapID, xytrap$trapID) ]
        z  <- as.matrix( seed[,seedNames] )
        
        if(length(snames) > 1)z <- rowSums(z,na.rm=T)
        
        w1 <- which(z > 0)
        w0 <- which(z == 0)
        if(length(w1) > 0){
          z <- z/seedMax
          z <- 5*sc*z*trapScale
          symbols(sx[w1], sy[w1], squares=z[w1], inches = F,
                  xlab='', ylab='',bg = .getColor('black',.3), 
                  fg = .getColor('black',.5), add=add,
                  xaxt='n', yaxt='n', xlim = xlim, ylim = ylim)
          for(i in 1:4)axis(i, labels=F, tck=.01)
          add <- T
        }
        if(add == F){
          plot(NULL, xlim = xlim, ylim = ylim, xlab='', ylab='',
               axes = F)
          for(i in 1:4)axis(i, labels=F, tck=.01)
        }
        if(length(w0) > 0)points(sx[w0], sy[w0], pch=0, cex=.5,
                                 col=.getColor('black',.6))
        add <- T
      }
      
      if(WO & !WP){ #est fecundity
        
        treek <-  fmat[fmat$year == jyr[k],]
        
        if(nrow(treek) == 0)next
        
        sx <- xytree$x[ match(treek$treeID, xytree$treeID) ]
        sy <- xytree$y[ match(treek$treeID, xytree$treeID) ]
        
        z  <- treek$fecEstMu
        if(!all(is.na(z))){
          mmm <- max(z, na.rm=T)
          if(mmm == 0)mmm <- 1
          z <- 1*sc*z/mmm*treeScale
          ic <- match(treek$species, specNames)
          
          symbols(sx, sy, circles=z*1.3, inches = F, add=add,
                  fg = .getColor('white',.5), 
                  bg =.getColor('white',.5), xaxt='n', yaxt='n')
          
          symbols(sx, sy, circles=z, inches = F, add=add,
                  fg = .getColor(specCol[ specNames[ic] ],.6), 
                  bg=.getColor(specCol[ specNames[ic] ],.3), xaxt='n', yaxt='n')
          if(!add)for(i in 1:4)axis(i, labels=F, tck=.01)
        }
      }
      
      if(!add)next
      
      cex <- sum(mfrow)^(-.1)
      .plotLabel(jyr[k], 'topright',cex = cex, above=F) 
      .plotLabel(mapPlot[j],'topleft',cex = cex,above=F)
    }
  }
  if(PREDICT)colList <- colList[ !duplicated(names(colList)) ]
  
  if(SCALEBAR)scaleBar('m', value = scaleValue, yadj=.07, cex=.8)
  if(LEGEND){
    cornerLegend('bottomright', names(specAll),
           text.col = specCol[match(names(specAll),specNames)],
           cex=.7, bty='n')
  }
  
  if(COLORSCALE & PREDICT){
    
    # use last plot (bottom or right side)
    cols <- colList[[length(colList)]]
    nss   <- length(cols$species)
    endLabels <- NULL
    
    for(k in 1:nss){
      
      if(k == nss)endLabels <- c(0, signif(seedMax,1) )
      
      ck <- .getColor(specCol[cols$species[k]], cols$colorLevels)
      
      clist <- list( kplot=k, ytick=NULL, text.col = 'black',
                     cols=ck, labside='right', text.col=col,
                     bg='grey', endLabels=endLabels) 
      cornerScale( clist )
  #    kk <- kk + 1
    }
  }
  invisible(add)
}

cornerLegend <- function(...) {
  suppressWarnings(
  opar <- par(fig=c(0, 1, 0, 1), oma=c(0, 0, 0, 0), 
              mar=c(0, 0, 0, 0), new=TRUE)
  )
  on.exit(par(opar))
  plot(0, 0, type='n', bty='n', xaxt='n', yaxt='n')
  legend(...)
}


cornerScale <- function( clist ) {
  
  opar <- par(fig=c(0, 1, 0, 1), oma=c(0, 0, 0, 0), 
              mar=c(0, 0, 0, 0), new=T)
  on.exit(par(opar))
  plot(0, 0, type='n', bty='n', xaxt='n', yaxt='n')
  
  .cornerLegendScale( clist )
}

.cornerLegendScale <- function( clist ){  
  
  # left and right corners: xx = (x1,x2), y = (y1,y2)
  # bg is color of border
  # cols  - matching color sequence
   kplot <- 1
  
  opar <- par(fig=c(0, 1, 0, 1), oma=c(0, 0, 0, 0), 
              mar=c(0, 0, 0, 0), new=TRUE)
  on.exit(par(opar))
  
  xx <- yy <- cols <- NULL
  ytick <- scale <- text.col <- text.col <- bg <- endLabels <- NULL
  labside <- 'right'
  
  for(k in 1:length(clist))assign( names(clist)[k], clist[[k]] ) 
  
  xx <- -1.08 + .1*(kplot - 1) + c(.1, .2)
  yy <- c(-1, -.75)
  
  nn <- length(cols)
  ys <- seq(yy[1],yy[2],length=nn)
  if(is.null(scale))scale <- ys
  
  for(j in 1:(length(scale)-1)){
    
    rect(xx[1],ys[j],xx[2],ys[j+1],col=cols[j],border=NA)
  }
  if(!is.null(bg))rect(xx[1],yy[1],xx[2],yy[2],border=bg,lwd=1)
  if(!is.null(ytick)){
    
    ys <- diff(yy)/diff(range(ytick))*ytick
    yt <- ys - min(ys) + yy[1]
    
    for(j in 1:length(yt)){
      lines(xx,yt[c(j,j)])
    }
  }
  if(!is.null(endLabels)){ 
    if(labside == 'right')text(diff(xx)+c(xx[2],xx[2]),yy,endLabels)
    if(labside == 'left')text(c(xx[1],xx[1]),yy,endLabels,pos=2)
  }
}

values2grid <- function(x, y, z, nx=NULL, ny=NULL, dx=NULL, dy=NULL,
                        ksearch = 4, MATFORMAT=T){
  
  xl <- range(x)
  yl <- range(y)
  
  xs <- seq(xl[1],xl[2],length=nx)
  ys <- seq(yl[1],yl[2],length=ny)
  
  grid <- as.matrix(expand.grid(xs,ys))
  
  tmp <- nn2(cbind(x,y),grid,k=ksearch)
  nn  <- tmp[[1]]
  wt  <- tmp[[2]]
  mn  <- min(wt[wt > 0])
  wt  <- 1/(wt + mn)
  
  zz  <- matrix(z[nn],nrow(nn),ncol(nn))
  zz  <- rowSums(zz*wt)/rowSums(wt)
  
  if(!MATFORMAT)return(  cbind(grid,zz) )
  
  zmat <- matrix(NA, nx, ny)
  ix  <- match(grid[,1],xs)
  iy  <- match(grid[,2],ys)
  zmat[ cbind(ix,iy) ] <- zz
  
  rownames(zmat) <- xs
  colnames(zmat) <- ys
  
  list( x = xs, y = ys, z = zmat )
}



.pmap <- function(specNames=NULL, xytree=NULL, plot=NULL, 
                  year=NULL, seedPredGrid=NULL,
                  treePredGrid=NULL, xlim, ylim, treeScale, trapScale,
                  sCol = 'blue', add=F){
  
  #multiple species for single plot-year
  
  ADD <- F
  if(add)ADD <- T
  
  pnames   <- paste(specNames,'_meanM2',sep='')
  nspec    <- length(specNames)
  predCols <- c('x','y',pnames)
  
  if(length(sCol) == 1 & nspec > 1)sCol <- rep(sCol, nspec)
  
  wtree <- which( treePredGrid$year %in% year & 
                    treePredGrid$plot %in% plot )
  tmat <- treePredGrid[wtree,]
  tmat$x <- xytree$x[ match(tmat$treeID, xytree$treeID) ]
  tmat$y <- xytree$y[ match(tmat$treeID, xytree$treeID) ]
  
  fec  <- tmat[,'fecEstMu']
  matr <- tmat[,'matrPred']
  
  wtrap <- which( seedPredGrid$year %in% year & 
                    seedPredGrid$plot %in% plot )
  
  if(length(wtrap) == 0){
    if(add)plot(NULL, xlim = xlim, ylim = ylim, axes= F, xlab='',
                ylab='')
    return()
  }
    
  smat   <- as.matrix( seedPredGrid[wtrap,predCols] )
  
  dens <- nrow(smat)/diff(xlim)/diff(ylim)
  nx <- ny <- ceiling(1/dens)
  
  fecMax  <- max(fec,na.rm=T)
  seedMax <- max(smat[,pnames],na.rm=T)
  
 # xseq <- sort(unique(smat[,'x']))
 # yseq <- sort(unique(smat[,'y']))
 # nx   <- length(xseq)
 # ny   <- length(yseq)
 # zmat <- matrix(NA, nx, ny)
 # ix  <- match(smat[,'x'],xseq)
 # iy  <- match(smat[,'y'],yseq)
  
  rr  <- apply( rbind(xlim, ylim), 2, range)
  sc  <- max( apply(rr,1,diff) )/20
  
  wspec <- which(colSums(smat[,pnames, drop=F]) > 0)
  
  sn <- 4
  if(seedMax > 1)sn <- 3
  if(seedMax > 10)sn <- 2
  
  q <- seq(0, 1, length=10)^.3
  q[1] <- .3
  
  levels <- signif( quantile(smat[,pnames], q ) , sn)
  levels <- c(levels, signif(max(seedMax)*1.3, sn) )
  levels <- sort(unique(levels))
  colorLevels <- seq(.001, .99, length=length(levels))^3
  
  for(k in wspec){
    
    
    tmp <- values2grid(x=smat[,'x'],y=smat[,'y'],z=smat[,pnames[k]],nx=nx,ny=ny)
    
    xseq <- tmp$x
    yseq <- tmp$y
    zmat <- tmp$z
    
   # zmat[ cbind(ix, iy) ] <- smat[,pnames[k]]
    
    col <- .getColor(sCol[specNames[k]],colorLevels)
    contour(xseq, yseq, zmat, levels=levels, add=add, 
            col=col, labcex = 1, frame.plot=F,
            drawlabels=F, axes = F)
    .filled.contour(xseq, yseq, zmat, levels=levels, col=col)
    if(!add){
      for(m in 1:4)axis(m, labels=F, tck=.03)
    }
    if(!ADD)add <- T
    
  }
  
  for(k in wspec){
    wk <- which( tmat$species == specNames[k] )
    if(length(wk) == 0) next
    
    z   <- 5*sc*fec[wk]*treeScale
    z[z > 0] <- z[z > 0]/fecMax
    .mapSpec(x = tmat[wk,'x'],y = tmat[wk,'y'], z*1.4,  add = add, 
             mapx = xlim, mapy = ylim,
             colVec='white', fill = 'white')
    .mapSpec(x = tmat[wk,'x'],y = tmat[wk,'y'], z,  add = T, 
             mapx = xlim, mapy = ylim,
             colVec=sCol[specNames[k]], fill = sCol[specNames[k]])
  }
  invisible( list(species = wspec, colorLevels = colorLevels, add = add) )
}

.colPaste <- function(...){
  
  # ... are vectors, all same length to be combined row-wise
  
  mm <- list(...)
  nn <- as.character(mm[[1]])
  for(i in 2:length(mm)){
    mi <- as.character(mm[[i]])
    nn <- apply( cbind(nn,mi),1,paste0,collapse='-')
  }
  nn
}
 
.updateCovariance <- function(SS, priorSS, n, df){
  
  SI   <- solveRcpp(SS + df*priorSS)
  sinv <- .rwish(n + df, SI)
  solveRcpp(sinv)
}

.updateAlphaRand <- function(yy, xfecA, zz, sg, reIndexA, reGroups,
                             xrandCols, Arand, priorVA, dfA, minmax = 2){
  
  ONEA <- F
  if(length(Arand) == 1)ONEA <- T
  
  wg <- xfecA[zz == 1, xrandCols, drop=F]
  alphaRand <- randEffectRcpp(reIndexA[zz == 1], reGroups, 
                              wg, yy[zz == 1], sg, solve(Arand))
  alphaRand[alphaRand < -minmax] <- -minmax
  alphaRand[alphaRand > minmax]  <- minmax
  rz    <- which( rowSums(alphaRand) != 0)
  
  if(length(rz) < 2)return(list(alphaRand = alphaRand, Arand = Arand))
  
  if(ONEA){
    Arand <- matrix(1/rgamma(1,1 + length(rz)/2, 1 + 1/2*sum(alphaRand[rz]^2)),1)
  }else{
    AA    <- crossprod(alphaRand[rz,])
    Arand <- .updateCovariance(AA, priorVA, length(rz), dfA)
  }
  
  list(alphaRand = alphaRand, Arand = Arand)
}

.rwish <- function(df,SS){
  z  <- matrix(rnorm(df*nrow(SS)),df,nrow(SS))%*%chol(SS)
  crossprod(z)
}

.riwish <- function(df,S){
  solveRcpp(.rwish(df,solveRcpp(S)))
}

print.mastif <- function(x, ...){
  
  rMu <- rSe <- usigma <- betaYrMu <- betaYrSe <- NULL
  
  cat("\nDIC:\n")
  print( round(x$fit$DIC,0) )
  
  cat("\nFecundity coefficients:\n")
  print( signif(x$parameters$betaFec, 3) )
  
  cat("\nMaturation coefficients:\n")
  print( signif(x$parameters$betaRep, 3) )
  
  if('betaYrMu' %in% names(x$parameters)){
    cat("\nYear effects:\n")
    print( signif(betaYrMu, 3) )
    print( signif(betaYrSe, 3) )
  }
  
  if('rgibbs' %in% names(x$chains)){
    cat("\nSpecies to seed type matrix R:\n")
    print( signif(rMu, 3) )
    print( signif(rSe, 3) )
  }
  
  cat("\nSigma, RMSPE:\n")
  usigma <- x$parameters$sigma
  print( signif(usigma, 3) )
  
  cat("\nDispersal parameter u (m^2):\n")
  upars <- x$parameters$upars
  print( signif(upars, 3) )
  
  cat("\nKernel mean distance (m):\n")
  dpars <- x$parameters$dpars
  print( signif(dpars, 3) )
  
}

summary.mastif <- function(object,...){ 
  
  betaFec   <- object$parameters$betaFec
  betaRep   <- object$parameters$betaRep
  rownames(betaFec) <- .replaceString(rownames(betaFec),'species','')
  rownames(betaRep) <- .replaceString(rownames(betaRep),'species','')
  seedNames <- object$inputs$seedNames
  specNames <- object$inputs$specNames
  tdata     <- object$inputs$treeData
  sdata     <- object$inputs$seedData
  plots     <- levels(tdata$plot)
  ntype     <- length(seedNames)
  nseed     <- nrow(sdata)
  
  out <- list()

  AR <- YR <- RANDOM <- SAMPR <- F
  
  RMSPE <- object$fit$RMSPE
  DIC   <- object$fit$DIC
  
  if("arList" %in% names(object$data))AR <- T
  
  #data summary
  wd <- which(!duplicated(tdata$treeID))
  trees <- table(tdata$species[wd],tdata$plot[wd])[,plots,drop=F]
  rownames(trees) <- paste('trees',rownames(trees))
  
  wd <- which(!duplicated(sdata$trapID))
  traps <- table(sdata$plot[wd])[plots]
  ntr   <- names(traps)
  traps <- matrix(traps,1)
  colnames(traps) <- ntr
  rownames(traps) <- 'traps'
  
  treeYears <- table(tdata$species,tdata$plot)[,plots,drop=F]
  rownames(treeYears) <- paste('tree-yrs',rownames(treeYears))
  trapYears <- table(sdata$plot)[plots]
  dataTab <- rbind(trees,treeYears,traps,trapYears)
  
  seedCount <- as.matrix(sdata[,seedNames])
  colnames(seedCount) <- seedNames
  
  wm <- match(as.character(sdata$plot),plots)
  ii <- rep(wm,ntype)
  jj <- rep(1:ntype,each=nseed)
  totalSeed <- t( .myBy(as.vector(seedCount),ii,jj,fun='sum') )
  colnames(totalSeed) <- plots
  rownames(totalSeed) <- paste('seeds',seedNames)
  dataTab <- rbind(dataTab,totalSeed)
  
  total   <- rowSums(dataTab)
  dataTab <- cbind(dataTab,total)
  
  cat('\nData summary:\n')
  print(dataTab)
  
  out$data <- dataTab
  
  cat("\nFecundity parameters:\n")
  print( signif(betaFec, 3) )
  
  out$betaFec <- signif(betaFec, 3)
  
  cat("\nMaturation parameters:\n")
  print( signif(betaRep, 3) )
  
  out$betaRep <- signif(betaRep, 3)
  
  if('betaYrMu' %in% names(object$parameters)){
    YR <- T
    byr <- t( rbind(object$parameters$betaYrMu, object$parameters$betaYrSe) )
    colnames(byr) <- c('estimate','se')
    cat("\nYear effects, only mature individuals:\n")
    print( signif(byr, 3) )
    
    out$betaYr <- signif(byr, 3)
  }
  
  if('betaYrRand' %in% names(object$parameters)){
    cat("\nYear effects, random group means:\n")
    print( object$parameters$betaYrRand )
    cat("\nYear effects, standard deviation between groups:\n")
    print( object$parameters$betaYrSe )
    
    out$betaYrRand <- object$parameters$betaYrRand
    out$betaYrGroupSE <- object$parameters$betaYrSe
  }
  
  
  if('aMu' %in% names(object$parameters)){
    RANDOM <- T
    amu <- diag( object$parameters$aMu)
    ase <- diag( object$parameters$aSe )
    
    arand <- signif(cbind(amu,ase),3)
    rownames(arand) <- .replaceString(rownames(arand),'species','')
    colnames(arand) <- c('estimate','std err')
    
    cat("\nDiagonal elements of random effects covariance matrix (aMu):\n")
    print( arand )
    
    out$arand <- arand
  }
  
  if('rgibbs' %in% names(object$chains)){
    SAMPR <- T
    cat("\nSpecies to seed type matrix R:\n")
    print( signif(object$parameters$rMu, 3) )
  
    cat("\nStandard errors for R:\n")
    print( signif(object$parameters$rSe, 3) )
    
    out$rMu <- signif(object$parameters$rMu, 3)
    out$rSe <- signif(object$parameters$rSe, 3)
  }
  
  cat("\nSigma, RMSPE:\n")
  usigma <- object$parameters$sigma
  print( signif(usigma, 3) )
  
  out$sigma <- signif(usigma, 3)
  
  upars <- object$parameters$upars
  dpars <- object$parameters$dpars
  utab  <- signif(rbind(upars,dpars),3)
  group <- rownames(utab)
  
  parameter <- c(rep('u (m^2)',nrow(upars)), 
                 rep('d (m)',nrow(dpars))) 
  parameter[nrow(upars)] <- 'u (m^2)'
  parameter[nrow(upars)*2] <- 'd (m)'
  rownames(utab) <- NULL
  utab <- data.frame( parameter = as.factor(parameter), utab) 
  if(!is.null(group))utab <- cbind(group,utab)
  rownames(utab) <- NULL
  
  cat("\nKernel estimates:\n")
  print(utab)
  
  out$u <- utab
  
  if(AR){
    out$eigenMu <- signif(object$parameters$eigenMu, 3)
    out$eigenSe <- signif(object$parameters$eigenSe, 3)
  }
  
  
  ntreeYr  <- nrow(object$data$setupData$tdata)
  ntrapYr  <- nrow(object$data$setupData$sdata)
  years    <- object$data$setupData$years
  ntree    <- ncol(object$data$setupData$distall)
  ntrap    <- nrow(object$data$setupData$distall)
  plots    <- object$data$setupData$plots
  nplot    <- length(plots)
  nyr      <- length(years)
  
  
  words <- paste("Sample contains ", ntreeYr, " tree-years and ",
                 ntrapYr, " trap-years on ", ntree, " individuals and ", 
                 ntrap, " seed traps on ", nplot, " plots, over ", nyr, 
                 " years. The RMSPE is ",
                 signif(object$fit$RMSPE,3),
                 ", and the DIC is ",round(object$fit$DIC),".", sep="")
  cat("\n",words)
  
  if(RANDOM){
    
    morewords <- paste(" Random effects were fitted on",
                            length(object$data$setupRandom$rnGroups), 
                            " individuals\n")
    cat("\n",morewords)
  }
  
  fit <- c(round(object$fit$DIC), object$fit$RMSPE)
  names(fit) <- c('DIC','RMSPE')
  
  out$fit   <- fit
  out$words <- words
 
  class(out) <- "summary.mastif"
  invisible(out) 
}

.getMuCv <- function(seed, tran){
  
  mu  <- .myBy(as.vector(seed),tran[,1],tran[,2],fun='mean')
  s1  <- .myBy(as.vector(seed)*0 + 1,tran[,1],tran[,2],fun='sum')
  v1  <- .myBy(as.vector(seed),tran[,1],tran[,2],fun='sum')
  v2  <- .myBy(as.vector(seed^2),tran[,1],tran[,2],fun='sum')
  vs  <- v2/s1 - mu^2
  cv  <- sqrt(vs)/mu
  list(mu = mu, sd = sqrt(vs), cv = cv)
}

.getTran <- function(grid, new){
  newx <- seq(min(grid[,'x']),max(grid[,'x']),by=new)
  newy <- seq(min(grid[,'y']),max(grid[,'y']),by=new)
  tx   <- findInterval(grid[,'x'],newx)
  ty   <- findInterval(grid[,'y'],newy)
  cbind(tx,ty)
}

mastSim <- function(sim){       # setup and simulate fecundity data
  
  if(!'seedNames' %in% names(sim))
    stop('sim must include seedNames')
  if(!'specNames' %in% names(sim))
    stop('sim must include specNames')
  
  .dsim( sim ) 
}
 
.dsim <- function(sim){
  
  nyr <- 5; ntree  <-  10; ntrap <- 20; plotWide  <-  100; nplot  <-  3
  meanDist  <-  20; Q  <-  2
  minDiam <- 10
  maxDiam <- 40
  yearEffect <- NULL
  facLevels <- character(0)
  seedNames <- specNames <- c('piceaGlauca','piceaMariana','piceaUNKN')
  SPECS <- SAMPR <- AR <- F
  maxF <- 1e+7
  
  for(k in 1:length(sim))assign( names(sim)[k], sim[[k]] )
  S         <- length(specNames)
  ntreePlot <- rpois(nplot, ntree) + 1 # trees per plot
  ntrapPlot <- rpois(nplot, ntrap) + 4        # traps per plot
  nyrPlot   <- rpois(nplot, nyr) + 1
  plotNames <- paste('p',1:nplot,sep='')
  yearNames <- 2017 - c(max(nyrPlot):1)
  if(length(specNames) > 1)SPECS <- T
  if(length(seedNames) > 1)SAMPR <- T
  nyr <- max(nyrPlot)
  
  upar <- (2*meanDist/pi)^2
  
  year <- plot <- tree <- trap <- yrsd <- plsd <- numeric(0)
  
  for(j in 1:nplot){
    
    tree <- c( tree, rep(1:ntreePlot[j], each=nyrPlot[j]) )
    year <- c( year, rep(1:nyrPlot[j], ntreePlot[j]) ) 
    plot <- c( plot, rep(j, (nyrPlot[j]*ntreePlot[j]) ) )
    
    trap <- c( trap, rep(1:ntrapPlot[j], each=nyrPlot[j]) )
    yrsd <- c( yrsd, rep(1:nyrPlot[j], ntrapPlot[j]) ) 
    plsd <- c( plsd, rep(j, (nyrPlot[j]*ntrapPlot[j]) ) )
  }
  
  tree <- data.frame(plot = plotNames[plot], year = yearNames[year], tree )
  tree <- tree[order(tree$plot, tree$tree, tree$year),]
  
  tree$treeID  <- columnPaste( tree$plot, tree$tree )
  treeID       <- as.character( tree$treeID ) 
  
  id           <- unique(as.character(tree$treeID))
  species      <- sample(specNames, length(id), replace=T) 
  tree$species <- factor(species[match(tree$treeID,id)], levels = specNames)
  
  tree$dcol <- match(as.character(tree$treeID),id)
  
  tree$plotYr <- columnPaste(tree$plot, tree$year )
  plotyr      <- unique(tree$plotYr)
  tree$plotyr <- match(tree$plotYr, plotyr )
  
  years <- sort(unique(tree$year))
  
  trap <- data.frame(plot = plotNames[plsd], year = yearNames[yrsd], trap )
  trap <- trap[order(trap$plot, trap$trap, trap$year),]
  
  trap$trapID  <- columnPaste( trap$plot, trap$trap )
  trapid       <- as.character( trap$trapID )
  drow         <- unique(trapid)
  trap$drow    <- match(trapid,drow)

  trap$plotYr <- columnPaste(trap$plot, trap$year )
  plotyr      <- unique(trap$plotYr)
  trap$plotyr <- match(trap$plotYr, plotyr )
  
  n <- nrow(tree)
  
  xfec <- round( matrix( .tnorm(n*(Q-1), 5, 50, 35, 5), n, (Q-1) ), 3)
  xnames <- paste('x',1:(Q-1),sep='')
  xnames[1] <- 'diam'
  colnames(xfec) <-  xnames
 
  xdata   <- data.frame( species = tree$species, xfec)
  
  formulaFec <- as.formula( '~ I(log(diam))'  )
  formulaRep <- as.formula( '~ I(log(diam))' )
  
  if(SPECS){
    formulaFec <- as.formula( '~ species*I(log(diam))'  ) 
    formulaRep <- as.formula( ' ~ species*I(log(diam))' )
  }
  xfec    <- model.matrix(formulaFec, xdata)
  Qf      <- ncol(xfec)
  xrep    <- model.matrix(formulaRep, xdata)
  Qr      <- ncol(xrep)
  
  if(!SPECS){
    ss     <- paste('species',specNames,sep='')
    xnames <- paste(ss,colnames(xfec),sep=':')
    xnames[1] <- .replaceString(xnames[1],':(Intercept)','')
    colnames(xfec) <- xnames
  }
    
  xytree <- xytrap <- distall <- numeric(0)
  
  for(j in 1:nplot){
    
    xy1 <- matrix( runif(2*ntreePlot[j],0,plotWide), ntreePlot[j],2)  
    xy2 <- matrix( runif(2*ntrapPlot[j],0,plotWide), ntrapPlot[j],2)
    xy1 <- round(xy1,1)
    xy2 <- round(xy2,1)
    
    rownames(xy1) <- columnPaste(rep(plotNames[j],ntreePlot[j]), 
                                 c(1:ntreePlot[j]), '-')
    rownames(xy2) <- columnPaste(rep(plotNames[j],ntrapPlot[j]), 
                                 c(1:ntrapPlot[j]), '-')
    xytree  <- rbind(xytree, xy1)
    xytrap  <- rbind(xytrap, xy2)
  }
  
  xdata <- xdata[,!colnames(xdata) %in% colnames(tree),drop=F]
  
  treeData <- cbind(tree,xdata)
  count    <- matrix(0, nrow(trap), length(seedNames))
  colnames(count) <- seedNames
  seedData <- data.frame(trap, count)
  
  colnames(xytree) <- colnames(xytrap) <- c('x','y')
  xy <- columnSplit(rownames(xytree),'-')
  colnames(xy) <- c('plot','tree')
  xytree <- data.frame( xy, xytree )
  
  xy <- columnSplit(rownames(xytrap),'-')
  colnames(xy) <- c('plot','trap')
  xytrap <- data.frame( xy, xytrap )
  
 # formulaFec <- .specFormula(formulaFec)
 # formulaRep <- .specFormula(formulaRep)
  
  seedData$active <- 1
  seedData$area   <- 1
  
  ntree    <- nrow(xytree)
  nyr      <- max(nyrPlot)
  dmat     <- matrix(runif(ntree*nyr, .2, .5), ntree, nyr )
  
  dmat[,1] <- .tnorm(ntree, 1, 60, 20, 10) 
  small    <- sample(ntree,round(ntree/2))
  dmat[small,1] <- .tnorm(length(small), 1, 20, 8, 20)
  if(nyr > 1)dmat <- round(t( apply(dmat,1,cumsum) ), 2)
  
  tyindex  <- cbind(treeData$dcol, match(treeData$year,years))
  treeData$diam <- dmat[ tyindex ]
  
  tmp <- .setupData(formulaFec, formulaRep, treeData, seedData, 
                    facLevels, xytree, xytrap, specNames, seedNames, 
                    AR=F, YR=F, yearEffect, minDiam, maxDiam)
  treeData  <- tmp$tdata
  seedData  <- tmp$sdata
  distall   <- tmp$distall
  xytree    <- tmp$xytree
  xytrap    <- tmp$xytrap
  plotNames <- tmp$plotNames
  plots     <- tmp$plots
  years     <- tmp$years
  xfec      <- tmp$xfec
  xrep      <- tmp$xrep
  nseed     <- nrow(seedData)
  scode     <- tmp$scode
  nplot     <- length(plotNames)
  n         <- nrow(xfec)
  ntobs     <- table(treeData$plot) 
  nsobs     <- table(seedData$plot)
  ttab      <- table(treeData$plot, treeData$year)
  wtab      <- which(ttab > 0, arr.ind=T) 
  xfecMiss <- tmp$xfecMiss
  xrepMiss <- tmp$xrepMiss
  xfecCols <- tmp$xfecCols
  xrepCols <- tmp$xrepCols
  ntree    <- nrow(xytree)
  ntrap    <- nrow(xytrap)
  xfecU    <- tmp$xfecU; xfecT <- tmp$xfecT
  xrepU    <- tmp$xrepU; xrepT <- tmp$xrepT
  xfecs2u  <- tmp$xfecs2u
  xfecu2s  <- tmp$xfecu2s
  xreps2u  <- tmp$xreps2u
  xrepu2s  <- tmp$xrepu2s
  
  nspec    <- length(specNames)
  

  betaFec <- matrix( runif(ncol(xfec), -10, -8), ncol=1)
  rownames(betaFec) <- colnames(xfec)
  wcol  <- grep('diam',colnames(xfec)) 
  slope <- (.9*log(maxF) + -betaFec[-wcol])/max(log(dmat))
  betaFec[wcol,1] <- slope
  
  fec <- xfec%*%betaFec
  
  lo  <- hi <- 0*fec
  hi[fec > 0] <- log(maxF)
  lo[fec <= 0] <- -log(maxF)/2
  fec <- .tnorm(nrow(xfec), lo, hi, xfec%*%betaFec, .01)

  d0 <- -betaFec[-wcol,]/betaFec[wcol,]   # zero
  z  <- fec*0
  
  betaRep <- matrix(0,ncol(xrep),1)
  
  
  int <- -8
  slope <- 3
  
  betaRep[wcol,1] <- slope
  betaRep[-wcol,1] <- int
  
  w <- rnorm(length(fec), xrep%*%betaRep)
  z <- w*0
  z[ w > 0 ] <- 1
  betaRep <- solve(crossprod(xrep))%*%crossprod(xrep,w)
  
  lo  <- hi <- 0*fec
  hi[z > 0] <- log(maxF)
  lo[z <= 0] <- -log(maxF)/2
  fec <- .tnorm(nrow(xfec), lo, hi, xfec%*%betaFec, .01)
  betaFec <- solve( crossprod(xfec) )%*%crossprod(xfec,fec)
  
  zmat <- matrix(0, ntree, nyr )
  zmat[ tyindex ] <- z
  if(nyr > 1)zmat <- round(t( apply(zmat,1,cumsum) ), 2)
  zmat[zmat > 1] <- 1
  
  z <- zmat[ tyindex ]
  lo <- hi <- z*0
  lo[z == 0] <- -Inf
  hi[z == 1] <- Inf
  w <- .tnorm(length(w),lo,hi,xrep%*%betaRep,1)
  betaRep <- solve(crossprod(xrep))%*%crossprod(xrep,w)
  
  ztrue <- z
  zmat[ sample(n,n/20) ] <- NA
  
  treeData$repr <- zmat[ tyindex ]

  tmp <- .getRepro(zmat, dmat, minDiam, maxDiam)
  last0first1 <- tmp$last0first1
  
  q <- which(ztrue == 1)
  betaFec <- solve( crossprod(xfec[q,]) )%*%crossprod(xfec[q,],fec[q])
  hi[ztrue > 0] <- log(maxF)
  lo[ztrue <= 0] <- -log(maxF)/2
  fec <- .tnorm(nrow(xfec), lo, hi, xfec%*%betaFec, .01)
  betaFec <- solve( crossprod(xfec[q,]) )%*%crossprod(xfec[q,],fec[q])
 
  attr(distall,'group') <- rep(1, ncol(distall) )
  
  seedData$active <- seedData$area <- 1
  
  fec <- exp(fec)
  
  # R matrix
  
  fill <- 0
  if(length(seedNames) == 1)fill <- 1
  
  R <- matrix(fill, length(plots)*nspec, length(seedNames))
  colnames(R) <- seedNames
  rr <- as.vector( outer(specNames, plots, paste, sep='-') )
  rownames(R) <- rr
    
  wun <- grep('UNKN', seedNames)
  if(length(wun) > 0){
    kk <- c(1:length(seedNames))[-wun]
    for(k in kk){
      wsp <- grep(seedNames[k],rownames(R))
      R[wsp,seedNames[k]] <- 1
    }
    R[,wun] <- 2
    R <- sweep(R, 1, rowSums(R), '/')
  }
  tmp <- columnSplit(rownames(R),'-')
  attr(R, 'species') <- tmp[,1]
  attr(R, 'plot')    <- tmp[,2]
  
  treeData$specPlot <- columnPaste(treeData$species,treeData$plot)
   
  lambda <- .getLambda(treeData, seedData, 
                       seedData$active*seedData$area, upar, fec, z, R,
                       SAMPR, distall, years, PERAREA=F) 
  lambda <- lambda + 1e-8
  ss     <- matrix(rpois(length(lambda),lambda),nrow(lambda), ncol(lambda))
  seedData[,seedNames] <-   ss
  seedData$active <- 1
  seedData$area   <- 1
  
  seedData$active <- 1
  seedData$area   <- 1
  
  stab <- with( seedData, table(plot, year) )
  ttab <- with( treeData, table(plot, year) )
  sc   <- colSums(stab)
  stab <- stab[,sc > 0, drop=F]
  ttab <- ttab[,sc > 0, drop=F]
  
  form <- as.character( formulaFec )
  form <- .replaceString( form, 'species *', '')
  formulaFec <- as.formula(form)
  
  form <- as.character( formulaRep )
  form <- .replaceString( form, 'species *', '')
  formulaRep <- as.formula(form)
  
  
  bfecSave <- xfecs2u%*%betaFec      # unstandardize
  brepSave <- xreps2u%*%betaRep
  
  xx <- cbind(1, treeData$diam)
  
  trueValues <- list(fec = fec, repr = ztrue, betaFec = bfecSave, 
                     betaRep = brepSave,upar = upar, R = R)
  
  
  treeData <- treeData[,c('plot','tree','year','species','diam','repr','repMu')]
  seedData <- seedData[,c('plot','trap','year','area','active',seedNames)]
  xytree   <- xytree[,c('plot','tree','x','y')]
  xytrap   <- xytrap[,c('plot','trap','x','y')]
  
  out <- list(trueValues = trueValues, treeData = treeData, seedData = seedData, 
       distall = distall, xytree = xytree, xytrap = xytrap, formulaFec = formulaFec,
       formulaRep = formulaRep, plots = plots, years = years,
       sim = sim, xfec = xfec, xrep = xfec, seedNames = seedNames, 
       specNames = specNames, R = R,
       sample = list(seeds = stab, trees = ttab))
  orr <- order(names(out))
  out <- out[orr]
  out
}
     
.seedFormat <- function(sfile, lfile, trapFile = NULL, seedNames = NULL, 
                        genusName = NULL, omitNames = NULL, plot, trapID ) {
  
  if(is.null(trapFile))trapFile <- "/Users/jimclark/makeMast/dataFiles/seedTrapArea.txt"
  
  
  midCD <- c( 273890.6, 3938623.3 )  #plot center for (x,y) at GSNP_CD
  
  loc  <- read.table(lfile, header=T)
  if(!'x' %in% colnames(loc)){
    xy <- loc[,c('UTMx','UTMy')]
    xy <- round(sweep(xy,2,colMeans(xy),'-'),1)
    loc$x <- xy[,1]
    loc$y <- xy[,2]
  }
  loc  <- loc[is.finite(loc[,'UTMx']) & is.finite(loc[,'UTMy']),]
  pcol <- rep(plot, nrow(loc))
  id   <- apply( cbind(plot, loc[,trapID]), 1, paste0, collapse='-')
  loc  <- data.frame(trapID = id, trap = loc[,trapID], 
                     loc[,!colnames(loc) == trapID])
  loc$plot <- pcol
  loc$UTMx <- loc[,'UTMx']
  loc$UTMy <- loc[,'UTMy']
  
  if(plot == "GSNP_CD"){
    loc$x <- loc$UTMx - midCD[1]
    loc$y <- loc$UTMy - midCD[2]
  }
  
  counts <- read.table(sfile, header=T)
  counts[is.na(counts$month),'month'] <- 3
  counts[counts[,'month'] < 9,'year'] <- counts[counts[,'month'] < 9,'year'] - 1
  
  if(is.null(genusName)){
    sj <- counts[,colnames(counts) %in% seedNames, drop=F]
  }else{
    sj <- counts[,grep(genusName, colnames(counts)), drop=F]
    seedNames <- colnames(sj)
  }
  if(!is.null(omitNames)){
    sj <- sj[,!colnames(sj) %in% omitNames, drop=F]
    seedNames <- colnames(sj)
  }
  
  if(length(sj) == 0){
    sj <- matrix(0,nrow(counts),1)
    seedNames <- colnames(sj) <- paste(genusName,'UNKN',sep='')
  }
  
  yr <- sort(unique(counts[,'year']))
  tn <- sort(unique(counts[,trapID]))
  jj <- match(counts[,'year'],yr)
  ii <- match(counts[,trapID],tn)
  smat <- matrix(0, length(tn), length(yr))
  
  seedj <- numeric(0)
  
  for(k in 1:ncol(sj)){
    
    # seed counts
    ck <- sj[,k]
    wk <- which(is.finite(ck))
    ck <- ck[wk]
    ik <- ii[wk]
    jk <- jj[wk]
    ky <- .myBy(ck, ik, jk, summat=smat, fun='sum')
    sy <- .myBy(ck*0 + 1, ik, jk, summat=smat, fun='sum')
    colnames(ky) <- yr
    rownames(ky) <- tn
    seedj <- cbind(seedj, as.vector(ky))
    if(k == 1)active <- sy/matrix( apply(sy,2,max), nrow(sy), ncol(sy), byrow=T)
  }
  colnames(seedj) <- colnames(sj)
  seed <- matrix(0, nrow(seedj), length(seedNames))
  colnames(seed) <- seedNames
  
  active <- as.vector(active)
  area <- read.table(trapFile,header=T)[,c('plot','seedArea')]
  area <- area[area$plot == as.character(plot[1]),'seedArea']
  
  seed[,colnames(seedj)] <- seedj
  
  year <- rep(yr, each = length(tn))
  trap <- rep(tn, length(yr) )
  plot <- rep(plot, length(trap))
  tr <- apply(cbind(plot, trap), 1, paste0, collapse='-')
  sd   <- data.frame( plot=plot, trapID=tr, trap=trap, year=year,
                      active=active, area = area) 
  seed <- cbind(sd, seed)
  rownames(seed) <- 
    apply( cbind(plot, trap, year), 1, paste0, collapse='-')
  
  seed <- seed[seed$trapID %in% loc$trapID,]
  
  
  list(counts = seed, xy = loc, active = active, seedNames = seedNames )
}

.getRepro <- function(rmat, dmat, minDiam, maxDiam){
  
  nc <- ncol(rmat)
  nr <- nrow(rmat)
  
  zknown <- rmat
  
  zeros <- ones <- ind <- matrix(1:nc, nr, nc, byrow=T)
  ones  <- ones*rmat
  ones[ones == 0] <- NA
  first1  <- suppressWarnings(apply(ones, 1, min, na.rm=T)) # this is first 1
  ones  <- matrix( first1, nr, nc)
  ones[ones <= ind] <- 1
  ones[ones > ind] <- 0
  first1[first1 == Inf] <- ncol(rmat) + 1
  all1 <- which(rowSums(ones) == ncol(ones))
  
  zeros <- zeros*(1 - rmat)
  zeros[ ones == 1 ] <- 0
  last0  <- suppressWarnings(apply(zeros, 1, max, na.rm=T)) # this is last 0
  zeros  <- matrix( last0, nr, nc)
  zeros[ zeros == 0] <- NA
  zeros[ zeros == -Inf ] <- NA
  zeros[zeros < ind] <- 0
  zeros[zeros >= ind] <- 1
  rmat[ ones == 1 ] <- 1
  rmat[ zeros == 1 ] <- 0
  last0[last0 < 0] <- 0
  all0 <- which(rowSums(zeros) == ncol(zeros))
  
  last0first1 <- cbind(last0, first1)
  
  mm <- matrix(0, length(last0), 2)
  colnames(mm) <- c('all0','all1')
  mm[all0,1] <- 1
  mm[all1,2] <- 1
  last0first1 <- cbind(last0first1,mm)
  
  rownames(last0first1) <- paste('drow',1:nrow(last0first1),sep='_')
  
  mmat <- matrix(1:nc, nr, nc, byrow=T)
  
  mid <- (maxDiam + minDiam)/2
  ss  <- (mid - minDiam)/2
  pd  <- pnorm(dmat,mid,ss)
  pd  <- rowMeans(pd, na.rm=T)
  rr  <- rbinom(length(pd),1,pd)
  rr[all1] <- 1
  rr[all0] <- 0
  
  wm   <- which(rr == 1)
  wf   <- which(first1 <= nc)
  ones <- zeros <- mmat*0
  ones[ cbind(1:nr, first1)[wf,] ] <- 1
  ones  <- t(apply(ones,1,cumsum))
  
  wl  <- which(last0 > 0)
  zeros[ cbind(1:nr, last0)[wl,] ] <- 1
  zeros[,nc] <- 0
  
  zeros <- t(apply(zeros,1,rev))
  zeros <- t(apply(zeros,1,cumsum))
  zeros <- t(apply(zeros,1,rev))
  zeros[all0,] <- 1
  
  mmat[rr == 1,] <- 1
  mmat[ones == 1] <- 1
  mmat[zeros > 0] <- 0
  
  mmat[mmat == 0] <- 1000
  
  myr <- apply(mmat, 1, which.min)
  myr[all0] <- nc + 1
  
  mmat[mmat == 1000] <- 0
  mmat[mmat > 0] <- 1
  
  
  list(rmat = mmat, matYr = myr, last0first1 = last0first1)
}

.treeFormat <- function( tfile, specNames = NULL, genusName = NULL, 
                         changeNames = NULL, plot, years, yrCols  ){
  
  
  midCD <- c( 273890.6, 3938623.3 )  #plot center for (x,y) at GSNP_CD
  
  region <- unlist(strsplit(plot,'_'))[1]
  trees  <- read.table( tfile, header=T)
  
  if(plot == "GSNP_CD"){
    trees$x <- trees$UTMx - midCD[1]
    trees$y <- trees$UTMy - midCD[2]
  }

  
  if( is.null(genusName) ){
    wj     <- which(trees$species %in% specNames & is.finite(trees$x) &
                      is.finite(trees$y))
  }else{
    ww <- grep(genusName,as.character(trees$species) )
    wj <- ww[which(is.finite(trees$x[ww]) & is.finite(trees$y[ww]))]
    specNames <- sort(unique(trees$species[wj]))
  }
  if(length(wj) == 0)return( list(yrDat = numeric(0)) )
  
  trees  <- trees[drop=F,wj,]
  wchange <- which( as.character(trees$species) %in% changeNames[,1])
  if(length(wchange) > 0){
    ts <- as.character(trees$species)
    wm <- match(ts[wchange],changeNames[,1])
    ts[wchange] <- changeNames[wm,2]
    trees$species <- as.factor(ts)
  }
  
  for(k in 1:length(yrCols)){                      # yrCol2017
    
    tk <- as.matrix( trees[,grep(yrCols[k],colnames(trees))] )
    wmin <- suppressWarnings( apply(tk,1,min, na.rm=T) )
    wna <- which(! is.finite( wmin ) )
    if(length(wna) > 0)trees <- trees[-wna,]
  }
  if(nrow(trees) == 0)return( list(yrDat = numeric(0)) )
  
  trees$species <- droplevels(trees$species) 
  
  # trees with multiple stems
 # im      <- trees[,'ID']
 # ID      <- floor(im)
  ID      <- as.character(trees$ID)
  dup     <- which(duplicated(ID))
  diamCol <- grep('diam',colnames(trees))
  
  if(length(dup) > 0){
    
    omit <- numeric(0)
    
    for(j in dup){
      
      idup <- which(ID == ID[j])
      tj   <- apply(trees[idup,diamCol], 1, max, na.rm=T)
      tj   <- which.max(tj)
      omit <- c(omit, idup[-tj])
    }
    omit <- unique(omit)
    trees <- trees[-omit,]
    trees$ID <- floor(trees$ID)
  }
  
  scols <- trees[,grep('sex',colnames(trees))]
  
  ccol <- matrix( unlist( strsplit(colnames(scols),'sex') ), ncol=2,byrow=2)[,2]
  colnames(scols) <- ccol
  
  snew <- matrix(NA, nrow(trees), length(years))
  colnames(snew) <- years
  
  snot <- which(scols == 'N', arr.ind=T)
  srep <- which(scols == 'R' | scols == 'F', arr.ind=T)
  sfem <- which(scols == 'F', arr.ind=T)
  smal <- which(scols == 'M', arr.ind=T)
  
  smal <- smal[!smal[,1] %in% sfem[,1],]  # female overrides male
  if(!is.matrix(smal))smal <- matrix(smal,1)
  
  matr <- repr <- matrix(NA, nrow(trees), ncol(scols))
  colnames(matr) <- ccol
  repr[smal[,1],] <- 0
  repr[sfem[,1],] <- 1
  
  matr[snot] <- 0
  matr[srep] <- 1
  
  #########
  
  snew[,colnames(matr)] <- matr
  
  #  repStatus <- tmp$repStatus
  
  yrDat  <- matrix(NA, length(years)*nrow(trees), length(yrCols))
  id     <- sort( unique(as.character(trees[,'ID'])) )
  tindex <- rep(id, each=length(years) )
  yindex <- rep(years, length(id) )
  
  elev <- NA
  if('elev' %in% colnames(trees))elev <- trees[,'elev']
  
  xytree <- trees[,c('x','y','UTMx','UTMy')]
  id     <- apply( cbind(plot, as.character(trees[,'ID'])), 1, paste0, 
                   collapse='-')
  xytree <- data.frame(treeID = id, plot = plot, tree = trees[,'ID'], xytree)
  xytree$elev <- elev
  
  index <- apply(cbind(tindex,yindex),1,paste0,collapse='-')
  rownames(yrDat) <- index
  colnames(yrDat) <- yrCols
  
  for(k in 1:length(yrCols)){                      # yrCol2017
    
    tk <- as.matrix( trees[,grep(yrCols[k],colnames(trees))] )
    kk <- matrix( unlist(strsplit(colnames(tk),yrCols[k])), ncol=2,byrow=T)[,2]
    yk <- as.numeric(kk)
    
    syr  <- trees[,'censinyr']                             # left censored
    imat <- matrix(years, nrow(trees), length(years), byrow=T) 
    wmin <- which(is.finite(tk),arr.ind=T)                 # measurements exist
    tmp  <- .myBy( years[wmin[,2]], wmin[,1], wmin[,1]*0+1, fun='min')
    syr  <- pmin(syr, tmp, na.rm=T)
    
    syr  <- matrix(syr, nrow(trees), length(years))           # start
    eyr  <- suppressWarnings( apply(trees[,c('deathyr','censoryr')],1,min, na.rm=T) )
    eyr  <- matrix(eyr, nrow(trees), length(years))               # end
    eyr[eyr == Inf] <- max(years) 
    imat[imat < syr | imat > eyr] <- NA
    tmat <- matrix(trees[,'ID'], nrow(tk), length(years))
    
    dmat <- matrix(NA, nrow(tk), length(years))
    icol <- match(yk, years)
    irow <- rep(1:nrow(tk),each=length(icol))
    icol <- rep(icol, nrow(tk))
    jcol <- rep(1:length(yk),nrow(tk))
    dmat[ cbind(irow,icol) ] <- tk[ cbind(irow, jcol) ]   # interpolate here
    
    start <- match( suppressWarnings( apply(imat,1,min,na.rm=T)), years)
    end   <- match( suppressWarnings( apply(imat,1,max,na.rm=T)), years)
    
    wf <- which(is.finite(start) & is.finite(end) & end > start)
    
    minVal <- 0
    maxVal <- 100
    if(yrCols[k] == 'canopy'){
      minVal <- 0
      maxVal <- 2
    } 
    
    tmp <- .interpRows(dmat[drop=F,wf,],startIndex=start[wf],endIndex=end[wf],
                       INCREASING=F,minVal=minVal,maxVal=maxVal,
                       defaultValue=NULL,tinySlope=.001)
    dmat[wf,] <- tmp
    
    if(yrCols[k] == 'canopy'){
      maxVal <- suppressWarnings( apply(dmat[drop=F,wf,], 1, max,na.rm=T) )
      wmm    <- which(is.finite(maxVal))
      mmat   <- matrix(maxVal[wmm],length(wmm),ncol(tmp))
      ww     <- which(tmp[wmm,] > mmat)
      tmp[wmm,][ww] <- mmat[ww]
    }
    
    dvec <- dmat[is.finite(imat)]
    
    it   <- tmat[is.finite(imat)]
    iy   <- imat[is.finite(imat)]
    inn  <- apply(cbind(it,iy),1,paste0,collapse='-')
    yrDat[match(inn, index),k] <- dvec
    
    if(yrCols[k] == 'diam'){
      dinc <- t(apply(dmat,1,diff,na.rm=T))
      dinc <- cbind(dinc[,1],dinc)
      dinc <- dinc[is.finite(imat)]
      growth  <- rep(0,nrow(yrDat))
      growth[match(inn, index)] <- dinc
    }
  }
  
  if('diam' %in% yrCols){
    yrDat <- cbind(yrDat,growth)
  }
  
  
    
  #  tmp <- .getRepro(snew, diam = trees[,diamCol], minDiam = 2, maxDiam = 80)
  #  snew <- tmp$rmat
  
  
  
  spec <- trees[match(tindex,trees[,'ID']),'species']
  repr <- rep(NA, nrow(yrDat))
  repr[match(inn, index)] <- snew[is.finite(imat)]
  
  
  plot   <- rep(plot, length(tindex))
  treeID <- apply(cbind(plot, tindex),1,paste0, collapse='-')
  yrDat <- data.frame(treeID = treeID, tree = tindex, species = spec, 
                      year = yindex, region = region,
                      plot = plot, repr = repr, yrDat)
  yrDat <- yrDat[is.finite(yrDat$diam) & yrDat$diam > 1,]
  xytree <- xytree[xytree$treeID %in% yrDat$treeID,]
  
  rownames(yrDat) <- NULL
  
  list(yrDat = yrDat, xytree = xytree)
}

.fac2num <- function(xx){ 
  
  dims <- dn <- NULL
  
  if(!is.null(ncol(xx))){
    dims <- dim(xx)
    dn   <- dimnames(xx)
  }
  xx <- if(is.list(xx))unlist(xx)
  xx <- as.numeric(as.character(xx)) 
  if(!is.null(dims))xx <- matrix(xx, dims[1], dims[2], 
                                 dimnames = dn)
  xx
}

.replaceString <- function(xx,now='_',new=' '){  #replace now string in vector with new
  
  ww <- grep(now,xx,fixed=T)
  if(length(ww) == 0)return(xx)
  
  for(k in ww){
    s  <- unlist( strsplit(xx[k],now,fixed=T) )
    ss <- s[1]
    if(length(s) == 1)ss <- paste( ss,new,sep='')
    if(length(s) > 1)for(kk in 2:length(s)) ss <- paste( ss,s[kk],sep=new)
    xx[k] <- ss
  }
  xx
}

.Iformat2Var <- function(iname){
  
  tt <- .replaceString(iname, 'I(','')
  tt <- .replaceString(tt, 'log(','')
  tt <- .replaceString(tt, 'sqrt(','')
  tt <- .replaceString(tt, '^2','')
  tt <- .replaceString(tt, ')','')
  tt <- .replaceString(tt, ' ','')
  tt
}

.getDesign <- function(formula, data){
  
  # one set of columns for each tree species, retain NAs
  
  specNames <- attr(data$species,'levels')
  nspec     <- length(specNames)
  
  attr(data$species,'contrasts') <- contrasts(data$species, contrasts=F)
  
  tmp1 <- model.frame(formula, data, na.action=NULL )
  tn1  <- attr( terms(tmp1), 'dataClasses' )
  sp1  <- names(tn1)[tn1 == 'numeric' | tn1 == 'nmatrix.1']
  sp1  <- .Iformat2Var(sp1)
  miss <- which(is.na(data[,sp1]),arr.ind=T)
  
  if(length(miss) > 0){
    xmean <- colMeans(data[,sp1], na.rm=T)
    data[,sp1][miss] <- 1e+8
  }
  
  x  <- model.matrix(formula, data)
  if(nspec > 1){
    ws <- grep('species',colnames(x))
    x  <- x[,ws]
  }
  missx <- which(x > 1e+5,arr.ind=T)
  
  if(length(missx) > 0){
    data[,sp1][miss] <- xmean[miss[,2]]
  }
  x  <- model.matrix(formula, data)
  if(nspec > 1){
    ws <- grep('species',colnames(x))
    x  <- x[,ws]
  }
  
  # specNames <- attr(data$species,'levels')
  specCols  <- numeric(0)
  if(nspec > 1){
    for(j in 1:length(specNames)){
      specCols <- rbind(specCols, grep( paste('species',specNames[j],sep=''),
                                        colnames(x)))
    }
    rownames(specCols) <- specNames
  }
  
  list(x = x, missx = missx, specCols = specCols)
}

columnSplit <- function(vec, sep='_', ASFACTOR = F, ASNUMERIC=F,
                        LASTONLY=F){
  
  vec <- as.character(vec)
  nc  <- length( strsplit(vec[1], sep, fixed=T)[[1]] )
  
  mat <- matrix( unlist( strsplit(vec, sep) ), ncol=nc, byrow=T )
  if(LASTONLY & ncol(mat) > 2){
    rnn <- mat[,1]
    for(k in 2:(ncol(mat)-1)){
      rnn <- columnPaste(rnn,mat[,k])
    }
    mat <- cbind(rnn,mat[,ncol(mat)])
  }
  if(ASNUMERIC){
    mat <- matrix( as.numeric(mat), ncol=nc )
  }
  if(ASFACTOR){
    mat <- data.frame(mat)
  }
  if(LASTONLY)mat <- mat[,2]
  mat
}

columnPaste <- function(c1, c2, sep='-'){
  
  FACT <- T
  if(!is.factor(c1))FACT <- F
  c1    <- as.character(c1)
  c2    <- as.character(c2)
  c12   <- apply( cbind(c1, c2) , 1, paste0, collapse=sep)
  c12   <- .replaceString(c12, ' ', '')
  if(FACT) c12 <- as.factor(c12)
  c12
}

.setupData <- function(formulaFec, formulaRep, tdata, sdata, facLevels,
                       xytree, xytrap, specNames, seedNames, AR, YR, 
                       yearEffect, minDiam, maxDiam){
  
  
  # formulas have 'species *' already
  arList <- numeric(0)
  
  if(!'active' %in% colnames(sdata)){
    cat('"\nactive" column added to sdata\n')
    sdata$active <- 1
  }
  if(!'area' %in% colnames(sdata)){
    cat('\n"area" column added to sdata\n')
    sdata$area <- 1
  }
  wna <- which(is.na(sdata$active) | sdata$active == 0)
  if(length(wna) > 0){
    cat('\nsome sdata$active values undefined\n')
    sdata$active[wna] <- .1
  }
  wna <- which(is.na(sdata$area) | sdata$area == 0)
  if(length(wna) > 0){
    cat('\nsome sdata$area values undefined or zero\n')
    sdata$area[wna] <- 1
  }
  
  #####################
  
  wna <- which(!sdata$trap %in% xytrap$trap)
  if(length(wna) > 0){
    mm <- paste0(sdata$trap[wna], collapse=', ')
    stop(paste('missing traps in xytrap: ',mm,sep=''))
  }
  
  tree <- tdata
  w    <- which(tree$species %in% specNames)
  tree <- tree[w,]
  tree$species <- droplevels(tree$species)
  
  tree          <- tree[order(as.character(tree$plot), tree$tree, tree$year),]
  tree$treeID   <- columnPaste( tree$plot, tree$tree )
  xytree$treeID <- columnPaste(xytree$plot, xytree$tree)
  
  sdata        <- sdata[order(as.character(sdata$plot), 
                                    sdata$trap, sdata$year),]
  sdata$trapID  <- columnPaste( sdata$plot, sdata$trap )
  xytrap$trapID <- columnPaste( xytrap$plot, xytrap$trap )
  xytrap        <- xytrap[xytrap$trapID %in% sdata$trapID,]
  
  plotNames <- sort(unique( as.character(xytrap$plot) ))
  years     <- as.numeric( sort(unique( as.character(sdata$year) )) )
  nplot     <- length(plotNames)
  
  tid     <- tree[!duplicated(tree$treeID),]
  specTab <- table(tid$species)
  wna     <- which(specTab < 8)
  
  if(length(wna) > 0){               # species too rare
    bad  <- names(specTab)[wna]
    ww   <- which(!as.character(tree$species) %in% bad)
    tree <- tree[ww,]
    rn   <- paste0(names(specTab)[wna],collapse=', ')
    specNames <- specNames[!specNames %in% names(specTab)[wna]]
    nspec <- length(specNames)
    tree$species <- droplevels(tree$species)
    cat('\ntoo rare, removed:\n')
    print( rn )
    
    if(bad %in% colnames(sdata)){
      wb <- which(colnames(sdata) == bad)
      sdata <- sdata[,-wb]
      seedNames <- seedNames[seedNames != bad]
    }
    
    
    ptab <- table(tree$plot)
    
    w0   <- which(ptab == 0)
    
    if(length(w0) > 0){
      pmiss <- names(ptab)[w0]
      wm <- which(as.character(xytree$plot) %in% pmiss)
      if(length(wm) > 0)xytree <- xytree[-wm,]
      
      wm <- which(as.character(sdata$plot) %in% pmiss)
      if(length(wm) > 0)sdata <- sdata[-wm,]
      
      wm <- which(as.character(xytrap$plot) %in% pmiss)
      if(length(wm) > 0)xytrap <- xytrap[-wm,]
      #    plots  <- plots[!plots %in% pmiss]
      plotNames  <- plotNames[!plotNames %in% pmiss]
    }
    tree   <- cleanFactors(tree)
    sdata  <- cleanFactors(sdata)
    xytree <- cleanFactors(xytree)
    xytrap <- cleanFactors(xytrap)
  }
  
  wna <- which(is.na(tree$diam))
  if(length(wna) > 0){
    mm <- paste0(tree$tree[wna], collapse=', ')
    tree <- tree[-wna,]
    cat(paste('\nremoved trees with missing diam:\n ',mm,sep=''))
    tree$species <- droplevels(tree$species)
  }
  wna <- which(!sdata$trapID %in% xytrap$trapID)
  if(length(wna) > 0){
    mm <- paste0(sdata$trapID[wna], collapse=', ')
    stop(paste('missing traps in xytrap: ',mm,sep=''))
  }
  wna <- which(!tree$treeID %in% xytree$treeID)
  if(length(wna) > 0){
    mm <- paste0(tree$treeID[wna], collapse=', ')
    stop(paste('missing trees in xytree: ',mm,sep=''))
  }
  
  # remove seed data where there are no trees
  plotRm <- character(0)
  
  for(j in plotNames){
    tj <- which(as.character(xytree$plot) == j)
    sj <- which(as.character(xytrap$plot) == j)
    t2 <- which(as.character(tdata$plot) == j)
    s2 <- which(as.character(sdata$plot) == j)
    
    if(length(tj) == 0 | length(t2) == 0){
      plotRm <- c(plotRm,j)
      cat(paste('\nplot ',j,' is absent from xytree or treeData\n'))
      next
    }
  }
  
  if(length(plotRm) > 0){
    plotNames <- plotNames[!plotNames %in% plotRm]
    wr        <- which(as.character(sdata$plot) %in% plotRm)
    sdata  <- sdata[-wr,]
    wr        <- which(as.character(xytrap$plot) %in% plotRm)
    xytrap    <- xytrap[-wr,]
    wr        <- which(as.character(tdata$plot) %in% plotRm)
    tdata  <- tdata[-wr,]
    wr        <- which(as.character(xytree$plot) %in% plotRm)
    xytee    <- xytree[-wr,]
    
    xytrap <- cleanFactors(xytrap)
    sdata  <- cleanFactors(sdata)
    tdata  <- cleanFactors(tdata)
    xytree <- cleanFactors(xytree)
  }
  
  # retain trees in plot-years that have seed traps
  py1  <- columnPaste( tree$plot, tree$year )
  py2  <- columnPaste( sdata$plot, sdata$year )
  plotYears <- sort(unique(as.character((py2))))
  keep <- which( as.character(py1) %in% plotYears )
  tree <- tree[keep,]
  py1  <- py1[keep]

  tree$plotYr <- py1
  tree$plotyr <- match( as.character(py1), plotYears )
  
  sdata$plotYr <- py2
  sdata$plotyr <- match( as.character(py2), plotYears )
  
  tree      <- cleanFactors(tree)
  xytree    <- xytree[xytree$treeID %in% tree$treeID,]
  specNames <- attr(tree$species,'levels')
  nspec     <- length(specNames)
  
  years  <- sort(unique(tree$year))
  
  print(dim(tree))
  
  vtypes <- getVarType(colnames(tree), tree, i=tree$treeID, j = tree$year)
  
  if('repr' %in% names(tree))vtypes$repr = 'ij'
  if('repMu' %in% names(tree))vtypes$repMu = 'ij'
  if('repSd' %in% names(tree))vtypes$repSd = 'ij'
  
  wnull <- which(is.null(vtypes))
  if(length(wnull) > 0){
    for(k in 1:length(wnull)) vtypes[[k]] <- 'ij'
  }
  
  # when tree census not in seed year
  tree$treeID <- as.factor(tree$treeID)
  tree$treeID <- droplevels(tree$treeID)
  treeIDs     <- attr(tree$treeID, 'levels')
  tree$i      <- match(as.character(tree$treeID), treeIDs)
  ww <- which(!as.character(sdata$plotYr) %in% 
                as.character(tree$plotYr) ) #trees not measured in seed year
  
  if(length(ww) > 0){
    
    ijIndex   <- cbind( match(as.character(tree$treeID), treeIDs),
                        match(tree$year, years) )
    
    py <- sdata$year[ww]
    pl <- as.character(sdata$plot[ww])
    pa <- sort(unique(pl))
    
    ijFull <- ijIndex
    
    for(k in 1:length(pa)){
      
      ij <- numeric(0)
      
      wl <- which(as.character(tree$plot) %in% pl[k] & !tree$year %in% py  &
                    tree$year %in% (py+1) )
      wh <- which(as.character(tree$plot) %in% pl[k] & !tree$year %in% py  &
                    tree$year %in% (py-1) ) 
      if(length(wl) > 0){
        ijk   <- cbind( match(as.character(tree$treeID[wl]), treeIDs),
                        match(tree$year[wl]-1, years) )
        ij <- rbind(ij,ijk)
      }
      if(length(wh) > 0){
        ijk   <- cbind( match(as.character(tree$treeID[wh]), treeIDs),
                        match(tree$year[wh]+1, years) )
        ij <- rbind(ij,ijk)
      }
      ijFull <- rbind(ijFull, ij)
    }
  
    tree <- fillMissing(vtypes, tree, icol='i', 
                       jcol='year', jtimes=years, ijIndex = ijIndex,
                       ijFull = ijFull)$data

    tree$plotYr <- as.factor( columnPaste( as.character(tree$plot),tree$year ) )
  }
  
  allYears    <- sort(unique(tree$year))
  tree$plotyr <- match( as.character(tree$plotYr), plotYears )
  tree$times  <- match(tree$year,allYears)
  
  treePlotYr <- columnPaste(tree$treeID,tree$year,sep='-')
  wdup <- which(duplicated(treePlotYr))
  if(length(wdup) > 0){
    tree <- tree[-wdup,]
  }
  
  if(!'repr' %in% colnames(tree))tree$repr <- NA
  if(!'repMu' %in% colnames(tree))tree$repMu <- NA
  
  tree$repMu[is.na(tree$repr) & tree$diam > minDiam] <- .5   # prior mean
  tree$repMu[is.na(tree$repr) & tree$diam < minDiam] <- 0
  tree$repr[is.na(tree$repr) & tree$diam < minDiam] <- 0
  tree$repr[is.na(tree$repr) & tree$diam > maxDiam] <- 1
  tree$repMu[is.na(tree$repr) & tree$diam > maxDiam] <- 1
  
  
  if(AR){
    plag   <- yearEffect$p
    gnames <- character(0)
    if('specGroups' %in% names(yearEffect))gnames <- yearEffect$specGroups
    if('plotGroups' %in% names(yearEffect)){
      if(length(gnames) == 0){
        gnames <- yearEffect$plotGroups
      }else{
        gnames <- c(gnames,yearEffect$plotGroups)
      }
    }
    
    group <- tree[,gnames]
    if(length(gnames) > 1)group <- apply( tree[,gnames], 1, paste0, collapse='-')
    
    tree$group <- as.factor(group)
    tree$g     <- match(as.character(tree$group),attr(tree$group,'levels'))
    
    allYears <- c( (min(years) - plag):(max(years) + plag) )
    tree$times <- match(tree$year,allYears)
    
    tmp   <- msarSetup(tree, variables=vtypes, plag, icol='treeID', 
                       jcol = 'times', gcol='group')
    tree       <- tmp$xdata
    times      <- tmp$times
    groupByInd <- tmp$groupByInd
    betaYr     <- tmp$betaYr
    yeGr       <- rownames(betaYr)
    ngroup     <- length(yeGr)
    yrIndex    <- tree[,c('g','j')]
    colnames(yrIndex) <- c('group','year')
    
    specPlot  <- as.character( columnPaste(tree$spec, tree$plot) )
    specPlots <- sort(unique(specPlot))
    specPlot  <- match(specPlot,specPlots)
    yrIndex   <- cbind(yrIndex,specPlot)
    
    tree$year <- allYears[ tree$times ]
    # reorder plot years
    tree$plotYr  <- as.factor( columnPaste(tree$plot,tree$year) )
    plotYears    <- attr(tree$plotYr, 'levels')
    tree$plotyr  <- match(tree$plotYr, plotYears)
    
    wkeep <- which( sdata$plotYr %in% plotYears ) # some trees have been removed
    sdata <- sdata[wkeep,]
    sdata$plotyr <- match(sdata$plotYr, plotYears)
    
    tree$plotTreeYr <- columnPaste(tree$treeID, tree$year)
    
    sdata <- cleanFactors(sdata)
    tree  <- cleanFactors(tree)
    
    nocol <- c('i','j','g')
    tree <- tree[,!colnames(tree) %in% nocol]
 
    lagMat <- msarLagTemplate(plag, tree, icol='treeID', jcol='times', 
                                gcol='group', ocol='obs')
    rownames(lagMat[[1]]) <- tree$plotTreeYr[lagMat[[1]][,1]]

    arList <- list(times = times, groupByInd = groupByInd, betaYr = betaYr,
                   yeGr = rownames(betaYr), ngroup = length(yeGr),
                   yrIndex = yrIndex, lagMatrix = lagMat$matrix, 
                   lagGroup = lagMat$group)
  }else{
    tree$obs <- 1
  }
 
  tree$plotTreeYr <- columnPaste(tree$treeID, tree$year )
  
  seedSummary <- with(sdata, table(plot, year) )
  treeSummary <- with(tree, table(plot, year) )
  plotNames   <- rownames(seedSummary)
  
  
  if(nspec == 1){
    formulaFec <- as.formula( .replaceString( as.character(formulaFec), 'species *','') )
    formulaRep <- as.formula( .replaceString( as.character(formulaRep), 'species *','') )
  }
  
  tmp   <- model.frame(formulaFec, tree, na.action=NULL)
  if(length(tmp) == 0)tmp <- numeric(0)
  tmp1  <- model.frame(formulaRep, tree, na.action=NULL)
  
  if(length(tmp1) > 0){
    
    xfecMiss <- which(is.na(tmp1),arr.ind=T)
    wnew <- which(!colnames(tmp1) %in% colnames(tmp))
    
    if(length(wnew) > 0){                               # all unique columns
      tmp <- cbind(tmp,tmp1[,wnew])
    }
  }
  xallNames <- colnames(tmp)
  
  scode <- names(tmp[ which(sapply( tmp, is.factor )) ])
  if(length(scode) > 0){
    for(j in 1:length(scode)) tree[,scode[j]] <- droplevels(tree[,scode[j]])
  }
  
  specNames <- attr(tree$species,'levels')
  
  standX <- character(0)
  xmean <- xsd <- numeric(0)
  
  wstand <- which(!colnames(tmp) %in% scode)
  
  
  if(length(wstand) > 0){
    
    standX    <- colnames(tmp)[wstand]
    
      wlog <- grep( "log(", standX, fixed = T)
      if(length(wlog) > 0)standX <- standX[ -wlog ]
      wlog <- grep( "sqrt(", standX, fixed = T)
      if(length(wlog) > 0)standX <- standX[ -wlog ]
    
    if(length(standX) > 0){
      
      treeUnstand <- tree[,standX, drop=F]
      
      xmean <- colMeans(tree[,standX, drop=F],na.rm=T)
      xsd   <- apply(tree[,standX, drop=F],2, sd, na.rm=T)
      xss <- t( (t(tree[,standX, drop=F]) - xmean)/xsd )
      tree[,colnames(treeUnstand)] <- xss
    }
  }
  
  # relevel factors
  tmp1 <- model.frame(formulaFec,tree)
  tmp2 <- model.frame(formulaRep,tree)
  tmp  <- unique( c(names(tmp1), names(tmp2)) )
  
  wf <- which( tmp %in% names(facLevels) )
  if(length(wf) > 0){
    for(k in wf){
      wk <- which(names(facLevels) == tmp[k])
      wt <- which(names(tree) == tmp[wf])
      levels(tree[[wt]]) <- facLevels[[wk]]
    }
  }
  
  tmp  <- .getDesign(formulaFec, tree)
  xfec <- tmp$x
  if(nspec > 1)xfec <- xfec[,grep('species',colnames(xfec)),drop=F]
  xfecMiss <- tmp$missx
  xfecCols <- tmp$specCols
  
  tmp  <- .getDesign(formulaRep, tree)
  xrep <- tmp$x
  if(nspec > 1)xrep <- xrep[,grep('species',colnames(xrep)),drop=F]
  xrepMiss <- tmp$missx
  xrepCols <- tmp$specCols
  
  rank <- qr(xfec)$rank
  if(rank < ncol(xfec))stop('fecundity design not full rank')
  rank <- qr(xrep)$rank
  if(rank < ncol(xrep))stop('maturation design not full rank')
  
  xfecU <- xfec
  xrepU <- xrep
  xfecT <- xrepT <- NULL
  
  xfecs2u <- diag(ncol(xfec))
  xreps2u <- diag(ncol(xrep))
  colnames(xfecs2u) <- rownames(xfecs2u) <- colnames(xfec)
  colnames(xreps2u) <- rownames(xreps2u) <- colnames(xrep)
  
  xfecu2s <- xfecs2u
  xrepu2s <- xreps2u
  
  if(length(xmean) > 0){   # unstandardized
    
    tree[,standX] <- treeUnstand
    tmp <- .unstandBeta(formula = formulaFec, 
                        xdata = tree, xnow = xfec, xmean = xmean)
    xfecU   <- tmp$x
    xfecs2u <- tmp$unstand    #unstandardize
    xfecu2s <- tmp$stand
    
    tmp <- .unstandBeta(formula = formulaRep, 
                        xdata = tree, xnow = xrep, xmean = xmean)
    xrepU   <- tmp$x
    xreps2u <- tmp$unstand
    xrepu2s <- tmp$stand
    
    xmean[xmean < 1e-10] <- 0
  }
  
  tree  <- cleanFactors(tree)
  sdata <- cleanFactors(sdata)
  
  treeIDs <- attr(tree$treeID, 'levels')
  trapIDs <- attr(sdata$trapID, 'levels')
  
  distall  <- numeric(0)
  treeid   <- trapid <- character(0)
  
  for(j in plotNames){
    
    tj <- which(xytree$plot == j)
    sj <- which(xytrap$plot == j)
    
    if(length(tj) == 0){
      plotRm <- c(plotRm,j)
      message(paste('plot ',j,' is absent from xytree'))
      next
    }
    if(length(sj) == 0)stop( paste('plot', j ,'has no traps in xytrap') )
    
    xy1     <- xytree[tj,]
    xy2     <- xytrap[sj,]
    treeid  <- c(treeid,xytree$treeID[tj])
    trapid  <- c(trapid,xytrap$trapID[sj])
    da      <- .distmat(xy1[,'x'],xy1[,'y'],xy2[,'x'],xy2[,'y']) 
    da      <- round(da, 1)
    colnames(da) <- xy1$treeID
    rownames(da) <- xy2$trapID
    distall <- .blockDiag(distall,da)
  }
  
  distall <- distall[,treeIDs]
  distall <- distall[trapIDs,]
  
  dcol  <- match(as.character(tree$treeID), treeIDs)
  drow  <- match(as.character(sdata$trapID), trapIDs)
  tree$dcol <- dcol
  sdata$drow <- drow
  
  distall[distall == 0] <- 10000
  
  tree$species <- droplevels(tree$species)
  
  seedNames <- seedNames[seedNames %in% colnames(sdata)]
  
  keepCol <- c('trapID','plot','year','trap','plotYr','plotyr','drow',
               'area','active',seedNames)
  sdata <- sdata[,keepCol]

  if(!AR){
    ngroup  <- 1
    yrIndex <- matrix(1,nrow(tree),4)
    colnames(yrIndex) <- c('dcol','group','year','groupYr')
    yrIndex[,'year']  <- match(tree$year,years)
    yrIndex[,'dcol']  <- tree$dcol
  }
  
  if(YR){
    pgroup <- numeric(0)
    if( 'plotGroups' %in% names(yearEffect) ){
      pgroup <- tree[,yearEffect$plotGroups]
    }
    if('specGroups' %in% names(yearEffect)){
      if(length(pgroup) == 0){
        pgroup <- tree[,yearEffect$specGroups]
      }else{
        pgroup <- columnPaste( as.character(pgroup), 
                               as.character(tree[,yearEffect$specGroups]), '-')
      }
    }
    
    pall <- unique(pgroup)
    group <- match(pgroup,pall)
    yrIndex[,'group'] <- group
    gy     <- apply( yrIndex, 1, paste0, collapse='-')
    gyall  <- sort(unique(gy))
    groupYr  <- match(gy,gyall)
    yrIndex[,'groupYr']  <- groupYr
  }
  
  specPlot  <- as.character( columnPaste(tree$spec, tree$plot) )
  tree$specPlot <- specPlot
  specPlots <- sort(unique(specPlot))
  specPlot  <- match(specPlot,specPlots)
  yrIndex   <- cbind(yrIndex,specPlot)
  
  ic <- sapply(tree,is.character)
  wc <- which(ic)
  if(length(wc) > 0){
    for(j in wc){ tree[[j]] <- as.factor(tree[[j]]) }
  }
  
  list(tdata = tree, sdata = sdata, distall = distall, seedNames = seedNames,
       specNames = specNames, arList = arList, plotYears = plotYears,
       xytree = xytree, xytrap = xytrap, plotNames = plotNames,
       plots = plotNames, years = years, xfec = xfec, xrep = xrep, scode = scode,
       xfecMiss = xfecMiss, 
       xrepMiss = xrepMiss, xfecCols = xfecCols, xrepCols = xrepCols,
       xmean = xmean, xsd = xsd, xfecU = xfecU, 
       xfecs2u = xfecs2u, xfecu2s = xfecu2s, xreps2u = xreps2u, xrepu2s = xrepu2s,
       xrepU = xrepU, xrepT = xrepT, specPlots = specPlots, yrIndex = yrIndex)
}

.stripI <- function( vname ){
  
  vname <- .replaceString( vname, "I(log(", "")
  vname <- .replaceString( vname, "I(sqrt(", "")
  vname <- .replaceString( vname, "I(", "")
  vname <- .replaceString( vname, "^2)", "")
  vname <- .replaceString( vname, "))", "")
  vname
}

.unstandBeta <- function(formula, xdata, xnow, xmean=NULL){
  
  # xnow  - current standardized matrix
  # xdata - data.frame, unstandardized variables
  
  tmp <- model.frame(formula, xdata, na.action=NULL)
  
  xterm <- names( tmp )
  st    <- grep('species',xterm)
  if(length(st) > 0)xterm <- xterm[-st]
  
  if(length(xmean) > 0){
    xterm <- xterm[xterm %in% names(xmean)]
    if(length(xterm) == 0)return(list(x = xnow, trans = NULL) )
  }
  
  xfu  <- .getDesign(formula, xdata)$x    # unstandardized
  
  if(length(st) > 0 & length(xterm) > 1)xfu  <- xfu[,grep('species',colnames(xfu))]
  
  XX <- crossprod(xfu)
  diag(XX) <- diag(XX) + .00000001
  unstand <- solve(XX)%*%crossprod(xfu,xnow) 
  unstand[abs(unstand) < 1e-10] <- 0
  
  UU <- crossprod(xnow)
  stand <- solve(UU)%*%crossprod(xnow,xfu)
  stand[abs(stand) < 1e-10] <- 0
  
  list(x = xfu, unstand = unstand, stand = stand)
}

.blockDiag <- function(mat1,mat2){
  
  #creates block diagional
  
  if(length(mat1) == 0)return(mat2)
  
  namesc <- c(colnames(mat1),colnames(mat2))
  namesr <- c(rownames(mat1),rownames(mat2))
  
  nr1 <- nrow(mat1)
  nr2 <- nrow(mat2)
  nc1 <- ncol(mat1)
  nc2 <- ncol(mat2)
  nr  <- nr1 + nr2
  nc  <- nc1 + nc2
  
  new <- matrix(0,nr,nc)
  new[ 1:nr1, 1:nc1 ] <- mat1
  new[ (nr1+1):nr, (nc1+1):nc ] <- mat2
  colnames(new) <- namesc
  rownames(new) <- namesr
  new
}

.getLambda <- function(tdat, sdat, AA, ug, ff, zz, R, SAMPR, 
                       distance, yrs, PERAREA=F, SPECPRED = F){
  
  # tdata needs species, year, dcol
  # seedData needs year, drow, active, area
  
  # PERAREA - from per-trap to per-area
  # if length(AA === 1) then it must equal 1
  # SPECPRED - predict species rather than seed types
  
  nf <- length(ff)
  fz <- ff*zz
  
  if(SAMPR | length(R) > 1){
    if(SPECPRED){
      ff <- matrix(0, nf, nrow(R))
      jj <- match(as.character(tdat$specPlot),rownames(R)) # predict species
      ff[ cbind(1:nf,jj) ] <- fz
    }else{
      ff <- matrix(fz,length(ff),ncol=ncol(R))*
            R[drop=F,as.character(tdat$specPlot),] # predict types
    }
  }else{
    ff <- matrix(fz,ncol=1)
  }
  
 # lindex <- c(1:ncol(ff)) - 1
  
  uvec <- ug[ attr(distance,'group') ]
  dmat <- t(uvec/pi/(uvec + t(distance)^2)^2)
  dmat[dmat < 1e-8] <- 0

  
  lambda <- kernYrRcpp(dmat, ff, yrs, seedyear = sdat[,'year'],
                    treeyear = tdat[,'year'], seedrow = sdat[,'drow'],
                    treecol = tdat[,'dcol'])
  
  if(SPECPRED){
    colnames(lambda) <- rownames(R)
    
    sname <- sort(unique(attr(R,'species')))
    ii <- rep( c(1:nrow(lambda)), ncol(lambda) )
    jj <- match(attr(R,'species'),sname)
    jj <- rep(jj, each=nrow(lambda))
    
    lambda <- .myBy(as.vector(lambda), ii, jj, fun='sum')
    colnames(lambda) <- sname
    
  }else{
    colnames(lambda) <- colnames(R)
  }
  
  if(PERAREA | length(AA) == 1) return(lambda)    # per area
  
  lambda*matrix( AA, nrow(lambda), ncol(lambda) )  # per trap
}

.tnorm <- function(n,lo,hi,mu,sig, tiny=0){   
  
  #normal truncated lo and hi
  
  if(length(lo) == 1 & length(mu) > 1)lo <- rep(lo,length(mu))
  if(length(hi) == 1 & length(mu) > 1)hi <- rep(hi,length(mu))
  
  q1 <- pnorm(lo,mu,sig)
  q2 <- pnorm(hi,mu,sig) 
  
  z <- runif(n,q1,q2)
  z <- qnorm(z,mu,sig)
  
  z[z == Inf]  <- lo[z == Inf] + tiny
  z[z == -Inf] <- hi[z == -Inf] + tiny
  z
}

.getPlotLayout <- function(np){
  
  # np - no. plots
  
  if(np == 1)return( list( mfrow=c(1,1), left=1, bottom=c(1,2) ) )
  if(np == 2)return( list( mfrow=c(1,2), left = 1, bottom = c(1,2) ) )
  if(np == 3)return( list( mfrow=c(1,3), left = 1, bottom = c(1:3) ) )
  if(np <= 4)return( list( mfrow=c(2,2), left = c(1, 3), bottom = c(3:4) ) )
  if(np <= 6)return( list( mfrow=c(2,3), left = c(1, 4), bottom = c(4:6) ) )
  if(np <= 9)return( list( mfrow=c(3,3), left = c(1, 4, 7), bottom = c(7:9) ) )
  if(np <= 12)return( list( mfrow=c(3,4), left = c(1, 5, 9), bottom = c(9:12) ) )
  if(np <= 16)return( list( mfrow=c(4,4), left = c(1, 5, 9, 13), 
                            bottom = c(13:16) ) )
  if(np <= 20)return( list( mfrow=c(4,5), left = c(1, 6, 11, 15), 
                            bottom = c(15:20) ) )
  if(np <= 25)return( list( mfrow=c(5,5), left = c(1, 6, 11, 15, 20), 
                            bottom = c(20:25) ) )
  if(np <= 25)return( list( mfrow=c(5,6), left = c(1, 6, 11, 15, 20, 25), 
                            bottom = c(25:30) ) )
  return( list( mfrow=c(6,6), left = c(1, 7, 13, 19, 25, 31), bottom = c(31:36) ) )
}

.seedProb <- function(tspec, tyear, tdcol, ug, ff, distall, sdata, 
                      seedNames, zz, R, SAMPR, years){
  
  tt <- data.frame(specPlot = tspec, year = tyear, dcol = tdcol)
  
  lambda <- .getLambda(tt, sdata, sdata$active*sdata$area, ug, ff, zz, R, SAMPR, 
                       distall, years, PERAREA=F)
  lambda <- lambda + 1e-8
  ss     <- as.matrix(sdata[,seedNames])
  dpois(ss, lambda, log=T)
}
  
.myBy <- function(x, i, j, summat=matrix(0,max(i),max(j)), 
                    totmat=summat, fun='mean'){  
  
  nn <- length(x)
  if( nn != length(i) | nn != length(j) )
    stop('vectors unequal in byFunctionRcpp')
  if( nrow(summat) < max(i) | ncol(summat) < max(j) )
    stop('matrix too small')
  
  ww <- which(is.na(x))
  if(length(ww) > 0){
    x <- x[-ww]
    i <- i[-ww]
    j <- j[-ww]
  }
  
  frommat <- cbind(i,j,x)
  
  nr  <- nrow(frommat)
  
  maxmat <- summat*0 - Inf
  minmat <- summat*0 + Inf
  
  tmp <- byRcpp(nr, frommat, totmat, summat, minmat, maxmat)
  
  if(fun == 'sum')return(tmp$sum)
  if(fun == 'mean'){
    mu <- tmp$sum/tmp$total
    mu[is.na(mu)] <- 0
    return(mu)
  }
  if(fun == 'min'){
    return( tmp$min )
  }
  tmp$max
}

.tnormMVNmatrix <- function(avec, muvec, smat, 
                            lo=matrix(-1000,nrow(muvec),ncol(muvec)), 
                            hi=matrix(1000,nrow(muvec),ncol(muvec)),
                            whichSample = c(1:nrow(smat))){
  
  # lo, hi must be same dimensions as muvec,avec
  # each sample is a row
  
  lo[lo < -1000] <- -1000
  hi[hi > 1000]  <- 1000
  
  if(max(whichSample) > length(muvec))
    stop('whichSample outside length(muvec)')
  
  nd <- dim(avec)
  
  r <- avec
  a <- trMVNmatrixRcpp(avec, muvec, smat, 
                       lo, hi, whichSample, 
                       idxALL = c(0:(nrow(smat)-1)) ) 
  r[,whichSample] <- a[,whichSample]
  r
}

.initEM <- function(last0first1, yeGr, distall, ug, priorU, minDiam, minU, maxU,
                    tdata, sdata, specNames, seedNames, R, SAMPR, years, 
                    plotYears, z, xfec, zobs, nsim=100){
  
  tiny <- 1e-3
  id <- sort(unique(as.character(tdata$treeID)))
  tt <- tdata[ match(id,as.character(tdata$treeID)), ]  #unique trees
  
  tdata$fecMin[zobs == 0 & tdata$fecMin > .0001] <- .0001
  tdata$fecMin[zobs == 1 & tdata$fecMin < 1] <- 1
  tdata$fecMax[zobs == 0 & tdata$fecMax > 1] <- 1
  tdata$fecMax[zobs == 1 & tdata$fecMax < 1.001] <- 1.001
  
  
  ngroup <- max( attr(distall, 'group') )
  nspec  <- length(specNames)
  
  dcol <- numeric(0)
  
  for(j in 1:ngroup){
    
    if(ngroup > 1){
      wj <- which(tt$group == yeGr[j])
      tj <- tt[drop=F,wj,]
    }else{
      wj <- 1:nrow(tt)
      tj <- tt
    }
    
    dtmp <- distall[,tt$dcol[wj],drop=F]
    dtmp[dtmp > 1000] <- NA
    
    dmin  <- apply(dtmp,2,min, na.rm=T)
    qdiam <- quantile(tj$diam,.5)
    
    ww  <- which(dmin < 30 & tj$diam > qdiam)
    dcol <- c(dcol, tj$dcol[ww])
  }
  wtree <- which(tdata$dcol %in% dcol & zobs != 0)
  
  fg <- rep(.5,nrow(tdata))
  
  ss   <- sdata[,seedNames,drop=F]
  ff <- sapply(ss,is.factor)
  if(!all(!ff))ss <- .fac2num( ss )
  ss <- rowSums(ss, na.rm=T)
  
  for(j in 1:length(plotYears)){
    
    i  <- which(as.character(tdata$plotYr) == plotYears[j] )
    m  <- i[i %in% wtree]
    k  <- which(sdata$plotYr == plotYears[j])
    sj <- sum(ss[k])
    
    if(length(k) == 0)next
    if(length(i) == 0)next
    if(length(m) == 0 & sj > 0)m <- i
    if(length(m) == 0)next
    
    d <- unique(tdata$dcol[m])
    dj <- d[d %in% dcol]
    ij <- m[d %in% dcol]
    if(length(dj) < 1){
      ij <- which(tdata$plotYr == plotYears[j] & tdata$dcol %in% d)
      dj <- tdata$dcol[ij]
    }
 
    dk     <- distall[ sdata[k,'drow'], dj, drop=F ]
    kern   <- priorU/pi/(priorU + dk^2)^2
    fg[ij] <- .getF( kern, gg = ss[k]/sdata$area[k]/
                       (.1 + sdata$active[k]) )
  }
  fg[!is.finite(fg)] <- .1
 # fg[zobs == 1 & fg < 1] <- 1.00001
 # fg[zobs == 0 & fg >= 1] <- .5
  
  nn <- length(fg)
  
  lo <- tdata$fecMin
  hi <- tdata$fecMax
 # lo[z == 1] <- 1
#  hi[z == 0] <- 1
  
  fg[fg < lo] <- lo[fg < lo]
  fg[fg > hi] <- hi[fg > hi]
  
   fg <- .tnorm(nn, lo, hi, fg, .1)
   
   propF <- fg/20
   propU <- .1
  # au <- 0
   
   # by plot-yr
   ii <- tdata[,'plotyr']     #consider bnow[znow == 0] = 0
   sm <- matrix(0, max(c(tdata$plotyr, sdata$plotyr)), 1)
   
   pcheck <- seq(1,nsim,by=20)
   
   cat("\ninitializing\n")
   
   pbar <- txtProgressBar(min=1,max=nsim,style=1)
   
   
   ug <- rep(priorU, ngroup )
   nn <- nrow(tdata)
   
  pall <- -1e+10
  count <- 0
  
  for(g in 1:nsim){
    
    fnew <- .tnorm(nn, lo, hi, fg, rexp(nn,1/propF))
    fnew[fnew < 0] <- 0
  #  unew <- exp(.tnorm(length(ug),log(minU),log(maxU),log(ug),propU))

    pnow <- .seedProb(tdata$specPlot, tdata$year, tdata$dcol,
                      ug, fg, distall, sdata, 
                      seedNames, z, R, SAMPR, years)
    pnew <- .seedProb(tdata$specPlot, tdata$year, tdata$dcol,
                      ug, fnew, distall, sdata, 
                      seedNames, z, R, SAMPR, years)
    
    pnow[pnow < -1e+8] <- -1e+8   # intensity parameter is zero
    pnew[pnew < -1e+8] <- -1e+8
    
    # by plot-yr
    ii <- sdata[,'plotyr']
    ii <- rep(ii, length(seedNames))
    
    pnow <- .myBy(as.vector(pnow), ii, ii*0 + 1, summat = sm*0, fun='sum')
    pnew <- .myBy(as.vector(pnew), ii, ii*0 + 1, summat = sm*0, fun='sum')
    
    if(g == 1)accept <- pnow*0
    
    a  <- exp(pnew - pnow)        #wt on seed data
    az  <- runif(length(a),0,1)
    aw  <- which(az < a)
    
    if(length(aw) > 0){
      wa <- which(tdata[,'plotyr'] %in% aw)
      fg[ wa ] <- fnew[ wa ]
      accept[aw] <- accept[aw] + 1
    }
    
  #  a <- exp(sum(pnew) - sum(pnow))
  #  if(a > runif(1,0,1)){
  #    ug <- unew
  #    au <- au + 1
  #  }
    if(g %in% pcheck){
      whi <- which(accept > g/2)
      if(length(whi) > 0)propF[whi] <- propF[whi]*2
      wlo <- which(accept < g/5)
      if(length(wlo) > 0)propF[wlo] <- propF[wlo]/2
   #   if(au > g/2)propU <- propU*2
   #   if(au < g/5)propU <- propU/2
      
      pq   <- sum(pnow)
      dl   <- pq - pall
      pall <- pq
      
      if(dl < 0){
        count <- count + 1
        if(count > 4)break
      }
      
    }
    setTxtProgressBar(pbar,g)
  }
  fg[fg < tiny] <- tiny
  
 # fstart <- rep(1,nn)
 # fstart[wtree] <- fg[wtree]
  
  XX <- crossprod(xfec[wtree,])
  diag(XX) <- diag(XX) + .000001
  
  bf <- solve(XX)%*%crossprod(xfec[wtree,],log(fg[wtree]))
  mu <- exp(xfec%*%bf)
  mu[wtree] <- fg[wtree]
  
  fstart <- .tnorm(nn,lo,hi,mu,.1)
  fstart[wtree] <- fg[wtree]
  fstart[fstart >= .95 & z == 0] <- .95
  fstart[fstart < 1 & z == 1] <- 1.01
  fstart
}

.ranEffVar <- function(yy, xf, Arand, sg, xrandCols){   # rewrite in cpp
  
  # marginalized variance from random effects
  # xr = xfec[,xrandCols]
  
  xx    <- xf[,xrandCols]%*%Arand           # w'A
  revar <- sg + rowSums(xx*xf[,xrandCols])  # sigma + w'Aw
  XR    <- t(xf/revar)                      
  VI    <- XR%*%xf                          # Q x Q
  V     <- solve(VI)
  v     <- XR%*%yy
  XX    <- XR*t(xf)                         # (x'x/s)_it
  XY    <- XR*yy                            # (x'y/s)_it
  u     <- rowSums(1/XX*XY)
  list(v = v, V = V)
}

.wrapperBeta <- function( priorB, priorIVB, minU, maxU, 
                          priorU, priorVU, SAMPR, obsRows, obsYr, obsRowSeed,
                          sdata, tdata, seedNames, xfecCols,
                          xrepCols, last0first1, ntree, nyr, 
                          betaPrior, years, distall, YR, AR, yrIndex,
                          RANDOM, reIndex, xrandCols, RANDYR,
                          tau1, tau2){
  
  function(pars, xfec, xrep, R, propU, propF, w, 
           z, zmat, matYr, muyr){
    
    fg        <- pars$fg
    ug        <- pars$ug
    umean     <- pars$umean
    uvar      <- pars$uvar
    sg        <- pars$sg
    bgFec     <- pars$bgFec
    bgRep     <- pars$bgRep
    betaYrF   <- pars$betaYrF
    betaYrR   <- pars$betaYrR
    alphaRand <- pars$alphaRand
    Arand     <- pars$Arand
    ngroup    <- length(ug)
    
    qf <- ncol(xfec)
    qr <- ncol(xrep)
    
    accept <- 0
    nspec  <- nrow(R)

    ONEF <- ONER <- ONEA <- F
    if(ncol(xfec) == 1)ONEF <- T
    if(ncol(xrep) == 1)ONER <- T
    if(length(Arand) == 1)ONEA <- T
    
    nxx <- length(obsRows)
    yg  <- log(fg)[obsRows]
    
    yeffect <- reffect <- 0
    
    xfz <- xfec[drop=F,obsRows,]
    bgf <- bgFec 
    
    w0  <- which( colSums(xfz) == 0 )  
    if(length(w0) > 0){
      xfz <- xfz[,-w0, drop=F]
      bgf <- bgf[drop=F,-w0,]
    }
    
    if(YR){                            # yeffect in mean
      yg <- yg - betaYrF[yrIndex[obsRows,'year']] 
      if(RANDYR) yg <- yg - betaYrR[yrIndex[obsRows,c('group','year')]]
    }
    if(AR){
      yg <- yg - muyr[obsRows]
    }
    
    if(RANDOM){                        
      reffect <- xfec[,xrandCols]*alphaRand[reIndex,]
      if(!ONEA)reffect <- rowSums( reffect )
      yg <- yg - reffect[obsRows]
    }
    
    if(ONEF){
      
      V <- 1/(sum(z[obsRows])/sg + .1)
      v <- sum(yg[z[obsRows] == 1])/sg
      bgf <- matrix( rnorm(1,V*v,sqrt(V)), 1)
      
    }else{
      
      XX    <- 1/sg*crossprod(xfz[z[obsRows] == 1,]) + diag(.1, qf)
      testv <- try( chol(XX) ,T)
      if( inherits(testv,'try-error') ){
        diag(XX)  <- diag(XX) + .001
        testv <- try(chol(XX, pivot=T),T)
      }
      V  <- chol2inv(testv)
      v  <- 1/sg*crossprod(xfz[z[obsRows] == 1,],yg[z[obsRows] == 1]) 
      if(is.null(betaPrior)){
        bgf  <- t( rmvnormRcpp(1, V%*%v, V) )
        ww   <- which( abs(bgf) > 50 )
        if(length(ww) > 0){
          bgf <- t(.tnormMVNmatrix( avec=t(bgf), muvec=t(bgf), smat=V,
                             lo=matrix(-20,1, length(bgf)), 
                             hi=matrix(20, 1, length(bgf)),
                             whichSample = ww) )
        }
      }else{
        diag(V) <- diag(V) + .001
        lims <- betaPrior$fec
        bgf  <- t(.tnormMVNmatrix( avec=t(V%*%v), muvec=t(V%*%v), smat=V,
                                   lo=matrix(lims[,1],1), 
                                   hi=matrix(lims[,2],1)))
      }
    }
    
    if(length(w0) > 0){      # no mature individuals
      bgFec <- bgFec*0
      bgFec[-w0,] <- bgf
      bgf <- bgFec
    }
    
    bgFec <- bgf
    
    # maturation
    if(ONER){
      
      V <- 1/(nxx/sg + .1)
      v <- sum(w[obsRows])/sg
      bgRep <- rnorm(1,V*v,sqrt(V))
      
    }else{
      V  <- solve(1/sg*crossprod(xrep[obsRows,])  + diag(.1, qr))
      v  <- 1/sg*crossprod(xrep[obsRows,],w[obsRows]) 
      if(is.null(betaPrior)){
        bgRep <-t( rmvnormRcpp(1, V%*%v, V) )
      }else{
        lims <- betaPrior$rep
        bgRep <- t(.tnormMVNmatrix( avec=matrix(bgRep,1), muvec=t(V%*%v), 
                                    smat=V,
                                    lo=matrix(lims[,1],1), 
                                    hi=matrix(lims[,2],1)))
      }
    }
    rownames(bgFec) <- colnames(xfec) 
    rownames(bgRep) <- colnames(xrep) 
    
    
    unew <- .tnorm(ngroup, minU, maxU, ug, rexp(ngroup, 1/propU) )
    names(unew) <- names(ug)
    
    if(!RANDYR){
      umean <- priorU
      uvar  <- priorVU
    }
    
    
    pnow <- .seedProb(tdata$specPlot[obsRows], tdata$year[obsRows], tdata$dcol[obsRows],
                      ug, fg[obsRows], distall, sdata[obsRowSeed,], seedNames, 
                      z[obsRows], R, SAMPR, obsYr)
    pnew <- .seedProb(tdata$specPlot[obsRows], tdata$year[obsRows], tdata$dcol[obsRows],
                      unew, fg[obsRows], distall, sdata[obsRowSeed,], seedNames, 
                      z[obsRows], R, SAMPR, obsYr) 
    
    pnow[pnow < -1e+9] <- -1e+9   # intensity parameter is zero
    pnew[pnew < -1e+9] <- -1e+9
    pnow <- sum(pnow) + sum( dnorm(ug, umean, sqrt(uvar),log=T) )
    pnew <- sum(pnew) + sum( dnorm(unew, umean, sqrt(uvar),log=T) )
    
    pdif <- pnew - pnow
    
    #############
    
  #  AA   <- sdata$active*sdata$area
  #  lnow <- .lambda(ug, fg[obsRows], z[obsRows], R, 
  #                  tdata[,c('specPlot','year','dcol')], 
  #                  sdata[,c('drow','year')], 
  #                  obsRows, obsRowSeed, obsYr, distall, AA=AA, SAMPR, PERAREA=F)
  #  lnew <- .lambda(unew, fg[obsRows], z[obsRows], R, 
  #                  tdata[,c('specPlot','year','dcol')], 
  #                  sdata[,c('drow','year')],
 #                  obsRows, obsRowSeed, obsYr, distall, AA=AA, SAMPR, PERAREA=F)
  #  lnow <- lnow + 1e-5
  #  lnew <- lnew + 1e-5
  #  pdif2 <- as.matrix(sdata[,seedNames])*(log(lnew) - log(lnow)) - 
  #          AA*(lnew - lnow)
  #  qdif <- 1/2/uvar*sum(unew^2 - ug^2 + 2*umean*(ug - unew))
  #  pdif2 <- sum(pdif2) + qdif
    ############
    
    a <- exp(pdif)
    if(is.finite(a)){
      if( runif(1,0,1) < a){
        ug   <- unew
        propU <- propU*2
      }else{
        propU <- propU*.9
      }
    }
    
    if(RANDYR){  # prior
      
      V <- 1/(ngroup/uvar + 1/priorVU)
      v <- 1/uvar*sum(ug) + priorU/priorVU
      umean <- .tnorm(1, minU, maxU, V*v, sqrt(V)) 
      uvar  <- 1/rgamma(1, tau1 + ngroup/2, tau2 + .5*sum( (ug - umean)^2 ) )
    }
    
    list(ug = ug, umean = umean, uvar = uvar, 
         bgFec = bgFec, bgRep = bgRep, propU = propU)
  }
}

.lambda <- function(ug, ff, zz, R, tdata, sdata, obsRows, obsRowSeed, obsYr, 
                    distall, AA, SAMPR, PERAREA=F, SPECPRED = F){
  
  # PERAREA - from per-trap to per-area
  # if length(AA === 1) then it must equal 1
  # SPECPRED - predict species rather than seed types
  #AA   - active*area, vector or one number
  
  nf <- length(ff)
  fz <- ff*zz
  
  if( SAMPR | length(R) > 1 ){
    if(SPECPRED){
      fk <- matrix(0, nf, nrow(R))
      jj <- match(as.character(tdata$specPlot[obsRows]),rownames(R)) 
      fk[ cbind(1:nf,jj) ] <- fz
      fz <- fk
      colnames(fz) <- rownames(R)
    }else{
      fz <- matrix(fz,length(ff),ncol=ncol(R))*
        R[drop=F,as.character(tdata$specPlot)[obsRows],] 
    }
  }else{
    fz <- matrix(fz,ncol=1)
  }
  
  uvec <- ug[ attr(distall,'group') ]
  dmat <- t(uvec/pi/(uvec + t(distall)^2)^2)
  dmat[dmat < 1e-8] <- 0
  
  lambda <- kernYrRcpp(dmat, fz, years = obsYr, 
                       seedyear = sdata$year[obsRowSeed],
                       treeyear = tdata$year[obsRows], 
                       seedrow = sdata$drow[obsRowSeed],
                       treecol = tdata$dcol[obsRows])
  if(SPECPRED){
    colnames(lambda) <- rownames(R)
    
    sname <- sort(unique(attr(R,'species')))
    ii <- rep( c(1:nrow(lambda)), ncol(lambda) )
    jj <- match(attr(R,'species'),sname)
    jj <- rep(jj, each=nrow(lambda))
    
    lambda <- .myBy(as.vector(lambda), ii, jj, fun='sum')
    colnames(lambda) <- sname
    
  }else{
    colnames(lambda) <- colnames(R)
  }
  
  if( PERAREA | length(AA) == 1 ) return(as.matrix(lambda))   # per area
  
  if( length(AA) > 1 )AA <- AA[obsRowSeed]
  
  as.matrix( lambda*matrix( AA, nrow(lambda), ncol(lambda) ) )  # per trap
}

.wrapperStates <- function( maxF, SAMPR, RANDOM, obsTimes, plotYears,
                            sdata, tdat, seedNames, last0first1, distall, 
                            YR, AR, obsRows, obsYr, predYr, obsRowSeed, 
                            ntree, years, nyr, xrandCols, reIndex, yrIndex, 
                            plag, groupByInd, RANDYR, updateProp){
  
  function(g, pars, xfec, xrep, propF, z, zmat, matYr, muyr, 
           epsilon = .00001, pHMC = .1){
    
    # tdat    - species, dcol, year, plotyr
    # pHMC    - fraction of steps that are Hamiltonian
    
    fg        <- pars$fg
    fecMin    <- pars$fecMin
    fecMax    <- pars$fecMax
    ug        <- pars$ug
    sg        <- pars$sg
    bgFec     <- pars$bgFec
    bgRep     <- pars$bgRep
    betaYrF   <- pars$betaYrF
    betaYrR   <- pars$betaYrR
    alphaRand <- pars$alphaRand
    Arand     <- pars$Arand
    R         <- pars$R
    nxx       <- length(fg)
    ngroup    <- nrow(betaYrR)
    bottom    <- -15
    
    accept  <- 0
      
    nspec  <- nrow(R)
    
    ONEF <- ONER <- ONEA <- F
    if(ncol(xfec) == 1)ONEF <- T
    if(ncol(xrep) == 1)ONER <- T
    if(length(Arand) == 1)ONEA <- T
    
    if(AR)lindex <- 1:plag
    
    fg[fg > fecMax] <- fecMax[fg > fecMax]
    fg[fg < fecMin] <- fecMin[fg < fecMin]
    propF[propF > .1*fg] <- .1*fg[propF > .1*fg]
    
    yg <- log(fg)
    yg[yg < bottom] <- bottom
    
    yeffect <- reffect <- 0
    
    xfz <- xfec
    bgf <- bgFec 
    
    w0  <- which( colSums(xfz) == 0 )  
    if(length(w0) > 0){
      xfz <- xfz[,-w0]
      bgf <- bgf[drop=F,-w0,]
    }
    
    if(YR & !AR){                            # yeffect in mean
      yeffect <- betaYrF[yrIndex[,'year']]
      if(RANDYR)yeffect <- yeffect + betaYrR[yrIndex[,c('group','year')]]
    }
    
    if(RANDOM){                        
      reffect <- xfec[,xrandCols]*alphaRand[reIndex,]
      if(!ONEA)reffect <- rowSums( reffect )
  #    yg <- yg - reffect
    }
    
    lmu           <- xfec%*%bgFec
    if(YR)lmu     <- lmu + yeffect
    if(RANDOM)lmu <- lmu + reffect
    nall <- length(fg)
    
    tt <- rbinom(1,1,pHMC)
    
    if(tt == 1){
      
      fg[fg < 1e-6] <- 1e-6
      tmp <- HMC(ff = fg[obsRows], fecMin[obsRows], fecMax[obsRows], 
                 ep = epsilon[obsRows],   
                 L = 4, tdat[obsRows,], sdata[obsRowSeed,], ug, 
                 mu = lmu[obsRows], sg, zz = z[obsRows], R, SAMPR, 
                 distance = distall, obsRows, obsYr, seedNames)
      fg[obsRows] <- tmp$fg
      epsilon[obsRows] <- tmp$epsilon
      
      fg[fg > fecMax] <- fecMax[fg > fecMax]
      fg[fg < fecMin] <- fecMin[fg < fecMin]
      
      fg[fg < 1e-6] <- 1e-6
      
      return( list(fg = fg, fecMin = fecMin, fecMax = fecMax,
                   z = z, zmat = zmat, matYr = matYr, propF = propF,
                   accept = tmp$accept, epsilon = epsilon) )
    }
    
    # maturation
    pr  <- pnorm( xrep%*%bgRep )
    iss <- sdata$plotyr
    ii  <- rep(iss, length(seedNames))
    
    if( !AR ){               
      
      tmp      <- .propZ(zmat, last0first1, matYr)
      zmatNew  <- tmp$zmat
      znew     <- zmatNew[ yrIndex[,c('dcol','year')] ] 
      matYrNew <- tmp$matYr 
      mnow <- z*log(pr) + (1 - z)*log(1 - pr)
      mnew <- znew*log(pr) + (1 - znew)*log(1 - pr)
      
      dz <- znew - z
      lo <- fecMin
      hi <- fecMax
      lo[dz == 1] <- 1
      hi[dz == 1] <- maxF
      lo[dz == -1] <- 1e-6
      hi[dz == -1] <- 1
      
      fnew <- .tnorm(nall, lo, hi, fg, rexp(nxx,1/propF), .001)
      fnew[fnew < 1e-6] <- 1e-6
      ynew <- log(fnew)
      
      # fecundity
      bnow <- bnew <- fg*0
      bnow[z == 1] <- dnorm(yg[z == 1], lmu[z == 1], sqrt(sg), log=T) 
      bnew[znew == 1] <- dnorm(ynew[znew == 1], lmu[znew == 1], sqrt(sg), log=T) 
      w0 <- which(z == 0 | znew == 0)
      bnew[w0] <- bnow[w0] <- 0
      
      # seed
      pnow <- pnew <- matrix(0,nrow(sdata),length(seedNames))
      pnow[obsRowSeed,] <- .seedProb(tdat$specPlot[obsRows], tdat$year[obsRows], 
                                    tdat$dcol[obsRows], ug, fg[obsRows], distall, 
                                    sdata[obsRowSeed,], seedNames, z[obsRows], R, 
                                    SAMPR, obsYr)
      pnew[obsRowSeed,] <- .seedProb(tdat$specPlot[obsRows], tdat$year[obsRows], 
                                    tdat$dcol[obsRows], ug, fnew[obsRows], distall, 
                                    sdata[obsRowSeed,], seedNames, z[obsRows], R, 
                                    SAMPR, obsYr)
      pnow[pnow < -1e+8] <- -1e+8   # intensity parameter is zero
      pnew[pnew < -1e+8] <- -1e+8
      
      # by plot-yr
      ii <- tdat$plotyr     #consider bnow[znow == 0] = 0
      sm <- matrix(0, max(c(tdat$plotyr, sdata$plotyr)), 1)
      
      mdif <- .myBy(mnew - mnow, ii, ii*0+1,summat = sm, fun='sum')
  #    mnew <- .myBy(mnew, ii, ii*0+1,summat = sm, fun='sum')
      
      bdif <- .myBy(bnew - bnow, i = ii, j = ii*0 + 1, summat = sm*0, fun='sum')
    #  bnew <- .myBy(bnew, i = ii, j = ii*0 + 1, summat = sm*0, fun='sum')
      
      ii <- sdata$plotyr
      ii <- rep(ii, length(seedNames))
      
      pdif <- .myBy(as.vector(pnew - pnow), ii, ii*0 + 1, summat = sm*0, fun='sum')
   #   pnew <- .myBy(as.vector(pnew), ii, ii*0 + 1, summat = sm*0, fun='sum')
      
      
  #    p <- pnew - pnow
  #    b <- bnew - bnow
  #    m <- mnew - mnow
      
      a  <- exp( pdif + bdif + mdif )        
      az  <- runif(length(a),0,1)
      aw  <- which(az < a)
      
      accept <- accept + length(aw)
      
      propF <- propF/2
      
      if(length(aw) > 0){
        
        wa <- which(tdat$plotyr %in% aw)
        
        yg[ wa ] <- ynew[ wa ]
        z[ wa ]  <- znew[ wa ]
        fecMin[wa] <- lo[wa]
        fecMax[wa] <- hi[wa]
        zmat[ yrIndex[,c('dcol','year')] ] <- z  
        
        tmp <- apply(zmat,1,which.max)
        tmp[rowSums(zmat) == 0] <- ncol(zmat)
        matYr <- tmp
        
        
        if(g %in% updateProp){
          propF[wa] <- propF[wa]*5
  #        propF[propF > maxProp] <- maxProp
        }
        
      }else{
        propF <- propF*.5
      }
      
    }else{         # AR
      
      #independence sampler
      
      p2s <- rep(0,length(plotYears))
      
      yy <- mu <- matrix(0, ntree, nyr)
      yy[ yrIndex[,c('dcol','year')] ] <- yg
      mu[ yrIndex[,c('dcol','year')] ] <- lmu
      
      #prior for backcast based on obsRows
      
      yprior <- .myBy(yg[obsRows], yrIndex[obsRows,'dcol'], 
                      obsRows*0+1, fun='mean')
      
      for(t in 1:length(predYr)){        
        
        tii <- which(yrIndex[,'year'] == t)  # row in tdat
        yii <- yrIndex[tii,'dcol']           # location in ytmp
        nii <- length(yii)
        fmin <- log(fecMin[tii])
        fmax <- log(fecMax[tii])
        
        zt <- zmat[yii,t]
        zp <- rbinom(length(zt),1,.5)
        
        #new and current
        ln <- lo <- fmin
        hn <- hi <- fmax
   #     lo <- ln <- zt*0 + bottom
   #     hi <- hn <- zt*0 + log(maxF)
        
        wp <- last0first1[yii,'first1'] <= t
        wn <- last0first1[yii,'last0'] >= t
        if(t > 1)  wp <- wp | zmat[yii,t-1] == 1
        if(t < nyr)wn <- wn | zmat[yii,t+1] == 0
        zp[wp] <- 1
        zp[wn] <- 0
        
   #     hi[zt == 0] <- 0
   #     lo[zt == 1] <- 1e-6
        dz <- zp - zt
        ln[dz == 1] <- 0  #from zero to one
        hn[dz == 1] <- log(maxF)
        ln[dz == -1] <- -log(maxF)  #from one to zero
        hn[dz == -1] <- 0
        
        mnow <- zt*log(pr[tii]) + (1 - zt)*log(1 - pr[tii])
        mnew <- zp*log(pr[tii]) + (1 - zp)*log(1 - pr[tii])
        
        a  <- exp( mnew - mnow )        
        az  <- runif(nii,0,1)
        aw  <- which(az < a)
        
        accept <- accept + length(aw)
        
        pindex <- t - (1:plag)
        w0     <- which(pindex > 0)
        mt     <- mu[,t]
        VI     <- rep(1,ntree)       # prior
        if(t > 1){                   # m_t and VI
          pindex <- pindex[w0]
          byr <- matrix(betaYrF[w0],ntree,length(w0),byrow=T)
          if(RANDYR)byr <- byr + betaYrR[groupByInd,w0,drop=F]
          mt <- mt + rowSums( byr*yy[,pindex] )
          VI <- VI + rowSums( byr^2 )
        }
        V <- sg/VI
        
        if(t > max(obsTimes)){    ###### predict forward
          
          if(length(aw) > 0){
            z[ tii[aw] ] <- zp[ aw ]
       #     lo[aw] <- ln[aw]
       #     hi[aw] <- hn[aw]
            zmat[ yii[aw],t]   <- zp[ aw ]
            fecMin[ tii[aw] ]  <- exp(ln[aw])
            fecMax[ tii[aw] ]  <- exp(hn[aw])
          }
          yy[yii,t] <- .tnorm(nii,lo,hi,mt[yii],sqrt(sg))
          next
        }
        
        vt <- mt            # v_t
        
        for(k in 1:plag){

          tindex <- t + k - lindex
          wt  <- which(lindex != k & tindex > 0)
          
          if(length(wt) == 0)next
          
          byr1 <- matrix(betaYrF[wt],ntree,length(wt),byrow=T) 
          byr2 <- betaYrF[k]
          if(RANDYR){
            byr1 <- byr1 + betaYrR[groupByInd,wt,drop=F]
            byr2 <- byr2 + betaYrR[groupByInd,k]
          }
          
          ntl <- yy[,t+k] - mu[,t+k] - rowSums( byr1*yy[,tindex[wt]] )
          vt  <- vt + ntl*byr2
        }
        vt <- vt/sg 
        
        if(t <= plag){       #### imput backward
          
          if(length(aw) > 0){
            z[ tii[aw] ] <- zp[ aw ]
       #     lo[aw] <- ln[aw]
       #     hi[aw] <- hn[aw]
            zmat[ yii[aw],t] <- zp[ aw ]
            fecMin[ tii[aw] ]  <- exp(ln[aw])
            fecMax[ tii[aw] ]  <- exp(hn[aw])
            
          }
          V[yii]  <- 1/(1/V[yii] + 1/.5)
          vt[yii] <- vt[yii] + yprior[yii]/.5
          yy[yii,t] <- .tnorm(nii,lo,hi,(V*vt)[yii],sqrt(V[yii]))
          
          next
        }
        
        
        # seed data
        ynew <- .tnorm(nii,ln,hn,(V*vt)[yii],sqrt(V[yii])) # from conditional
        
        sii  <- which(sdata$year == years[t]) # year row in seedData
        
        tree <- tdat[tii,]
        seed <- sdata[sii,]
        
        wtt <- which(!tree$plotyr %in% seed$plotyr) 
        
        if(length(wtt) > 0){         #year before seed data, draw from conditional
          if(length(aw) > 0){
            aww <- aw[aw %in% wtt]   #year before and update
            z[ tii[aww] ] <- zp[ aww ]
      #      lo[aww] <- ln[aww]
      #      hi[aww] <- hn[aww]
            zmat[ yii[aww],t] <- zp[ aww ]
            fecMin[ tii[aww] ]  <- exp(ln[aww])
            fecMax[ tii[aww] ]  <- exp(hn[aww])
            
          }
          yy[yii[wtt],t] <- .tnorm(length(wtt),lo[wtt],hi[wtt],
                                   (V*vt)[yii[wtt]],sqrt(V[yii[wtt]]))
          wtk <- which(tree$plotyr %in% seed$plotyr)
          tii <- tii[ wtk ]
          yii <- yii[ wtk ]
          
          tree <- tdat[tii,]
      #    fecMin[tii] <- lo
      #    fecMax[tii] <- hi
          mnow  <- mnow[wtk]
          mnew  <- mnew[wtk]
          ynew  <- ynew[wtk]
          zp    <- zp[wtk]
        }
        
        pnow <- .seedProb(tree$specPlot, tree$year, tree$dcol,
                          ug, fg[tii], distall, seed,
                          seedNames, zz=z[tii], R, SAMPR, years[t])
        pnew <- .seedProb(tree$specPlot, tree$year, tree$dcol,
                          ug, exp(ynew), distall, seed,
                          seedNames, zz=z[tii], R, SAMPR, years[t])
        pnow[pnow < -1e+9] <- -1e+9   # intensity parameter is zero
        pnew[pnew < -1e+9] <- -1e+9
        
        iy <- match(seed$plotYr,plotYears)
        iy <- rep(iy,length(seedNames))
        
        pnow <- .myBy(as.vector(pnow), iy, iy*0 + 1, fun='sum')
        pnew <- .myBy(as.vector(pnew), iy, iy*0 + 1, fun='sum')
        
        pyID <- unique(iy)       # plot yr for pnew/pnow
        p2s[pyID] <- p2s[pyID] + pnew[pyID] - pnow[pyID]
        
        ip <- match(tree$plotYr,plotYears)
        sm <- matrix(0,max(ip),1)
        
        mnow <- .myBy(mnow, ip, ip*0+1,summat = sm, fun='sum')
        mnew <- .myBy(mnew, ip, ip*0+1,summat = sm, fun='sum')
        
        myID <- unique(ip)
        p2s[myID] <- p2s[myID] +  mnew[myID] - mnow[myID]
        
        a  <- exp( p2s )        
        az  <- runif(length(a),0,1)
        aw  <- which(az < a)
        
        accept <- length(aw)  # no. plot-years
        
        if(length(aw) > 0){
          wa <- which(tree$plotyr %in% aw)  #rows in tdat[tii,]
          yy[ yii[wa],t] <- ynew[wa]
          z[ tii[wa] ]  <- zp[ wa ]
          fecMin[ tii[wa] ]  <- exp(ln[wa])
          fecMax[ tii[wa] ]  <- exp(hn[wa])
    #      zmat[wa,t] <- z[ tii[wa] ]
        }
      }
      yg <- yy[ cbind(tdat$dcol, tdat$times) ]
      tmp <- apply(zmat,1,which.max)
      tmp[rowSums(zmat) == 0] <- ncol(zmat)
      matYr <- tmp
    }
    
    fg <- exp(yg)
    
    wf <- which(propF < fg/100)
    
    if(length(wf) > 0)propF[wf] <- fg[wf]/100
    
    fecMin[fecMin < 1e-6] <- 1e-6
    fecMax[fecMax < 1] <- 1
    
    fg[fg > fecMax] <- fecMax[fg > fecMax]
    fg[fg < fecMin] <- fecMin[fg < fecMin]
    fg[fg < 1e-6] <- 1e-6
    fg[fg > maxF] <- maxF
            
    list(fg = fg, fecMin = fecMin, fecMax = fecMax, 
         z = z, zmat = zmat, matYr = matYr, propF = propF, 
         epsilon = epsilon, accept = accept) 
  } 
}

.getF <- function(kern, gg ){
  
  tiny <- .0001
  
  fec <- rep(0, ncol(kern))
  
  kk <- kern
  K   <- crossprod(kk)
  K   <- K + diag(tiny*diag(K), nrow(K), nrow(K))
  fec <- solve(K)%*%crossprod(kk, gg)
  
  fec[fec < tiny] <- tiny
  fec
}

.specFormula <- function(formula, NOINTERCEPT=F){
  
  form <- paste0( as.character(formula), collapse=' ')
  form <- .replaceString(form, '~', '~ species*')
  form <- .replaceString(form, ' + ', '+ species*')
  form <- .replaceString(form, '* 1','')
  form <- .replaceString(form, '*1','')
  if(NOINTERCEPT) form <- paste(form, '-1')
  as.formula(form)
}

.getBetaPrior <- function(betaPrior, bgFec, bgRep, specNames){
  
  fecHi <- bgFec*0 + 10
  fecLo <- bgFec*0 - 10
  repHi <- bgRep*0 + 10
  repLo <- bgRep*0 - 10
  
  if('pos' %in% names(betaPrior)){
    for(j in 1:length(betaPrior$pos)){
      jn <- paste('species',specNames,sep='')
      if(betaPrior$pos[j] != 'intercept')jn <- paste(jn,':',betaPrior$pos[j],sep='')
      fecLo[rownames(fecLo) %in% jn] <- 0
      repLo[rownames(repLo) %in% jn] <- 0
    }
  }
  if('neg' %in% names(betaPrior)){
    for(j in 1:length(betaPrior$neg)){
      jn <- paste('species',specNames,sep='')
      if(betaPrior$neg[j] != 'intercept')jn <- paste(jn,':',betaPrior$neg[j],sep='')
      fecHi[rownames(fecHi) %in% jn] <- 0
      repHi[rownames(repHi) %in% jn] <- 0
    }
  }
  list(fec = cbind(fecLo, fecHi), rep = cbind(repLo, repHi) )
}

.updateBetaYr <- function(yg, z, sg, sgYr, betaYrF, betaYrR, yrIndex, yeGr,
                          RANDYR, obs){
  
  #fixed effects
  
  wz <- which(z == 1 & obs == 1)
  nk <- max(yrIndex[,'year'])              # no. years
  yk <- yrIndex[wz,c('group','year')]      # year groups, years
  G  <- length(yeGr)
  
  yfix <- yg[wz] 
  if(RANDYR)yfix <- yfix - betaYrR[yrIndex[wz,c('group','year')]]
  
  ygroup <- .myBy(yfix, yk[,2]*0+1, yk[,2], 
                  summat=matrix(0, 1, nk), fun='sum')
  ngr  <- .myBy(yfix*0+1, yk[,2]*0+1, yk[,2], 
                  summat=matrix(0, 1, nk), fun='sum')
  v <- ygroup/sg
  V <- 1/(ngr /sg + .1)
  bf <- matrix( .tnorm(length(v), -1.5, 1.5, V*v, sqrt(V)), 1, nk)
  bf <- bf - mean(bf)                                   # sum to zero
  
  if(!RANDYR)return( list(betaYrF = bf, betaYrR = bf*0, sgYr = sgYr, 
                 wfinite = 1:nk) )
  
  # random effects
  
  yfix <- yg[wz] - bf[yrIndex[wz,'year']]
  
  ygroup <- .myBy(yfix, yk[,1], yk[,2], 
                  summat=matrix(0, G, nk), fun='sum')
  ngr  <- .myBy(yfix*0+1, yk[,1], yk[,2], 
                  summat=matrix(0, G, nk), fun='sum')
  v  <- ygroup/sg 
  V  <- 1/(ngr /sg + matrix(1/sgYr, G, nk, byrow=T))
  nc <- ngr 
  nc[nc > 1] <- 1
  ns <- colSums(nc)
  nc[,ns == 1] <- 0     # no group effect if only one group
  
  br <- matrix( .tnorm(length(v), -1.5, 1.5, V*v, sqrt(V)), G, nk )
  br <- br*nc
  
  rs <- rowSums(nc)
  ws <- which(rs > 0)
  
  br[ws,] <- sweep(br[drop=F,ws,], 1, rowSums(br[drop=F,ws,])/rs[ws], '-')*nc[drop=F,ws,]
  
  sgYr <- 1/rgamma(nk, 2 + ns/2, 1 + .5* colSums(br^2))
  
  list(betaYrF = bf, betaYrR = br, sgYr = sgYr, wfinite = which(nc > 0))
}

.multivarChainNames <- function(rowNames,colNames){
  as.vector( t(outer(colNames,rowNames,paste,sep='_')) )
}

.rdirichlet <- function(pmat, n = nrow(pmat) ){
  
  # pmat - rows are parameter vectors
  # if pmat is a vector, then n is the number of random vectors
  
  if(!is.matrix(pmat)){
    pmat <- matrix(pmat,1)
    pmat <- pmat[rep(1,n),]
  }
  pmat <- matrix( rgamma(n*ncol(pmat),pmat,1), n, ncol(pmat))
  sweep(pmat, 1, rowSums(pmat), '/')
}

.updateR <- function(ug, ff, zz, SAMPR, distall, sdata, seedNames, 
                     tdat, R, priorR, priorRwt, years, posR, plots){
  
 # a <- .tnorm(length(posR), 0, 2, 1, 1)
 # b <- a/R[posR]
  
  mnew <- R
  mnew[posR] <- .tnorm(length(posR), 0, 1, R[posR], .02)
  
  mnew <- sweep(mnew, 1, rowSums(mnew,na.rm=T), '/')
  mnew[-posR] <- 0

  qnow <- 2*priorRwt*log(R)
  qnow[-posR] <- 0
  qnew <- 2*priorRwt*log(mnew)
  qnew[-posR] <- 0
  
  jj <- rep(attr(R,'plot'), ncol(R))
  qnow <- tapply(as.vector(qnow), jj, sum, na.rm=T)
  qnew <- tapply(as.vector(qnew), jj, sum, na.rm=T)
  
  tnow <- .seedProb(tdat$specPlot, tdat$year, tdat$dcol,
                    ug, ff, distall, sdata, seedNames,
                        zz, R, SAMPR, years)
  tnew <- .seedProb(tdat$specPlot, tdat$year, tdat$dcol,
                    ug, ff, distall, sdata, seedNames,
                    zz, mnew, SAMPR, years)
  tnow[!is.finite(tnow)] <- -8   # intensity parameter is zero
  tnew[!is.finite(tnew)] <- -8
  
  ii <- match(sdata$plot, plots)
  pnow <- .myBy( rowSums(tnow), ii, ii*0+1, fun='sum')[,1]
  pnew <- .myBy( rowSums(tnew), ii, ii*0+1, fun='sum')[,1]
  
  a <- exp(pnew + qnew - pnow - qnow) 
  wa <- which(a > runif(length(a),0,1))
  if(length(wa) > 0){
    R[attr(R,'plot') %in% plots[wa],] <- mnew[attr(R,'plot') %in% plots[wa],]
  }
  R
}

.distmat <- function(x1,y1,x2,y2){
    xd <- outer(x1,x2,function(x1,x2) (x1 - x2)^2)
    yd <- outer(y1,y2,function(y1,y2) (y1 - y2)^2)
    t(sqrt(xd + yd)) 
}

.updateVariance <- function(y,mu,s1=1,s2=1){
  
  u1 <- s1 + length(y)/2
  u2 <- s2 + .5*sum( (y - mu)^2 )
  1/rgamma(1,u1,u2) 
  
}

sqrtSeq <- function(maxval){ #labels for sqrt scale
  
  # maxval on sqrt scale
  
  by <- .2
  if(maxval > 1)    by <- 2
  if(maxval >= 3)   by <- 4
  if(maxval >= 4)   by <- 8
  if(maxval >= 5)   by <- 10 
  if(maxval >= 10)  by <- 20 
  if(maxval >= 15)  by <- 50
  if(maxval >= 20)  by <- 100 
  if(maxval >= 30)  by <- 200 
  if(maxval >= 50)  by <- 500
  if(maxval >= 70)  by <- 1000
  if(maxval >= 100) by <- 2000
  if(maxval >= 200) by <- 10000
  if(maxval >= 500) by <- 50000
  if(maxval >= 700) by <- 100000
  if(maxval >= 1000)by <- 200000
  if(maxval >= 1500)by <- 400000
  if(maxval >= 2000)by <- 1000000
  if(maxval >= 3000)by <- 2000000
  if(maxval >= 4000)by <- 5000000
  
  labs <- seq(0, maxval^2, by = by)
  at   <- sqrt(labs)

  list(at = at, labs = labs)
  
}

.plotObsPred <- function(obs, yMean, ySE=NULL, opt = NULL){
  
  nbin <- nPerBin <- xlimit <- ylimit <- NULL
  add <- log <- SQRT <- F
  xlabel <- 'Observed'
  ylabel <- 'Predicted'
  trans <- .4
  col <- 'black'
  bins <- NULL
  atx <- aty <- labx <- laby <- NULL
  
  for(k in 1:length(opt))assign( names(opt)[k], opt[[k]] )
  
  if(!is.null(bins))nbin <- length(bins)
  
  if(log & SQRT)stop('cannot have both log and SQRT scale')
  
  yMean <- as.matrix(yMean)
  obs   <- as.matrix(obs)
  
  if(SQRT){
    xlim <- sqrt(xlimit)
    ylim <- sqrt(ylimit)
    obs   <- as.vector(sqrt(obs))
    yMean <- as.vector(sqrt(yMean))
    if(!is.null(bins))bins <- sqrt(bins)
    xlimit <- sqrt(range(obs,na.rm=T))
    xlimit[2] <- xlimit[2]*2
    ylimit <- sqrt(range(yMean,na.rm=T))
    ylimit[2] <- 1.2*ylimit[2]
 
    maxy <- max(yMean,na.rm=T)
    maxx   <- max(obs,na.rm=T)
    maxval <- max( c(maxx, maxy) )
    
    tt   <- sqrtSeq(1.2*maxx)
    if(is.null(atx))atx   <- tt$at
    if(is.null(labx))labx <- tt$labs
    
    
    if(ylimit[2] < xlimit[2]) ylimit[2] <- xlimit[2]
    if(xlimit[2] < xlim[2])   xlimit[2] <- xlim[2]
    if(ylimit[2] < ylim[2])   ylimit[2] <- ylim[2]
    
    tt   <- sqrtSeq(1.2*ylimit[2])
    if(is.null(aty))aty   <- tt$at
    if(is.null(laby))laby <- tt$labs

  }
    
  if(is.null(xlimit))xlimit <- range(obs)
  if(is.null(ylimit) & !add){                      # can only happen if !SQRT
    if(!log){
      plot(obs,yMean,col=.getColor('black',.2),cex=.3, xlim=xlimit,
           xlab=xlabel,ylab=ylabel)
      if(log) suppressWarnings( plot(obs,yMean,col=.getColor('black',.2),cex=.3,
                                     xlim=xlimit,xlab=xlabel,ylab=ylabel,log='xy') )
    }
  }
    
  if(!is.null(ylimit)){
    if(!log & !add){
      if(!SQRT){
        plot(obs,yMean,col=.getColor('black',trans),cex=.3,
                 xlim=xlimit,xlab=xlabel,ylab=ylabel,ylim=ylimit)
      }else{
        plot(obs,yMean,col=.getColor('black',trans),cex=.3,
             xlim=xlimit,xlab=xlabel,ylab=ylabel,ylim=ylimit,
             xaxt='n',yaxt='n')
        
        axis(1, at = atx, labels = labx)
        axis(2, at = aty, labels = laby, las=2)
      }
    }
    if(log & !add) plot(obs,yMean,col=.getColor('black',trans),cex=.3,
                 xlim=xlimit,xlab=xlabel,ylab=ylabel,log='xy',ylim=ylimit)
  }
  if(!is.null(ySE)){
    ylo <- yMean - 1.96*ySE
    yhi <- yMean + 1.96*ySE
    for(i in 1:length(obs))lines(c(obs[i],obs[i]),c(ylo[i],yhi[i]),
                                 col='grey',lwd=2)
  }
  
  if( !is.null(nbin) | !is.null(nPerBin) ){
    
    if(is.null(bins)){
      nbin <- 20
      bins <- seq(min(obs,na.rm=T),max(obs,na.rm=T),length=nbin)
    }else{
      nbin <- length(bins)
    }
    
    if(!is.null(nPerBin)){
      nbb <- nPerBin/length(obs)
      nbb <- seq(0,1,by=nbb)
      if(max(nbb) < 1)nbb <- c(nbb,1)
      bins <- quantile(obs,nbb,na.rm=T)
      bins <- bins[!duplicated(bins)]
      nbin <- length(bins)
    }
    
    xxk <- findInterval(obs,bins)
    
    if(SQRT){
      opos <- obs[obs > 0]
      qq <- seq(0, 1, length=15)
      bins <- quantile(opos, qq)
      bins <- bins[!duplicated(bins)]
      
  #    bins <- bins[-2]
      
      nbin <- length(bins)
      
      xxk <- findInterval(obs,bins)
      xxk[xxk == max(xxk)] <- max(xxk) - 1
    }
    xxk[xxk == nbin] <- nbin - 1
    
    wide <- diff(bins)/2
    db   <- 1
    for(k in 2:(nbin-1)){
      
      qk <- which(is.finite(yMean) & xxk == k)
      q  <- quantile(yMean[qk],c(.5,.025,.158,.841,.975),na.rm=T)
      
      if(!is.finite(q[1]))next
      if(q[1] == q[2])next
      
      ym <- mean(yMean[qk])
      xx <- mean(bins[k:(k+1)])
      rwide <- wide[k]
      
      if(k > 1)db <- bins[k] - bins[k-1]
      
      if( xx > (bins[k] + db) ){
        xx <- bins[k] + db
        rwide <- wide[ max(c(1,k-1)) ]
      }
      
      suppressWarnings(
        arrows(xx, q[2], xx, q[5], lwd=2, angle=90, code=3, col=.getColor(col,.8),
               length=.05)
      )
      lines(c(xx-.5*rwide,xx+.5*rwide),q[c(1,1)],lwd=2, 
            col=.getColor(col,.8))
      rect(xx-.4*rwide,q[3],xx+.4*rwide,q[4], col=.getColor(col,.5))
    }
  }
  invisible( list(atx = atx, labx = labx, aty = aty, laby = laby) )
}

.getKern <- function(u,dij){
  
  uvec <- u[ attr(dij,'group') ]
  kk <- t(uvec/pi/(uvec + t(dij)^2)^2)
  
 # kk <- u/pi/(u + dij^2)^2
  kk[is.na(kk)] <- 0
  kk
}

.mapSpec <- function(x, y, z, mapx=range(x), mapy=range(y), scale=0,
                     add=F, sym='circles',
                     colVec=rep(1,length(x)), fill=F){
  
  fillCol <- NA
  if(is.logical(fill))fillCol <- colVec
  if(is.character(fill))fillCol <- fill
  
  opin <- par()$pin
  
  if(scale > 0).mapSetup(mapx,mapy,scale)
  if(!add){
    plot(NA, xlim=mapx, ylim=mapy, axes = F, xlab='', ylab='')
    Axis(side=1, labels=FALSE)
    Axis(side=2, labels=FALSE)
    add <- T
  }
  
  if(sym == 'circles'){

    symbols(x,y,circles=z/10,inches=F,
                              xlim=mapx,ylim=mapy,fg=colVec,bg=fillCol,
                              lwd=2,add=add)
  }
  if(sym == 'squares'){
    symbols(x,y,squares=z/10,inches=F,
                              xlim=mapx,ylim=mapy,fg=colVec,bg=fillCol,
                              lwd=2,add=add)
  }
  par(pin = opin)
}

scaleBar = function(label, value = 1, fromLeft = .5, yadj = .1, 
                    lwd = 3, cex = 1) {
  
  xl <- par("usr")[1:2]
  yl <- par("usr")[3:4]
  
  xm <- xl[1] + fromLeft*diff(xl)
  x1 <- xm - value/2
  x2 <- xm + value/2
  
  y  <- yl[1] + .05*diff(yl)
  ym <- y + yadj*diff(yl)
    
  lines(c(x1,x2),c(y,y), lwd=lwd + 2, col='white')
  lines(c(x1,x2),c(y,y), lwd=lwd)
  
  lab <- paste(value, label)
  text(xm, ym, lab, cex = cex)
}

.mapSetup<- function(xlim,ylim,scale){  #scale is m per inch

  px   <- diff(xlim)/scale
  py   <- diff(ylim)/scale
  pin  <- c(px, py)
  par(pin = pin)
  invisible( pin )
}

.getColor <- function(col,trans){
  
  # trans - transparency fraction [0, 1]
  
  tmp <- col2rgb(col)
  rgb(tmp[1,], tmp[2,], tmp[3,], maxColorValue = 255, 
      alpha = 255*trans, names = paste(col,trans,sep='_'))
}

.interp <- function(y,INCREASING=F,minVal=-Inf,maxVal=Inf,defaultValue=NULL,
                   tinySlope=NULL){  #interpolate vector x
  
  if(is.null(defaultValue))defaultValue <- NA
  
  tiny <- .00001
  if(!is.null(tinySlope))tiny <- tinySlope
  
  y[y < minVal] <- minVal
  y[y > maxVal] <- maxVal
  
  n  <- length(y)
  wi <- which(is.finite(y))
  
  if(length(wi) == 0)return(rep(defaultValue,n))
  if(length(wi) == 1)ss <- tiny
  
  xx  <- c(1:n)
  z  <- y
  
  if(wi[1] != 1) wi <- c(1,wi)
  if(max(wi) < n)wi <- c(wi,n)
  
  ss <- diff(z[wi])/diff(xx[wi])
  
  ss[is.na(ss)] <- 0
  
  if(length(ss) > 1){
    if(length(ss) > 2)ss[1] <- ss[2]
    ss[length(ss)] <- ss[length(ss)-1]
  }
  if(INCREASING)ss[ss < tiny] <- tiny
  
  if(is.na(y[1]))  z[1] <- z[wi[2]] - xx[wi[2]]*ss[1]
  if(z[1] < minVal)z[1] <- minVal
  if(z[1] > maxVal)z[1] <- maxVal
  
  for(k in 2:length(wi)){
    
    ki <- c(wi[k-1]:wi[k])
    yk <- z[wi[k-1]] + (xx[ki] - xx[wi[k-1]])*ss[k-1]
    yk[yk < minVal] <- minVal
    yk[yk > maxVal] <- maxVal
    z[ki] <- yk
  }
  z
}

.interpRows <- function(x, startIndex=rep(1,nrow(x)), endIndex=rep(ncol(x),nrow(x)),
                       INCREASING=F, minVal=-Inf, maxVal=Inf,
                       defaultValue=NULL,tinySlope=.001){  
  #interpolate rows of x subject to increasing
  
  nn  <- nrow(x)
  p  <- ncol(x)
  xx <- c(1:p)
  
  if(length(minVal) == 1)minVal <- rep(minVal,nn)
  if(length(maxVal) == 1)maxVal <- rep(maxVal,nn)
  
  ni   <- rep(NA,nn)
  flag <- numeric(0)
  
  z <- x
  
  for(i in 1:nn){
    if(startIndex[i] == endIndex[i]){
      z[i,-startIndex[i]] <- NA
      next
    }
    z[i,startIndex[i]:endIndex[i]] <- .interp(x[i,startIndex[i]:endIndex[i]],
                                             INCREASING,minVal[i],maxVal[i],
                                             defaultValue,tinySlope)
  }
  
  z
}

.shadeInterval <- function(xvalues,loHi,col='grey',PLOT=T,add=T,
                           xlab=' ',ylab=' ', xlim = NULL, ylim = NULL, 
                           LOG=F, trans = .5){
  
  #draw shaded interval
  
  tmp <- smooth.na(xvalues,loHi)

  xvalues <- tmp[,1]
  loHi <- tmp[,-1]
  
  xbound <- c(xvalues,rev(xvalues))
  ybound <- c(loHi[,1],rev(loHi[,2]))
  if(is.null(ylim))ylim <- range(as.numeric(loHi))
  if(is.null(xlim))xlim <- range(xvalues)
  
  if(!add){
    if(!LOG)plot(NULL, xlim = xlim, ylim=ylim, 
                 xlab=xlab, ylab=ylab)
    if(LOG)suppressWarnings( plot(NULL,  xlim = xlim, ylim=ylim, 
                xlab=xlab, ylab=ylab, log='y') )
  }
 
  
  if(PLOT)polygon(xbound,ybound, border=NA,col=.getColor(col, trans))
  
  invisible(cbind(xbound,ybound))
  
}

smooth.na <- function(x,y){   
  
  #remove missing values
  #x is the index
  #y is a matrix with rows indexed by x
  
  if(!is.matrix(y))y <- matrix(y,ncol=1)
  
  wy <- which(!is.finite(y),arr.ind =T)
  if(length(wy) == 0)return(cbind(x,y))
  wy <- unique(wy[,1])
  ynew <- y[-wy,]
  xnew <- x[-wy]
  
  return(cbind(xnew,ynew))
}

.appendMatrix <- function(m1,m2,fill=NA,SORT=F,asNumbers=F){  
  
  # matches matrices by column names
  # asNumbers: if column heads are numbers and SORT, then sort numerically
  
  if(length(m1) == 0){
    if(is.matrix(m2)){
      m3 <- m2
    } else {
      m3 <- matrix(m2,nrow=1)
    }
    if( !is.null(names(m2)) )colnames(m3) <- names(m2)
    return(m3)
  }
  if(length(m2) == 0){
    if(!is.matrix(m1))m1 <- matrix(m1,nrow=1)
    return(m1)
  }
  if( is.vector(m1) | (length(m1) > 0 & !is.matrix(m1)) ){
    nn <- names(m1)
    if(is.null(nn))message('cannot append matrix without names')
    m1 <- matrix(m1,1)
    colnames(m1) <- nn
  }  
  if( is.vector(m2) | (length(m2) > 0 & !is.matrix(m2)) ){
    nn <- names(m2)
    if(is.null(nn))message('cannot append matrix without names')
    m2 <- matrix(m2,1)
    colnames(m2) <- nn
  }
  
  c1 <- colnames(m1)
  c2 <- colnames(m2)
  r1 <- rownames(m1)
  r2 <- rownames(m2)
  n1 <- nrow(m1)
  n2 <- nrow(m2)
  
  allc <-  unique( c(c1,c2) ) 
  if(SORT & !asNumbers)allc <- sort(allc)
  if(SORT & asNumbers){
    ac <- as.numeric(allc)
    allc <- as.character( sort(ac) )
  }
  
  nr <- n1 + n2
  nc <- length(allc)
  
  if(is.null(r1))r1 <- paste('r',c(1:n1),sep='-')
  if(is.null(r2))r2 <- paste('r',c((n1+1):nr),sep='-')
  new <- c(r1,r2)
  
  mat1 <- match(c1,allc)
  mat2 <- match(c2,allc)
  
  out <- matrix(fill,nr,nc)
  colnames(out) <- allc
  rownames(out) <- new
  
  out[1:n1,mat1] <- m1
  out[(n1+1):nr,mat2] <- m2
  out
}

.myBoxPlot <- function(mat, tnam, snames, specColor, label){
  
  # tnam is columns of mat, with values of snames used to match specColor
  
  ord <- order(colMeans(mat),decreasing=F)
  mat  <- mat[,ord]
  tnam <- tnam[ord]
  bb   <- specColor[ match(tnam, snames) ]
  ry   <- range(mat)
  ymin <- min(mat) - diff(ry)*.15
  ymax <- max(mat) + diff(ry)*.15
  bx   <- .getColor(bb,.4)
  
  tmp <- .boxplotQuant( mat,xaxt='n',outline=F,ylim=c(ymin,ymax),
                        col=bx, border=bb, xaxt='n',lty=1)
  abline(h=0,lwd=2,col='grey',lty=2)
  
  dy <- .05*diff(par()$yaxp[1:2])
  
  cext <- .fitText2Fig(tnam,fraction=1)
  text((1:length(ord)) - .1,dy + tmp$stats[5,],tnam,srt=70,pos=4,
       col=bb, cex=cext)
  
  pl    <- par('usr')
  xtext <- pl[1]
  ytext <- pl[3] + diff(pl[3:4])*.85
  .plotLabel(label,location='topleft', cex=1.0)
}

.boxplotQuant <- function( xx, ..., boxfill=NULL ){
  
  tmp <- boxplot( xx, ..., plot=F)
  ss  <- apply( xx, 2, quantile, pnorm(c(-1.96,-1,0,1,1.96)), na.rm=T ) 
  tmp$stats <- ss
  
  pars <- list(...)
  if( 'col' %in% names(pars) )boxfill <- pars$col
  
  bxp( tmp, ..., boxfill = boxfill )
  
  invisible( tmp )
}

.fitText2Fig <- function(xx, width=T, fraction=1, cex.max=1){
  
  # returns cex to fit xx within fraction of the current plotting device
  # width - horizontal labels stacked vertically
  #!width - vertical labels plotted horizontally
  
  px <- par('pin')[1]
  py <- par('pin')[2]
  cl <- max( strwidth(xx, units='inches') )
  ch <- strheight(xx, units='inches')[1]*length(xx)  # ht of stacked vector
  
  if(width){              #horizontal labels stacked vertically
    xf <- fraction*px/cl
    yf <- fraction*py/ch
  } else {                #vertical labels plotted horizontally
    xf <- fraction*px/ch
    yf <- fraction*py/cl
  }
  
  cexx <- min(c(xf,yf))
  if(cexx > cex.max)cexx <- cex.max
  cexx
}

deltaScore <- function(mu, vr, scale, score){
  
#  mu <- colMeans(mu/scale, na.rm=T)
#  vr <- colMeans(vr/scale, na.rm=T)
  
  mu <- matrix( mu[nrow(mu),], nrow(mu), ncol(mu), byrow=T)
  vr <- matrix( vr[nrow(vr),], nrow(mu), ncol(mu), byrow=T)
  
  muScore <- log(mu) + 1/2*( log(scale) - log(vr))
  score - muScore
}



