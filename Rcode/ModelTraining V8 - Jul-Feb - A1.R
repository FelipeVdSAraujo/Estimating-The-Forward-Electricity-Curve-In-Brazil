


## Local settings
set.seed(2020)

windowsFonts(A = windowsFont("Times New Roman")) # Set plot font

## Global Functions
createECPGObject <- function(dataset, discountVector, forwardPrice, buyerParam, sellerParam) {
  VaR0 = apply(dataset, 1, quantile, 1)
  VaR1 = apply(dataset, 1, quantile, 0)
  meanPLD = apply(dataset, 1, mean)
  
  ECPGObject = list(
    dataset=dataset,
    discountVector=discountVector,
    forwardPrice=forwardPrice,
    buyerParam=buyerParam,
    sellerParam=sellerParam,
    VaR0=VaR0,
    VaR1=VaR1,
    meanPLD=meanPLD
  )
  return(ECPGObject)
}

ECPG.phi <- function(parameters, dataset, auxList, buyer=TRUE) {
  if (any(parameters < 0) || any(parameters > 1)) return(NA)
  sumLambda = parameters['Lambda1'] + parameters['Lambda2']
  if (sumLambda > 1) return(NA)
  Lambda0 = 1 - sumLambda
  if (Lambda0 < 0) return(NA)
  
  if (buyer) {
    mean.pld = auxList$meanPLD
    VaR0 = auxList$VaR0
    VaR1 = auxList$VaR1
  } else {
    dataset = -dataset
    mean.pld = -auxList$meanPLD
    VaR0 = auxList$VaR1
    VaR1 = auxList$VaR0
  }
  
  VaRAlpha1 = apply(dataset, 1, quantile, (1-parameters['Alpha1']))
  VaRAlpha2 = apply(dataset, 1, quantile, (1-parameters['Alpha2']))
  
  CVaR1 = CVaR2 = NA
  
  CVaR1.data = dataset
  CVaR1.data[dataset > VaRAlpha1] <- NA
  
  CVaR2.data = dataset
  CVaR2.data[dataset > VaRAlpha2] <- NA
  
  CVaR1 = apply(CVaR1.data, 1, mean, na.rm=TRUE)
  CVaR2 = apply(CVaR2.data, 1, mean, na.rm=TRUE)
  
  Sum1 = parameters['Lambda1']*VaRAlpha1 + parameters['Lambda2']*VaRAlpha2
  Sum2 = parameters['Lambda1']*(CVaR1 - VaRAlpha1) + parameters['Lambda2']*(CVaR2 - VaRAlpha2)
  
  Limit1 = Lambda0*VaR0 + Sum1
  Limit2 = Lambda0*VaRAlpha1 + Sum1
  Limit3 = Lambda0*VaRAlpha2 + Sum1
  Limit4 = Lambda0*VaR1 + Sum1
  
  ECPG = Lambda0*mean.pld + parameters['Lambda1']*CVaR1 + parameters['Lambda2']*CVaR2
  
  nIter = nrow(dataset)
  Lambda1 = rep(parameters['Lambda1'], nIter)
  Lambda2 = rep(parameters['Lambda2'], nIter)
  
  if (any(ECPG>=Limit2)) Lambda1[ECPG>=Limit2] = 0
  if (any(ECPG>=Limit3)) Lambda2[ECPG>=Limit3] = 0
  
  Q = Lambda0 + Lambda1/(1-parameters['Alpha1']) + Lambda2/(1-parameters['Alpha2'])
  
  A = Sum2 + Lambda0*mean.pld + 
    Lambda1*VaRAlpha1/(1-parameters['Alpha1']) +
    Lambda2*VaRAlpha2/(1-parameters['Alpha2'])
  
  Eq = A/Q
  
  if (any(is.na(Eq))) return(NA)
  return(Eq)
}

yearly.mean <- function(x) {
  len = length(x)
  if (len/12 > len %/%12) stop("Error in yearly.mean function: Length of data is not divisible by 12!")
  
  mean.vec=NULL
  
  iter = len/12
  for (i in 1:iter) {
    mean.vec[i] = mean(x[(1+(i-1)*12):(i*12)])
  }
  return(mean.vec)
}

ECPG.equilibrium <- function (ECPGobject) {
  auxList = list(
    VaR0 = ECPGobject$VaR0,
    VaR1 = ECPGobject$VaR1,
    meanPLD = ECPGobject$meanPLD
  )
  
  buyer.phi = ECPG.phi(ECPGobject$buyerParam, ECPGobject$dataset, auxList)
  seller.phi = -ECPG.phi(ECPGobject$sellerParam, ECPGobject$dataset, auxList, buyer=FALSE)
  
  equilibrium.prices = (buyer.phi + seller.phi)/2
  
  adjusted.prices = equilibrium.prices/ECPGobject$discountVector
  
  averagePrice = yearly.mean(adjusted.prices)
  if (length(averagePrice) != length(ECPGobject$forwardPrice)) stop("Error in ECPG.equilibrium function: Length of fowardPrice vector is different from the length of data divided by 12!")
  
  ECPGobject$buyerPhi = buyer.phi
  ECPGobject$sellerPhi = seller.phi
  ECPGobject$equilibriumPrices = equilibrium.prices
  ECPGobject$adjustedPrices=adjusted.prices
  ECPGobject$averagePrice=averagePrice
  
  return(ECPGobject)
}

optim.foo <- function(parameters, ECPGobject) {
  param_b = parameters[1:4]
  param_s = parameters[5:8]
  names(param_b) = names(param_s) = c('Lambda1', 'Lambda2', 'Alpha1', 'Alpha2')
  
  ECPGobject$buyerParam = param_b
  ECPGobject$sellerParam = param_s
  
  averagePrice = ECPG.equilibrium(ECPGobject)$averagePrice
  
  return(mean((averagePrice - ECPGobject$forwardPrice)^2))
}

getDiscountVector <- function(index, modifier, anualRate) {
  yearAdj = (11+modifier*12)
  monthlyRate = (1+anualRate)^(1/12)-1
  completeDiscountVector <- rep(NA, yearAdj)
  for (i in 1:yearAdj){
    completeDiscountVector[i] <- (1+monthlyRate)^(i-1)
  }
  
  getDiscountBounds <- function(index, modifier, completeDiscountVector) {
    lower.limit = 1+12*modifier-index
    upper.limit = 12+12*modifier-index
    
    return(completeDiscountVector[lower.limit:upper.limit])
  }
  
  discountVector = as.numeric(sapply(index, getDiscountBounds, modifier=modifier, completeDiscountVector=completeDiscountVector))
  return(discountVector)
}

series.prettyPlot <- function(ECPGobject, ylim=c(0,350), lwd = 2, legendPos="bottomright", inset=0.05, cex=1, legCex=0.7, horiz=FALSE) {
  
  averagePrice = NULL
  forwardPrice = NULL
  
  iter = length(ECPGobject$forwardPrice)
  for (i in 1:iter){
    forwardPrice[(1+(i-1)*12):(12*i)] = ECPGobject$forwardPrice[i]
    averagePrice[(1+(i-1)*12):(12*i)] = ECPGobject$averagePrice[i]
  }
  
  plot.matrix = cbind(ECPGobject$meanPLD, forwardPrice, ECPGobject$adjustedPrices, averagePrice)
  colnames(plot.matrix) = c("Average PLD", "Forward Price", "Model Price", "Model Average")
  rownames(plot.matrix) = rownames(ECPGobject$dataset)
  
  Nseries = ncol(plot.matrix)
  color.pallette = c("#56B4E9", "#F0E442", "#999999", "#E69F00") # Selected for easy visualization for colorblind
  line.types = c(2, 1, 2, 1)
  
  Nrows = 1:nrow(plot.matrix)
  labels=row.names(plot.matrix)
  
  par(family="serif", cex=cex)
  matplot(Nrows, plot.matrix, type='l', xlab="", ylab='Electricity forward price(R$/MWh)', ylim=ylim, lty=line.types, lwd=lwd, col=color.pallette, axes=F)
  axis(2)
  axis(side=1,at=Nrows,labels=labels)
  grid (NULL,NULL, lty = 6, col = "grey")
  legend(x=legendPos, inset=inset, legend=colnames(plot.matrix),col=color.pallette, lty=line.types, cex=legCex, horiz=horiz, lwd=lwd)
}

#### MODEL TRAINING WITH A1 DATA FROM FIRST EIGHT MONTHS ####

## Set Variables
buyer.initial.param = c(Lambda1 = 0, Lambda2 = 0, Alpha1 = 0.5, Alpha2 = 0.95) # Initial Parameters
seller.initial.param = c(Lambda1 = 0, Lambda2 = 0, Alpha1 = 0.5, Alpha2 = 0.95) # Initial Parameter

anualRate = 0.05

data.index = rbind(
  c(7, 2019),
  c(8, 2019),
  c(9, 2019),
  c(10, 2019),
  c(11, 2019),
  c(12, 2019),
  c(1, 2020),
  c(2, 2020)
)

seriesA = 1

## Load Data
forwardPrice.data = read.csv(file='../Data/ForwardPrices.csv')

data.abbrev = paste(month.abb[data.index[, 1]], data.index[,2], sep="_")

price.data = list()
for (n in 1:8) {
  price.data[[data.abbrev[n]]] <- t(read.csv(file=paste('../Data/Data', data.abbrev[n], paste0('A', seriesA, '.csv'), sep="_"), col.names=paste(month.abb, data.index[n,2]+seriesA, sep="-")))
}

## Allocate and verify data

training.discountVector = getDiscountVector(data.index[,1], seriesA, anualRate)

len = length(training.discountVector)
if (len/12 > len %/%12) stop("Error in discountVector allocation: Length of data is not divisible by 12!")

training.forwardPrice = as.numeric(forwardPrice.data[seriesA, 2:9])
if(len/12 != length(training.forwardPrice)) stop("Error in forwardPrice: Length of data is not compatible with discountVector!")


training.data = Reduce(rbind, price.data)
if (nrow(training.data) != len) stop("Error in trainingDataA1: the number of rows is not compatible with the discountVector!")

training.ECPG_object = createECPGObject(dataset = training.data, buyerParam = buyer.initial.param, sellerParam = seller.initial.param, discountVector = training.discountVector, forwardPrice = training.forwardPrice)

## Test functions with initial data (optimization will fail if the result from initial parameters is NA)

training.initialResults = ECPG.equilibrium(training.ECPG_object)

## Print results onscreen
print(training.initialResults$adjustedPrices)

## Run the optimization - uses proc.time to measure the execution time
ptm = proc.time()
optim.results = optim(c(buyer.initial.param, seller.initial.param), optim.foo, ECPGobject=training.ECPG_object, method="SANN")
timeSpent = proc.time() - ptm

## Get ECPG parameters from optimization results

training.ECPG_object$buyerParam = optim.results$par[1:4]
training.ECPG_object$sellerParam = optim.results$par[5:8]

## Calculate results from optimized parameters

training.ECPG_object = ECPG.equilibrium(training.ECPG_object)

## Print results onscreen
training.ECPG_object$buyerParam
training.ECPG_object$sellerParam

## Plot price series
#series.prettyPlot(training.ECPG_object, cex=1.2, legendPos="bottom", horiz = TRUE)

