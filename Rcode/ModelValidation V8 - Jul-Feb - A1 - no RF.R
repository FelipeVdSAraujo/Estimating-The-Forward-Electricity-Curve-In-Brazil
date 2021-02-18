


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

ECPG.equilibrium <- function (ECPGobject) {
  auxList = list(
    VaR0 = ECPGobject$VaR0,
    VaR1 = ECPGobject$VaR1,
    meanPLD = ECPGobject$meanPLD
  )
  
  buyer.phi = ECPG.phi(ECPGobject$buyerParam, ECPGobject$dataset, auxList)
  seller.phi = -ECPG.phi(ECPGobject$sellerParam, ECPGobject$dataset, auxList, buyer=FALSE)
  
  equilibrium.prices = (buyer.phi + seller.phi)/2
  
  adjusted.prices = equilibrium.prices/ECPGobject$discount
  
  averagePrice = mean(adjusted.prices)
  
  ECPGobject$buyerPhi = buyer.phi
  ECPGobject$sellerPhi = seller.phi
  ECPGobject$equilibriumPrices = equilibrium.prices
  ECPGobject$adjustedPrices=adjusted.prices
  ECPGobject$averagePrice=averagePrice
  
  return(ECPGobject)
}

getCompleteDiscountVector <- function(index, modifier, anualRate) {
  yearAdj = (11+modifier*12)
  monthlyRate = (1+anualRate)^(1/12)-1
  completeDiscountVector <- rep(NA, yearAdj)
  for (i in 1:yearAdj){
    completeDiscountVector[i] <- (1+monthlyRate)^(i-1)
  }
  
  return(completeDiscountVector)
}

getDiscountBounds <- function(index, modifier, completeDiscountVector) {
  lower.limit = 1+12*modifier-index
  upper.limit = 12+12*modifier-index
  
  return(completeDiscountVector[lower.limit:upper.limit])
}



series.prettyPlot <- function(ECPGobject, ylim=c(0,350), lwd = 2, legendPos="bottomright", inset=0.05, cex=1, legCex=0.7, horiz=FALSE) {
  plot.matrix = cbind(ECPGobject$meanPLD, ECPGobject$forwardPrice, ECPGobject$adjustedPrices, ECPGobject$averagePrice)
  colnames(plot.matrix) = c("Monthly PLD", "Forward Price", "Model Price", "Model Average")
  rownames(plot.matrix) = rownames(ECPGobject$dataset)
  
  Nseries = ncol(plot.matrix)
  color.pallette = c("#56B4E9", "#F0E442", "#999999", "#E69F00") # Selected for easy visualization for colorblind
  line.types = c(2, 1, 2, 1)
  
  Nrows = 1:nrow(plot.matrix)
  labels=row.names(plot.matrix)
  
  par(family="serif", cex=cex)
  matplot(Nrows, plot.matrix, type='l', xlab="", ylab='Electricity forward price (R$/MWh)', ylim=ylim, lty=line.types, lwd=lwd, col=color.pallette, axes=F)
  axis(2)
  axis(side=1,at=Nrows,labels=labels)
  grid (NULL,NULL, lty = 6, col = "grey", lwd = 1)
  legend(x=legendPos, inset=inset, legend=colnames(plot.matrix),col=color.pallette, lty=line.types, cex=legCex, horiz=horiz, lwd=lwd)
}

data.prettyPlot <- function(DataObject, ylim=c(0,350), lwd = 2, legendPos="bottomright", inset=0.05, cex=1, legCex=0.7, horiz=FALSE) {
  plot.matrix = cbind(DataObject$meanPLD, DataObject$forwardPrice, DataObject$averagePrice)
  colnames(plot.matrix) = c("NEWAVE Spot", "DCide Forward", "Model Forward")
  rownames(plot.matrix) = DataObject$labels
  
  Nseries = ncol(plot.matrix)
  color.pallette = c("#56B4E9", "#F0E442", "#E69F00") # Selected for easy visualization for colorblind
  line.types = c(2, 1, 1)
  
  Nrows = 1:nrow(plot.matrix)
  labels=row.names(plot.matrix)
  
  LastInSample = which(data.index[,3]=="OUTSAMPLE")[1] - 1
  
  par(family="serif", cex=cex)
  matplot(Nrows, plot.matrix, type='l', xlab="", ylab='Electricity forward price (R$/MWh)', ylim=ylim, lty=line.types, lwd=lwd, col=color.pallette, axes=F)
  axis(2)
  axis(side=1,at=Nrows,labels=labels)
  grid (NULL,NULL, lty = 6, col = "grey", lwd = 1)
  abline(v=LastInSample + 0.5, lty = 2)
  text(x=LastInSample/2 + 0.5, y=25, "In sample")
  text(x=(nrow(plot.matrix) - LastInSample)/2 + LastInSample + 0.5, y=25, "Out-of-sample")
  legend(x=legendPos, inset=inset, legend=colnames(plot.matrix),col=color.pallette, lty=line.types, cex=legCex, horiz=horiz, lwd=lwd)
}

## Obtain parameters from training 

buyer.param = c(Lambda1 = 0.3063046, Lambda2 = 0.1453824, Alpha1 = 0.2892819, Alpha2 = 0.7380821)
seller.param = c(Lambda1 = 0.4820085, Lambda2 = 0.3398702, Alpha1 = 0.7671451, Alpha2 = 0.9881327)


anualRate = 0.0

data.index = rbind(
  c(7, 2019, "insample"),
  c(8, 2019, "insample"),
  c(9, 2019, "insample"),
  c(10, 2019, "insample"),
  c(11, 2019, "insample"),
  c(12, 2019, "insample"),
  c(1, 2020, "insample"),
  c(2, 2020, "insample"),
  c(3, 2020, "OUTSAMPLE"),
  c(4, 2020, "OUTSAMPLE"),
  c(5, 2020, "OUTSAMPLE"),
  c(6, 2020, "OUTSAMPLE"),
  c(7, 2020, "OUTSAMPLE"),
  c(8, 2020, "OUTSAMPLE"),
  c(9, 2020, "OUTSAMPLE"),
  c(10, 2020, "OUTSAMPLE")
)

seriesA = 1

loadPlots = FALSE
savePDF = FALSE
writeData = TRUE

## Load Data

forwardPrice.data = read.csv(file='../Data/ForwardPrices.csv')
colnames(forwardPrice.data) = c("Period", paste(c(rep(2019,6), rep(2020, 10)), c(7:12,1:10), sep="_"))

selectedData = paste(data.index[,2], data.index[,1], sep="_")
forwardPrices = as.numeric(forwardPrice.data[seriesA,selectedData])

completeDiscountVector = getCompleteDiscountVector(as.numeric(data.index[,1]), seriesA, anualRate)
data.abbrev = paste(month.abb[as.numeric(data.index[,1])], data.index[,2], sep="_")

all.data = list()
averagePrices = c()
spotPrices = c()

for (n in 1:nrow(data.index)) {
  
  ## Set Variables
  discountVector = getDiscountBounds(as.numeric(data.index[n,1]), seriesA, completeDiscountVector)
  
  forwardPrice = forwardPrices[n]
  
  ## Load Data
  original.data = read.csv(file=paste('../Data/Data', data.abbrev[n], paste0('A', seriesA, '.csv'), sep="_"), col.names=paste(month.abb, as.numeric(data.index[n,2])+seriesA, sep="-"))
  
  spotPrices[n] = mean(colMeans(original.data))
  
  ECPG_object = createECPGObject(dataset = t(original.data), buyerParam = buyer.param, sellerParam = seller.param, discountVector = discountVector, forwardPrice = forwardPrice)
  
  ## Calculate results from optimized parameters
  ECPG_object = ECPG.equilibrium(ECPG_object)
  
  ## Print results onscreen
  #print(ECPG_object$adjustedPrices)
  
  all.data[[data.abbrev[n]]] <- ECPG_object
  averagePrices[n] = all.data[[n]]$averagePrice

  ## Plot price series  
  if (loadPlots) {
    if (savePDF) pdf(file=paste0("../Plots/A", seriesA,"_", selectedData[n], "_", data.index[n,3], ".pdf"), width=8, height=5)
    series.prettyPlot(all.data[[n]], lwd=3, cex=1.2, horiz=TRUE, legCex=0.7, legendPos="top")
    if (savePDF) dev.off()
  }
}

modelError = forwardPrices - averagePrices

resultMatrix = rbind(spotPrices, forwardPrices, averagePrices, modelError)
colnames(resultMatrix) = paste(selectedData, t(data.index[,3]), sep="_")

if (writeData) {
  write.csv(resultMatrix, file='ModelResults.csv')
  write.csv2(resultMatrix, file='ModelResults2.csv')
}

errorStats = matrix(c(
  mean(modelError[data.index[,3]=="insample"]),
  mean(modelError[data.index[,3]=="OUTSAMPLE"]),
  sd(modelError[data.index[,3]=="insample"]),
  sd(modelError[data.index[,3]=="OUTSAMPLE"])),
  nrow=2, byrow=TRUE, dimnames=list(c("mean", "sd"), c("insample", "OUTSAMPLE")))

consolidated.data <- list()
consolidated.data$meanPLD = spotPrices
consolidated.data$forwardPrice = forwardPrices
consolidated.data$averagePrice = averagePrices
consolidated.data$labels = data.abbrev

data.prettyPlot(consolidated.data, lwd=3, cex=1.2, horiz=TRUE, legCex=1, legendPos="top")

