


## Local settings
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

series.prettyPlot <- function(ECPGobject, ylim=c(0,350), lwd = 2, legendPos="bottomright", inset=0.05, cex=1, legCex=0.7, horiz=FALSE) {
  plot.matrix = cbind(ECPGobject$meanPLD, ECPGobject$forwardPrice, ECPGobject$adjustedPrices, ECPGobject$averagePrice)
  colnames(plot.matrix) = c("Average PLD", "Forward Price", "Model Price", "Model Average")
  rownames(plot.matrix) = rownames(ECPGobject$dataset)
  
  Nseries = ncol(plot.matrix)
  color.pallette = c("#56B4E9", "#F0E442", "#999999", "#E69F00") # Selected for easy visualization for colorblind
  line.types = c(2, 1, 2, 1)
  
  Nrows = 1:nrow(plot.matrix)
  labels=row.names(plot.matrix)
  
  par(family="A", cex=cex)
  matplot(Nrows, plot.matrix, family="A", type='l', xlab="", ylab='Energy Price (R$/MWh)', ylim=ylim, lty=line.types, lwd=lwd, col=color.pallette, axes=F)
  axis(2)
  axis(side=1,at=Nrows,labels=labels)
  grid (NULL,NULL, lty = 6, col = "grey")
  legend(x=legendPos, inset=inset, legend=colnames(plot.matrix),col=color.pallette, lty=line.types, cex=legCex, horiz=horiz, lwd=lwd)
}


## Obtain parameters from training 
buyer.param = c(Lambda1 = 0.5712411, Lambda2 = 0.1400453, Alpha1 = 0.5845306, Alpha2 = 0.5995087)
seller.param = c(Lambda1 = 0.3483552, Lambda2 = 0.5599081, Alpha1 = 0.7256726, Alpha2 = 0.9737610)

# Set risk free rates
anualRate = 0.05
monthlyRate = (1+anualRate)^(1/12)-1

#### MODEL VALIDATION WITH DATA FROM AUGUST A1 ####  

## Set Variables
temp.discount.vector <- rep(NA, 16)
for (i in 1:16){
  temp.discount.vector[i] <- (1+monthlyRate)^(i-1)
}
augustA1.discountVector = temp.discount.vector[5:16]

augustA1.forwardPrice = 192.75

## Load Data
augustA1.original.data = read.csv(file='../Data/DataAugustA1.csv', col.names=paste(month.abb, 2020, sep="-"))
augustA1.data = t(augustA1.original.data)

augustA1.ECPG_object = createECPGObject(dataset = augustA1.data, buyerParam = buyer.param, sellerParam = seller.param, discountVector = augustA1.discountVector, forwardPrice = augustA1.forwardPrice)

## Calculate results from optimized parameters
augustA1.ECPG_object = ECPG.equilibrium(augustA1.ECPG_object)

## Print results onscreen
print(augustA1.ECPG_object$adjustedPrices)

## Plot price series
series.prettyPlot(augustA1.ECPG_object, lwd=3, cex=1.3, horiz=TRUE, legCex=0.8, legendPos="bottom")

#### MODEL VALIDATION WITH DATA FROM JANUARY A2 ####  

## Set Variables
temp.discount.vector <- rep(NA, 35)
for (i in 1:35){
  temp.discount.vector[i] <- (1+monthlyRate)^(i-1)
}
januaryA2.discountVector = temp.discount.vector[24:35]

januaryA2.forwardPrice = 179.26

## Load Data
januaryA2.original.data = read.csv(file='../Data/DataJanuaryA2.csv', col.names=paste(month.abb, 2021, sep="-"))
januaryA2.data = t(januaryA2.original.data)

januaryA2.ECPG_object = createECPGObject(dataset = januaryA2.data, buyerParam = buyer.param, sellerParam = seller.param, discountVector = januaryA2.discountVector, forwardPrice = januaryA2.forwardPrice)

## Calculate results from optimized parameters
januaryA2.ECPG_object = ECPG.equilibrium(januaryA2.ECPG_object)

## Print results onscreen
print(januaryA2.ECPG_object$adjustedPrices)

## Plot price series
series.prettyPlot(januaryA2.ECPG_object, lwd=3, cex=1.3, horiz=TRUE, legCex=0.8, legendPos="top")