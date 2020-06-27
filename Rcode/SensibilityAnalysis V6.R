


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

sensibility.foo <- function(i, arg, output, combination.matrix, ECPGobject) {
  
  ECPGobject[[arg]][TRUE] = combination.matrix[i, ]
  
  return(ECPG.equilibrium(ECPGobject)[[output]])
}

sensibility.prettyPlot <- function(plot.matrix, xlabel="", ylim=c(0,250), lwd = 2, legendTitle, legendPos="bottomright", inset=0.05, ncol=1, cex=1, legCex=0.9) {
  Nseries = ncol(plot.matrix)
  viridis.pallette = viridis(Nseries)
  line.types = rep(1, Nseries)
  
  Nrows = 1:nrow(plot.matrix)
  labels=row.names(plot.matrix)
  
  par(family="A", cex=cex)
  matplot(Nrows, plot.matrix, type='l', xlab=xlabel, ylab='Energy Price (R$/MWh)', ylim=ylim, lty=line.types, col=viridis.pallette, lwd=lwd, axes=F)
  axis(2)
  axis(side=1,at=Nrows,labels=labels)
  grid (NULL,NULL, lty = 6, col = "grey")
  legend(x=legendPos, inset=inset, legend=colnames(plot.matrix), horiz=FALSE, lty=line.types, col=viridis.pallette, lwd=lwd, cex=legCex, ncol=ncol, title=legendTitle)
}


## Set Variables
buyer.param = c(Lambda1 = 0.5712411, Lambda2 = 0.1400453, Alpha1 = 0.5845306, Alpha2 = 0.5995087)
seller.param = c(Lambda1 = 0.3483552, Lambda2 = 0.5599081, Alpha1 = 0.7256726, Alpha2 = 0.9737610)

anualRate = 0.05
monthlyRate = (1+anualRate)^(1/12)-1

temp.discount.vector <- rep(NA, 23)
for (i in 1:23){
  temp.discount.vector[i] <- (1+monthlyRate)^(i-1)
}
januaryA1.discountVector = temp.discount.vector[12:23]

januaryA1.forwardPrice = 199.32

## Load Data
januaryA1.original.data = read.csv(file='../Data/DataJanuaryA1.csv', col.names=paste(month.abb, 2020, sep="-"))
januaryA1.data = t(januaryA1.original.data)

januaryA1.ECPG_object = createECPGObject(dataset = januaryA1.data, buyerParam = buyer.param, sellerParam = seller.param, discountVector = januaryA1.discountVector, forwardPrice = januaryA1.forwardPrice)

# Get the original results for comparison
januaryA1.Results = ECPG.equilibrium(januaryA1.ECPG_object)

# Set sensibility ranges

lambdas.vec = seq(0, 0.45, by=0.05)
lambdas.matrix = as.matrix(expand.grid(lambdas.vec, lambdas.vec))
colnames(lambdas.matrix) = c('Lambda1', 'Lambda2')
nLambdas = nrow(lambdas.matrix)

buyer.lambdas.matrix = cbind(lambdas.matrix, Alpha1=rep(buyer.param['Alpha1'], nLambdas), Alpha2=rep(buyer.param['Alpha2'], nLambdas))

lambdaBuyer.sensibility = sapply(1:nLambdas,
                                 sensibility.foo,
                                 arg = "buyerParam",
                                 output = "averagePrice",
                                 combination.matrix = buyer.lambdas.matrix,
                                 januaryA1.ECPG_object)

lambdaBuyer.sensiMatrix = matrix(lambdaBuyer.sensibility,
                                 ncol=length(lambdas.vec),
                                 byrow = TRUE,
                                 dimnames = list(lambdas.vec, lambdas.vec))

seller.lambdas.matrix = cbind(lambdas.matrix, Alpha1=rep(seller.param['Alpha1'], nLambdas), Alpha2=rep(seller.param['Alpha2'], nLambdas))

lambdaSeller.sensibility = sapply(1:nLambdas,
                                  sensibility.foo,
                                  arg = "sellerParam",
                                  output = "averagePrice",
                                 combination.matrix = seller.lambdas.matrix,
                                 januaryA1.ECPG_object)

lambdaSeller.sensiMatrix = matrix(lambdaSeller.sensibility,
                                 ncol=length(lambdas.vec),
                                 byrow = TRUE,
                                 dimnames = list(lambdas.vec, lambdas.vec))

# Plot sensibility results

require(plotly)

ax.xy <- list(ticketmode = 'array', ticktext = as.character(lambdas.vec), tickvals=0:9)
axz <- list(title = "Power Price (R$/MWh)")


fig.Buyer <- plot_ly(width=700, height=500) %>%
  add_surface(z = lambdaBuyer.sensiMatrix) %>%
  layout(scene = list(xaxis=c(ax.xy, title = "Buyer Lambda1"),yaxis=c(ax.xy, title = "Buyer Lambda2"),zaxis=axz))

#fig.Buyer

fig.Seller <- plot_ly(width=700, height=500) %>%
  add_surface(z = lambdaSeller.sensiMatrix) %>%
  layout(scene = list(xaxis=c(ax.xy, title = "Seller Lambda1"),yaxis=c(ax.xy, title = "Seller Lambda2"),zaxis=axz))

#fig.Seller


## Sensibility for Lambda and Alpha on ECP

lambdaAlpha.vec = seq(0, 0.9, by=0.1)
lambdaAlpha.matrix = as.matrix(expand.grid(lambdaAlpha.vec, lambdaAlpha.vec))
colnames(lambdaAlpha.matrix) = c('Lambda1', 'Alpha1')
nLambdaAlpha = nrow(lambdaAlpha.matrix)

buyer.lambdaAlpha.matrix = cbind(lambdaAlpha.matrix, Lambda2=rep(0, nLambdaAlpha), Alpha2=rep(buyer.param['Alpha2'], nLambdaAlpha))
buyer.lambdaAlpha.matrix = buyer.lambdaAlpha.matrix[, c('Lambda1', 'Lambda2', 'Alpha1', 'Alpha2')]

lambdaAlphaBuyer.sensibility = sapply(1:nLambdaAlpha,
                                 sensibility.foo,
                                 arg="buyerParam",
                                 output = "averagePrice",
                                 combination.matrix = buyer.lambdaAlpha.matrix,
                                 within(januaryA1.ECPG_object, buyerParam[TRUE] <- rep(0,4)))

lambdaAlphaBuyer.sensiMatrix = matrix(lambdaAlphaBuyer.sensibility,
                                 ncol=length(lambdaAlpha.vec),
                                 byrow = TRUE,
                                 dimnames = list(lambdaAlpha.vec, lambdaAlpha.vec))

seller.lambdaAlpha.matrix = cbind(lambdaAlpha.matrix, Lambda2=rep(0, nLambdaAlpha), Alpha2=rep(seller.param['Alpha2'], nLambdaAlpha))
seller.lambdaAlpha.matrix = seller.lambdaAlpha.matrix[, c('Lambda1', 'Lambda2', 'Alpha1', 'Alpha2')]

lambdaAlphaSeller.sensibility = sapply(1:nLambdaAlpha,
                                  sensibility.foo,arg="sellerParam",
                                  output = "averagePrice",
                                  combination.matrix = seller.lambdaAlpha.matrix,
                                  within(januaryA1.ECPG_object, sellerParam[TRUE] <- rep(0,4)))

lambdaAlphaSeller.sensiMatrix = matrix(lambdaAlphaSeller.sensibility,
                                  ncol=length(lambdaAlpha.vec),
                                  byrow = TRUE,
                                  dimnames = list(lambdaAlpha.vec, lambdaAlpha.vec))


# Plot sensibility for Lambda and Alpha
#matplot(lambdaAlphaBuyer.sensiMatrix, type='l')
#matplot(lambdaAlphaSeller.sensiMatrix, type='l')

require(viridis)

sensibility.prettyPlot(lambdaAlphaBuyer.sensiMatrix, xlabel=expression(paste("Buyer ", alpha, 1)), ylim=c(170, 205), ncol=2, cex=1.4, lwd=3, legendTitle=expression(paste("Buyer ", lambda, 1)))

sensibility.prettyPlot(lambdaAlphaSeller.sensiMatrix, xlabel=expression(paste("Seller ", alpha, 1)), ylim=c(0, 180), legendTitle=expression(paste("Seller ", lambda, 1)), legendPos = "topleft", ncol=5, inset=0.02, cex=1.4, lwd=3)


#### Sensibility for Risk-Free Rate on ECP WITH DATA FROM JANUARY A2 ####  

## Set Variables
januaryA2.discountVector = rep(1, 12)
januaryA2.forwardPrice = 179.26

## Load Data
januaryA2.original.data = read.csv(file='../Data/DataJanuaryA2.csv', col.names=paste(month.abb, 2021, sep="-"))
januaryA2.data = t(januaryA2.original.data)

januaryA2.ECPG_object = createECPGObject(dataset = januaryA2.data, buyerParam = buyer.param, sellerParam = seller.param, discountVector = januaryA2.discountVector, forwardPrice = januaryA2.forwardPrice)

## Prepare variables for the sensibility analysis

annualRates.vec = seq(0, 0.09, by=0.01)
monthlyRates.vec = (1+annualRates.vec)^(1/12)-1
nRates = length(monthlyRates.vec)

discount.matrix <- matrix(NA, ncol=35, nrow=nRates)
for (ii in 1:nRates) {
  for (i in 1:35){
    discount.matrix[ii, i] <- (1+monthlyRates.vec[ii])^(i-1)
  }
}
discount.matrix = discount.matrix[, 24:35]
rownames(discount.matrix) = annualRates.vec

rates.sensibility = sapply(1:nRates,
                            sensibility.foo,
                            arg="discountVector",
                           output = "averagePrice", 
                            combination.matrix = discount.matrix,
                            januaryA2.ECPG_object)
names(rates.sensibility) = annualRates.vec

rates.adjustedPrices = sapply(1:nRates,
                           sensibility.foo,
                           arg="discountVector",
                           output = "adjustedPrices", 
                           combination.matrix = discount.matrix,
                           januaryA2.ECPG_object)

rates.sensibility.matrix = matrix(rates.adjustedPrices,
                                  ncol = nRates,
                                  byrow = FALSE,
                                  dimnames = list(rownames(januaryA2.ECPG_object$dataset), annualRates.vec))

# Check results
januaryA2.Results = ECPG.equilibrium(within(januaryA2.ECPG_object, discountVector <- discount.matrix["0.05", ]))
januaryA2.Results$averagePrice == rates.sensibility["0.05"]

# Plot
plot(annualRates.vec, rates.sensibility, ylim=c(0,200), ylab="Energy Price (R$/MWh)", xlab="Annual Risk-Free Rate", col="steelblue")

sensibility.prettyPlot(rates.sensibility.matrix, ylim=c(0, 205), ncol=5, cex=1.4, lwd=3, legendTitle="Risk-Free Rate aa.")

# Plotly Charts must be run outside of source
# fig.Buyer
# fig.Seller
