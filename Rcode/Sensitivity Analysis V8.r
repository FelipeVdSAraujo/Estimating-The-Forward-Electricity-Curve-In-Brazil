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

sensitivity.foo <- function(i, arg, output, combination.matrix, ECPGobject) {
  
  ECPGobject[[arg]][TRUE] = combination.matrix[i, ]
  
  return(ECPG.equilibrium(ECPGobject)[[output]])
}

sensitivity.prettyPlot <- function(plot.matrix, xlabel="", ylim=c(0,250), lwd = 2, legendTitle, legendPos="bottomright", inset=0.05, ncol=1, cex=1, legCex=0.9) {
  Nseries = ncol(plot.matrix)
  viridis.pallette = viridis(Nseries)
  line.types = rep(1, Nseries)
  
  Nrows = 1:nrow(plot.matrix)
  labels=row.names(plot.matrix)
  
  par(family="A", cex=cex)
  matplot(Nrows, plot.matrix, type='l', xlab=xlabel, ylab='Electricity forward price (R$/MWh)', ylim=ylim, lty=line.types, col=viridis.pallette, lwd=lwd, axes=F)
  axis(2)
  axis(side=1,at=Nrows,labels=labels)
  grid (NULL,NULL, lty = 6, col = "grey")
  legend(x=legendPos, inset=inset, legend=colnames(plot.matrix), horiz=FALSE, lty=line.types, col=viridis.pallette, lwd=lwd, cex=legCex, ncol=ncol, title=legendTitle)
}


## Set Variables
buyer.param = c(Lambda1 = 0.25159927, Lambda2 = 0.09390924, Alpha1 = 0.32334330, Alpha2 = 0.69293115)
seller.param = c(Lambda1 = 0.1712721, Lambda2 = 0.7466865, Alpha1 = 0.7641934, Alpha2 = 0.9912755)


anualRate = 0.05
monthlyRate = (1+anualRate)^(1/12)-1

temp.discount.vector <- rep(NA, 23)
for (i in 1:23){
  temp.discount.vector[i] <- (1+monthlyRate)^(i-1)
}
julyA1.discountVector = temp.discount.vector[5:16]

julyA1.forwardPrice = 205.43

## Load Data
julyA1.original.data = read.csv(file='../Data/Data_Jul_2019_A1.csv', col.names=paste(month.abb, 2020, sep="-"))
julyA1.data = t(julyA1.original.data)

julyA1.ECPG_object = createECPGObject(dataset = julyA1.data, buyerParam = buyer.param, sellerParam = seller.param, discountVector = julyA1.discountVector, forwardPrice = julyA1.forwardPrice)

# Get the original results for comparison
julyA1.Results = ECPG.equilibrium(julyA1.ECPG_object)

# Set sensitivity ranges

lambdas.vec = seq(0, 0.45, by=0.05)
lambdas.matrix = as.matrix(expand.grid(lambdas.vec, lambdas.vec))
colnames(lambdas.matrix) = c('Lambda1', 'Lambda2')
nLambdas = nrow(lambdas.matrix)

buyer.lambdas.matrix = cbind(lambdas.matrix, Alpha1=rep(buyer.param['Alpha1'], nLambdas), Alpha2=rep(buyer.param['Alpha2'], nLambdas))

lambdaBuyer.sensitivity = sapply(1:nLambdas,
                                 sensitivity.foo,
                                 arg = "buyerParam",
                                 output = "averagePrice",
                                 combination.matrix = buyer.lambdas.matrix,
                                 julyA1.ECPG_object)

lambdaBuyer.sensiMatrix = matrix(lambdaBuyer.sensitivity,
                                 ncol=length(lambdas.vec),
                                 byrow = TRUE,
                                 dimnames = list(lambdas.vec, lambdas.vec))

seller.lambdas.matrix = cbind(lambdas.matrix, Alpha1=rep(seller.param['Alpha1'], nLambdas), Alpha2=rep(seller.param['Alpha2'], nLambdas))

lambdaSeller.sensitivity = sapply(1:nLambdas,
                                  sensitivity.foo,
                                  arg = "sellerParam",
                                  output = "averagePrice",
                                  combination.matrix = seller.lambdas.matrix,
                                  julyA1.ECPG_object)

lambdaSeller.sensiMatrix = matrix(lambdaSeller.sensitivity,
                                  ncol=length(lambdas.vec),
                                  byrow = TRUE,
                                  dimnames = list(lambdas.vec, lambdas.vec))

# Plot sensitivity results

require(plotly)

ax.xy <- list(ticketmode = 'array', ticktext = as.character(lambdas.vec), tickvals=0:9)
axz <- list(title = "Electricity forward price (R$/MWh)")


fig.Buyer <- plot_ly(width=1050, height=750) %>%
  add_surface(z = lambdaBuyer.sensiMatrix) %>%
  layout(scene = list(xaxis=c(ax.xy, title = "Buyer Lambda1"),yaxis=c(ax.xy, title = "Buyer Lambda2"),zaxis=axz))

#fig.Buyer

fig.Seller <- plot_ly(width=1050, height=750) %>%
  add_surface(z = lambdaSeller.sensiMatrix) %>%
  layout(scene = list(xaxis=c(ax.xy, title = "Seller Lambda1"),yaxis=c(ax.xy, title = "Seller Lambda2"),zaxis=axz))

#fig.Seller


## Sensitivity for Lambda and Alpha on ECP

lambdaAlpha.vec = seq(0, 0.9, by=0.1)
lambdaAlpha.matrix = as.matrix(expand.grid(lambdaAlpha.vec, lambdaAlpha.vec))
colnames(lambdaAlpha.matrix) = c('Lambda1', 'Alpha1')
nLambdaAlpha = nrow(lambdaAlpha.matrix)

buyer.lambdaAlpha.matrix = cbind(lambdaAlpha.matrix, Lambda2=rep(0, nLambdaAlpha), Alpha2=rep(buyer.param['Alpha2'], nLambdaAlpha))
buyer.lambdaAlpha.matrix = buyer.lambdaAlpha.matrix[, c('Lambda1', 'Lambda2', 'Alpha1', 'Alpha2')]

lambdaAlphaBuyer.sensitivity = sapply(1:nLambdaAlpha,
                                      sensitivity.foo,
                                      arg="buyerParam",
                                      output = "averagePrice",
                                      combination.matrix = buyer.lambdaAlpha.matrix,
                                      within(julyA1.ECPG_object, buyerParam[TRUE] <- rep(0,4)))

lambdaAlphaBuyer.sensiMatrix = matrix(lambdaAlphaBuyer.sensitivity,
                                      ncol=length(lambdaAlpha.vec),
                                      byrow = TRUE,
                                      dimnames = list(lambdaAlpha.vec, lambdaAlpha.vec))

seller.lambdaAlpha.matrix = cbind(lambdaAlpha.matrix, Lambda2=rep(0, nLambdaAlpha), Alpha2=rep(seller.param['Alpha2'], nLambdaAlpha))
seller.lambdaAlpha.matrix = seller.lambdaAlpha.matrix[, c('Lambda1', 'Lambda2', 'Alpha1', 'Alpha2')]

lambdaAlphaSeller.sensitivity = sapply(1:nLambdaAlpha,
                                       sensitivity.foo,arg="sellerParam",
                                       output = "averagePrice",
                                       combination.matrix = seller.lambdaAlpha.matrix,
                                       within(julyA1.ECPG_object, sellerParam[TRUE] <- rep(0,4)))

lambdaAlphaSeller.sensiMatrix = matrix(lambdaAlphaSeller.sensitivity,
                                       ncol=length(lambdaAlpha.vec),
                                       byrow = TRUE,
                                       dimnames = list(lambdaAlpha.vec, lambdaAlpha.vec))


# Plot sensitivity for Lambda and Alpha
#matplot(lambdaAlphaBuyer.sensiMatrix, type='l')
#matplot(lambdaAlphaSeller.sensiMatrix, type='l')

require(viridis)

sensitivity.prettyPlot(lambdaAlphaBuyer.sensiMatrix, xlabel=expression(paste("Buyer ", alpha, 1)), ylim=c(170, 250), ncol=2, cex=1.4, lwd=3, legendTitle=expression(paste("Buyer ", lambda, 1)))

sensitivity.prettyPlot(lambdaAlphaSeller.sensiMatrix, xlabel=expression(paste("Seller ", alpha, 1)), ylim=c(0, 180), legendTitle=expression(paste("Seller ", lambda, 1)), legendPos = "topleft", ncol=5, inset=0.02, cex=1.4, lwd=3)

## Prepare variables for the sensitivity analysis

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

rates.sensitivity = sapply(1:nRates,
                           sensitivity.foo,
                           arg="discountVector",
                           output = "averagePrice", 
                           combination.matrix = discount.matrix,
                           julyA1.ECPG_object)
names(rates.sensitivity) = annualRates.vec

rates.adjustedPrices = sapply(1:nRates,
                              sensitivity.foo,
                              arg="discountVector",
                              output = "adjustedPrices", 
                              combination.matrix = discount.matrix,
                              julyA1.ECPG_object)

rates.sensitivity.matrix = matrix(rates.adjustedPrices,
                                  ncol = nRates,
                                  byrow = FALSE,
                                  dimnames = list(rownames(julyA1.ECPG_object$dataset), annualRates.vec))

# Check results
julyA1.Results = ECPG.equilibrium(within(julyA1.ECPG_object, discountVector <- discount.matrix["0.05", ]))
julyA1.Results$averagePrice == rates.sensitivity["0.05"]

# Plot
plot(annualRates.vec, rates.sensitivity, ylim=c(0,250), ylab="Electricity forward price (R$/MWh)", xlab="Annual Risk-Free Rate", col="steelblue")

sensitivity.prettyPlot(rates.sensitivity.matrix, ylim=c(0, 300), ncol=5, cex=1.4, lwd=3, legendTitle="Risk-Free Rate aa.")


# Plotly contour plot


fig.BuyerLambdaFlat <- plot_ly(
  x = lambdas.vec, 
  y = lambdas.vec, 
  z = lambdaBuyer.sensiMatrix, 
  type = "contour" 
)

fig.BuyerLambdaFlat <- fig.BuyerLambdaFlat %>% colorbar(
  len=0.85,
  title = "Electricity\nForward\nPrice\n(R$/MWh)")

fig.BuyerLambdaFlat <- fig.BuyerLambdaFlat %>% layout(
                      xaxis = list(title = "Buyer's Lambda 1"),
                      yaxis = list(title = "Buyer's Lambda 2"),
                      font=list(size=24)
                      )


fig.SellerLambdaFlat <- plot_ly(
  x = lambdas.vec, 
  y = lambdas.vec, 
  z = lambdaSeller.sensiMatrix, 
  type = "contour" 
)

fig.SellerLambdaFlat <- fig.SellerLambdaFlat %>% colorbar(
  len=0.85,
  title = "Electricity\nForward\nPrice\n(R$/MWh)")

fig.SellerLambdaFlat <- fig.SellerLambdaFlat %>% layout(
  xaxis = list(title = "Seller's Lambda 1"),
  yaxis = list(title = "Seller's Lambda 2"),
  font=list(size=24)
)

# Plotly Charts must be run outside of source
# fig.Buyer
# fig.Seller
# fig.SellerLambdaFlat
# fig.BuyerLambdaFlat
