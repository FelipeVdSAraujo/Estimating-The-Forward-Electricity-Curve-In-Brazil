# Estimating The Forward Electricity Curve In Brazil With A Model Of Two Agents Using Contracts By Difference And ECP_G Function

Authors: Felipe Van de Sande Araujo, Cristina Spineti Luz, Leonardo Lima Gomes, Luís Eduardo Teixeira Brandão  

Abstract: The development of simple and effective mechanisms to estimate the value of the forward curve of power could enable market participants to better price hedging or speculative positions. This could in turn provide transparency in future price definition to all market participants and lead to more safety and liquidity in the market for electricity futures and power derivatives. This work presents a model for two market participants, a buyer and a seller of a contract for difference on the future spot price of electricity in southwest Brazil. It is shown that this model is representative of all market participants that have exposure to the future price of power. Each participant’s utility function is modeled using a Generalized Extended CVaR Preference (ECP_G) and the market equilibrium is obtained through the minimization of the quadratic difference between the certainty equivalent of both agents. The results are compared with prediction of the future spot price of power made by market specialists and found to yield reasonable results when using out of sample data.

The files provided here are the calculations for the paper. 
They are organized in the following way:
Data - contain datasets used by the authors and loaded by the R scripts.
Rcode - contain the scripts that can be executed provided the data is organized in the same structure as this repository.

In the main level the calculations are saved in HTML files created by R Notebooks which show the executed commands and their results, along with plots.
The sensitivity analysis uses an interactive *plotly* widget that can take a while to fully load.
