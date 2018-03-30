library(PortfolioAnalytics)
library(quantmod)
library(PerformanceAnalytics)
library(zoo)
library(DEoptim)
library(ROI)
require(ROI.plugin.glpk)
require(ROI.plugin.quadprog)

# # Get data
# getSymbols(c("MSFT", "SBUX", "IBM", "AAPL", "^GSPC", "AMZN"))
# Get data
fromDate <- '2015-01-01'
HK1313 <- getSymbols("1313.HK", from = fromDate, auto.assign = FALSE)
HK1800 <- getSymbols("1800.HK", from = fromDate, auto.assign = FALSE)
HK1972 <- getSymbols("1972.HK", from = fromDate, auto.assign = FALSE)
HK0002 <- getSymbols("0002.HK", from = fromDate, auto.assign = FALSE)
USEEMV <- getSymbols("EEMV", from = fromDate, auto.assign = FALSE)
USFM <- getSymbols("FM", from = fromDate, auto.assign = FALSE)
# USFLOT <- getSymbols("FLOT", from = fromDate, auto.assign = FALSE)

# Assign to dataframe
# Get adjusted prices
# prices.data <- merge.zoo(HK1313[,6], HK1800[,6], HK1972[,6], HK0002[,6], USEEMV[,6], USFLOT[,6], USFM[,6])
prices.data <- merge.zoo(HK1313[,6], HK1800[,6], HK1972[,6], HK0002[,6],  USEEMV[,6], USFM[,6])

# Calculate returns
returns.data <- CalculateReturns(prices.data)
returns.data <- na.omit(returns.data)

# Set names
colnames(returns.data) <- c("HK.1313","HK.1800","HK.1972","HK.0002","US.EEMV","US.FM")
funds <- colnames(returns.data)

# Save mean return vector and sample covariance matrix
meanReturns <- colMeans(returns.data)
covMat <- cov(returns.data)

# Start with the names of the assets
port <- portfolio.spec(assets = funds)

# Add the full investment constraint that specifies the weights must sum to 1.
port <- add.constraint(portfolio=port, type="weight_sum", min_sum=1, max_sum=1)

# Box
# port <- add.constraint(port, type = "box", min = 0.05,max = 0.25)
port <- add.constraint(portfolio=port,type="box",
                       min=c(0.10, 0.10, 0.10, 0.10, 0.15, 0.15),
                       max=c(0.15, 0.15, 0.15, 0.15, 0.3, 0.3))
# Leverage
port <- add.constraint(portfolio = port, type = "full_investment")

# Generate random portfolios
rportfolios <- random_portfolios(port, permutations = 5000, rp_method = "sample")

# return objective to maximize the portfolio mean return
maxRet.portfolio <- add.objective(portfolio=port, type='return', name='mean')
maxRet.opt <- optimize.portfolio(R=returns.data, portfolio=maxRet.portfolio, optimize_method="ROI", trace=TRUE)
result_maxRetplot <- plot(maxRet.opt, risk.col="StdDev", return.col="mean",
                    main="Maximum Return Optimization", chart.assets=TRUE)

minSD.portfolio <- add.objective(portfolio=port, type="risk", name="var")
minSD.opt <- optimize.portfolio(R = returns.data, portfolio = minSD.portfolio, optimize_method = "ROI", trace = TRUE)
result_minSDplot <- plot(minSD.opt, risk.col="StdDev",return.col="mean",
                    main="Minimum Variance Optimization",  chart.assets=TRUE,xlim=c(0.005,0.027))

# Generate vector of returns
minret <- 0.00/100
maxret <- maxRet.opt$weights %*% meanReturns

vec <- seq(minret, maxret, length.out = 100)

eff.frontier <- data.frame(Risk = rep(NA, length(vec)),Return = rep(NA, length(vec)),SharpeRatio = rep(NA, length(vec)))
frontier.weights <- mat.or.vec(nr = length(vec), nc = ncol(returns.data))
colnames(frontier.weights) <- colnames(returns.data)

for(i in 1:length(vec)){
  eff.port <- add.constraint(port, type = "return", name = "mean", return_target = vec[i])
  eff.port <- add.objective(eff.port, type = "risk", name = "var")
  eff.port <- optimize.portfolio(returns.data, eff.port, optimize_method = "ROI")
  eff.frontier$Risk[i] <- sqrt(t(eff.port$weights) %*% covMat %*% eff.port$weights)
  eff.frontier$Return[i] <- eff.port$weights %*% meanReturns
  eff.frontier$Sharperatio[i] <- eff.port$Return[i] / eff.port$Risk[i]
  frontier.weights[i,] = eff.port$weights
  print(paste(round(i/length(vec) * 100, 0), "% done..."))
}

feasible.sd <- apply(rportfolios, 1, function(x){
  return(sqrt(matrix(x, nrow = 1) %*% covMat %*% matrix(x, ncol = 1)))
})

feasible.means <- apply(rportfolios, 1, function(x){
  return(x %*% meanReturns)
})
feasible.sr <- feasible.means / feasible.sd