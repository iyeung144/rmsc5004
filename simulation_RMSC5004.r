######################
# Stock price simulation
#
#
######################

library("plyr")
library("ggplot2")
library("ggthemes")

tickers <- read.csv("D:\\OneDrive\\Documents\\RMSC\\Course materials\\RMSC5004\\RMSC5004 Group Project\\simulation_data.csv",header=T,stringsAsFactors = FALSE)
sim_price <- c(rep(0,nrow(tickers)))
attach(tickers)
cbind(Ticker,sim_price)

p <- 1000
N <- 52
mu <- 0.0282 #10-Year Treasury Constant Maturity Rate (DGS10)

t <- (0:N)/N
dt <- 1/N

mydata <- list()

for (k in 1:nrow(tickers)) {
  x <- matrix(rep(0,p*(N+1)),nrow=(N+1))
  y <- matrix(rep(0,p*(N+1)),nrow=(N+1))
  for (j in 1:p) {
    
    S0 <- close_price[k]
    sigma <- volatility_260D[k]
    nu <- mu -sigma*sigma*0.5

    z <- rnorm(N,0,1)
    y[1,j] <- S0
    for (i in 1:N) {
      x[i+1,j] <- sqrt(dt)*sum(z[1:i])
      y[i+1,j] <- S0*exp(nu*t[i+1]+sigma*x[i+1,j])
    }
    mydata[[k]] <- y
  }
  sim_price[k] <- mean(y[N+1,])
}

par(mfrow=c(4,2))
par(bg = 'cornsilk')
for (m in 1:nrow(tickers)) {
  matplot(t,mydata[[m]],type="l",xlab="",ylab="Simulated stock price",main = Ticker[m])
}
# png(filename="D:\\OneDrive\\Documents\\RMSC\\Course materials\\RMSC5004\\RMSC5004 Group Project\\simulation.png", 
#     type="cairo",
#     units="mm",
#     width=210,
#     height=297, 
#     pointsize=12, 
#     res=300)
# dev.off()
tickers <- cbind(tickers,sim_price)
ret<-formatC(c((tickers[,4]/tickers[,3]-1)*100),digits=2, format="f")
tickers <- cbind(tickers,ret)
colnames(tickers) <- c("Ticker","Volatility 260 Days", "Last close price","Simulated price", "Return (%)")
print(tickers)


