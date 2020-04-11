## Shark bites analysis
## Accompanies the article: "BRADSHAW, CJA, P MEAGHER, MJ THIELE, CA SIMPFENDORFER, RG HARCOURT, C HUVENEERS. In review. Predicting
##  potential future reduction in shark attack. People and Nature. In review
## April 2020
## Corey J. A. Bradshaw (corey.bradshaw@flinders.edu.au)
## Global Ecology, Flinders University
## GlobalEcologyFlinders.com

## remove everything
rm(list = ls())

# required R libraries
library(TSA)
library(stats)
library(zoo)
library(lme4)
library(boot)
library(Hmisc)
library(ggplot2)
library(igraph)
library(zoo)
library(TTR)
library(ggplot2)
library(plotly)
library(biwavelet)
library(pracma)
library(MuMIn)
library(nlme)

# source files
source("new_lmer_AIC_tables3.R") 
source("r.squared.R")

# Set functions
# beta distribution shape parameter estimator function
estBetaParams <- function(mu, var) {
  alpha <- ((1 - mu) / var - 1 / mu) * mu ^ 2
  beta <- alpha * (1 / mu - 1)
  return(params = list(alpha = alpha, beta = beta))
}

AICc <- function(...) {
  models <- list(...)
  num.mod <- length(models)
  AICcs <- numeric(num.mod)
  ns <- numeric(num.mod)
  ks <- numeric(num.mod)
  AICc.vec <- rep(0,num.mod)
  for (i in 1:num.mod) {
    if (length(models[[i]]$df.residual) == 0) n <- models[[i]]$dims$N else n <- length(models[[i]]$residuals)
    if (length(models[[i]]$df.residual) == 0) k <- sum(models[[i]]$dims$ncol) else k <- (length(models[[i]]$coeff))+1
    AICcs[i] <- (-2*logLik(models[[i]])) + ((2*k*n)/(n-k-1))
    ns[i] <- n
    ks[i] <- k
    AICc.vec[i] <- AICcs[i]
  }
  return(AICc.vec)
}

delta.AIC <- function(x) x - min(x) ## where x is a vector of AIC
weight.AIC <- function(x) (exp(-0.5*x))/sum(exp(-0.5*x)) ## Where x is a vector of dAIC
ch.dev <- function(x) ((( as.numeric(x$null.deviance) - as.numeric(x$deviance) )/ as.numeric(x$null.deviance))*100) ## % change in deviance, where x is glm object

linreg.ER <- function(x,y) { # where x and y are vectors of the same length; calls AICc, delta.AIC, weight.AIC functions
  fit.full <- lm(y ~ x); fit.null <- lm(y ~ 1)
  AIC.vec <- c(AICc(fit.full),AICc(fit.null))
  dAIC.vec <- delta.AIC(AIC.vec); wAIC.vec <- weight.AIC(dAIC.vec)
  ER <- wAIC.vec[1]/wAIC.vec[2]
  r.sq.adj <- as.numeric(summary(fit.full)[9])
  return(c(ER,r.sq.adj))
}


# set wd
setwd("~/.../data/") ## data files (request from authors)

## import
bites <- read.csv("sharkbite.exp.csv", header=T)
dim(bites)

# choose which subset
bites2 <- subset(bites, Year>1899)
#bites2 <- subset(bites, Year>1899 & Recovery.Status == "FATAL") # fatal only
#bites2 <- subset(bites, Year>1899 & (Provoked == "Unprovoked" | Provoked == "UNPROVOKED")) # unprovoked only only
#bites2 <- subset(bites, Year>1899 & (Shark.Common.Name == "White shark" | Shark.Common.Name == "White Shark" | Shark.Common.Name == "WHITE SHARK" | Shark.Common.Name == "White Shark " | Shark.Common.Name == "GREAT WHITE SHARK"))  # white sharks only only
#bites2 <- subset(bites, Year>1899 & (Shark.Common.Name == "White shark" | Shark.Common.Name == "White Shark" | Shark.Common.Name == "WHITE SHARK" | Shark.Common.Name == "White Shark " | Shark.Common.Name == "Bull Shark" | Shark.Common.Name == "BULL SHARK" | Shark.Common.Name == "TIGER" | Shark.Common.Name == "Tiger Shark" | Shark.Common.Name == "Tiger Shark " | Shark.Common.Name == "TIGER SHARK"))  # white, bull & tiger sharks only only

# tables
bites.yr <- as.numeric(table(bites2$Year))
yr.lab <- as.numeric(attr(table(bites2$Year), "names"))
plot(yr.lab, bites.yr, type="c", lty=3, xlab="", ylab="bites/yr")

# by state
table(bites2$State)
nsw.bites <- subset(bites2, State=="NSW")
nsw.bites.yr <- as.numeric(table(nsw.bites$Year))
nsw.yr.lab <- as.numeric(attr(table(nsw.bites$Year), "names"))
lines(nsw.yr.lab, nsw.bites.yr, type="c", col="grey", lwd=2)

qld.bites <- subset(bites2, State=="QLD")
qld.bites.yr <- as.numeric(table(qld.bites$Year))
qld.yr.lab <- as.numeric(attr(table(qld.bites$Year), "names"))
lines(qld.yr.lab, qld.bites.yr, type="l", col="red", lwd=2)

wa.bites <- subset(bites2, State=="WA")
wa.bites.yr <- as.numeric(table(wa.bites$Year))
wa.yr.lab <- as.numeric(attr(table(wa.bites$Year), "names"))
lines(wa.yr.lab, wa.bites.yr, type="l", col="blue", lwd=2)

sa.bites <- subset(bites2, State=="SA")
sa.bites.yr <- as.numeric(table(sa.bites$Year))
sa.yr.lab <- as.numeric(attr(table(sa.bites$Year), "names"))
lines(sa.yr.lab, sa.bites.yr, type="l", col="green", lwd=2)

vic.bites <- subset(bites2, State=="VIC")
vic.bites.yr <- as.numeric(table(vic.bites$Year))
vic.yr.lab <- as.numeric(attr(table(vic.bites$Year), "names"))
lines(vic.yr.lab, vic.bites.yr, type="l", col="yellow", lwd=2)

tas.bites <- subset(bites2, State=="TAS")
tas.bites.yr <- as.numeric(table(tas.bites$Year))
tas.yr.lab <- as.numeric(attr(table(tas.bites$Year), "names"))
lines(tas.yr.lab, tas.bites.yr, type="l", col="purple", lwd=2)

nt.bites <- subset(bites2, State=="NT")
nt.bites.yr <- as.numeric(table(nt.bites$Year))
nt.yr.lab <- as.numeric(attr(table(nt.bites$Year), "names"))
lines(nt.yr.lab, nt.bites.yr, type="l", col="pink", lwd=2)

# plot as matrix
par(mfrow=c(2,4),pty="m")
plot(yr.lab, bites.yr, type="c", lty=1, lwd=1, col="black",xlab="",ylab="",main="ALL", ylim=c(0,max(bites.yr)))
plot(nsw.yr.lab, nsw.bites.yr, type="c", col="black", lwd=1,xlab="",ylab="",main="NSW", ylim=c(0,max(bites.yr)))
plot(qld.yr.lab, qld.bites.yr, type="c", col="black", lwd=1,xlab="",ylab="",main="QLD", ylim=c(0,max(bites.yr)))
plot(wa.yr.lab, wa.bites.yr, type="c", col="black", lwd=1,xlab="",ylab="",main="WA", ylim=c(0,max(bites.yr)))
plot(sa.yr.lab, sa.bites.yr, type="c", col="black", lwd=1,xlab="",ylab="",main="SA", ylim=c(0,max(bites.yr)))
plot(vic.yr.lab, vic.bites.yr, type="c", col="black", lwd=1,xlab="",ylab="",main="VIC", ylim=c(0,max(bites.yr)))
plot(tas.yr.lab, tas.bites.yr, type="c", col="black", lwd=1,xlab="",ylab="",main="TAS", ylim=c(0,max(bites.yr)))
plot(nt.yr.lab, nt.bites.yr, type="c", col="black", lwd=1,xlab="",ylab="",main="NT", ylim=c(0,max(bites.yr)))
par(mfrow=c(1,1))


# plot as matrix (log10)
all.lfit <- lm(log10(bites.yr) ~ yr.lab)
linreg.ER(yr.lab, log10(bites.yr))
nsw.lfit <- lm(log10(nsw.bites.yr) ~ nsw.yr.lab)
linreg.ER(nsw.yr.lab, log10(nsw.bites.yr))
qld.lfit <- lm(log10(qld.bites.yr) ~ qld.yr.lab)
linreg.ER(qld.yr.lab, log10(qld.bites.yr))
wa.lfit <- lm(log10(wa.bites.yr) ~ wa.yr.lab)
linreg.ER(wa.yr.lab, log10(wa.bites.yr))

par(mfrow=c(2,4),pty="m")
plot(yr.lab, log10(bites.yr), type="c", lty=1, lwd=1, col="black",xlab="",ylab="",main="ALL", ylim=c(0,max(log10(bites.yr))))
abline(all.lfit,lty=2, lwd=1)
plot(nsw.yr.lab, log10(nsw.bites.yr), type="c", col="black", lwd=1,xlab="",ylab="",main="NSW", ylim=c(0,max(log10(bites.yr))))
abline(nsw.lfit,lty=2, lwd=1)
plot(qld.yr.lab, log10(qld.bites.yr), type="c", col="black", lwd=1,xlab="",ylab="",main="QLD", ylim=c(0,max(log10(bites.yr))))
abline(qld.lfit,lty=2, lwd=1)
plot(wa.yr.lab, log10(wa.bites.yr), type="c", col="black", lwd=1,xlab="",ylab="",main="WA", ylim=c(0,max(log10(bites.yr))))
abline(wa.lfit,lty=2, lwd=1)
plot(sa.yr.lab, log10(sa.bites.yr), type="c", col="black", lwd=1,xlab="",ylab="",main="SA", ylim=c(0,max(log10(bites.yr))))
plot(vic.yr.lab, log10(vic.bites.yr), type="c", col="black", lwd=1,xlab="",ylab="",main="VIC", ylim=c(0,max(log10(bites.yr))))
plot(tas.yr.lab, log10(tas.bites.yr), type="c", col="black", lwd=1,xlab="",ylab="",main="TAS", ylim=c(0,max(log10(bites.yr))))
plot(nt.yr.lab, log10(nt.bites.yr), type="c", col="black", lwd=1,xlab="",ylab="",main="NT", ylim=c(0,max(log10(bites.yr))))
par(mfrow=c(1,1))


# proportional
pbites.yr <- bites.yr/max(bites.yr)
nsw.pbites.yr <- nsw.bites.yr/max(nsw.bites.yr)
qld.pbites.yr <- qld.bites.yr/max(qld.bites.yr)
wa.pbites.yr <- wa.bites.yr/max(wa.bites.yr)
sa.pbites.yr <- sa.bites.yr/max(sa.bites.yr)
vic.pbites.yr <- vic.bites.yr/max(vic.bites.yr)
tas.pbites.yr <- tas.bites.yr/max(tas.bites.yr)
nt.pbites.yr <- nt.bites.yr/max(nt.bites.yr)

plot(yr.lab, pbites.yr, type="c", lty=1, lwd=2, col="black")
lines(nsw.yr.lab, nsw.pbites.yr, type="c", col="grey", lwd=1)
lines(qld.yr.lab, qld.pbites.yr, type="c", col="red", lwd=1)
lines(wa.yr.lab, wa.pbites.yr, type="c", col="blue", lwd=1)
lines(sa.yr.lab, sa.pbites.yr, type="c", col="green", lwd=1)
lines(vic.yr.lab, vic.pbites.yr, type="c", col="yellow", lwd=1)
lines(tas.yr.lab, tas.pbites.yr, type="c", col="purple", lwd=1)
lines(nt.yr.lab, nt.pbites.yr, type="c", col="pink", lwd=1)

# plot as matrix
par(mfrow=c(2,4),pty="m")
plot(yr.lab, pbites.yr, type="c", lty=1, lwd=1, col="black",xlab="",ylab="",main="ALL")
plot(nsw.yr.lab, nsw.pbites.yr, type="c", col="black", lwd=1,xlab="",ylab="",main="NSW")
plot(qld.yr.lab, qld.pbites.yr, type="c", col="black", lwd=1,xlab="",ylab="",main="QLD")
plot(wa.yr.lab, wa.pbites.yr, type="c", col="black", lwd=1,xlab="",ylab="",main="WA")
plot(sa.yr.lab, sa.pbites.yr, type="c", col="black", lwd=1,xlab="",ylab="",main="SA")
plot(vic.yr.lab, vic.pbites.yr, type="c", col="black", lwd=1,xlab="",ylab="",main="VIC")
plot(tas.yr.lab, tas.pbites.yr, type="c", col="black", lwd=1,xlab="",ylab="",main="TAS")
plot(nt.yr.lab, nt.pbites.yr, type="c", col="black", lwd=1,xlab="",ylab="",main="NT")
par(mfrow=c(1,1))

# spectral analysis (to identify cycles)
# add 2 missing years to time series (1908 & 1970)
fbites.yr <- c(bites.yr[1:8],0,bites.yr[9:69],0,bites.yr[70:118])
acf(coredata(fbites.yr), main="", xlab="lag", ylab="autocorrelation function", demean=T)
fbites.yr.acf <- acf(coredata(fbites.yr),lag.max = 120, main="", xlab="lag", ylab="autocorrelation function")
str(fbites.yr.acf)
pacf((fbites.yr))

fbites.spec <- spectrum(fbites.yr, log="no", spans=c(2,2), plot=FALSE)
delta <- 1/10
specx <- fbites.spec$freq/delta
specy <- 2*fbites.spec$spec
plot(specx, specy, xlab="period (decades)", ylab="spectral density", type="l")

spectrum(fbites.yr, method="ar")
periodogram(fbites.yr, log="no")
spec.pgram(fbites.yr, detrend=T, demean=T, main="")
fbites.yr.pgram <- spec.pgram(fbites.yr, detrend=T, demean=T, main="")

# human population size
pop <- read.csv("auspop.csv", header=T)
dim(pop)
head(pop)

plot(pop$year, pop$AUS, type="l", lwd=1, xlab="", ylab="N", xlim=c(1900,2066))
lines(pop$year, pop$NSW, lty=2, lwd=1, col="grey")
lines(pop$year, pop$QLD, lty=2, lwd=1, col="grey")
lines(pop$year, pop$VIC, lty=2, lwd=1, col="grey")
lines(pop$year, pop$WA, lty=2, lwd=1, col="grey")
lines(pop$year, pop$SA, lty=2, lwd=1, col="grey")
lines(pop$year, pop$TAS, lty=2, lwd=1, col="grey")
lines(pop$year, pop$NT, lty=2, lwd=1, col="grey")
lines(pop$year, pop$ACT, lty=2, lwd=1, col="grey")
abline(v=2020,col="dark grey", lty=3, lwd=2)

# bites per capita
pop1900 <- subset(pop, year > 1899)

all.bites.dat <- data.frame(yr.lab,bites.yr)
colnames(all.bites.dat) <- c("year","bites")
all.pop <- pop1900[,c(1,10)]
all.bites.cap <- merge(all.bites.dat, all.pop, by="year")
all.bites.cap$bites.cap <- all.bites.cap[,2] / all.bites.cap[,3]

nsw.bites.dat <- data.frame(nsw.yr.lab,nsw.bites.yr)
colnames(nsw.bites.dat) <- c("year","bites")
nsw.pop <- pop1900[,1:2]
nsw.bites.cap <- merge(nsw.bites.dat, nsw.pop, by="year")
nsw.bites.cap$bites.cap <- nsw.bites.cap[,2] / nsw.bites.cap[,3]

qld.bites.dat <- data.frame(qld.yr.lab,qld.bites.yr)
colnames(qld.bites.dat) <- c("year","bites")
qld.pop <- pop1900[,c(1,4)]
qld.bites.cap <- merge(qld.bites.dat, qld.pop, by="year")
qld.bites.cap$bites.cap <- qld.bites.cap[,2] / qld.bites.cap[,3]

wa.bites.dat <- data.frame(wa.yr.lab,wa.bites.yr)
colnames(wa.bites.dat) <- c("year","bites")
wa.pop <- pop1900[,c(1,6)]
wa.bites.cap <- merge(wa.bites.dat, wa.pop, by="year")
wa.bites.cap$bites.cap <- wa.bites.cap[,2] / wa.bites.cap[,3]

sa.bites.dat <- data.frame(sa.yr.lab,sa.bites.yr)
colnames(sa.bites.dat) <- c("year","bites")
sa.pop <- pop1900[,c(1,5)]
sa.bites.cap <- merge(sa.bites.dat, sa.pop, by="year")
sa.bites.cap$bites.cap <- sa.bites.cap[,2] / sa.bites.cap[,3]

vic.bites.dat <- data.frame(vic.yr.lab,vic.bites.yr)
colnames(vic.bites.dat) <- c("year","bites")
vic.pop <- pop1900[,c(1,3)]
vic.bites.cap <- merge(vic.bites.dat, vic.pop, by="year")
vic.bites.cap$bites.cap <- vic.bites.cap[,2] / vic.bites.cap[,3]

tas.bites.dat <- data.frame(tas.yr.lab,tas.bites.yr)
colnames(tas.bites.dat) <- c("year","bites")
tas.pop <- pop1900[,c(1,7)]
tas.bites.cap <- merge(tas.bites.dat, tas.pop, by="year")
tas.bites.cap$bites.cap <- tas.bites.cap[,2] / tas.bites.cap[,3]

nt.bites.dat <- data.frame(nt.yr.lab,nt.bites.yr)
colnames(nt.bites.dat) <- c("year","bites")
nt.pop <- pop1900[,c(1,8)]
nt.bites.cap <- merge(nt.bites.dat, nt.pop, by="year")
nt.bites.cap$bites.cap <- nt.bites.cap[,2] / nt.bites.cap[,3]

# plot as matrix
par(mfrow=c(2,4),pty="m")
plot(all.bites.cap$year, all.bites.cap$bites.cap, type="c", lty=1, lwd=1, col="black",xlab="",ylab="",main="ALL")
plot(nsw.bites.cap$year, nsw.bites.cap$bites.cap, type="c", lty=1, lwd=1, col="black",xlab="",ylab="",main="NSW")
plot(qld.bites.cap$year, qld.bites.cap$bites.cap, type="c", lty=1, lwd=1, col="black",xlab="",ylab="",main="QLD")
plot(wa.bites.cap$year, wa.bites.cap$bites.cap, type="c", lty=1, lwd=1, col="black",xlab="",ylab="",main="WA")
plot(sa.bites.cap$year, sa.bites.cap$bites.cap, type="c", lty=1, lwd=1, col="black",xlab="",ylab="",main="SA")
plot(vic.bites.cap$year, vic.bites.cap$bites.cap, type="c", lty=1, lwd=1, col="black",xlab="",ylab="",main="VIC")
plot(tas.bites.cap$year, tas.bites.cap$bites.cap, type="c", lty=1, lwd=1, col="black",xlab="",ylab="",main="TAS")
plot(nt.bites.cap$year, nt.bites.cap$bites.cap, type="c", lty=1, lwd=1, col="black",xlab="",ylab="",main="NT")
par(mfrow=c(1,1))

# plot as matrix (with secondary y-axis with pop trend by region)
par(mar=c(5,5,2,5))
plot(all.bites.cap$year, all.bites.cap$bites.cap, type="c", lty=1, lwd=1, col="black", xlab="",ylab="bites/person",main="ALL")
par(new=T)
with(all.bites.cap, plot(all.bites.cap$year, log10(all.bites.cap[,3]), axes=F, type="l", ylab=NA, xlab=NA, lty=2, col="grey"))
axis(side=4)
mtext(side=4, line=3, "log N", col="grey")

par(mar=c(5,5,2,5))
plot(nsw.bites.cap$year, nsw.bites.cap$bites.cap, type="c", lty=1, lwd=1, col="black",xlab="",ylab="bites/person",main="NSW")
par(new=T)
with(nsw.bites.cap, plot(nsw.bites.cap$year, log10(nsw.bites.cap[,3]), axes=F, type="l", ylab=NA, xlab=NA, lty=2, col="grey", ylim=c(min(log10(nt.bites.cap[,3]), na.rm=T),max(log10(nsw.bites.cap[,3])))))
axis(side=4)
mtext(side=4, line=3, "log N", col="grey")

par(mar=c(5,5,2,5))
plot(qld.bites.cap$year, qld.bites.cap$bites.cap, type="c", lty=1, lwd=1, col="black",xlab="",ylab="",main="QLD")
par(new=T)
with(qld.bites.cap, plot(qld.bites.cap$year, log10(qld.bites.cap[,3]), axes=F, type="l", ylab=NA, xlab=NA, lty=2, col="grey", ylim=c(min(log10(nt.bites.cap[,3]), na.rm=T),max(log10(nsw.bites.cap[,3])))))
axis(side=4)
mtext(side=4, line=3, "log N", col="grey")

plot(wa.bites.cap$year, wa.bites.cap$bites.cap, type="c", lty=1, lwd=1, col="black",xlab="",ylab="",main="WA")
par(new=T)
with(wa.bites.cap, plot(wa.bites.cap$year, log10(wa.bites.cap[,3]), axes=F, type="l", ylab=NA, xlab=NA, lty=2, col="grey", ylim=c(min(log10(nt.bites.cap[,3]), na.rm=T),max(log10(nsw.bites.cap[,3])))))
axis(side=4)
mtext(side=4, line=3, "log N", col="grey")

plot(sa.bites.cap$year, sa.bites.cap$bites.cap, type="c", lty=1, lwd=1, col="black",xlab="",ylab="bites/person",main="SA")
par(new=T)
with(sa.bites.cap, plot(sa.bites.cap$year, log10(sa.bites.cap[,3]), axes=F, type="l", ylab=NA, xlab=NA, lty=2, col="grey", ylim=c(min(log10(nt.bites.cap[,3]), na.rm=T),max(log10(nsw.bites.cap[,3])))))
axis(side=4)
mtext(side=4, line=3, "log N", col="grey")

plot(vic.bites.cap$year, vic.bites.cap$bites.cap, type="c", lty=1, lwd=1, col="black",xlab="",ylab="",main="VIC")
par(new=T)
with(vic.bites.cap, plot(vic.bites.cap$year, log10(vic.bites.cap[,3]), axes=F, type="l", ylab=NA, xlab=NA, lty=2, col="grey", ylim=c(min(log10(nt.bites.cap[,3]), na.rm=T),max(log10(nsw.bites.cap[,3])))))
axis(side=4)
mtext(side=4, line=3, "log N", col="grey")

plot(tas.bites.cap$year, tas.bites.cap$bites.cap, type="c", lty=1, lwd=1, col="black",xlab="",ylab="",main="TAS")
par(new=T)
with(tas.bites.cap, plot(tas.bites.cap$year, log10(tas.bites.cap[,3]), axes=F, type="l", ylab=NA, xlab=NA, lty=2, col="grey", ylim=c(min(log10(nt.bites.cap[,3]), na.rm=T),max(log10(nsw.bites.cap[,3])))))
axis(side=4)
mtext(side=4, line=3, "log N", col="grey")

plot(nt.bites.cap$year, nt.bites.cap$bites.cap, type="c", lty=1, lwd=1, col="black",xlab="",ylab="",main="NT")
par(new=T)
with(nt.bites.cap, plot(nt.bites.cap$year, log10(nt.bites.cap[,3]), axes=F, type="l", ylab=NA, xlab=NA, lty=2, col="grey", ylim=c(min(log10(nt.bites.cap[,3]), na.rm=T),max(log10(nsw.bites.cap[,3])))))
axis(side=4)
mtext(side=4, line=3, "log N", col="grey")

par(mfrow=c(1,1))




# scale & trends
all.bites.cap$bites.cap.sc <- scale(all.bites.cap$bites.cap, scale=T, center=F)
nsw.bites.cap$bites.cap.sc <- scale(nsw.bites.cap$bites.cap, scale=T, center=F)
qld.bites.cap$bites.cap.sc <- scale(qld.bites.cap$bites.cap, scale=T, center=F)
wa.bites.cap$bites.cap.sc <- scale(wa.bites.cap$bites.cap, scale=T, center=F)
sa.bites.cap$bites.cap.sc <- scale(sa.bites.cap$bites.cap, scale=T, center=F)
vic.bites.cap$bites.cap.sc <- scale(vic.bites.cap$bites.cap, scale=T, center=F)
tas.bites.cap$bites.cap.sc <- scale(tas.bites.cap$bites.cap, scale=T, center=F)
nt.bites.cap$bites.cap.sc <- scale(nt.bites.cap$bites.cap, scale=T, center=F)

all.bcs <- lm(all.bites.cap$bites.cap.sc ~ all.bites.cap$year)
nsw.bcs <- lm(nsw.bites.cap$bites.cap.sc ~ nsw.bites.cap$year)
qld.bcs <- lm(qld.bites.cap$bites.cap.sc ~ qld.bites.cap$year)
wa.bcs <- lm(wa.bites.cap$bites.cap.sc ~ wa.bites.cap$year)
sa.bcs <- lm(sa.bites.cap$bites.cap.sc ~ sa.bites.cap$year)
vic.bcs <- lm(vic.bites.cap$bites.cap.sc ~ vic.bites.cap$year)
tas.bcs <- lm(tas.bites.cap$bites.cap.sc ~ tas.bites.cap$year)
nt.bcs <- lm(nt.bites.cap$bites.cap.sc ~ nt.bites.cap$year)

linreg.ER(all.bites.cap$year, all.bites.cap$bites.cap.sc)
linreg.ER(nsw.bites.cap$year, nsw.bites.cap$bites.cap.sc)
linreg.ER(qld.bites.cap$year, qld.bites.cap$bites.cap.sc)
linreg.ER(wa.bites.cap$year, wa.bites.cap$bites.cap.sc)
linreg.ER(sa.bites.cap$year, sa.bites.cap$bites.cap.sc)
linreg.ER(vic.bites.cap$year, vic.bites.cap$bites.cap.sc)
linreg.ER(tas.bites.cap$year, tas.bites.cap$bites.cap.sc)
linreg.ER(nt.bites.cap$year, nt.bites.cap$bites.cap.sc)

par(mfrow=c(2,4),pty="m")
plot(all.bites.cap$year, all.bites.cap$bites.cap.sc, type="c", lty=1, lwd=1, col="black",xlab="",ylab="",main="ALL")
abline(all.bcs,lty=2)
plot(nsw.bites.cap$year, nsw.bites.cap$bites.cap.sc, type="c", lty=1, lwd=1, col="black",xlab="",ylab="",main="NSW")
abline(nsw.bcs,lty=2)
plot(qld.bites.cap$year, qld.bites.cap$bites.cap.sc, type="c", lty=1, lwd=1, col="black",xlab="",ylab="",main="QLD")
abline(qld.bcs,lty=2)
plot(wa.bites.cap$year, wa.bites.cap$bites.cap.sc, type="c", lty=1, lwd=1, col="black",xlab="",ylab="",main="WA")
abline(wa.bcs,lty=2)
plot(sa.bites.cap$year, sa.bites.cap$bites.cap.sc, type="c", lty=1, lwd=1, col="black",xlab="",ylab="",main="SA")
abline(sa.bcs,lty=2)
plot(vic.bites.cap$year, vic.bites.cap$bites.cap.sc, type="c", lty=1, lwd=1, col="black",xlab="",ylab="",main="VIC")
abline(vic.bcs,lty=2)
plot(tas.bites.cap$year, tas.bites.cap$bites.cap.sc, type="c", lty=1, lwd=1, col="black",xlab="",ylab="",main="TAS")
abline(tas.bcs,lty=2)
plot(nt.bites.cap$year, nt.bites.cap$bites.cap.sc, type="c", lty=1, lwd=1, col="black",xlab="",ylab="",main="NT")
abline(nt.bcs,lty=2)
par(mfrow=c(1,1))

# multi-model comparisons
# intercept-only vs. linear vs. quadratic
# response: scaled bites/capita
# GLM 

# define list
dat.list <- list(all.bites.cap, nsw.bites.cap, qld.bites.cap, wa.bites.cap, sa.bites.cap, vic.bites.cap, tas.bites.cap, nt.bites.cap)

par(mfrow=c(2,4), pty="m")
for (i in 1:length(dat.list)) {

  # set data
  data.use <- dat.list[[i]]
  data.use$year2 <- data.use$year^2
  data.use$year3 <- data.use$year^3
  
  # model set
  m1 <- "bites.cap.sc ~ year + year2 + year3"
  m2 <- "bites.cap.sc ~ year + year2"
  m3 <- "bites.cap.sc ~ year"
  m4 <- "bites.cap.sc ~ 1"
  
  ## Make model vector
  mod.vec <- c(m1,m2,m3,m4)
  
  ## Define n.mod
  n.mod <- length(mod.vec)
  
  # Model fitting and logLik output loop
  Modnum <- length(mod.vec)
  LL.vec <- SaveCount <- AICc.vec <- BIC.vec <- pc.dev.vec <- k.vec <- rep(0,Modnum)
  mod.list <- list()
  mod.num <- seq(1,Modnum,1)
  
  for(j in 1:Modnum) {
    fit <- glm(as.formula(mod.vec[j]),family=Gamma(link="identity"), data=data.use, na.action=na.omit)
    assign(paste("fit",j,sep=""), fit)
    mod.list[[j]] <- fit
    print(j)
  }
  
  sumtable <- aicW(mod.list, finite = TRUE, null.model = NULL, order = F)
  row.names(sumtable) <- mod.vec
  summary.table <- sumtable[order(sumtable[,4],decreasing=F),1:9]
  summary.table
  
  mod.vec.lab <- c("cubic","quadratic","linear","intercept")
  top.mod <- mod.vec.lab[which(mod.vec == row.names(summary.table)[1])]
  
  plot(data.use$year, data.use$bites.cap.sc, type="c", lwd=0.7, pch=19, col="dark grey", cex=0.6, xlab="", ylab="", main=colnames(data.use)[3], sub=paste(top.mod,": ", "wAICc=",summary.table[1,5],"; %DE=",summary.table[1,9],sep=""))
  abline(mod.list[[4]], lty=2, lwd=1)
  abline(mod.list[[3]], lty=2, lwd=2, col="red")
  lin.pred <- coef(mod.list[[3]])[1] + (coef(mod.list[[3]])[2]*data.use$year)
  quad.pred <- coef(mod.list[[2]])[1] + (coef(mod.list[[2]])[2]*data.use$year) + (coef(mod.list[[2]])[3]*data.use$year2)
  lines(data.use$year, quad.pred, lty=2, lwd=2, col="green")
  cub.pred <- coef(mod.list[[1]])[1] + (coef(mod.list[[1]])[2]*data.use$year) + (coef(mod.list[[1]])[3]*data.use$year2) + (coef(mod.list[[1]])[4]*data.use$year3)
  lines(data.use$year, cub.pred, lty=2, lwd=2, col="purple")
  
}
par(mfrow=c(1,1), pty="m")

## NCEI PDO (Pacific Decadal Oscillation) from https://www.ncdc.noaa.gov/teleconnections/pdo/
pdo <- read.csv("pdo.csv", header=T)
dim(pdo)
head(pdo)
pdo$year.mo <- pdo$year + pdo$month/12
pdo1900 <- subset(pdo, year > 1899)
plot(pdo1900$year.mo, pdo1900$PDO, type="l", xlab="", ylab="PDO")


##############################################################################################
# select dataset (select relevant state for prediction of number of bites potentiall averted)
data.use <- all.bites.cap
#data.use <- nsw.bites.cap
#data.use <- qld.bites.cap
#data.use <- wa.bites.cap

data.use$year2 <- data.use$year^2
data.use$year3 <- data.use$year^3
  
  ## modify for hypothesised underreporting?
  underreport <- 0 # no
  #underreport <- 1 # yes
  
  ## incorporate a scenario of temporal change in reporting probability
  ## assume x% underreporting in early part of 20th Century, followed by rise through television era (1950-1960s), and maxing by internet age (1990s)
  ur.max <- 0.40 # early-Century underreporting %
  ur.min <- 0
  inflex <- mean(c((1-ur.max),(1-ur.min)))
  uri <- (0.2*(inflex - (1-ur.max))) + (1-ur.max)
  iur <- (0.85*((1-ur.min) - inflex)) + (inflex)
  
  urmax.yr <- 1900
  urinfl.yr <- 1950
  urmi.yr <- (0.55*(urinfl.yr - (urmax.yr))) + (urmax.yr)  
  urmin.yr <- 1990
  iurm.yr <- (0.55*(urmin.yr - (urinfl.yr))) + (urinfl.yr)
  
  rr.vec <- c(1-ur.max, uri, inflex, iur, 1-ur.min)
  rryr.vec <- c(urmax.yr, urmi.yr, urinfl.yr, iurm.yr, urmin.yr)
  plot(rryr.vec, rr.vec, type="b", pch=19)
  rr.dat <- data.frame(rryr.vec, rr.vec)
  colnames(rr.dat) <- c("rryr", "rr")
  rr.dat
  
  # fit function
  # y = γ + ((1-γ)/(1+exp(-α-(β*log(x))))) (dose-response log-logistic)
  param.init <- c(5.803e-01, -1.076e03, 1.42e02)
  fit.rr <- nls(rr ~ g + ((1-g)/(1+exp(-a-(b*log(rryr))))), 
                 data = rr.dat,
                 algorithm = "port",
                 start = c(g = param.init[1], a = param.init[2], b = param.init[3]),
                 trace = TRUE,      
                 nls.control(maxiter = 1000, tol = 1e-05, minFactor = 1/1024))
  fit.rr.summ <- summary(fit.rr)
  plot(rr.dat$rryr, rr.dat$rr, pch=19, xlab="year", ylab="reporting probability", xlim=c(1900,2020))
  rryr.cont <- seq(1900,2020,1)
  pred.rr.fx <- coef(fit.rr)[1] + ((1-coef(fit.rr)[1])/(1+exp(-coef(fit.rr)[2]-(coef(fit.rr)[3]*log(rryr.cont)))))
  lines(rryr.cont, pred.rr.fx,lty=2,lwd=3,col="red")
  
  g.rr <- coef(fit.rr)[1]
  a.rr <- coef(fit.rr)[2]
  b.rr <- coef(fit.rr)[3]
  
  ## assume x% underreporting in early part of 20th Century, followed by rise through television era (1950-1960s), and maxing by social media age (2000s)
  ur.max <- 0.60 # early-Century underreporting %
  ur.min <- 0
  inflex <- mean(c((1-ur.max),(1-ur.min)))
  uri <- (0.2*(inflex - (1-ur.max))) + (1-ur.max)
  iur <- (0.85*((1-ur.min) - inflex)) + (inflex)
  
  urmax.yr <- 1900
  urinfl.yr <- 1950
  urmi.yr <- (0.55*(urinfl.yr - (urmax.yr))) + (urmax.yr)  
  urmin.yr <- 2000
  iurm.yr <- (0.55*(urmin.yr - (urinfl.yr))) + (urinfl.yr)
  
  rr.vec <- c(1-ur.max, uri, inflex, iur, 1-ur.min)
  rryr.vec <- c(urmax.yr, urmi.yr, urinfl.yr, iurm.yr, urmin.yr)
  plot(rryr.vec, rr.vec, type="b", pch=19)
  rr.dat <- data.frame(rryr.vec, rr.vec)
  colnames(rr.dat) <- c("rryr", "rr")
  rr.dat
  
  # fit function
  # y = γ + ((1-γ)/(1+exp(-α-(β*log(x))))) (dose-response log-logistic)
  param.init <- c(5.803e-01, -1.076e03, 1.42e02)
  fit.rr <- nls(rr ~ g + ((1-g)/(1+exp(-a-(b*log(rryr))))), 
                data = rr.dat,
                algorithm = "port",
                start = c(g = param.init[1], a = param.init[2], b = param.init[3]),
                trace = TRUE,      
                nls.control(maxiter = 1000, tol = 1e-05, minFactor = 1/1024))
  fit.rr.summ <- summary(fit.rr)
  plot(rr.dat$rryr, rr.dat$rr, pch=19, xlab="year", ylab="reporting probability", xlim=c(1900,2020))
  rryr.cont <- seq(1900,2020,1)
  pred.rr.fx <- coef(fit.rr)[1] + ((1-coef(fit.rr)[1])/(1+exp(-coef(fit.rr)[2]-(coef(fit.rr)[3]*log(rryr.cont)))))
  lines(rryr.cont, pred.rr.fx,lty=2,lwd=3,col="red")
  
  rr.out <- data.frame(rryr.cont, pred.rr.fx)
  setwd("~/Documents/Papers/Fish/Sharks/Shark bites/results/")
  write.table(rr.out,file="rrcorrfunc2.csv", sep=",", row.names = F, col.names = T)
  
  g.rr <- coef(fit.rr)[1]
  a.rr <- coef(fit.rr)[2]
  b.rr <- coef(fit.rr)[3]
  
  if (underreport == 1) {
    rr.pred <- g.rr + ((1-g.rr)/(1+exp(-a.rr-(b.rr*log(data.use$year)))))
    bcsc.corr <- data.use$bites.cap.sc + ((1-rr.pred)*data.use$bites.cap.sc)
    par(mar = c(5, 5, 3, 4))
    plot(data.use$year, data.use$bites.cap.sc, type="c", lwd=2, pch=19, ylim=c(0,max(bcsc.corr)), col="black", cex=0.6, xlab="", ylab="scaled bites/capita", main=colnames(data.use)[3])
    lines(data.use$year, bcsc.corr, lty=2, col="grey")
    bcsc.corr.out <- data.frame(data.use$year, data.use$bites.cap.sc, bcsc.corr)
    colnames(bcsc.corr.out) <- c("year","bcsc.orig","bcsc.corr")
    setwd("~/Documents/Papers/Fish/Sharks/Shark bites/results/")
    write.table(bcsc.corr.out,file="bcsc.corr.csv", sep=",", row.names = F, col.names = T)
    data.use$bites.cap.sc <- bcsc.corr # replace bites.cap.sc
  }
  
  if (underreport == 0) {
    par(mar = c(5, 5, 3, 4))
    plot(data.use$year, data.use$bites.cap.sc, type="c", lwd=2, pch=19, ylim=c(0,max(data.use$bites.cap.sc)), col="black", cex=0.6, xlab="", ylab="scaled bites/capita", main=colnames(data.use)[3])
  }
  
  # autocorrelation
  acf(coredata(data.use$bites.cap.sc),lag.max = 120, main="", xlab="lag", ylab="autocorrelation function", demean=T)
  bitescapsc.acf <- acf(coredata(data.use$bites.cap.sc),lag.max = 120, main="", xlab="lag", ylab="autocorrelation function", demean=T)
  str(bitescapsc.acf)
  acf.out <- data.frame(bitescapsc.acf$lag, bitescapsc.acf$acf)
  alpha <- 0.95
  conf.lims <- c(-1,1)*qnorm((1 + alpha)/2)/sqrt(bitescapsc.acf$n.used)
  conf.lims
  
  ## need to interpolate for pacf() with approx()
  yr.int.vec <- seq(from=min(data.use$year), to=max(data.use$year), by=1)
  bites.cap.sc.intp <- approx(data.use$year, data.use$bites.cap.sc, xout=yr.int.vec, method="linear")
  dat.intp <- data.frame(yr.int.vec, bites.cap.sc.intp$y)
  plot(dat.intp$yr.int.vec, dat.intp$bites.cap.sc.intp.y, type="l")
  
  pacf(dat.intp$bites.cap.sc.intp, lag.max=120)
  bitescapsc.pacf <- pacf(dat.intp$bites.cap.sc.intp.y, lag.max=120)
  str(bitescapsc.pacf)
  pacf.out <- data.frame(bitescapsc.pacf$lag, bitescapsc.pacf$acf)
  conf.lims <- c(-1,1)*qnorm((1 + alpha)/2)/sqrt(bitescapsc.pacf$n.used)
  conf.lims
  
  # model set
  m1 <- "bites.cap.sc ~ year + year2 + year3"
  m2 <- "bites.cap.sc ~ year + year2"
  m3 <- "bites.cap.sc ~ year"
  m4 <- "bites.cap.sc ~ 1"
  
  ## Make model vector
  mod.vec <- c(m1,m2,m3,m4)
  
  ## Define n.mod
  n.mod <- length(mod.vec)
  
  # Model fitting and logLik output loop
  Modnum <- length(mod.vec)
  LL.vec <- SaveCount <- AICc.vec <- BIC.vec <- pc.dev.vec <- k.vec <- rep(0,Modnum)
  mod.list <- list()
  mod.num <- seq(1,Modnum,1)
  
  for(i in 1:Modnum) {
    fit <- glm(as.formula(mod.vec[i]),family=Gamma(link="identity"), data=data.use, na.action=na.omit)
    assign(paste("fit",i,sep=""), fit)
    mod.list[[i]] <- fit
    print(i)
  }
  
  sumtable <- aicW(mod.list, finite = TRUE, null.model = NULL, order = F)
  row.names(sumtable) <- mod.vec
  summary.table <- sumtable[order(sumtable[,4],decreasing=F),1:9]
  summary.table
  
  mod.vec.lab <- c("cubic","quadratic","linear","intercept")
  top.mod <- mod.vec.lab[which(mod.vec == row.names(summary.table)[1])]
  
  par(mar = c(5, 5, 3, 4))
  plot(data.use$year, data.use$bites.cap.sc, type="c", lwd=2, pch=19, ylim=c(-1.5,4), col="black", cex=0.6, xlab="", ylab="scaled bites/capita", main=colnames(data.use)[3], sub=paste(top.mod,": ", "wAICc=",summary.table[1,5],"; %DE=",summary.table[1,9],sep=""))
  cub.pred <- coef(mod.list[[1]])[1] + (coef(mod.list[[1]])[2]*data.use$year) + (coef(mod.list[[1]])[3]*data.use$year2) + (coef(mod.list[[1]])[4]*data.use$year3)
  lines(data.use$year, cub.pred, lty=2, lwd=2, col="purple")
  par(new = TRUE)
  plot(pdo1900$year.mo, pdo1900$PDO, type = "l", xaxt = "n", yaxt = "n", lwd=0.4, col="grey", ylab = "", xlab = "", ylim=c(-4,4))
  mtext("PDO", side = 4, line = 3)
  axis(side = 4)
  
  # Southern Oscillation Index (SOI) from http://www.bom.gov.au/climate/current/soihtm1.shtml
  soi <- read.csv("soi.csv", header=T)
  dim(soi)
  head(soi)
  soi.mat <- matrix(data=0,nrow = dim(soi)[1]*12, ncol = 2)
  st.seq <- seq(1, (dim(soi)[1]*12), 12)
  en.seq <- seq(12, (dim(soi)[1]*12), 12)
  
  for (j in 1:dim(soi)[1]) {
    soi.mat[st.seq[j]:en.seq[j], 1] <- rep(soi[j,1], 12) + seq(1,12,1)/12
    soi.mat[st.seq[j]:en.seq[j], 2] <- as.numeric(soi[j,2:13])
  }
  soi.dat <- as.data.frame(soi.mat)
  colnames(soi.dat) <- c("year.mo", "soi")
  
  soi1900 <- subset(soi.dat, year.mo >= 1900)
  plot(soi1900$year.mo, soi1900$soi, type="l", xlab="", ylab="SOI", cex=0.7, lwd=0.7, col="grey")
  rm.win <- 5 # years
  soi1900$soi.rm <- runMean(soi1900$soi, rm.win*12)
  lines(soi1900$year.mo, soi1900$soi.rm, lty=1, lwd=2, col="red")
  
  # PDO vs. SOI
  plot(pdo1900$PDO, soi1900$soi[-length(soi1900$soi)], pch=19, cex=0.5)
  fit.soi.pdo <- lm(soi1900$soi[-length(soi1900$soi)] ~ pdo1900$PDO)
  abline(fit.soi.pdo, col="red")
  linreg.ER(pdo1900$PDO, soi1900$soi[-length(soi1900$soi)])
  pdosoi.out <- data.frame(pdo1900$PDO, soi1900$soi[-length(soi1900$soi)], pdo1900$PDO, soi1900$soi.rm[-length(soi1900$soi)])
  colnames(pdosoi.out) <- c("PDO","SOI","SOIrm")
  
  par(mar = c(5, 5, 3, 4))
  plot(data.use$year, data.use$bites.cap.sc, type="c", lwd=2, pch=19, ylim=c(-1.5,3.7), col="black", cex=0.6, xlab="", ylab="scaled bites/capita", main=colnames(data.use)[3], sub=paste(top.mod,": ", "wAICc=",summary.table[1,5],"; %DE=",summary.table[1,9],sep=""))
  cub.pred <- coef(mod.list[[1]])[1] + (coef(mod.list[[1]])[2]*data.use$year) + (coef(mod.list[[1]])[3]*data.use$year2) + (coef(mod.list[[1]])[4]*data.use$year3)
  lines(data.use$year, cub.pred, lty=2, lwd=2, col="purple")
  par(new = TRUE)
  plot(soi1900$year.mo, soi1900$soi, type = "l", xaxt = "n", yaxt = "n", lwd=0.4, col="grey", ylab = "", xlab = "", ylim=c(-38,45))
  mtext("SOI", side = 4, line = 3)
  axis(side = 4)
  
  ## predict to 2066
  yr.proj <- 2020:2066; yr.proj2 <- yr.proj^2; yr.proj3 <- yr.proj^3
  cub.pred3 <- coef(mod.list[[1]])[1] + (coef(mod.list[[1]])[2]*yr.proj) + (coef(mod.list[[1]])[3]*yr.proj2) + (coef(mod.list[[1]])[4]*yr.proj3)
  
  # soi average
  soi.mn <- apply(soi[,2:13], 1, mean)
  soi.mn.dat <- data.frame(soi$Year, soi.mn)
  colnames(soi.mn.dat) <- c("year", "soi.mn")
  soi.mn.dat1900 <- subset(soi.mn.dat, year > 1899)
  
  plot(soi.mn.dat1900$year, soi.mn.dat1900$soi.mn, type="l", xlab="", ylab="SOI", cex=0.7, lwd=0.7, col="grey")
  rm.win <- 5 # years
  soi.mn.dat1900$soi.mn.rm <- runMean(soi.mn.dat1900$soi, rm.win)
  lines(soi.mn.dat1900$year, soi.mn.dat1900$soi.mn.rm, lty=1, lwd=2, col="red")
  
  
  par(mar = c(5, 5, 3, 4))
  plot(data.use$year, data.use$bites.cap.sc, type="c", lwd=2, pch=19, ylim=c(-1,3.7), col="black", cex=0.6, xlab="", ylab="scaled bites/capita", main=colnames(data.use)[3], sub=paste(top.mod,": ", "wAICc=",summary.table[1,5],"; %DE=",summary.table[1,9],sep=""))
  cub.pred <- coef(mod.list[[1]])[1] + (coef(mod.list[[1]])[2]*data.use$year) + (coef(mod.list[[1]])[3]*data.use$year2) + (coef(mod.list[[1]])[4]*data.use$year3)
  lines(data.use$year, cub.pred, lty=2, lwd=2, col="purple")
  par(new = TRUE)
  plot(soi.mn.dat1900$year, soi.mn.dat1900$soi.mn, type = "l", xaxt = "n", yaxt = "n", lwd=0.4, col="grey", ylab = "", xlab = "", ylim=c(-20,25))
  mtext("annual mean SOI", side = 4, line = 3)
  axis(side = 4)
  
  
  data2use <- soi.mn.dat1900
  data2use$year2 <- data2use$year^2
  data2use$year3 <- data2use$year^3
  data2use$soi.sc <- scale(data2use$soi.mn, scale=T, center=T)
  data2use$soi.sc.rm <- runMean(data2use$soi.sc, 5)
  
  # model set
  m1 <- "soi.sc.rm ~ year + year2 + year3"
  m2 <- "soi.sc.rm ~ year + year2"
  m3 <- "soi.sc.rm ~ year"
  m4 <- "soi.sc.rm ~ 1"
  
  ## Make model vector
  mod.vec <- c(m1,m2,m3,m4)
  
  ## Define n.mod
  n.mod <- length(mod.vec)
  
  # Model fitting and logLik output loop
  Modnum <- length(mod.vec)
  LL.vec <- SaveCount <- AICc.vec <- BIC.vec <- pc.dev.vec <- k.vec <- rep(0,Modnum)
  mod.list <- list()
  mod.num <- seq(1,Modnum,1)
  
  for(i in 1:Modnum) {
    fit <- glm(as.formula(mod.vec[i]),family=gaussian(link="identity"), data=data2use, na.action=na.omit)
    assign(paste("fit",i,sep=""), fit)
    mod.list[[i]] <- fit
    print(i)
  }
  
  sumtable <- aicW(mod.list, finite = TRUE, null.model = NULL, order = F)
  row.names(sumtable) <- mod.vec
  summary.table <- sumtable[order(sumtable[,4],decreasing=F),1:9]
  summary.table
  
  plot(data2use$year, data2use$soi.sc.rm, type="c", lwd=0.7, pch=19, col="dark grey", cex=0.6, xlab="", ylab="")
  cub.pred2 <- coef(mod.list[[1]])[1] + (coef(mod.list[[1]])[2]*data2use$year) + (coef(mod.list[[1]])[3]*data2use$year2) + (coef(mod.list[[1]])[4]*data2use$year3)
  lines(data2use$year, cub.pred2, lty=2, lwd=2, col="purple")
  
  par(mar = c(5, 5, 3, 4))
  plot(data.use$year, data.use$bites.cap.sc, type="c", lwd=2, pch=19, ylim=c(-0.7,4.2), col="black", cex=0.6, xlab="", ylab="scaled bites/capita", main=colnames(data.use)[3], sub=paste(top.mod,": ", "wAICc=",summary.table[1,5],"; %DE=",summary.table[1,9],sep=""))
  lines(data.use$year, cub.pred, lty=2, lwd=2, col="purple")
  par(new = TRUE)
  plot(data2use$year, data2use$soi.sc.rm, type = "l", xaxt = "n", yaxt = "n", lwd=0.7, col="grey", ylab = "", xlab = "", ylim=c(-1.2,1.8))
  mtext("5-yr running mean annual SOI", side = 4, line = 3)
  axis(side = 4)
  lines(data2use$year, cub.pred2, lty=2, lwd=2, col="grey")
  
  ## coherence (wavelets)
  # first time series
  len1 <- length(data.use$year)
  st.yr1 <- data.use$year[1]
  en.yr1 <- data.use$year[len1]
  yr.seq1 <- as.data.frame(st.yr1:en.yr1)
  colnames(yr.seq1) <- "year"
  t1 <- as.data.frame(cbind(data.use$year, data.use$bites.cap.sc))
  colnames(t1) <- c("year","cub.pred")
  t1.mrg <- merge(yr.seq1,t1,by="year",all.x = T)
  t1.interp <- as.matrix(na.approx(t1.mrg))
  t1.intp <- cbind(t1.interp[,1], t1.interp[,2])
  dim(t1.intp)
  str(t1.intp)
  
  # second time series
  len2 <- length(data2use$year)
  t2 <- as.matrix(cbind(data2use$year, data2use$soi.sc.rm))
  dim(t2)
  str(t2)
  
  # first year with non-NA values
  t1.1st <- na.omit(t1.intp)[1,1]
  t2.1st <- na.omit(t2)[1,1]
  lowest.1st <- max(t1.1st, t2.1st)
  lowest.t1 <- which(t1.intp[,1] == lowest.1st)
  lowest.t2 <- which(t2[,1] == lowest.1st)
  
  # detrend
  t1.det <- cbind(t1.intp[lowest.t1:length(t1.intp[,1]),1], detrend(t1.intp[,2], tt="linear")[lowest.t1:length(t1.intp[,1]),1])
  plot(t1.det[,1], t1.det[,2], type="l",xlab=NA,ylab=NA)
  par(new=T)
  t2.det <- cbind(t2[lowest.t2:length(t2[,1]),1], detrend(t2[,2], tt="linear")[lowest.t2:length(t2[,1]),1])
  plot(t2.det[,1], t2.det[,2], type="l",col="red",lty=2,xlab=NA,ylab=NA, axes=F)
  
  # estimate autoregression coefficients for each time series
  ar1 <- arima(t1.det[,2], method="ML")
  ar2 <- arima(t2.det[,2], method="ML")
  
  # plot wavelet correspondence
  t1t2.wtc <- wtc(t1.det, t2.det, mother="morlet", lag1=as.vector(c(coef(ar1),coef(ar2))), nrands=1000, sig.level=0.05, sig.test=1)
  par(mar = c(5, 5, 3, 10))
  plot(t1t2.wtc, plot.cb=T, plot.phase=F, ylab="period (years)", xlab=NA, cex.lab=1.5, cex.axis=1.5)
  t1t2.wtc$signif
  
  
  ## predict to 2066
  ## fit sinusoidal model
  dat.sin <- data.frame(data.use$year,cub.pred)
  colnames(dat.sin) <- c("year","biteppscpred")
  
  # fit function
  # sinusoidal
  # y = a + (b*cos(c*x + d))
  param.init <- c(8.13e-01, 4.64e-01, 5.78e-02, -1.06e01)
  fit.sin <- nls(biteppscpred ~ a1 + (b2*cos(c3*year + d4)), 
                  data = dat.sin,
                  algorithm = "default",
                  start = c(a1 = param.init[1], b2 = param.init[2], c3 = param.init[3], d4 = param.init[4]),
                  trace = TRUE,      
                  nls.control(maxiter = 1000, tol = 1e-05, minFactor = 1/1024))
  par(mar = c(5, 5, 3, 4))
  plot(data.use$year, data.use$bites.cap.sc, type="c", lwd=2, pch=19, ylim=c(-0.7,4.2), col="black", cex=0.6, xlim=c(1900,2066), xlab="", ylab="scaled bites/capita", main=colnames(data.use)[3], sub=paste(top.mod,": ", "wAICc=",summary.table[1,5],"; %DE=",summary.table[1,9],sep=""))
  sin.pred <- as.numeric(coef(fit.sin)[1]) + (as.numeric(coef(fit.sin)[2])*cos(as.numeric(coef(fit.sin)[3])*yr.proj + as.numeric(coef(fit.sin)[4])))
  lines(yr.proj, sin.pred, lty=2, lwd=1, col="red")
  lines(data.use$year, cub.pred, lty=2, lwd=2, col="purple")
  par(new = TRUE)
  plot(data2use$year, data2use$soi.sc.rm, type = "l", xaxt = "n", yaxt = "n", lwd=0.7, col="grey", ylab = "", xlab = "", xlim=c(1900, 2066), ylim=c(-1.2,1.8))
  mtext("5-yr running mean annual SOI", side = 4, line = 3)
  axis(side = 4)
  #lines(data2use$year, cub.pred2, lty=2, lwd=2, col="grey")
  
  ## calculate residuals from best-fit (cubic, quadratic)
  ## to predict stochastic departure from sinusoidal projection to 2066
  bfit.resid <- data.use$bites.cap.sc - cub.pred
  plot(cub.pred, bfit.resid, pch=19, cex=0.7)
  fit.resid.cube <- lm(bfit.resid ~ cub.pred)
  coef(fit.resid.cube)
  abline(fit.resid.cube, lty=2, col="blue")
  plot(data.use$year, bfit.resid, pch=19, cex=0.7)
  hist(bfit.resid, prob=T)
  bfit.resid.fun <- density(bfit.resid, bw=bw.SJ(bfit.resid), adjust=4)
  lines(bfit.resid.fun, col="blue", lwd=2)
  bfit.resid.pdf <- bfit.resid.fun$y/max(bfit.resid.fun$y)
  presid.trend <- coef(fit.resid.cube)[1] + coef(fit.resid.cube)[2]*bfit.resid.fun$x
  presid.range <- range(presid.trend)
  presid.range
  bfit.resid.range <- range(bfit.resid)
  max.shift <- (presid.range[2] - presid.range[1]) / (bfit.resid.range[2] - bfit.resid.range[1])
  shift.vec <- 1 - seq(0,max.shift,by=max.shift/length(bfit.resid.pdf))
  bfit.resid.out <- data.frame(bfit.resid.fun$x, bfit.resid.pdf, shift.vec[1:length(bfit.resid.pdf)])
  colnames(bfit.resid.out) <- c("presid", "prob", "shift")
  bfit.resid.out$prob.shift <- bfit.resid.out$prob * bfit.resid.out$shift
  head(bfit.resid.out)
  
  # separate by positive and negative residual groups
  # positive residuals
  bfit.resid.pos <- bfit.resid[which(bfit.resid > 0)]
  plot(cub.pred[which(bfit.resid > 0)], bfit.resid.pos, pch=19, cex=0.7)
  fit.resid.cube.pos <- lm(bfit.resid.pos ~ cub.pred[which(bfit.resid > 0)])
  coef(fit.resid.cube.pos)
  abline(fit.resid.cube.pos, lty=2, col="blue")
  plot(data.use$year[which(bfit.resid > 0)], bfit.resid.pos, pch=19, cex=0.7)
  hist(bfit.resid.pos, prob=T)
  bfit.resid.fun.pos <- density(bfit.resid.pos, bw=bw.SJ(bfit.resid.pos), adjust=4)
  lines(bfit.resid.fun.pos, col="blue", lwd=2)
  bfit.resid.pdf.pos <- bfit.resid.fun.pos$y/max(bfit.resid.fun.pos$y)
  presid.trend.pos <- coef(fit.resid.cube)[1] + coef(fit.resid.cube)[2]*bfit.resid.fun.pos$x
  presid.range.pos <- range(presid.trend.pos)
  presid.range.pos
  bfit.resid.range.pos <- range(bfit.resid.pos)
  max.shift.pos <- (presid.range.pos[2] - presid.range.pos[1]) / (bfit.resid.range.pos[2] - bfit.resid.range.pos[1])
  shift.pos.vec <- 1 - seq(0,max.shift.pos,by=max.shift.pos/length(bfit.resid.pdf.pos))
  bfit.resid.out.pos <- data.frame(bfit.resid.fun.pos$x, bfit.resid.pdf.pos, shift.pos.vec[1:length(bfit.resid.pdf.pos)])
  colnames(bfit.resid.out.pos) <- c("presid", "prob", "shift")
  bfit.resid.out.pos$prob.shift <- bfit.resid.out.pos$prob * bfit.resid.out.pos$shift
  head(bfit.resid.out.pos)
  
  # negative residuals
  bfit.resid.neg <- bfit.resid[which(bfit.resid < 0)]
  plot(cub.pred[which(bfit.resid < 0)], bfit.resid.neg, pch=19, cex=0.7)
  fit.resid.cube.neg <- lm(bfit.resid.neg ~ cub.pred[which(bfit.resid < 0)])
  coef(fit.resid.cube.neg)
  abline(fit.resid.cube.neg, lty=2, col="blue")
  plot(data.use$year[which(bfit.resid < 0)], bfit.resid.neg, pch=19, cex=0.7)
  hist(bfit.resid.neg, prob=T)
  bfit.resid.fun.neg <- density(bfit.resid.neg, bw=bw.SJ(bfit.resid.neg), adjust=4)
  lines(bfit.resid.fun.neg, col="blue", lwd=2)
  bfit.resid.pdf.neg <- bfit.resid.fun.neg$y/max(bfit.resid.fun.neg$y)
  presid.trend.neg <- coef(fit.resid.cube)[1] + coef(fit.resid.cube)[2]*bfit.resid.fun.neg$x
  presid.range.neg <- range(presid.trend.neg)
  presid.range.neg
  bfit.resid.range.neg <- range(bfit.resid.neg)
  max.shift.neg <- (presid.range.neg[2] - presid.range.neg[1]) / (bfit.resid.range.neg[2] - bfit.resid.range.neg[1])
  shift.neg.vec <- 1 - seq(0,max.shift.neg,by=max.shift.neg/length(bfit.resid.pdf.neg))
  bfit.resid.out.neg <- data.frame(bfit.resid.fun.neg$x, bfit.resid.pdf.neg, shift.neg.vec[1:length(bfit.resid.pdf.neg)])
  colnames(bfit.resid.out.neg) <- c("presid", "prob", "shift")
  bfit.resid.out.neg$prob.shift <- bfit.resid.out.neg$prob * bfit.resid.out.neg$shift
  head(bfit.resid.out.neg)
  
  ## iterate to calculate predicted future confidence intervals
  iter <- 10000
  out.mat <- matrix(data=0, nrow = iter, ncol = length(sin.pred))
  for (i in 1:iter) {
    resid.iter.neg <- sample(bfit.resid.out.neg$presid, size=length(sin.pred), prob=bfit.resid.out.neg$prob.shift, replace=T)
    resid.iter.pos <- sample(bfit.resid.out.pos$presid, size=length(sin.pred), prob=bfit.resid.out.pos$prob.shift, replace=T)
    resid.iter.neg.sub <- which(rbinom(length(resid.iter.neg),1,0.5) == 0)
    resid.iter.pos.sub <- seq(1:length(resid.iter.neg))[-resid.iter.neg.sub]
    resid.iter.neg.samp <- resid.iter.neg[resid.iter.neg.sub]
    resid.iter.pos.samp <- resid.iter.pos[resid.iter.pos.sub]
    subs.vec <- c(resid.iter.neg.sub, resid.iter.pos.sub)
    resids.vec <- c(resid.iter.neg.samp, resid.iter.pos.samp)
    resid.dat <- data.frame(subs.vec, resids.vec)
    colnames(resid.dat) <- c("sub","resid")
    resid.sort <- resid.dat[order(resid.dat[,1],decreasing=F),1:2]
    biteppsc.iter <- sin.pred + resid.sort$resid
    out.mat[i, ] <- ifelse(biteppsc.iter < 0, 0, biteppsc.iter)
  }
  bitppsc.pred.lo <- apply(out.mat, MARGIN=2, quantile, probs=0.025, na.rm=T)
  bitppsc.pred.up <- apply(out.mat, MARGIN=2, quantile, probs=0.975, na.rm=T)
  
  par(mar = c(5, 5, 3, 4))
  plot(data.use$year, data.use$bites.cap.sc, type="c", lwd=2, pch=19, ylim=c(-0.1,1.05*max(data.use$bites.cap.sc)), col="black", cex=0.6, xlim=c(1900,2066), xlab="", ylab="scaled bites/capita", main=colnames(data.use)[3], sub=paste(top.mod,": ", "wAICc=",summary.table[1,5],"; %DE=",summary.table[1,9],sep=""))
  sin.pred <- as.numeric(coef(fit.sin)[1]) + (as.numeric(coef(fit.sin)[2])*cos(as.numeric(coef(fit.sin)[3])*yr.proj + as.numeric(coef(fit.sin)[4])))
  lines(yr.proj, sin.pred, lty=2, lwd=1, col="red")
  lines(data.use$year, cub.pred, lty=2, lwd=2, col="purple")
  lines(yr.proj, bitppsc.pred.lo, lwd=2, lty=2, col="green")
  lines(yr.proj, bitppsc.pred.up, lwd=2, lty=2, col="green")
  lines(yr.proj, out.mat[runif(1,min=1, max=iter),], lwd=1, lty=3, col="green") # print one random selection of predicted bites
  
  ## how many people bitten? (median)
  ## return to normal scale
  bites.cap.usc <- as.numeric(data.use$bites.cap.sc) * attr(data.use$bites.cap.sc, 'scaled:scale')
  sin.pred.usc <- sin.pred * attr(data.use$bites.cap.sc, 'scaled:scale')
  region.use <- colnames(data.use)[3]
  proj.pop <- pop[which(pop$year > 2019),which(colnames(pop) == region.use)]
  sin.pred.pop <- sin.pred.usc * proj.pop
  plot(data.use$year, data.use$bites, xlim=c(1900,2066), type="l",ylim=c(min(data.use$bites), max(sin.pred.pop)), xlab=NA, ylab="bites/year")
  lines(yr.proj, sin.pred.pop, lty=2, lwd=2, col="green")
  
  ## reduction?
  ## assume 60% reduction in bites if wearing shark shield (Huveneers et al. 2018)
  ## assume different proportions wearing
  bite.det.red <- 0.60
  bite.det.red.se <- 0.05 # assumed
  bred.alpha <- estBetaParams(bite.det.red, bite.det.red.se^2)$alpha
  bred.beta <- estBetaParams(bite.det.red, bite.det.red.se^2)$beta
  prop.wearing <- seq(0.1,1,0.1)
  
  # stochastically resample bite rate reduction
  bites.fewer.arr <- bites.red.arr <- bites.arr <- array(data=NA, c(length(prop.wearing), length(sin.pred.pop), iter))
  for (i in 1:iter) {
    #bites.red.mat <- bites.mat <- matrix(data=0, nrow=length(prop.wearing), ncol = length(sin.pred.pop))
    for (j in 1:length(prop.wearing)) {
      bites.red.arr[j,,i]  <- sin.pred.pop * (1 - (rbeta(1,bred.alpha,bred.beta) * prop.wearing[j]))
      bites.arr[j,,i] <- sin.pred.pop
    }
  }
  bites.fewer.arr <- bites.arr - bites.red.arr
  bites.mn <- apply(bites.arr, c(1,2), mean, na.rm=T)
  bites.red.mn <- apply(bites.red.arr, c(1,2), mean, na.rm=T)
  bites.fewer.mn <- apply(bites.fewer.arr, c(1,2), mean, na.rm=T)
  
  bites.lo <- apply(bites.arr, c(1,2), quantile, probs=0.025, na.rm=T)
  bites.red.lo <- apply(bites.red.arr, c(1,2), quantile, probs=0.025, na.rm=T)
  bites.fewer.lo <- apply(bites.fewer.arr, c(1,2), quantile, probs=0.025, na.rm=T)
  
  bites.up <- apply(bites.arr, c(1,2), quantile, probs=0.975, na.rm=T)
  bites.red.up <- apply(bites.red.arr, c(1,2), quantile, probs=0.975, na.rm=T)
  bites.fewer.up <- apply(bites.fewer.arr, c(1,2), quantile, probs=0.975, na.rm=T)
  
  # contour plot2
  # all bites
  cpyconta <- plot_ly(z = ~bites.mn, autocontour=T, type="contour", line = list(smoothing = 0.85), contours = list(showlabels = TRUE, labelfont=list(
    size=14, color="white"))) %>%
    colorbar(title = "all bites/year") %>%
    layout(
      xaxis = list(title="year"),
      yaxis = list(title="proportion wearing", ticketmode='array', ticktext=as.character(seq(0.1,1,0.1)), tickvals=c(0,seq(1,10,1))))
  cpyconta
  cpycont1 <- plot_ly(z = ~bites.red.mn, autocontour=T, type="contour", line = list(smoothing = 0.85), contours = list(showlabels = TRUE, labelfont=list(
    size=14, color="white"))) %>%
    colorbar(title = "reduced bites/year") %>%
    layout(
      xaxis = list(title="year"),
      yaxis = list(title="proportion wearing", ticketmode='array', ticktext=as.character(seq(0.1,1,0.1)), tickvals=c(0,seq(1,10,1))))
  cpycont1
  cpycont2 <- plot_ly(z = ~bites.fewer.mn, autocontour=T, type="contour", line = list(smoothing = 0.85), contours = list(showlabels = TRUE, labelfont=list(
    size=14, color="white"))) %>%
    colorbar(title = "fewer bites/year") %>%
    layout(
      xaxis = list(title="year"),
      yaxis = list(title="proportion wearing", ticketmode='array', ticktext=as.character(seq(0.1,1,0.1)), tickvals=c(0,seq(1,10,1))))
  cpycont2

  # 3D surface plot with upper, lower, mean
  # fewer
  fb3d <- plot_ly(showscale = FALSE) %>% 
    add_surface(z = ~bites.fewer.mn) %>%
    add_surface(z = ~bites.fewer.lo, opacity = 0.55) %>%
    add_surface(z = ~bites.fewer.up, opacity = 0.55) %>%
    layout(scene = list(
      xaxis = list(title="year"),
      yaxis = list(title="proportion wearing", ticketmode='array', ticktext=as.character(seq(0.1,1,0.1)), tickvals=c(0,seq(1,10,1))),
      zaxis = list(title="fewer bites/year")))
  fb3d
  
  
  ##########################################
  # stochastic resampling for bites avoided
  # accounting for full uncertainty in projected bites/person
  
  bites.mat.usc <- out.mat * attr(data.use$bites.cap.sc, 'scaled:scale')
  iter2 <- dim(bites.mat.usc)[1]
  region.use <- colnames(data.use)[3]
  proj.pop <- pop[which(pop$year > 2019),which(colnames(pop) == region.use)]
  
  bites.mat.usc.pop <- matrix(data=NA, nrow=iter2, ncol=length(proj.pop))
  for (b in 1:iter2) {
    bites.mat.usc.pop[b,] <- bites.mat.usc[b,] * proj.pop
  }
  
  plot(data.use$year, data.use$bites, xlim=c(1900,2066), type="l",ylim=c(min(data.use$bites), max(sin.pred.pop)), xlab=NA, ylab="bites/year")
  lines(yr.proj, sin.pred.pop, lty=2, lwd=2, col="green")
  
  ## reduction?
  ## assume 60% reduction in bites if wearing shark shield (Huveneers et al. 2018)
  ## assume different proportions wearing
  bite.det.red <- 0.60
  bite.det.red.se <- 0.05 # assumed
  bred.alpha <- estBetaParams(bite.det.red, bite.det.red.se^2)$alpha
  bred.beta <- estBetaParams(bite.det.red, bite.det.red.se^2)$beta
  prop.wearing <- seq(0.1,1,0.1)
  
  # stochastically resample bite rate reduction
  bites.fewer.arr <- bites.red.arr <- bites.arr <- array(data=NA, c(length(prop.wearing), length(sin.pred.pop), iter2))
  for (i in 1:iter2) {
    for (j in 1:length(prop.wearing)) {
      bites.red.arr[j,,i]  <- bites.mat.usc.pop[i,] * (1 - (rbeta(1,bred.alpha,bred.beta) * prop.wearing[j]))
      bites.arr[j,,i] <- bites.mat.usc.pop[i,]
    }
  }
  bites.fewer.arr <- bites.arr - bites.red.arr
  bites.mn <- apply(bites.arr, c(1,2), mean, na.rm=T)
  bites.red.mn <- apply(bites.red.arr, c(1,2), mean, na.rm=T)
  bites.fewer.mn <- apply(bites.fewer.arr, c(1,2), mean, na.rm=T)
    
  bites.lo <- apply(bites.arr, c(1,2), quantile, probs=0.025, na.rm=T)
  bites.red.lo <- apply(bites.red.arr, c(1,2), quantile, probs=0.025, na.rm=T)
  bites.fewer.lo <- apply(bites.fewer.arr, c(1,2), quantile, probs=0.025, na.rm=T)
    
  bites.up <- apply(bites.arr, c(1,2), quantile, probs=0.975, na.rm=T)
  bites.red.up <- apply(bites.red.arr, c(1,2), quantile, probs=0.975, na.rm=T)
  bites.fewer.up <- apply(bites.fewer.arr, c(1,2), quantile, probs=0.975, na.rm=T)
  
  f1 <- list(
    family = "Avenir Light",
    size = 26,
    color = "black"
  )
  f2 <- list(
    family = "Avenir Light",
    size = 18,
    color = "black"
  )
  cpycont2 <- plot_ly(z = ~bites.fewer.mn, autocontour=F, type="contour", line = list(smoothing = 0.90), contours = list(start=0, end=max(bites.fewer.mn), size=2, showlabels = TRUE, labelfont=list(
    size=18, family="Avenir Light", face="bold", color="white"))) %>%
    colorbar(title = "fewer bites/year", titlefont=f2, tickfont=f2) %>%
    layout(
      xaxis = list(title="", titlefont=f1, tickfont=f2, ticketmode='array', ticktext=as.character(seq(2020,2066,10)), tickvals=c(0,seq(10,40,10))),
      yaxis = list(title="proportion wearing", titlefont=f1, tickfont=f2, ticketmode='array', ticktext=as.character(seq(0.1,1,0.1)), tickvals=c(0,seq(1,10,1))))
  cpycont2
  
  plot(prop.wearing, rowSums(bites.fewer.mn), pch=19, type="l", xlab="proportion wearing", ylab="total fewer incidents 2020-2066", ylim=c(min(rowSums(bites.fewer.lo)),max(rowSums(bites.fewer.up))))
  lines(prop.wearing, rowSums(bites.fewer.lo), lty=2)
  lines(prop.wearing, rowSums(bites.fewer.up), lty=2)
  
  round(rowSums(bites.fewer.lo), 0)
  round(rowSums(bites.fewer.mn), 0)
  round(rowSums(bites.fewer.up), 0)


# 3D surface plot with upper, lower, mean
# fewer
f3 <- list(
  family = "Avenir Light",
  size = 16,
  color = "black"
)
fb3d <- plot_ly(showscale = FALSE) %>% 
  add_surface(z = ~bites.fewer.mn) %>%
  add_surface(z = ~bites.fewer.lo, opacity = 0.55) %>%
  add_surface(z = ~bites.fewer.up, opacity = 0.55) %>%
  layout(scene = list(
    xaxis = list(title="", titlefont=f1, tickfont=f3, ticketmode='array', ticktext=as.character(seq(2020,2066,10)), tickvals=c(0,seq(10,40,10))),
    yaxis = list(title="proportion wearing", titlefont=f1, tickfont=f3, ticketmode='array', ticktext=as.character(seq(0.1,1,0.1)), tickvals=c(0,seq(1,10,1))),
    zaxis = list(title="fewer bites/year", tickfont=f3, titlefont=f1)))
fb3d
