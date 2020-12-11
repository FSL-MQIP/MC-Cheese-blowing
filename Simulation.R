## Setting working directory
setwd("~/Cornell/Projects/MC Cheese blowing/MC-Cheese-blowing/")

## Simulate distribution of initial spore conc-----------------------------------------
#Load file
Spore = read.csv("Initial spore concentration.csv", header = T)

#Load pacakge
library(tidyr)
library(fitdistrplus)

##Data cleaning
colnames(Spore)[1] = "ID"
Spore = gather(Spore,"Month","Count",-"ID")

##Create data frame
Spore$Count[Spore$Count == "" ] = NA
count = na.omit(Spore$Count)
count[count == "<18"] = "4.5"

# box-cox
library(MASS)
boxcox(count ~ 1)

## Check the correct distribution to fit
#count = as.numeric(count)
#fit = fitdist(count_t,"norm", "mme")
#summary(fit)
#fit2 = fitdist(count, "pois", method = "mme")
#summary(fit2)
#fit3 = fitdist(count, "nbinom", method = "mme")
#summary(fit3)

# Visualize simulated distribution
count_t = as.numeric(count)^0.1 
#qqnorm(count_t);qqline(count_t)
#par(mfrow=c(1,1))
hist(count_t,probability = T, xlim = c(1,2), breaks = 20, 
     main = "Histogram of spore count distribution",
     xlab = "Spore count^0.1 MPN/L",
     ylim = c(0,4))
tempx = seq(1,2,length.out = 100)
tempy = dnorm(tempx, mean(count_t), sd(count_t)); lines(tempx, tempy,col='red',lwd=2)


## Simulation set-up -----------------------------------------------------------
## Cheese vat and block identifications
n_sim = 10000
vat_id = rep(1:n_sim, each = 20)
block_id = rep(1:20, n_sim)

## Simulate var spore count as MPN/kg milk in cheese vat

vat_norm_count = rnorm(n_sim, mean(count_t), sd(count_t))
vat_sim_count = vat_norm_count^10
#hist(vat_sim_count, breaks = 50)

## Combine ids with initial spore count
vat_count = round(rep(vat_sim_count, each = 20))
block_count = vector()
for (i in 1:length(vat_count)){
  block_count[i] = rpois(1, lambda = vat_count[i])
}
data = data.frame(vat_id, vat_count,block_id,block_count)
#head(data, 20)
#mean(data$block_sim_count)

## Ripening parameters
temp = 14
pH = runif(n_sim, min = 5.2, max = 5.6)
data$pH = round(rep(pH, each = 20),2)


#Functions-------------------------------------------------------------------------
baranyi_log10N = function(t,lag,mumax,LOG10N0,LOG10Nmax) {
  ans <- LOG10Nmax + log10((-1 + exp(mumax * lag) + exp(mumax * t))/(exp(mumax * t) - 1 + exp(mumax * lag) * 10^(LOG10Nmax -LOG10N0)))
  return(ans)
}

final_conc = function(N0, temp = 37, pH = 5.8, aw = 1, aw_min, hour){
  #Default growth rate 
  mumax = 5.415
  
  #Gamma concenpt
  gamma_T = ((temp -10)/(37-10))^2
  gamma_pH = (pH-4.65)*(7.3-pH)/(5.8-4.65)/(7.3-5.8)
  gamma_aw = (aw-aw_min)/(1-aw_min)
  
  #Adjusted values
  mu = mumax*gamma_T*gamma_pH*gamma_aw
  mu = if (mu<0) 0 else mu
  print(mu)
  
  #Model growth
  N = baranyi_log10N(hour, lag=0, mumax=mu, LOG10N0 = log10(N0), LOG10Nmax = 12.2)
  10^N
}


aw_inhibit = function(pH){
  salt_max_h = -4.7953*pH^2+55.187*pH-154.18
  salt_max_l = -3.6257*pH^2+41.971*pH-118.42
  aw_min_l = 0.995-0.00721*salt_max_h
  aw_min_h = 0.995-0.00721*salt_max_l
  aw_min = vector()
  for (i in 1:length(pH)){
    aw_min[i] = runif(1, aw_min_l[i],aw_min_h[i])
  }
  
  round(aw_min,3)
}


day2nogrowth=function(aw){
}


aw_atday = function(days){
  0.000004*days^2-0.0009*days+0.9896
}

#------------------------------------------
## add time limit based on pH
data$aw_inhibit = round(aw_inhibit(data$pH),3)

final_count = final_conc(data$block_count, temp = 14, pH=data$pH, aw=aw_atday(1),aw_min=data$aw_inhibit, hour=1*24)
final_count = cbind(final_count, 
                    final_conc(final_count, temp = 14, pH=data$pH, aw=aw_atday(2),aw_min=data$aw_inhibit, hour=1*24))
for (i in 3:37){
  N0 = final_count[,i-1]
  new = final_conc(N0, temp = 14, pH=data$pH, aw=aw_atday(i),aw_min=data$aw_inhibit,hour=1*24)
  final_count = cbind(final_count, new)
}

var_name = vector()
for (i in 1:37){
  var_name[i] = paste("day",as.character(i),sep = "_")
}

colnames(final_count) = var_name
data = cbind(data, final_count)

hist(log10(data$day_25), freq = F)

sum(log10(data$day_25)>9.185)/(20*n_sim)
