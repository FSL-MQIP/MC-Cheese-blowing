## Setting working directory
setwd("~/Cornell/Projects/MC Cheese blowing/MC-Cheese-blowing/")

## Simulate distribution of initial spore conc-----------------------------------------
#Load file
Spore = read.csv("Initial spore concentration.csv", header = T)

#Load pacakge
library(tidyr)
library(dplyr)
library(fitdistrplus)
library(EnvStats)
library(ggplot2)

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
#hist(count_t,probability = T, xlim = c(1,2), breaks = 20, 
#main = "Histogram of spore count distribution",
#xlab = "Spore count^0.1 MPN/L",
#ylim = c(0,4))
#tempx = seq(1,2,length.out = 100)
#tempy = dnorm(tempx, mean(count_t), sd(count_t)); lines(tempx, tempy,col='red',lwd=2)

## Functions ----------------------------------------------------------------------
baranyi_log10N = function(t,lag,mumax,LOG10N0,LOG10Nmax){
  ans = vector()
  for (i in 1:length(LOG10N0)){
    ans[i] <- LOG10Nmax + log10((-1 + exp(mumax[i] * lag) + exp(mumax[i] * t))/(exp(mumax[i] * t) - 1 + exp(mumax[i] * lag) * 10^(LOG10Nmax -LOG10N0[i])))
  }
  return(ans)
}

final_conc = function(N0, temp = 37, pH = 5.8, aw = 1, aw_min, hour){
  #Default growth rate 
  mumax = 1.44
  
  #Gamma concenpt
  gamma_T = ((temp -10)/(37-10))^2
  gamma_pH = (pH-4.65)*(7.3-pH)/(5.8-4.65)/(7.3-5.8)
  gamma_aw = (aw-aw_min)/(1-aw_min)
  
  #Adjusted values
  mu = mumax*gamma_T*gamma_pH*gamma_aw
  for (i in 1:length(pH)){
    mu[i] = if (mu[i]<0) 0 else mu[i]
  }
  print(mu)
  
  #Model growth
  N = baranyi_log10N(hour, lag=0, mumax=mu, LOG10N0 = log10(N0), LOG10Nmax = 12.2)
  10^N
}

mumax = function(temp, pH, aw=1){
  #Default growth rate 
  mumax = 1.44
  
  #Gamma concenpt
  gamma_T = ((temp -10)/(37-10))^2
  gamma_pH = (pH-4.65)*(7.3-pH)/(5.8-4.65)/(7.3-5.8)
  gamma_aw = (aw-0.948)/(1-0.948)
  
  #Adjusted values
  mu = mumax*gamma_T*gamma_pH*gamma_aw
  for (i in 1:length(pH)){
    mu[i] = if (mu[i]<0) 0 else mu[i]
  }
  mu
}

dbtime = function(mu){
  log(2, 1+mu)
}

salt2aw = function(salt){
  aw = 0.995-0.00721*salt
  aw
}


## Simulation set-up -----------------------------------------------------------
## Cheese vat and block identifications
set.seed(1)
n_sim = 10000
vat_id = rep(1:n_sim, each = 20)
block_id = rep(1:20, n_sim)

## Simulate var spore count as MPN/kg milk in cheese vat

vat_norm_count = rnorm(n_sim, mean(count_t), sd(count_t))
vat_sim_count = vat_norm_count^10
#hist(vat_sim_count, breaks = 50)

## Combine ids with initial spore count
vat_count = round(rep(vat_sim_count, each = 20))

## Worst case N0
#vat_count = round(rep(vat_sim_count*10, each = 20))
## Best case N0
#vat_count = round(rep(vat_sim_count/10, each = 20))


block_count = vector()
for (i in 1:length(vat_count)){
  block_count[i] = rpois(1, lambda = vat_count[i])
}

# Consider concentration effect during whey draining
conc_factor = 6
block_count = conc_factor*block_count


# Combining data
data = data.frame(vat_id, vat_count,block_id,block_count)
#head(data, 20)
#mean(data$block_sim_count)

## Ripening parameters
temp = 14
# Worst temp
#temp = 16
# Best temp
#temp = 12

#pH = runif(n_sim, min = 5.2, max = 5.3)
pH= rtri(n_sim, min = 5.1, max = 5.63, mode = 5.4)
data$pH = round(rep(pH, each = 20),2)
data$mumax = round(mumax(temp,data$pH,aw=0.953),5)

#------------------------------------------
# Threshold level calculation
low_mumax = mumax(13,5.04,0.962)
high_mumax = mumax(13,5.01,0.962)
low_bound = baranyi_log10N(t=70*24, lag=0, mumax=low_mumax,  log10(250), log10(40000000))
high_bound = baranyi_log10N(t=60*24, lag=0, mumax=high_mumax,  log10(2500), log10(40000000))


# Model simulation
final_count = vector()
for (i in c(50,60,70,80,90,120)){
  new = baranyi_log10N(t=i*24, lag=0, mumax=data$mumax,log10(data$block_count), log10(40000000))
  final_count = cbind(final_count, new) %>% round(3)
}

name = c("day50","day60","day70","day80","day90","day120")
colnames(final_count) = name
final_count = as.data.frame(final_count)


## Check results
result = vector()
for (i in 1:ncol(final_count)){
  high = sum(final_count[i]>low_bound)/(20*n_sim)
  low = sum(final_count[i]>high_bound)/(20*n_sim)
  prob = c(high, low)
  result = cbind(result, prob)
}
colnames(result) = name
rownames(result) = c("High prob", "Low prob")
result


hist(final_count$day60, freq = F)
data = data %>% cbind(final_count)


##-----------------------------------------------------------------------------------------------------------
# Sensitivity analysis

# Threshold level
# 



















## Error analysis---------------------------------------------------------------------------------------------------

# lower likelihood
# late blown cheese
blow_low = vector()
for (i in 1:6){
  blow = data %>% filter(final_count[,i]>=low_bound)
  mean_pH = mean(blow$pH)
  sd_pH = sd(blow$pH)
  mean_block_count=mean(log10(blow$block_count))
  sd_block_count = sd(log10(blow$block_count))
  sum = cbind(mean_pH, sd_pH, mean_block_count, sd_block_count)
  blow_low = rbind(blow_low, sum)
  
}
rownames(blow_low) = name
blow_low

# not late blown cheese
noblow_low = vector()
for (i in 1:6){
  noblow = data %>% filter(final_count[,i]<low_bound)
  mean_pH = mean(blow$pH)
  sd_pH = sd(blow$pH)
  mean_block_count=mean(log10(blow$block_count))
  sd_block_count = sd(log10(blow$block_count))
  sum = cbind(mean_pH, sd_pH, mean_block_count, sd_block_count)
  noblow_low = rbind(noblow_low, sum)
  
}
rownames(noblow_low) = name
noblow_low


# Higher likelihood
blow_high = vector()
for (i in 1:6){
  blow = data %>% filter(final_count[,i]>high_bound)
  mean_pH = mean(blow$pH)
  sd_pH = sd(blow$pH)
  mean_block_count=mean(log10(blow$block_count))
  sd_block_count = sd(log10(blow$block_count))
  sum = cbind(mean_pH, sd_pH, mean_block_count, sd_block_count)
  blow_high = rbind(blow_high, sum)
  
}
rownames(blow_high) = name
blow_high


# Boundary line
bound_low = data %>% filter(final_count$day60> (low_bound-0.01),  final_count$day60< (low_bound+0.01))
plot(x=bound_low$pH, y=log10(bound_low$block_count), xlab="pH", ylab="log count (MPN/kg)",
     main = "Boundary line at day 50 for low threshold",
     xlim = c(5.1,5.63))

bound_high = data %>% filter(final_count$day60> (high_bound-0.01),  final_count$day60< (high_bound+0.01))
plot(x=bound_high$pH, y=log10(bound_high$block_count), xlab="pH", ylab="log count (MPN/kg)",
     main = "Boundary line at day 50 for high threshold",
     xlim = c(5.1,5.63))

