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
baranyi_log10N = function(t,lag,mumax,LOG10N0,LOG10Nmax) {
  ans <- LOG10Nmax + log10((-1 + exp(mumax * lag) + exp(mumax * t))/(exp(mumax * t) - 1 + exp(mumax * lag) * 10^(LOG10Nmax -LOG10N0)))
  return(ans)
}

final_conc = function(N0, temp = 37, pH = 5.8, aw = 1, aw_min, hour){
  #Default growth rate 
  mumax = 5.37
  
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


salt2aw = function(salt){
  aw = 0.995-0.00721*salt
  aw
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

aw_worst = function(pH){
  salt_max_h = -4.7953*pH^2+55.187*pH-154.18
  aw_min_l = 0.995-0.00721*salt_max_h
  round(aw_min_l,3)
}

aw_best = function(pH){
  salt_max_l = -3.6257*pH^2+41.971*pH-118.42
  aw_min_h = 0.995-0.00721*salt_max_l
  round(aw_min_h,3)
}

day2nogrowth=function(aw){
  (0.9907-aw)/0.0008
}


aw_atday = function(days){
  4e-06*days^2-0.0009*days+0.9896
}

aw_atday_best = function(days){
  -0.0003*days+0.9893
}

pH_atday = function(pH, days){
  for (i in length(pH)){
    pH[i]=pH[i]+0.4/40*days
  }
  round(pH)
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

#------------------------------------------
## calculate inhibitory aw
#data$aw_inhibit = round(aw_inhibit(data$pH),3)

# Worst inhibitory aw
#data$aw_inhibit = round(aw_worst(data$pH),3)
# Best inhibitory aw
#data$aw_inhibit = round(aw_inhibit(data$pH),3)

salt_values = c(2.4, 3.0, 4.2)
#salt_values = c(0, 1.2, 2.4, 3.0, 3.6, 4.2)
aw_values = salt2aw(salt_values)
data$aw_inhibit = round(sample(aw_values, size = 20*n_sim, replace = T),3)

# Last day for growth
lday = 33

## Iteratively simulate aw and counts at each day
final_count = final_conc(data$block_count, temp = temp, pH=data$pH, aw=aw_atday(1),aw_min=data$aw_inhibit, hour=1*24)
final_count = cbind(final_count, 
                    final_conc(final_count, temp = temp, pH=data$pH, aw=aw_atday(2),aw_min=data$aw_inhibit, hour=1*24))
for (i in 3:lday){
  N0 = final_count[,i-1]
  new = final_conc(N0, temp = temp, pH=data$pH, aw=aw_atday(i),aw_min=data$aw_inhibit,hour=1*24)
  final_count = cbind(final_count, new)
}
final_count = round(final_count)

## Create the column names for count at each day
var_name = vector()
for (i in 1:lday){
  var_name[i] = paste("day",as.character(i),sep = "_")
}
colnames(final_count) = var_name

## Combine data
data = cbind(data, final_count)


## Check results
hist(log10(data$day_15), freq = F)
sum(log10(data$day_15)>11.61)/(20*n_sim)


mean(c(11.26,11.96))
for (i in (1:lday)){
  print(i)
  print(sum(log10(data[,i+6])>11.61)/(20*n_sim))
}


## Error analysis
blowat15 = data %>% filter(log10(day_15)>11.61)
length(blowat15[,1])/(20*n_sim)
range(blowat15$pH)
mean(blowat15$pH); sd(blowat15$pH)
mean(log10(blowat15$block_count)); sd(log10(blowat15$block_count))
mean(blowat15$aw_inhibit)

blowat18 = data %>% filter(log10(day_15)<11.61, log10(day_18)>11.61 )
length(blowat18[,1])/(20*n_sim)
range(blowat18$pH)
mean(blowat18$pH); sd(blowat18$pH)
mean(log10(blowat18$block_count)); sd(log10(blowat18$block_count))
mean(blowat18$aw_inhibit)

blowat21 = data %>% filter(log10(day_18)<11.61, log10(day_21)>11.61 )
length(blowat21[,1])/(20*n_sim)
range(blowat21$pH)
mean(blowat21$pH); sd(blowat21$pH)
mean(log10(blowat21$block_count)); sd(log10(blowat21$block_count))
mean(blowat21$aw_inhibit)

blowat24 = data %>% filter(log10(day_21)<11.61, log10(day_24)>11.61 )
length(blowat24[,1])/(20*n_sim)
range(blowat24$pH)
mean(blowat24$pH); sd(blowat24$pH)
mean(log10(blowat24$block_count)); sd(log10(blowat24$block_count))
mean(blowat24$aw_inhibit)

blowat27 = data %>% filter(log10(day_24)<11.61, log10(day_27)>11.61 )
length(blowat27[,1])/(20*n_sim)
range(blowat27$pH)
mean(blowat27$pH); sd(blowat27$pH)
mean(log10(blowat27$block_count)); sd(log10(blowat27$block_count))
mean(blowat27$aw_inhibit)

blowat30 = data %>% filter(log10(day_27)<11.61, log10(day_30)>11.61 )
length(blowat30[,1])/(20*n_sim)
range(blowat30$pH)
mean(blowat30$pH); sd(blowat30$pH)
mean(log10(blowat30$block_count)); sd(log10(blowat30$block_count))
mean(blowat30$aw_inhibit)

blowat33 = data %>% filter(log10(day_30)<11.61, log10(day_33)>11.61 )
length(blowat33[,1])/(20*n_sim)
range(blowat33$pH)
mean(blowat33$pH); sd(blowat33$pH)
mean(log10(blowat33$block_count)); sd(log10(blowat33$block_count))
mean(blowat33$aw_inhibit)

noblow = data %>% filter(log10(day_33)<11.61)
length(noblow[,1])/(20*n_sim)
range(noblow$pH)
mean(noblow$pH); sd(noblow$pH)
mean(log10(noblow$block_count)); sd(log10(noblow$block_count))
mean(noblow$aw_inhibit)