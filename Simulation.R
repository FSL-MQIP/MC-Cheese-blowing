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
#boxcox(count ~ 1)

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


## Simulation set-up -----------------------------------------------------------
## Cheese vat and block identifications
n_sim = 10000
vat_id = rep(1:n_sim, each = 20)
block_id = rep(1:20, n_sim)

## Simulate var spore count as MPN/kg milk in cheese vat
set.seed(1)
vat_norm_count = rnorm(n_sim, mean(count_t), sd(count_t))
vat_sim_count = vat_norm_count^10
#hist(vat_sim_count, breaks = 50)

## Simulate block spore count 
mean_sd_ratio = mean(count_t)/sd(count_t)# simulate sd for block spore count
block_sim_mean = vat_sim_count
block_sim_sd = block_sim_mean/mean_sd_ratio

block_sim_count = vector()
set.seed(1)
for (i in 1:n_sim){
  block_sim_count = c(block_sim_count, 
                      rnorm(20, block_sim_mean[i], block_sim_sd[i]))
}


## Combine ids with initial spore count
vat_count = round(rep(vat_sim_count, each = 20),1)
block_count = round(block_sim_count,1)
data = data.frame(vat_id, vat_count,block_id, block_count)
#head(data, 20)
#mean(data$block_sim_count)

## Ripening parameters
temp = 14
pH = runif(20*n_sim, min = 5.2, max = 5.6)
data$pH = round(pH,2)

## Function to calculate final concentration 
final_conc = function(N0, temp = 37, pH = 5.8, nitrite = 0, hour){
  #Default growth rate 
  mumax = 0.12
  
  #Randomize mu*lambda
  prod = runif(1, min = 1.02, max = 4.8)
  
  #Calculate the lag phase
  lag = prod/mumax
  
  #Gamma concenpt
  gamma_T = ((temp -10)/(37-10))^2
  gamma_pH = (pH-4.6)*(7.5-pH)/(5.8-4.6)/(7.5-5.8)
  gamma_nitrite = 1-nitrite/75
  
  #Add interaction term
  theta_T = ((37-temp)/(37-10))^3
  theta_pH = ((5.8-pH)/(5.8-4.6))^3
  epi = theta_T/(2*(1-theta_pH))+theta_pH/(2*(1-theta_T))
  if(epi<=0.5){
    theta = 1
  }  else if (epi>0.5&epi<1){
    theta=2*(1-epi)
  }  else{
    theta=0
  }
  
  #Adjusted values
  mu = mumax*gamma_T*gamma_pH*gamma_nitrite*theta
  lag = 1/(1/lag*gamma_T*gamma_pH*gamma_nitrite*theta)
  
  #Model growth
  time = hour - lag
  if (time < 0){
    time = 0
  }
  Nmax = N0*(1+mu)^time
  Nmax
}


## Simulate final conc at day 70
data$day70 = NA
for (i in 1:length(data$day70)){
data$day70[i] = final_conc(data$block_count[i], 
                           temp = temp, 
                           pH = pH[i], 
                           hour = 70*24)
}
data$day70 = round(data$day70,1)
#head(data, 40)

## Proportions of spoilage at day 70
sum(data$day70>(10^4.23))/length(data$day70) # Threshold at log average


## Simulate final conc at day 90
data$day90 = NA
for (i in 1:length(data$day90)){
  data$day90[i] = final_conc(data$block_count[i], 
                             temp = temp, 
                             pH = pH[i], 
                             hour = 90*24)
}
data$day90 = round(data$day90,1)
#head(data, 40)

## Proportions of spoilage at day 90
sum(data$day90>(10^4.23))/length(data$day90) # Threshold at log average


