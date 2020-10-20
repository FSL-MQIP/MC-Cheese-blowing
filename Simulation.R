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

# Visualize simulated distribution
count_t = as.numeric(count)^0.1 
par(mfrow=c(1,1))
hist(count_t,probability = T, xlim = c(1,2), breaks = 20, 
     main = "Histogram of spore count distribution",
     xlab = "Spore count^0.1 MPN/L",
     ylim = c(0,4))
tempx = seq(1,2,length.out = 100)
tempy = dnorm(tempx, mean(count_t), sd(count_t)); lines(tempx, tempy,col='red',lwd=2)


## Simulation set-up -----------------------------------------------------------
## Cheese vat and block identifications
n_sim = 1000
vat_id = rep(1:n_sim, each = 20)
block_id = rep(1:20, n_sim)

## Simulate var spore count as MPN/kg milk in cheese vat
set.seed(87)
vat_norm_count = rnorm(n_sim, mean(count_t), sd(count_t))
vat_sim_count = vat_norm_count^10
#hist(vat_sim_count, breaks = 50)

## Simulate block spore count 
mean_sd_ratio = mean(count_t)/sd(count_t)# simulate sd for block spore count
block_sim_mean = vat_sim_count
block_sim_sd = block_sim_mean/mean_sd_ratio

block_sim_count = vector()
set.seed(87)
for (i in 1:n_sim){
  block_sim_count = c(block_sim_count, 
                      rnorm(20, block_sim_mean[i], block_sim_sd[i]))
}

## Combine ids with initial spore count
data = data.frame(vat_id, block_id, block_sim_count)
#head(data, 20)
#mean(data$block_sim_count)

## Ripening parameters
temp = 14
pH = runif(20*n_sim, min = 5.2, max = 5.6)
data$pH = round(pH,2)

## Function to calculate final concentration 
final_conc = function(N0, temp = 37, pH = 5.6, nitrite = 0, hour){
  #Default growth rate and lag time
  mumax = 0.12; lag = 1/mumax
  
  #Gamma concenpt
  gamma_T = ((temp -9)/(37-9))^2
  gamma_pH = (pH-4.4)*(6.8-pH)/(5.6-4.4)/(6.8-5.6)
  gamma_nitrite = 1-nitrite/75
  
  #Adjusted values
  mu = mumax*gamma_T*gamma_pH*gamma_nitrite
  lag = 1/(1/lag*gamma_T*gamma_pH*gamma_nitrite)
  
  #Model growth
  time = hour - lag
  Nmax = N0*(1+mu)^time
  Nmax
}


## Simulate final conc at day 70
data$day70 = NA
for (i in 1:length(data$day70)){
data$day70[i] = final_conc(data$block_sim_count[i], 
                           temp = temp, 
                           pH = pH[i], 
                           hour = 70*24)
}
#head(data, 40)

## Proportions of spoilage at day 70
sum(data$day70>140000)/length(data$day70) # Threshold at 1.4e5
sum(data$day70>1500000)/length(data$day70) # Threshold at 1.5e6
sum(data$day70>477773)/length(data$day70) # Threshold at log average

## Simulate final conc at day 90
data$day90 = NA
for (i in 1:length(data$day90)){
  data$day90[i] = final_conc(data$block_sim_count[i], 
                             temp = temp, 
                             pH = pH[i], 
                             hour = 90*24)
}
#head(data, 40)

## Proportions of spoilage at day 90
sum(data$day90>140000)/length(data$day90) # Threshold at 1.4e5
sum(data$day90>1500000)/length(data$day90) # Threshold at 1.5e6
sum(data$day90>477773)/length(data$day90) # Threshold at log average

