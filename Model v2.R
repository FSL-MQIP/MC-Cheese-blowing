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
library(reshape2)
library(plotly)

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
  mumax = 0.76
  
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
  mumax = 0.76
  
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

## Combine ids with initial spore count
vat_count = round(rep(vat_sim_count, each = 20))


block_count = vector()
for (i in 1:length(vat_count)){
  block_count[i] = rpois(1, lambda = vat_count[i])
}

# Consider concentration effect during whey draining
conc_factor = 6
block_count = conc_factor*block_count

# Combining data
data = data.frame(vat_id, vat_count,block_id,block_count)

## Ripening parameters
temp = 13
pH= rtri(n_sim, min = 5.1, max = 5.63, mode = 5.4)
data$pH = round(rep(pH, each = 20),2)
data$mumax = round(mumax(temp,data$pH,aw=salt2aw(4.9)),5)

#------------------------------------------
# Threshold level calculation
low_mumax = mumax(13,5.0,salt2aw(3.7))
high_mumax = mumax(13,5.0,salt2aw(3.7))
low_bound = baranyi_log10N(t=70*24, lag=0, mumax=low_mumax,  log10(1300), log10(40000000))
high_bound = baranyi_log10N(t=60*24, lag=0, mumax=high_mumax,  log10(2500), log10(40000000))

#
#data = data %>% filter(vat_count>=1000)

# Model simulation
final_count = vector()
for (i in c(30,60,70,80,90,120)){
  new = baranyi_log10N(t=i*24, lag=0, mumax=data$mumax,log10(data$block_count), log10(40000000))
  final_count = cbind(final_count, new) %>% round(3)
}

name = c("day30","day60","day70","day80","day90","day120")
colnames(final_count) = name
final_count = as.data.frame(final_count)


## Check results
result = vector()
for (i in 1:ncol(final_count)){
  high = mean(final_count[i]>low_bound)
  low = mean(final_count[i]>high_bound)
  prob = c(high, low)
  result = cbind(result, prob)
}
colnames(result) = name
rownames(result) = c("High_prob", "Low_prob")
result

## Append the final count to the data
data = data %>% cbind(final_count)


## Error analysis
data[final_count$day30 > high_bound, ]$vat_count %>% log10() %>% mean()  #average conc. for cheese LBD
data[final_count$day30 > low_bound, ]$vat_count %>% log10() %>% mean()
data[final_count$day30 > high_bound, ]$pH %>% mean()  #average pH for cheese LBD
data[final_count$day30 > low_bound, ]$pH %>% mean()

data[final_count$day30 < high_bound & vat_count!=0, ]$vat_count %>% log10() %>% mean()  #average conc. for cheese without LBD
data[final_count$day30 < low_bound & vat_count!=0, ]$vat_count %>% log10() %>% mean()
data[final_count$day30 < high_bound, ]$pH %>% mean()  #average pH for cheese without LBD
data[final_count$day30 < low_bound, ]$pH %>% mean()

final_count$day30[final_count$day30 != -Inf] %>% mean()
final_count$day60[final_count$day60 != -Inf] %>% mean()
final_count$day90[final_count$day90 != -Inf] %>% mean()
final_count$day120[final_count$day120 != -Inf] %>% mean()


## Stack histograms of distribution at differnet day 
lbd = gather(final_count, key="day", value="logcount")
lbd = lbd %>% filter(day ==c("day60", "day90","day120"))
  
ggplot(lbd, aes(x=logcount, color=day,fill=day)) +
  geom_histogram(aes(y=..density..),position = "identity", alpha=0.5, binwidth = 0.1) +
  geom_density(alpha=.2)+
  scale_x_continuous(breaks = scales::pretty_breaks(10))+
  theme_classic()+
  theme(legend.position="top")+
  geom_vline(aes(xintercept=low_bound),
             color="red", linetype="dashed", size=1)+
  geom_vline(aes(xintercept=high_bound),
             color="red", linetype="dashed", size=1)+
  labs(title="C.tyrobutyricum logcount at different ripening time",x="logcount(MPN/kg)", y = "Density")



##Sensitivity analysis-----------------------------------------------------------------------------------------------------------


##What-if scenario ---------------------------------------------------------------------------------

# Microfiltration
high_count = data%>% filter(vat_count>=1000)
mean(high_count$day90>low_bound)
mean(high_count$day90>high_bound)

low_count = high_count %>% mutate(vat_count = round(vat_count*0.02))
for (i in 1:nrow(low_count)){
  low_count$block_count[i] = rpois(1, lambda = low_count$vat_count[i])
}
low_count$block_count = low_count$block_count*conc_factor
low_count$day90 = baranyi_log10N(t=90*24, lag=0, mumax=low_count$mumax,log10(low_count$block_count), log10(40000000))
mean(low_count$day90>low_bound)
mean(low_count$day90>high_bound)

## Cumulative probability of late blowing -------------------------------------------------------------------------------------------------
cum_prob = as.data.frame(t(result))
cum_prob = cum_prob %>% mutate(mean_prob = (High_prob+Low_prob)/2)
ggplot(data =cum_prob, aes(x=c(30,60,70,80,90,120),y=(mean_prob)))+
  geom_line(aes(y=mean_prob)) +
  geom_ribbon(aes(ymin = Low_prob, ymax= High_prob), fill = "grey70", alpha =0.3)+
  labs(title="Cumulative proportion of LBD at different ripening time",
       x="Ripening time (days)",
       y="Cumulative proportion")+
  theme_classic()+
  scale_x_continuous(breaks = scales::pretty_breaks(10))






## Contour plot--------------------------------------------------------------------------------------------------------------
# day 60
bound_low = data %>% filter(final_count$day60> (low_bound-0.0001),  final_count$day60< (low_bound+0.0001))
bound_high = data %>% filter(final_count$day60> (high_bound-0.001),  final_count$day60< (high_bound+0.001))
bound_level = c(rep("low", nrow(bound_low)),rep("high", nrow(bound_high)))
bound = bind_rows(bound_low, bound_high) %>% cbind(bound_level)

bound %>% ggplot(aes(x=pH,y=log10(vat_count), col=bound_level)) +
              geom_point()+
              geom_hline(yintercept=2, size =5, alpha=0.3)+
              geom_vline(xintercept=5.4, size=5, alpha=0.3)+
              theme_classic()+
  labs(title="Defect boundary line at day 60",x="pH", y = "Raw milk logcount (MPN/L)")

# day 90
bound_low = data %>% filter(final_count$day90> (low_bound-0.001),  final_count$day90< (low_bound+0.001))
bound_high = data %>% filter(final_count$day90> (high_bound-0.001),  final_count$day90< (high_bound+0.001))
bound_level = c(rep("low", nrow(bound_low)),rep("high", nrow(bound_high)))
bound = bind_rows(bound_low, bound_high) %>% cbind(bound_level)
bound %>% ggplot(aes(x=pH,y=log10(vat_count), col=bound_level)) +
  geom_point()+
  geom_hline(yintercept=2, size =5, alpha=0.3)+
  geom_vline(xintercept=5.4, size=5, alpha=0.3)+
  theme_classic()+
  labs(title="Defect boundary line at day 90",x="pH", y = "Raw milk logcount (MPN/L)")


## 3D plot
ripen30 = data %>% filter(final_count$day30> (high_bound-0.001),  final_count$day30< (high_bound+0.001))
ripen60 = data %>% filter(final_count$day60> (high_bound-0.001),  final_count$day60< (high_bound+0.001))
ripen70 = data %>% filter(final_count$day70> (high_bound-0.001),  final_count$day70< (high_bound+0.001)) 
ripen80 = data %>% filter(final_count$day80> (high_bound-0.001),  final_count$day80< (high_bound+0.001)) 
ripen90 = data %>% filter(final_count$day90> (high_bound-0.001),  final_count$day90< (high_bound+0.001))
ripen120 = data %>% filter(final_count$day120> (high_bound-0.001),  final_count$day120< (high_bound+0.001))
time = c(rep(30, nrow(ripen30)),rep(60, nrow(ripen60)),rep(70, nrow(ripen70)),rep(80, nrow(ripen80)),rep(90, nrow(ripen90)),rep(120, nrow(ripen120)))
ripen_data = rbind(ripen30, ripen60, ripen70, ripen80, ripen90, ripen120) %>% cbind(time)
fig = plot_ly(data=ripen_data, x=~pH, y=~vat_count, z=~time) %>% add_mesh()
fig

## Error analysis---------------------------------------------------------------------------------------------------

# lower likelihood
# late blown cheese
blow_low = vector()
for (i in 1:6){
  blow = data %>% filter(final_count[,i]>=low_bound)
  mean_pH = mean(blow$pH)
  sd_pH = sd(blow$pH)
  min_pH= min(blow$pH)
  max_pH = max(blow$pH)
  mean_vat_count=mean(log10(blow$vat_count))
  sd_vat_count = sd(log10(blow$vat_count))
  min_vat_count = min(log10(blow$vat_count))
  sum = cbind(mean_pH, sd_pH, min_pH, max_pH,mean_vat_count, sd_vat_count, min_vat_count)
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
  min_pH= min(blow$pH)
  max_pH = max(blow$pH)
  mean_block_count=mean(log10(blow$block_count))
  sd_block_count = sd(log10(blow$block_count))
  sum = cbind(mean_pH, sd_pH, min_pH, max_pH,mean_block_count, sd_block_count)
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
  mean_vat_count=mean(log10(blow$vat_count))
  sd_vat_count = sd(log10(blow$vat_count))
  sum = cbind(mean_pH, sd_pH, mean_vat_count, sd_vat_count)
  blow_high = rbind(blow_high, sum)
  
}
rownames(blow_high) = name
blow_high


