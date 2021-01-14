## Functions ----------------------------------------------------------------------
baranyi_log10N = function(t,lag,mumax,LOG10N0,LOG10Nmax) {
  ans <- LOG10Nmax + log10((-1 + exp(mumax * lag) + exp(mumax * t))/(exp(mumax * t) - 1 + exp(mumax * lag) * 10^(LOG10Nmax -LOG10N0)))
  return(ans)
}

final_conc = function(N0, temp = 37, pH = 5.8, aw = 1, aw_min, hour){
  #Default growth rate 
  mumax = 3.49
  
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

pH_LST = function(days){
  9e-05*days^2-0.0058*days+5.05
}

pH_HST = function(days){
  9e-05*days^2-0.0053*days+5.01
}

aw_LST = function(days){
  -6e-05*days+0.9871
}

aw_HST = function(days){
  -2e-06*days^2+8e-05*days+0.985
}

salt2aw = function(salt){
  pH = 0.995-0.00721*salt
  pH
}

## Simulation set-up -----------------------------------------------------------
## Initial spore concentration
N0_LST = 1200
N0_HST = 2500
pH_LST = 5.04
pH_HST = 5.01

## Ripening parameters
temp = 13


# Last day of growth
lday = 70


## Iteratively simulate aw and counts at each day
final_count = vector()
final_count = final_conc(N0_LST, temp = temp, pH=5.01, aw=1,aw_min=salt2aw(3), hour=1*24)
final_count = cbind(final_count, 
                    final_conc(final_count, temp = temp, pH=5.01, aw=1,aw_min=salt2aw(3), hour=1*24))
for (i in 3:lday){
  N0 = final_count[,i-1]
  new = final_conc(N0, temp = temp, pH=5.01, aw=1,aw_min=salt2aw(3),hour=1*24)
  final_count = cbind(final_count, new)
}
final_count = round(final_count)

## Create the column names for count at each day
var_name = vector()
for (i in 1:lday){
  var_name[i] = paste("day",as.character(i),sep = "_")
}
colnames(final_count) = var_name



## Check results
log10(final_count[60])



##
final_conc(N0_LST, temp = temp, pH = 5.04, aw = 1, aw_min = salt2aw(3), hour = 60*24) %>% log10()
final_conc(N0_HST, temp = temp, pH = 5.01, aw = 1, aw_min = salt2aw(3), hour = 60*24) %>% log10()
