## What-if scenarios##
###################### 

## 1. Bactofugation/microfiltration
## effect = 98% reduction of spores in raw milk

vat_count = 0.02*vat_count



## 2. Addition of nitrate at 2.5g/ 100L milk
## Growth rate = 0 for raw milk with <1 spore/mL
  
index = which(data$vat_count < 1000)
data$mumax[index] = 0

  

## 3. Addition of lysozyme at 2.5g/ 100L milk
## Growth rate = 0 for raw milk with <0.3 spore/mL
index = which(data$vat_count < 300)
data$mumax[index] = 0
  

## 4. Bacteriocinogenic LAB at 0.3% level
## Lower the mumax by 1.86
data$mumax = data$mumax/1.86



## 5. Premium raw milk
## Truncate raw milk spore count over 25 MPN/L
index = which(data$vat_count<=25)
final_count = final_count[index,]


## 6. Lower ripening temperature
## Set temp = 12C
temp = 12

  
  
