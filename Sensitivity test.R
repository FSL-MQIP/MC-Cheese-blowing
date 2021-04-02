## Sensitivity analysis##
#########################

#Load package
library(epiR)


#PRCC
input = data.frame(data$vat_count, data$pH)
output = data$day60
df = data.frame(input, output)
epi.prcc(df)
