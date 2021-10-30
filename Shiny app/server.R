#Load file
Spore = read.csv("Initial spore concentration.csv", header = T)

#Load pacakge
library(dplyr)
library(tidyr)
library(EnvStats)
library(ggplot2)
library(shiny)

##Data cleaning
colnames(Spore)[1] = "ID"
Spore = gather(Spore,"Month","Count",-"ID")

##Create data frame
Spore$Count[Spore$Count == "" ] = NA
count = na.omit(Spore$Count)
count[count == "<18"] = "4.5"
count_t = as.numeric(count)^0.1

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

salt2aw = function(salt){
    aw = 0.995-0.00721*salt
    aw
}


##--------------------------------------------------------------------------------------------------------------
## Base model
set.seed(1)
n_sim = 10000
vat_id = rep(1:n_sim, each = 20)
block_id = rep(1:20, n_sim)

## Simulate var spore count as MPN/kg milk in cheese vat
mean_to_sd = mean(count_t)/sd(count_t)
mean_count = (10^1.85)^0.1
vat_norm_count = rnorm(n_sim, mean_count, mean_count/mean_to_sd)
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

# Threshold level calculation
low_mumax = mumax(13,5.0,salt2aw(3.7))
high_mumax = mumax(13,5.0,salt2aw(3.7))
low_bound = baranyi_log10N(t=70*24, lag=0, mumax=low_mumax,  log10(1300), log10(40000000))
high_bound = baranyi_log10N(t=60*24, lag=0, mumax=high_mumax,  log10(2500), log10(40000000))


# Model simulation
final_count = vector()
for (i in c(30,60,90,120)){
    new = baranyi_log10N(t=i*24, lag=0, mumax=data$mumax,log10(data$block_count), log10(40000000))
    final_count = cbind(final_count, new) %>% round(3)
}

name = c("day30","day60","day90","day120")
colnames(final_count) = name
final_count = as.data.frame(final_count)


## Check results
result = vector()
for (i in 1:ncol(final_count)){
    high = round(mean(final_count[i]>low_bound)*100,3)
    low = round(mean(final_count[i]>high_bound)*100,3)
    avg = round(mean(c(high,low)),3)
    prob = c(high, low, avg)
    result = cbind(result, prob)
}
colnames(result) = name
rownames(result) = c("Upper_percent", "Lower_percent", "Avg_percent")

## Base plot
cum_prob_base = as.data.frame(t(result))
cum_prob_base = cum_prob_base %>% mutate(mean_prob = (Upper_percent+Lower_percent)/2)
base_plot = ggplot(data =cum_prob_base, aes(x=c(30,60,90,120),y=(mean_prob)))+
    geom_line(aes(y=mean_prob)) +
    geom_ribbon(aes(ymin = Lower_percent, ymax= Upper_percent), fill = "grey70", alpha =0.3)+
    labs(title="Default mean spore concentraton (1.85log10 MPN/L milk); no intervention)",
         x="Ripening time (days)",
         y="Cumulative proportion (%)")+
    theme_classic()+
    ylim(0,100)+
    scale_x_continuous(breaks = scales::pretty_breaks(10))
#------------------------------------------


# Define server logic 
shinyServer(function(input, output){
    
    model_result = reactive({
    ## Cheese vat and block identifications
    set.seed(1)
    n_sim = 1000
    vat_id = rep(1:n_sim, each = 20)
    
    ## Simulate var spore count as MPN/kg milk in cheese vat
    mean_to_sd = mean(count_t)/sd(count_t)
    mean_count = (10^input$count)^0.1
    vat_norm_count = rnorm(n_sim, mean_count, mean_count/mean_to_sd)
    
    vat_sim_count = vat_norm_count^10
    
    ## Combine ids with initial spore count
    vat_count = round(rep(vat_sim_count, each = 20))
    if (input$mf){vat_count = 0.02*vat_count}
    if (input$bf){vat_count = 0.02*vat_count}
    
    block_id = rep(1:20, length(vat_id))
    block_count = vector()
    for (i in 1:length(vat_count)){
        block_count[i] = rpois(1, lambda = vat_count[i])
    }
    
    # Consider concentration effect during whey draining
    conc_factor = 6#+input$cf_change
    block_count = conc_factor*block_count
    
    # Combining data
    data = data.frame(vat_id, vat_count,block_id,block_count)
    
    ## Ripening parameters
    temp = input$temp

    
    pH= rtri(length(vat_id), min = 5.1, max = 5.63, mode = 5.4)
    data$pH = round(rep(pH, each = 20),2)
    data$pH = data$pH #+ input$pH_change
    
    data$mumax = round(mumax(temp,data$pH,aw=salt2aw(4.9)),5)
    if (input$nitrate) {index = which(data$vat_count < 1000)
    data$mumax[index] = 0} 
    if (input$lysozyme) {index = which(data$vat_count < 300)
    data$mumax[index] = 0}
    if (input$lab) {data$mumax = data$mumax/1.86}
    #data$mumax = (1+input$mu_change)*data$mumax

    
    # Predicted conc
    final_count = vector()
    for (i in c(30,60,90,120)){
        new = baranyi_log10N(t=i*24, lag=0, mumax=data$mumax,log10(data$block_count), log10(40000000))
        final_count = cbind(final_count, new) %>% round(3)
    }
    
    name = c("day30","day60","day90","day120")
    colnames(final_count) = name
    final_count = as.data.frame(final_count)
    
    #high_bound = high_bound+input$tl_change
    #low_bound = low_bound+input$tl_change
    ## Check results
    result = vector()
    for (i in 1:ncol(final_count)){
        high = round(mean(final_count[i]>low_bound)*100,3)
        low = round(mean(final_count[i]>high_bound)*100,3)
        avg = round(mean(c(high,low)),3)
        prob = c(high, low, avg)
        result = cbind(result, prob)
    }
    colnames(result) = name
    rownames(result) = c("Upper_percent", "Lower_percent", "Avg_percent")
    result
    })

    output$plot <- renderPlot({
        
        base = !(input$mf | input$bf | input$nitrate | input$lysozyme | input$lab | input$temp !=13 )
                     #input$tl_change!=0 | input$pH_change!=0 | input$mu_change!=0|input$cf_change!=0)
        if (input$count == 1.85 & base) {base_plot}
        else {
        cum_prob = as.data.frame(t(model_result()))
        cum_prob = cum_prob %>% mutate(mean_prob = (Upper_percent+Lower_percent)/2)
        ggplot(data =cum_prob, aes(x=c(30,60,90,120),y=(mean_prob)))+
            geom_line(aes(y=mean_prob)) +
            geom_ribbon(aes(ymin = Lower_percent, ymax= Upper_percent), fill = "grey70", alpha =0.3)+
            labs(title="What-if scenarios",
                 x="Ripening time (days)",
                 y="Cumulative proportion (%)")+
            theme_classic()+
            ylim(0,100)+
            scale_x_continuous(breaks = scales::pretty_breaks(10))
        }
    })
    # Predict the proportion with uncertainty
    output$prop <- renderPrint({
        
        base = !(input$mf | input$bf | input$nitrate | input$lysozyme | input$lab | input$temp !=13)
                     #input$mu_change!=0 |
                     #input$tl_change!=0 | input$pH_change!=0 |input$cf_change!=0)
        if (input$count == 1.85 & base) {result}
        else {model_result()}
    })
        
    output$base_plot <- renderPlot({
        base_plot
    })
    
    output$base_prop <- renderPrint({
        result
    })
    

})
