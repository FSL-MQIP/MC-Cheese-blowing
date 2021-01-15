## Sensitivity analysis

# Define the function to help sensitivity analysis
prob_atday60 = function(N0_factor = 0,
                        mumax_factor = 0,
                        aw_min_factor = 0,
                        threshold_factor = 0){
  set.seed(1)
  n_sim = 10000
  vat_id = rep(1:n_sim, each = 20)
  block_id = rep(1:20, n_sim)
  
  ## Simulate var spore count as MPN/kg milk in cheese vat
  vat_norm_count = rnorm(n_sim, mean(count_t), sd(count_t))
  vat_sim_count = vat_norm_count^10
  
  ## Combine ids with initial spore count
  vat_count = round(rep(vat_sim_count, each = 20))
  vat_count = ifelse(vat_count == 0, 0, (10^(log10(vat_count)*(1+N0_factor))))
  
  
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
  data$mumax = round(mumax(temp,data$pH,aw=salt2aw(4.9))*(1+mumax_factor),5)
  
  # Threshold level calculation
  low_mumax = mumax(13,5.04,salt2aw(3.7))
  high_mumax = mumax(13,5.01,salt2aw(3.7))
  low_bound = baranyi_log10N(t=70*24, lag=0, mumax=low_mumax,  log10(250), log10(40000000))*(1+threshold_factor)
  high_bound = baranyi_log10N(t=60*24, lag=0, mumax=high_mumax,  log10(2500), log10(40000000))*(1+threshold_factor)
  
  #
  #data = data %>% filter(vat_count>=1000)
  
  # Model simulation
  final= baranyi_log10N(t=60*24, lag=0, mumax=data$mumax,log10(data$block_count), log10(40000000))
  final = final %>% round(3)
  
  high = mean(final>low_bound)
  low = mean(final>high_bound)
  prob = c(high, low)
  prob

}
## Testing ------------------------------------------------------------------------------------------
## N0
prob_atday60(N0_factor = 0)
prob_atday60(N0_factor = 0.2)
prob_atday60(N0_factor = 0.4)
prob_atday60(N0_factor = -0.2)
prob_atday60(N0_factor = -0.4)

## mumax
prob_atday60(mumax_factor = 0.2)
prob_atday60(mumax_factor = 0.4)
prob_atday60(mumax_factor = -0.2)
prob_atday60(mumax_factor = -0.4)

## Threshold level
prob_atday60(threshold_factor = 0.2)
prob_atday60(threshold_factor = 0.4)
prob_atday60(threshold_factor = -0.2)
prob_atday60(threshold_factor = -0.4)
