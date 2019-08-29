modelFitting = function(thisTrialData, fileName, paraNames, model, modelName){
  # load experiment paras
  load('expParas.RData')
  
  # rstan parameters 
  nChain = 4  
  nIter = 100
  controlList = list(adapt_delta = 0.99, max_treedepth = 11)
  
  # duration of one time step (namely one temporal state) 
  stepSec = 1
  
  # prepare inputs for fitting the model
  condition = unique(thisTrialData$condition)
  ## since each participant went through two environments
  ## we use max(tMaxs) to calculate the  maximal number of steps
  nStepMax = max(tMaxs) / stepSec
  ## ensure timeWaited = scheduledWait on rewarded trials
  thisTrialData = within(thisTrialData, {timeWaited[trialEarnings!= 0] = scheduledWait[trialEarnings!= 0]})
  ## terminal state in each trial
  ## Noticeably, ceiling(timeWaited / stepSec) gives the num of steps
  ## and adding 1 upon it gives the terminal state, namely the state following the final action. 
  Ts = with(thisTrialData, {round(ceiling(timeWaited / stepSec) + 1)}) 
  ## orgianze inputs into a list
  inputs <- list(
    iti = iti,
    stepSec = stepSec,
    nStepMax = nStepMax,
    N = length(thisTrialData$blockNum), # number of trials
    Rs = thisTrialData$trialEarnings, # rewards on each trial
    Ts = Ts)
  if(modelName %in% c("QL1", "QL2")){
    ## in Q-learning, the initial value of the iti state is proportional to 
    ## the discounted total rewards averaged across two conditions
    ## discount factor for one step is 0.85
    VitiIni = 0.9 * mean(unlist(optimRewardRates) * stepSec / (1 - 0.85))
    inputs$VitiIni = VitiIni
  }else{
    ## in R-learning, the initial reward rate is proportional to
    ## the optimal reward rates averaged across two conditions
    reRateIni = 0.9 * mean(unlist(optimRewardRates)) * stepSec;
    inputs$reRateIni = reRateIni     
  }
  
  # extract subName from fileName
  subName = str_extract(fileName, pattern = "s[0-9]*")
  # fit the model
  withCallingHandlers({
    fit = sampling(object = model, data = inputs, cores = 1, chains = nChain,
                   iter = nIter, control = controlList) 
    print(sprintf("Finish %s !", subName))
    write(sprintf("Finish %s !", subName), sprintf("outputs/%s_log.txt", modelName), append = T, sep = "n")
  }, warning = function(w){
    print(sprintf("Finish %s !", subName))
    write(sprintf("Finish %s !", subName), sprintf("outputs/%s_log.txt", modelName), append = T, sep = "n")
    warnText = paste(modelName, subName, w)
    write(warnText, sprintf("outputs/%s_log.txt", modelName), append = T, sep = "n")
  })
  
  # extract posterior samples
  samples = fit %>%
    rstan::extract(permuted = F, pars = c(paraNames, "LL_all"))
  # save posterior samples
  samples = samples %>% adply(2, function(x) x) %>% dplyr::select(-chains) 
  write.table(matrix(unlist(samples), ncol = length(paraNames) + 1), file = sprintf("%s.txt", fileName), sep = ",",
              col.names = F, row.names=FALSE) 
  # calculate WAIC and Efficient approximate leave-one-out cross-validation (LOO)
  log_lik = extract_log_lik(fit) 
  WAIC = waic(log_lik)
  LOO = loo(log_lik)
  save("WAIC", "LOO", file = sprintf("%s_waic.RData", fileName))
  # summarise posterior parameters and LL_all
  fitSummary <- summary(fit, pars = c(paraNames, "LL_all"), use_cache = F)$summary
  write.table(matrix(fitSummary, nrow = length(paraNames) + 1), file = sprintf("%s_summary.txt", fileName),  sep = ",",
              col.names = F, row.names=FALSE)
}