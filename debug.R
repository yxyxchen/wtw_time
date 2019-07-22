

modelName = "QL1"
a = read.csv(sprintf("genData/expModelFittingCV/%s/s1_f6.txt", modelName),
             header = F)
# load data
id = ids[sIdx]
load(sprintf("genData/expModelFittingCV/split/s%d.RData", id))
thisTrialData = trialData[[id]]
# excluded some trials
excluedTrialsHP = which(thisTrialData$trialStartTime > (blockSecs - tMaxs[1]) &
                          thisTrialData$condition == "HP")
excluedTrialsLP = which(thisTrialData$trialStartTime > (blockSecs - tMaxs[2]) &
                          thisTrialData$condition == "LP")
excluedTrials = c(excluedTrialsHP, excluedTrialsLP)
thisTrialData = thisTrialData[!(1 : nrow(thisTrialData)) %in% excluedTrials,]

f = 7
thisTrialData = thisTrialData[1 : nrow(thisTrialData) %in% as.vector(partTable[-f,]),]

# prepare the data
nTrial = length(thisTrialData$trialEarnings)
cond = thisTrialData$cond
trialEarnings = thisTrialData$trialEarnings
timeWaited = pmin(thisTrialData$timeWaited, max(tMaxs))
Ts = round(ceiling(timeWaited / stepDuration) + 1)
scheduledWait = thisTrialData$scheduledWait

source("subFxs/modelComparisonFxs.R")
likFun = getLikFun(modelName)
paras = as.double(a[1, 1:4])
lik_ = likFun(paras, cond, trialEarnings, timeWaited)$lik_


sum(sapply(1 : length(thisTrialData$blockNum), function(i){
  if(trialEarnings[i] > 0){
    junk = log(lik_[1 : max(Ts[i]-1, 1), i])
    junk[is.infinite(junk)] = -10000
    sum(junk)
  }else{
    junk = c(log(lik_[1:max(Ts[i] - 2,1), i]), log(1-lik_[Ts[i] - 1, i]))
    junk[is.infinite(junk)] = -10000
    sum(junk)
  }
}))


c = fit %>% rstan::extract(permuted = F, pars = c("LL_all", "Qwaits[20,1]", "Vitis[1]", "prior", "gamma")) %>% 
  adply(2, function(x) x) %>% dplyr::select(-chains) 


fit %>% rstan::extract(permuted = F, pars = c("LL_all", "Qwaits[20,1]")) %>% 
  adply(2, function(x) x) %>% dplyr::select(-chains) 