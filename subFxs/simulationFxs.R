getSimFun = function(modelName){
  if(modelName == "QL1") simFun = QL1
  else if(modelName == "QL2") simFun = QL2
  else if(modelName == "RL1") simFun = RL1
  else if(modelName == "RL2") simFun = RL2
  else if(modelName == "BL") simFun = BL
  else{
    return("wrong model name!")
  }
  return(simFun)
}

simulateUnit = function(paras, nSim, modelName, cb){
  # get simFun
  simFun = getSimFun(modelName)
  
  # initialize outputs
  aucHP_ = vector(length = nSim)
  aucLP_ = vector(length = nSim)
  wtwHP_ = matrix(NA, nrow = length(tGrid), ncol = nSim)
  wtwLP_ = matrix(NA, nrow = length(tGrid), ncol = nSim)
  reRate_ = vector(length = nSim)
  # usually can not use foreach to fill a matrix or a vector
  for(i in 1 : nSim){
    set.seed(i)
    thisTrialData = simFun(paras, cb)
    thisTrialData$Qwaits = NULL
    thisTrialData = as.data.frame(thisTrialData)
    # HP
    kmscResults = kmsc(thisTrialData[thisTrialData$condition == "HP",], min(tMaxs), "", F, kmGrid)
    aucHP_[i] = kmscResults$auc
    wtwtsResults = wtwTS(thisTrialData[thisTrialData$condition == "HP",], tGrid, min(tMaxs), "", F )
    wtwHP_[,i] = wtwtsResults$timeWTW
    # LP
    kmscResults = kmsc(thisTrialData[thisTrialData$condition == "LP",], min(tMaxs), "", F, kmGrid)
    aucLP_[i] = kmscResults$auc
    wtwtsResults = wtwTS(thisTrialData[thisTrialData$condition == "LP",], tGrid, min(tMaxs), "", F )
    wtwLP_[,i] = wtwtsResults$timeWTW
    
    junk = nrow(thisTrialData)
    
    reRate_[i] = mean(thisTrialData$reRates[(junk - 10) : junk])
  }
  # summarise 
  outputs = list(aucHP = mean(aucHP_),
                 aucLP = mean(aucLP_),
                 aucHPSD = sd(aucHP_),
                 aucLPSD = sd(aucLP_),
                 wtwHP = apply(wtwHP_, MARGIN = 1, mean),
                 wtwLP = apply(wtwLP_, MARGIN = 1, mean),
                 wtwHPSD = apply(wtwHP_, MARGIN = 1, sd),
                 wtwLPSD = apply(wtwHP_, MARGIN = 1, sd),
                 reRate = mean(reRate_))
  return(outputs)
}
RL2 = function(paras, cb){
  # parse para
  phi = paras[1]; phiP = paras[2]; tau = paras[3]; prior = paras[4]
  beta = paras[5]; betaP = paras[6]
  
  # prepare inputs
  nTrialMax = blockSecs / iti *2
  tMax= max(tMaxs)
  nTimeStep = tMax / stepDuration
  
  # initialize action values
  subOptimalRatio = 0.9
  wIni = mean(as.double(optimRewardRates)) * stepDuration * subOptimalRatio
  Viti = 0; reRate = wIni 
  Qwait = prior*0.1 - 0.1*(0 : (nTimeStep - 1)) + Viti
  
  # recording variables
  reRates = vector(length = nTrialMax); reRates[1] = reRate
  Vitis = vector(length = nTrialMax); Vitis[1] = Viti
  Qwaits = matrix(NA, nrow = nTimeStep, ncol = nTrialMax); Qwaits[,1] = Qwait 
  
  # initialize outputs 
  trialEarnings = rep(0, nTrialMax); timeWaited = rep(0, nTrialMax); sellTime = rep(0, nTrialMax);
  condition = rep(0, nTrialMax); scheduledWait = rep(0, nTrialMax);
  
  
  # loop over blocks
  tIdx = 1
  for(bkIdx in 1 : 2){
    # determine distrib
    seq = c()
    distrib = ifelse(cb[bkIdx] == "HP", "unif16", "exp32")
    elapsedTime = 0
    # loop over trials
    while(elapsedTime <= blockSecs) {
      # determine scheduledWait
      junk = drawSample(distrib, seq)
      thisScheduledWait = junk[['delay']]
      seq = junk[['seq']]
      # loop over timesteps
      t = 1
      while(t <= nTimeStep){
        # determine At
        pWait =  1 / sum(1  + exp((Viti -reRate - Qwait[t])* tau))
        action = ifelse(runif(1) < pWait, 'wait', 'quit')
        # observe St+1 and Rt+1
        rewardOccur = thisScheduledWait <= (t * stepDuration) && thisScheduledWait > ((t-1) * stepDuration)
        getReward = (action == 'wait' && rewardOccur);
        nextReward = ifelse(getReward, tokenValue, 0) 
        # dertime whether St+1 is the terminal state
        nextStateTerminal = (getReward || action == "quit")
        if(nextStateTerminal){
          elapsedTime = elapsedTime + ifelse(getReward, thisScheduledWait, t * stepDuration) + iti
          # only record values before the end of the block
          if(elapsedTime<= blockSecs){
            T = t+1
            trialEarnings[tIdx] = ifelse(nextReward == tokenValue, tokenValue, 0);
            timeWaited[tIdx] = ifelse(getReward, thisScheduledWait, t * stepDuration)
            sellTime[tIdx] = elapsedTime - iti
            condition[tIdx] = cb[bkIdx]
            scheduledWait[tIdx] = thisScheduledWait
          }
          break
        }else{
          t = t + 1
        }
      }# end of the loop over timesteps
      # update action values before the end of the block
      if(elapsedTime <= blockSecs){
        returns = sapply(1 : (T-1), function(t) nextReward - reRate * (T-t) + Viti)
        if(getReward){
          Qwait[1 : (T-1)] = Qwait[1 : (T-1)] + phi*(returns[1 : (T-1)] - Qwait[1 : (T-1)])
        }else{
          if(T > 2){
            Qwait[1 : (T-2)] = Qwait[1 : (T-2)] + phiP*(returns[1 : (T-2)] - Qwait[1 : (T-2)])
          }
        }
        # update Viti
        delta = (returns[1] - reRate * (iti / stepDuration) - Viti)
        Viti = ifelse(nextReward > 0, Viti + phi * delta, Viti + phiP* delta)
        # update reRate 
        reRate = ifelse(nextReward > 0, reRate + beta * delta, reRate + betaP* delta)   
        
        # record variables
        reRates[tIdx + 1] = reRate
        Vitis[tIdx + 1] = Viti
        Qwaits[,tIdx + 1] = Qwait
        # update tIdx and go to the next trial
        tIdx = tIdx + 1
      }
    } # end of the loop over trials
  }# end of the loop over blocks
  
  # cut off the last trial and return outputs
  trialNum = c(1 : sum(condition == cb[1]), 1 : sum(condition == cb[2]))
  outputs = list( 
    "trialNum" = trialNum, "trialEarnings" = trialEarnings[1 : (tIdx - 1)],
    "timeWaited" = timeWaited[1 : (tIdx - 1)], "sellTime" = sellTime[1 : (tIdx - 1)],
    "scheduledWait" = scheduledWait[1 : (tIdx - 1)], "condition" = condition[1 : (tIdx - 1)],
    "reRates" = reRates[1: (tIdx - 1)], "Vitis" = Vitis[1 : (tIdx - 1)],
    "Qwaits" = Qwaits[,1 : (tIdx - 1)]
  )
  return(outputs)
}
