# differences from the cluster version can be found by searching # needed in local computers
getMatData = function(){
  # create outfiles
  dir.create("genData") 
  dir.create("genData/expModelFitting") 
  dir.create(sprintf("genData/expModelFitting/%s", modelName))
  
  #load libraries
  library('plyr'); library(dplyr); library(ggplot2);library('tidyr');
  library("loo")
  library("coda") 
  source('subFxs/modelFittingFxs.R') # for fitting each single participant
  source('subFxs/loadFxs.R') # for load data
  source("subFxs/helpFxs.R") # for getparaNames
  source("subFxs/analysisFxs.R") # for block2session
  load("wtwSettings.RData")
  library("R.matlab")
  
  # load expData
  allData = loadAllData()
  hdrData = allData$hdrData           
  trialData = allData$trialData       
  idList = hdrData$ID                   
  n = length(idList)                    
  
  # loop over subjects
  dir.create("matData")
  for(i in 1 : n){
    thisID = idList[[i]]
    thisTrialData = trialData[[thisID]]
    # excluded some trials
    excluedTrials1 = which(thisTrialData$trialStartTime > (blockSecs - tMaxs[1]) &
                             thisTrialData$condition == conditions[1])
    excluedTrials2 = which(thisTrialData$trialStartTime > (blockSecs - tMaxs[2]) &
                             thisTrialData$condition == conditions[2])
    excluedTrials = c(excluedTrials1, excluedTrials2)
    thisTrialData = thisTrialData[!(1 : nrow(thisTrialData)) %in% excluedTrials,]
    # 
    timeWaited = thisTrialData$timeWaited;
    scheduledWait = thisTrialData$scheduledWait
    trialEarnings = thisTrialData$trialEarnings
    timeWaited[trialEarnings != 0] = scheduledWait[trialEarnings != 0]
    # initialize output
    subData = list("timeWaited" =  timeWaited,
                   "trialEarnings" = thisTrialData$trialEarnings,
                   "id" = thisID)
    writeMat(sprintf("matData/sub%d.mat", i), subData = subData)
  }
}
