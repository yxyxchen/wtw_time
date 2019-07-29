# this script is used to demonstrate the effect of 

library('ggplot2')
library('plyr')
library('dplyr')
library('tidyr')
load("wtwSettings.RData")
source('subFxs/repetitionFxs.R') # called by simulate 
source("subFxs/helpFxs.R") # getParas
source("subFxs/loadFxs.R") # load scheduledWait from empirical data
source("subFxs/analysisFxs.R") 
source("subFxs/plotThemes.R")
source("subFxs/simulationFxs.R")
source("subFxs/taskFxs.R")
# loop over participants 
library("doMC")
library("foreach")
nCore = parallel::detectCores() -1 # only for the local computer
registerDoMC(nCore)

# load expData
allData = loadAllData()
hdrData = allData$hdrData           
trialData = allData$trialData       
ids = hdrData$ID          
nSub = length(ids)   

# load expPara
modelName = "RL2"
paraNames = getParaNames(modelName)
nPara = length(paraNames)
expPara = loadExpPara(paraNames, dirName = sprintf("genData/expModelFitting/%sdb", modelName))
useID = getUseID(expPara, paraNames)
if(length(useID) != nSub){
  cat("The model doesn't converge completely!")
}

# determine paraSamples
nCut = 10
paraSamples = matrix(NA, nrow = nCut, ncol = nPara)
for(i in 1 : nPara){
  junk = expPara[,paraNames[i]]
  paraSamples[,i] = seq(min(junk), max(junk), length.out = nCut)
}

# median paras
medianParas = sapply(1 : nPara, function(i) median(expPara[,i]))

# simulate 
# determine simFun
simFun = getSimFun("RL2")
nSim = 5
# initialize 
auc_ = array(NA, dim = c(nSim, nCut, nPara, nCond))
wtw_ = array(NA, dim = c(length(tGrid), nSim, nCut, nPara, nCond))
aucSD_ = array(NA, dim = c(nSim, nCut, nPara, nCond))
wtwSD_ = array(NA, dim = c(length(tGrid), nSim, nCut, nPara, nCond))
set.seed(123)
for(cutIdx in 1 : nCut){
  paras = medianParas
  paras[pIdx] = paraSamples[cutIdx,pIdx]
  for(smIdx in 1 : nSim){

  }
  wtwHP_[[cutIdx]] = wtwMatrixHP
  wtwLP_[[cutIdx]] = wtwMatrixLP
}

modelName = "RL2"
cb = c("HP", "LP")
simulateUnit = function(paras, nSim, modelName, cb){
  # get simFun
  simFun = getSimFun(modelName)

  # initialize outputs
  aucHP_ = vector(length = nSim)
  aucLP_ = vector(length = nSim)
  wtwHP_ = matrix(NA, nrow = length(tGrid), ncol = nSim)
  wtwLP_ = matrix(NA, nrow = length(tGrid), ncol = nSim)
  
  # usually can not use foreach to fill a matrix or a vector
  for(i in 1 : nSim){
    set.seed(i)
    thisTrialData = simFun(paras, cb)
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
  }
  # summarise 
  outputs = list(aucHP = mean(aucHP_),
                 aucLP = mean(aucLP_),
                 aucHPSD = sd(aucHP_),
                 aucLPSD = sd(aucLP_),
                 wtwHP = apply(wtwHP_, MARGIN = 1, mean),
                 wtwLP = apply(wtwLP_, MARGIN = 1, mean),
                 wtwHPSD = apply(wtwHP_, MARGIN = 1, sd),
                 wtwLPSD = apply(wtwHP_, MARGIN = 1, sd))
  return(outputs)
}

# HP

# LP
kmscResults = kmsc(thisTrialData[thisTrialData$condition == "LP",], min(tMaxs), "", F, kmGrid)
aucLP_[smIdx, cutIdx] = kmscResults$auc
wtwtsResults = wtwTS(thisTrialData[thisTrialData$condition == "LP",], tGrid, min(tMaxs), "", F)
wtwMatrixLP[, smIdx] = wtwtsResults$timeWTW
# effects of a single parameter 

