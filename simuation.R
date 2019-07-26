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
expPara = loadExpPara(paraNames, dirName = sprintf("genData/expModelFitting/%sdb"), modelName)
useID = getUseID(expPara, paraNames)
if(length(useID) != nSub){
  cat("The model doesn't converge completely!")
}

# determine paraSamples
nSample = 10
paraSamples = matrix(NA, nrow = nSample, ncol = nPara)
for(i in 1 : nPara){
  junk = expPara[,paraNames[i]]
  paraSamples[,i] = seq(min(junk), max(junk), length.out = nSample)
}

# observe the effect of 