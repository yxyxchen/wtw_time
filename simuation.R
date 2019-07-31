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

# median paras
medianParas = sapply(1 : nPara, function(i) median(expPara[,i]))
meanParas = sapply(1 : nPara, function(i) mean(expPara[,i]))
                
meanParas = meanParas[1:5]
meanParas[2] = meanParas[2] / meanParas[1] 

medianParas = medianParas[1:5]
medianParas[2] = medianParas[2] / medianParas[1] 
# simulate for the original paras
nSim = 10
modelName = "RL2"

cb = c("HP", "LP")
tempt = simulateUnit(meanParas, nSim, modelName, cb)


# determine paraSamples
nCut = 5
paraSamples = matrix(NA, nrow = nCut, ncol = (nPara-1))
for(i in 1 : (nPara - 1)){
  if(i == 2){
    junk = expPara[,2] / expPara[,1]
  }else{
    junk = expPara[,i]
  } 
  qt2 = quantile(junk, 0.2)
  qt8 = quantile(junk, 0.8)
  paraSamples[,i] = seq(qt2, qt8, length.out = nCut)
}



# simulate 
# determine simFun
modelName = "RL2"
cb = c("HP", "LP")
nSim = 5
nCond = length(conditions)
# initialize 
auc_ = array(NA, dim = c(nCut, nPara - 1, nCond))
wtw_ = array(NA, dim = c(length(tGrid), nCut, nPara - 1, nCond))
aucSD_ = array(NA, dim = c(nCut, nPara - 1, nCond))
wtwSD_ = array(NA, dim = c(length(tGrid), nCut, nPara - 1, nCond))
reRate_ = array(NA, dim = c(nCut, nPara - 1))
for(pIdx in 1 : (nPara-1)){
  for(cutIdx in 1 : nCut){
    paras = meanParas
    paras[pIdx] = paraSamples[cutIdx,pIdx]
    tempt = simulateUnit(paras, nSim, modelName, cb)
    auc_[cutIdx, pIdx, 1] = tempt$aucHP
    auc_[cutIdx, pIdx, 2] = tempt$aucLP
    aucSD_[cutIdx, pIdx, 1] = tempt$aucHPSD
    aucSD_[cutIdx, pIdx, 2] = tempt$aucLPSD
    wtw_[ , cutIdx, pIdx, 1] = tempt$wtwHP
    wtw_[ , cutIdx, pIdx, 2] = tempt$wtwLP
    wtwSD_[ , cutIdx, pIdx, 1] = tempt$wtwHPSD
    wtwSD_[ , cutIdx, pIdx, 2] = tempt$wtwLPSD
    
    reRate_[cutIdx, pIdx] = tempt$reRate 
  }
}

# reorganizd the data 
junk = data.frame(rbind(auc_[,,1], auc_[,,2])) 
names(junk) = paraNames[1:5]
junk %>% as_tibble() %>% 
  mutate("condition" = rep(conditions, each = nCut),
         "percentile" = rep(1:nCut, 2)) %>%
  gather(key = para, value = auc, -condition, -percentile) %>%
  mutate(para = factor(para, levels = paraNames)) %>%
  ggplot(aes(percentile, auc)) + geom_point() + facet_grid(condition~para)

# interactions between paras
cb = c("LP", "HP")
for(p1 in 1 : (nPara - 1)){
  for(p2 in (p1 + 2) : (nPara - 1)){
    auc_ = array(NA, dim = c(nCut, nCut, nCond))
    wtw_ = array(NA, dim = c(length(tGrid), nCut, nCut, nCond))
    aucSD_ = array(NA, dim = c(nCut, nCut, nCond))
    wtwSD_ = array(NA, dim = c(length(tGrid), nCut, nCut, nCond))
    reRate_ = array(NA, dim = c(nCut, nCut))
    for(c1 in 1 : nCut){
      for(c2 in 1 : nCut){
        paras = meanParas
        paras[p1] = paraSamples[c1, p1]
        paras[p2] = paraSamples[c2, p2]
        tempt = simulateUnit(paras, nSim, modelName, cb)
        auc_[c1, c2, 1] = tempt$aucHP
        auc_[c1, c2, 2] = tempt$aucLP
        aucSD_[c1, c2, 1] = tempt$aucHPSD
        aucSD_[c1, c2, 2] = tempt$aucLPSD
        wtw_[ , c1, c2, 1] = tempt$wtwHP
        wtw_[ , c1, c2, 2] = tempt$wtwLP
        wtwSD_[ , c1, c2, 1] = tempt$wtwHPSD
        wtwSD_[ , c1, c2, 2] = tempt$wtwLPSD
        reRate_[c1, c2] = tempt$reRate 
      }
    }
      
  }
}

junk = data.frame(rbind(auc_[,,1], auc_[,,2])) 
names(junk) = 1:nCut
junk %>% as_tibble() %>% 
  mutate("condition" = rep(conditions, each = nCut),
         "p1" = rep(1:nCut, 2))  %>%
  gather(key = p2, value = auc, -condition, -p1) %>%
  ggplot(aes(p1, auc, color = p2)) + geom_point() + facet_grid(~condition)
