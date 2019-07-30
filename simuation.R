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
modelName = "RL2"
cb = c("LP", "HP")
nSim = 5
nCond = length(conditions)
# initialize 
auc_ = array(NA, dim = c(nCut, nPara, nCond))
wtw_ = array(NA, dim = c(length(tGrid), nCut, nPara, nCond))
aucSD_ = array(NA, dim = c(nCut, nPara, nCond))
wtwSD_ = array(NA, dim = c(length(tGrid), nCut, nPara, nCond))
reRate_ = array(NA, dim = c(nCut, nPara))
for(pIdx in 1 : nPara){
  for(cutIdx in 1 : nCut){
    paras = medianParas
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
names(junk) = paraNames
junk %>% as_tibble() %>% 
  mutate("condition" = rep(conditions, each = nCut),
         "percentile" = rep(1:nCut, 2)) %>%
  gather(key = para, value = auc, -condition, -percentile) %>%
  mutate(para = factor(para, levels = paraNames)) %>%
  ggplot(aes(percentile, auc)) + geom_point() + facet_grid(condition~para)


# reorganizd the data 
junk = data.frame(rbind(wtw_[,,6,1], wtw_[,,6,2])) 
names(junk) = 1:nCut
junk %>% as_tibble() %>%
  gather(key = cut, value = wtw) %>%
  mutate(t = rep(1 : (length(tGrid) * 2), nCut)) %>% 
ggplot(aes(t, wtw)) + geom_point() + facet_grid(~cut)

a[a$cut == 3,] %>% ggplot(aes(t, wtw)) + geom_point() 

junk = data.frame(reRate_) 
names(junk) = paraNames
junk %>% as_tibble() %>% 
  mutate("percentile" = 1:nCut) %>%
  gather(key = para, value = auc, -percentile) %>%
  mutate(para = factor(para, levels = paraNames)) %>%
  ggplot(aes(percentile, auc)) + geom_point() + facet_grid(~para)


cb = c("LP", "HP")
simFun = getSimFun("RL2")
set.seed(123)
tempt1 =  simFun(medianParas, cb)
set.seed(123)
tempt2 = simFun(medianParas, cb)
t = 170
data.frame(Qwait = c(tempt1$Qwaits[,t], tempt2$Qwaits[,t]),
           Qquit = c(rep(c(tempt1$Qwaits[1,t]) - tempt1$reRates[t]*3, 32),
                     rep(c(tempt2$Qwaits[1,t]) - tempt2$reRates[t]*3, 32)),
           t = rep(1:32, 2), betaP = rep(c("L", "H"), each = 32)) %>%
  gather(key = action, value = value, -t, - betaP) %>%
  ggplot(aes(t, value, color = action)) + geom_point() +
  facet_grid(~betaP)

mean(tempt1$reRates[(length(tempt1$reRates)-20) : length(tempt1$reRates)])
mean(tempt2$reRates[(length(tempt2$reRates)-20) : length(tempt2$reRates)])
