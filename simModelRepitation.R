# libraries and scripts
library("stringr")
library("ggplot2")
library("dplyr")
library("tidyr")
source("subFxs/helpFxs.R")
source("subFxs/loadFxs.R")
source("subFxs/modelComparisonFxs.R")
source("subFxs/plotThemes.R")
load("wtwSettings.RData")
# load model names
allData = loadAllData()
hdrData = allData$hdrData           
trialData = allData$trialData       
allIDs = hdrData$ID                   # column of subject IDs
n = length(allIDs) 
load("genData/expDataAnalysis/blockData.RData")
# select common useID
idList = hdrData$ID


#
modelName = "RL2"
paraNames = getParaNames(modelName)
nPara = length(paraNames)
expPara = loadExpPara(paraNames, sprintf("genData/expModelFitting/%sdb", "RL2"))
simPara = loadExpPara(paraNames, sprintf("genData/simModelFitting/%sdb", "RL2"))
expID = getUseID(expPara, paraNames)
simID = getUseID(simPara, paraNames)
useID = idList[idList %in% expID & idList %in% simID]


tempt  = cbind(exp = as.vector(as.matrix(expPara[,1:nPara])), sim = as.vector(as.matrix(simPara[,1 : nPara])))
data.frame(tempt, para = factor(rep(paraNames, each = nrow(simPara)), levels = paraNames),
           id = rep(simPara$id,  nPara)) %>% filter(id %in% useID) %>% 
  ggplot(aes(exp, sim)) + geom_point() + facet_wrap(.~para, scales = "free") +
  geom_abline(aes(slope = 1, intercept = 0), color = "red", linetype = "dashed") + 
  myTheme


tempt  = cbind(exp = as.vector(as.matrix(expPara[,1:nPara])), sim = as.vector(as.matrix(simPara[,1 : nPara])))
data.frame(tempt, para = factor(rep(paraNames, each = nrow(simPara)), levels = paraNames),
           id = rep(simPara$id,  nPara)) %>% filter(id %in% useID) %>% 
  ggplot(aes(exp, sim)) + geom_point() + facet_wrap(.~para, scales = "free") +
  geom_abline(aes(slope = 1, intercept = 0), color = "red", linetype = "dashed") + 
  myTheme
dir.create("figures/paraRecovery")
ggsave("figures/paraRecovery/RL2.png", width = 6, height = 4)
