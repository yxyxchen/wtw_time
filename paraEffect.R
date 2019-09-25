# this script is used to demonstrate the effect of 
library('ggplot2')
source("subFxs/plotThemes.R")
library('plyr')
library('dplyr')
library('tidyr')
load("expParas.RData")
source("subFxs/helpFxs.R") # getParas
source("subFxs/loadFxs.R") # load scheduledWait from empirical data
source("subFxs/analysisFxs.R") 


# loop over participants 
library("doMC")
library("foreach")
nCore = parallel::detectCores() -1 # only for the local computer
registerDoMC(nCore)

# load expData
allData = loadAllData()
hdrData = allData$hdrData           
trialData = allData$trialData       
ids = hdrData$ids        
nSub = length(ids)   

# load expPara
modelName = "QL2"
paraNames = getParaNames(modelName)
nPara = length(paraNames)
expPara = loadExpPara(paraNames, dirName = sprintf("genData/expModelFit/%s", modelName))

# determine paraSamples
nCut = 3
paraSamples = matrix(NA, nrow = nCut, ncol = (nPara))
for(i in 1 : nPara){
  qt2 = quantile(expPara[,i], 0.2)
  qt8 = quantile(expPara[,i], 0.8)
  paraSamples[,i] = seq(qt2, qt8, length.out = nCut)
}

paraSamples[,4] = paraSamples[,4]  - 22
# median paras
medianParas = paraSamples[2,]

# 

# simulate for a single condition
# determine simFun
modelName = "RL2"
blockDuration = 20
cb = c("HP")
nSim = 10
tGrid = seq(0, blockDuration * 60, by = 5)
# initialize 
auc_ = array(NA, dim = c(nCut, nPara))
wtw_ = array(NA, dim = c(length(tGrid) * length(cb), nCut, nPara))
aucSD_ = array(NA, dim = c(nCut, nPara))
wtwSD_ = array(NA, dim = c(length(tGrid) * length(cb), nCut, nPara))
reRate_ = array(NA, dim = c(blockDuration*60 / iti * length(cb), nCut, nPara))
auc1_ = array(NA, dim = c(nCut, nPara))
auc2_ = array(NA, dim = c(nCut, nPara))
auc3_ = array(NA, dim = c(nCut, nPara))
aucSD1_ = array(NA, dim = c(nCut, nPara))
aucSD2_ = array(NA, dim = c(nCut, nPara))
aucSD3_ = array(NA, dim = c(nCut, nPara))
for(pIdx in 1 : (nPara)){
  for(cutIdx in 1 : nCut){
    paras = medianParas
    paras[pIdx] = paraSamples[cutIdx,pIdx]
    tempt = simulateUnitSingle(paras, nSim, modelName, cb, blockDuration)
    auc_[cutIdx, pIdx] = tempt$auc
    aucSD_[cutIdx, pIdx] = tempt$aucSD
    wtw_[ , cutIdx, pIdx] = tempt$wtw
    wtwSD_[ , cutIdx, pIdx] = tempt$wtwSD
    auc1_[cutIdx, pIdx] = tempt$auc1
    auc2_[cutIdx, pIdx] = tempt$auc2
    auc3_[cutIdx, pIdx] = tempt$auc3
    aucSD1_[cutIdx, pIdx] = tempt$aucSD1
    aucSD2_[cutIdx, pIdx] = tempt$aucSD2
    aucSD3_[cutIdx, pIdx] = tempt$aucSD3    
    reRate_[,cutIdx, pIdx] = tempt$reRate 
  }
}

# reorganizd the data 
colnames(auc_) =  paraNames
colnames(aucSD_) =  paraNames
junk1 = auc_ %>% as_tibble() %>% 
  mutate("percentile" = rep(1:nCut, 1)) %>%
  gather(key = para, value = auc, -percentile) %>%
  mutate(para = factor(para, levels = paraNames))
junk2 =  aucSD_ %>% as_tibble() %>% 
  mutate("percentile" = rep(1:nCut, 1)) %>%
  gather(key = para, value = aucSD, -percentile) %>%
  mutate(para = factor(para, levels = paraNames))
data.frame(junk1, aucSD = junk2$aucSD) %>%
  mutate(min = auc - aucSD, max = auc + aucSD) %>%
  ggplot(aes(percentile, auc)) + geom_point() +
  geom_errorbar(aes(ymin = min, ymax = max)) + facet_grid(~para) +
  ylab("AUC (s)") + xlab("Value (a.u.)") + ggtitle(cb[1]) + 
  myTheme
dir.create("figures/simualtion")
ggsave(file = sprintf("figures/simualtion/auc_%s.png", cb[1]), width = 6, height = 4)

# I would like to know how much the results change
dimnames(wtw_) = list(tGrid + 1, 1:nCut, paraNames)
library("driver")
gather_array(wtw_, wtw, "t", "paraValue", "para") %>%
  mutate(paraValue = as.factor(paraValue),
         para = factor(para,labels = junk)) %>%
  ggplot(aes(t, wtw,  color = paraValue)) + geom_line() +
  facet_wrap(para~.)

# plot the change the wtw1 in the first, mindle and last 20 seconds
wtw1 = apply(wtw_, MARGIN = c(2,3), FUN = function(x) mean(x[1:60]))
wtw2 = apply(wtw_, MARGIN = c(2,3), FUN = function(x) mean(x[61:120]))
wtw3 = apply(wtw_, MARGIN = c(2,3), FUN = function(x) mean(x[121: 181]))
wtwSD1 = apply(wtwSD_, MARGIN = c(2,3), FUN = function(x) mean(x[1:60]))
wtwSD2 = apply(wtwSD_, MARGIN = c(2,3), FUN = function(x) mean(x[61:120]))
wtwSD3 = apply(wtwSD_, MARGIN = c(2,3), FUN = function(x) mean(x[121: 181]))
data.frame(wtw = c(as.vector(wtw1),as.vector(wtw2),as.vector(wtw3)),
           wtwSD = c(as.vector(wtwSD1), as.vector(wtwSD2), as.vector(wtwSD3)),
           paraValue = rep(rep(1 : nCut, nPara), 3),
           para = rep(rep(paraNames[1:5], each = nCut), 3),
           time = rep(1:3, each = nCut * (nPara))) %>% 
  mutate(para = factor(para, levels = paraNames[1:5]),
         paraValue = factor(paraValue),
         min = wtw - wtwSD, max = wtw + wtwSD) %>%
  ggplot(aes(time, wtw, color = paraValue)) + geom_line() +
  geom_point() + 
  facet_grid(~para) 


data.frame(auc = c(as.vector(auc1_),as.vector(auc2_),as.vector(auc3_)),
           aucSD = c(as.vector(aucSD1_), as.vector(aucSD2_), as.vector(aucSD3_)),
           paraValue = rep(rep(1 : nCut, nPara), 3),
           para = rep(rep(paraNames[1:5], each = nCut), 3),
           time = rep(1:3, each = nCut * (nPara))) %>% 
  mutate(para = factor(para, levels = paraNames[1:5]),
         paraValue = factor(paraValue),
         min = auc - aucSD, max = auc + aucSD) %>%
  ggplot(aes(time, auc, color = paraValue)) + geom_line() +
  geom_point() + 
  facet_grid(~para) 

# I want to know the results of reRates




# simulate for both conditions
# determine simFun
modelName = "RL2"
cb = c("HP")
nSim = 10
nCond = length(conditions)
# initialize 
auc_ = array(NA, dim = c(nCut, nPara, nCond))
wtw_ = array(NA, dim = c(length(tGrid), nCut, nPara, nCond))
aucSD_ = array(NA, dim = c(nCut, nPara, nCond))
wtwSD_ = array(NA, dim = c(length(tGrid), nCut, nPara, nCond))
reRate_ = array(NA, dim = c(nCut, nPara))
for(pIdx in 1 : (nPara)){
  for(cutIdx in 1 : nCut){
    paras = medianParas
    paras[pIdx] = paraSamples[cutIdx,pIdx]
    tempt = simulateUnitSingle(paras, nSim, modelName, cb, 20)
    auc_[cutIdx, pIdx, 1] = tempt$aucHP
    aucSD_[cutIdx, pIdx, 1] = tempt$aucHPSD
    wtw_[ , cutIdx, pIdx, 1] = tempt$wtwHP
    wtwSD_[ , cutIdx, pIdx, 1] = tempt$wtwHPSD
    
    auc_[cutIdx, pIdx, 2] = tempt$aucLP
    aucSD_[cutIdx, pIdx, 2] = tempt$aucLPSD
    wtw_[ , cutIdx, pIdx, 2] = tempt$wtwLP
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
for(p1 in 1 : (nPara)){
  for(p2 in (p1 + 2) : (nPara)){
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