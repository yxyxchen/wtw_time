library("ggplot2")
library("dplyr")
library("tidyr")
library("Hmisc")
library("coin")
source("subFxs/plotThemes.R")
source("subFxs/loadFxs.R") # load blockData and expPara
source("subFxs/helpFxs.R") # getParaNames
source("subFxs/analysisFxs.R") # plotCorrelation and getCorrelation
load("wtwSettings.RData")
# load trialData since we need scheduledWait 
allData = loadAllData()
hdrData = allData$hdrData 
trialData = allData$trialData       
allIDs = hdrData$ID 

modelName = "BL"

# create output directories
dir.create("figures/expParaAnalysis")
saveDir = sprintf("figures/expParaAnalysis/%s", modelName)
dir.create(saveDir)

# load expPara
paraNames = getParaNames(modelName)
nPara = length(paraNames)
parentDir = "genData/expModelFitting"
dirName = sprintf("%s/%sdb",parentDir, modelName)
expPara= loadExpPara(paraNames, dirName)
useID = getUseID(expPara, paraNames)


# plot hist
expPara %>% filter(id %in% useID) %>% select(c(paraNames)) %>%
  gather(key = "para", value = "value") %>%
  mutate(para = factor(para, levels = paraNames, labels = paraNames ))%>%
  ggplot(aes(value)) + geom_histogram(bins = 8) +
  facet_grid(~ para, scales = "free", labeller = label_parsed) + 
  myTheme + xlab(" ") + ylab(" ")
fileName = sprintf("%s/%s/hist.pdf", "figures/expParaAnalysis", modelName)
ggsave(fileName, width = 8, height = 3)

expPara %>% filter(id %in% useID) %>% select(c(paraNames)) %>%
  gather(key = "para", value = "value") %>%
  mutate(para = factor(para, levels = paraNames, labels = paraNames ))%>% 
  group_by(para) %>% summarise(mu = mean(value), median = median(value))


# explore relationships between para and blockData
load("genData/expDataAnalysis/blockData.RData")
expParaRL2 = loadExpPara(getParaNames("RL2"), "genData/expModelFitting/RL2db")
expParaQL2 = loadExpPara(getParaNames('QL2'), 'genData/expModelFitting/QL2db')
expParaPlus = data.frame(AUCHP = blockData$AUC[blockData$condition == "HP"],
                   AUCLP = blockData$AUC[blockData$condition == "LP"],
                   id = blockData$id[blockData$condition == "LP"],
                   gamma = expParaQL2$gamma)
expParaPlus$LL_all_RQ = expParaRL2$LL_all - expParaQL2$LL_all
expParaPlus$CV_RQ =logEvidence[,4] - logEvidence[,2] 
expParaPlus %>% filter(CV_RQ > - 1000) %>%
  ggplot(aes(gamma, CV_RQ)) + geom_point()


