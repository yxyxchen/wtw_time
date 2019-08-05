# in this dataset, only trials within the 7 mins will be kept. Therefore, we don't need to delete any data
# determine whether use truncated data
isTrun = T
# load libraries
source('subFxs/loadFxs.R') # for loading data 
source('subFxs/analysisFxs.R') # for analysis 
source("subFxs/plotThemes.R")
library("ggplot2")
library('dplyr')
dir.create("genData")
dir.create("genData/expDataAnalysis")

# load setting parameters 
load("wtwSettings.RData")

# load all data
allData = loadAllData()
hdrData = allData$hdrData           
trialData = allData$trialData       
allIDs = hdrData$ID                   # column of subject IDs
n = length(allIDs)                    # n
cat('Analyzing data for',n,'subjects.\n')

# define nBlock
nBlock = 2

# control which individual-level plots to generate
plotTrialwiseData = F
plotKMSC = F
plotWTW = F

# parameter for longtermR
window = 2 * 60
stepLen = 2 * 60
nWindow = (blockSecs - window) / stepLen + 1


# initialize outputs, organised by block
AUC = numeric(length =n * nBlock)
totalEarnings =  numeric(length =n * nBlock)
nExclude =  numeric(length =n * nBlock)
nAction = numeric(length =n * nBlock)
wtwEarly = numeric(length =n * nBlock)
timeWTW_ = vector(mode = "list", length = n * nBlock)
trialWTW_ = vector(mode = "list", length = n * nBlock)
kmOnGrid_ = vector(mode = "list", length = n * nBlock)
trialEndTime_ = vector(mode = "list", length = n * nBlock)
trialReRate_ = vector(mode = "list", length = n * nBlock)
longtermR_ = matrix(NA, nrow = nWindow, ncol = n * nBlock)
shorttermR_ = matrix(NA, nrow = nWindow, ncol = n * nBlock)
nQuit = numeric(length =n * nBlock)
nTrial = numeric(length =n * nBlock)
stdWd = numeric(length =n * nBlock)

# descriptive statistics for individual subjects and blocks
for (sIdx in 1 : n) {
  thisID = allIDs[sIdx]
  #if(blockData[blockData$id == thisID, "AUC"] > 20 & blockData$condition[blockData$id == thisID] == "LP"){
  for (bkIdx in 1: nBlock){
    noIdx = (sIdx - 1) * nBlock + bkIdx # 
    # select data 
    thisTrialData = trialData[[thisID]]
    thisBlockIdx = (thisTrialData$blockNum == bkIdx)
    thisTrialData = thisTrialData[thisBlockIdx,]
    cond = unique(thisTrialData$condition)
    cIdx = ifelse(cond == "HP", 1, 2)
    # truncate the last min(tMaxs) seconds, since here thisTrialData is not a data.frame
    if(isTrun){
      excluedTrials = which(thisTrialData$trialStartTime > (blockSecs - tMaxs[cIdx]))
      nExclude[[noIdx]] = length(excluedTrials)
      thisTrialData = thisTrialData[! (1 : length(thisTrialData$blockNum)) %in% excluedTrials, ]
    }
    # generate arguments for later analysis 
    label = sprintf('Subject %s, Cond %s, %s',thisID, unique(thisTrialData$condition), hdrData$stress[sIdx])
    tMax = min(tMaxs)
    
    # calcualte totalEarnings
    totalEarnings[noIdx] =  sum(thisTrialData$trialEarnings)
    timeWaited = thisTrialData$timeWaited
    trialEarnings = thisTrialData$trialEarnings
    scheduledWait = thisTrialData$scheduledWait
    timeWaited[trialEarnings > loseValue] = scheduledWait[trialEarnings > loseValue]
    nAction[noIdx] = sum(round(ifelse(trialEarnings > loseValue, ceiling(timeWaited / stepDuration), floor(timeWaited / stepDuration) + 1)))
    nTrial[noIdx] = length(timeWaited)
   
     # calculate varQuitTime
    nQuit[noIdx] = sum(trialEarnings == 0)
      
    # plot trial-by-trial data
    if (plotTrialwiseData) {
      trialPlots(thisTrialData,label)
      readline(prompt = paste('subject',thisID, "block", bkIdx, '(hit ENTER to continue)'))
      graphics.off()
    }
    
    # survival analysis
    kmscResults = kmsc(thisTrialData,min(tMaxs),label,plotKMSC, kmGrid)
    AUC[noIdx] = kmscResults[['auc']]
    kmOnGrid_[[noIdx]] = kmscResults$kmOnGrid
    stdWd[noIdx] = kmscResults$stdWd
    
    if (plotKMSC) {
      readline(prompt = paste('subject',thisID, "block", bkIdx, '(hit ENTER to continue)'))
      graphics.off()
    }

    # WTW time series
    wtwCeiling = min(tMaxs)
    wtwtsResults = wtwTS(thisTrialData, tGrid, wtwCeiling, label, plotWTW)
    timeWTW_[[noIdx]] = wtwtsResults$timeWTW
    trialWTW_[[noIdx]] = wtwtsResults$trialWTW
    wtwEarly[noIdx] =   wtwtsResults$trialWTW[1]
    # wait for input before continuing, if individual plots were requested
    if (plotWTW) {
      readline(prompt = paste('subject',thisID, "block", bkIdx, '(hit ENTER to continue)'))
      graphics.off()
    }
    
    # # calculate rewardRates
    # trialReRate_[[noIdx]] = trialEarnings / (timeWaited + iti)
    # trialEndTime_[[noIdx]] = thisTrialData$sellTime
    # 
    # # calculate longtermR
    # longtermR =     sapply(1 : nWindow, function(i) {
    #   startTime = (i-1) * stepLen
    #   endTime = (i-1) * stepLen + window
    #   junk = which(thisTrialData$trialStartTime >= startTime)
    #   if(length(junk) == 0){
    #     NA
    #   }else{
    #     startIdx = min(junk)
    #     endIdx = max(which(thisTrialData$sellTime < endTime))
    #     realStartTime = thisTrialData$trialStartTime[startIdx]
    #     realEndTime = thisTrialData$sellTime[endIdx]
    #     sum(thisTrialData$trialEarnings[startIdx : endIdx]) / (realEndTime - realStartTime)        
    #   }
    # }) 
    # longtermR_[,noIdx] = longtermR
    # 
    # # calculate shorttermR
    # # here timeWaited not includes reaction time
    # # not includes iti
    # shorttermR =   sapply(1 : nWindow, function(i) {
    #   startTime = (i-1) * stepLen
    #   endTime = (i-1) * stepLen + window
    #   junk = which(thisTrialData$trialStartTime >= startTime)
    #   if(length(junk) == 0){
    #     NA
    #   }else{
    #     startIdx = min(which(thisTrialData$trialStartTime >= startTime))
    #     endIdx = max(which(thisTrialData$sellTime < endTime))
    #     longtermR[i] = mean(thisTrialData$trialEarnings[startIdx : endIdx] / timeWaited[startIdx : endIdx])
    #   }
    # }) 
    # shorttermR_[,noIdx] = shorttermR
    
  } # loop over blocks
}
blockData = data.frame(id = rep(allIDs, each = nBlock), blockNum = rep( t(1 : nBlock), n),
                       condition = factor(rep(c("LP", "HP"), n), levels = c("HP", "LP")),
                       AUC = AUC, wtwEarly = wtwEarly,totalEarnings = totalEarnings,
                       nAction = nAction, nQuit = nQuit, nTrial = nTrial, stdWd = stdWd,
                       nExclude = nExclude)
# lastEndTime = sapply(1 : (nBlock * n), function(i) max(trialEndTime_[[i]]))
# hist(lastEndTime)
# range(lastEndTime)

# save data
save(kmOnGrid_, file = 'genData/expDataAnalysis/kmOnGridBlock.RData')
save(blockData, file = 'genData/expDataAnalysis/blockData.RData')

# descriptive statistics for individual subjects and blocks
# for (sIdx in 1 : n) {
#   thisID = allIDs[sIdx]
#   thisTrialData = trialData[[thisID]]
#   label = sprintf('Subject %s',thisID)
#   trialPlots(block2session(thisTrialData),label)
#   readline(prompt = paste('subject',thisID, '(hit ENTER to continue)'))
#   graphics.off()
# }

# plot AUC in two conditions
library("ggpubr")
blockData %>% ggplot(aes(condition, AUC)) + geom_boxplot() +
  geom_dotplot(binaxis='y', stackdir='center', aes(fill = condition)) +
  scale_fill_manual(values = conditionColors) + 
  xlab("") + ylab("Average WTW(s)") + myTheme +
  stat_compare_means(comparisons = list(c("HP", "LP")),
                     aes(label = ..p.signif..), label.x = 1.5, symnum.args= symnum.args,
                     bracket.size = 1, size = 6) + ylim(c(0, 20))
dir.create("figures")
dir.create("figures/expDataAnalysis")
ggsave(sprintf("figures/expDataAnalysis/zTruc_AUC.png"), width = 4, height = 3)

# plot wtw 
plotData = data.frame(wtw = unlist(timeWTW_), time = rep(tGrid, n),
           condition = rep(blockData$condition, each = length(tGrid))) %>% group_by(condition, time) %>%
  summarise(mean = mean(wtw), se = sd(wtw) / sqrt(length(wtw)), min = mean - se, max = mean + se) 

policy = data.frame(condition = c("HP", "LP"), wt = unlist(optimWaitTimes))
plotData %>% ggplot(aes(time, mean, color = condition,  fill = condition)) +
  geom_ribbon(aes(ymin=min, ymax=max), colour=NA, alpha = 0.3) +
  geom_line(size = 1) + facet_wrap(~condition, scales = "free") +
  scale_color_manual(values = conditionColors) + scale_fill_manual(values = conditionColors) +
  xlab("Cumulative task time (min)") +
  scale_x_continuous(breaks = seq(0, max(tGrid), by = 60 * 3),
                     label = seq(0, max(tGrid), by = 60 * 3) / 60) + 
  ylab("Willingness to wait (s)") +
  geom_hline(data = policy, aes(yintercept = wt, color = condition), linetype = "dotted", size = 1.5) +
  myTheme + ylim(c(0, 18))
ggsave("figures/expDataAnalysis/zTruc_wtw_timecourse.png", width = 6, height = 3)

# plot survival curve
data.frame(kmsc = unlist(kmOnGrid_), time = rep(kmGrid, n * nBlock),
                      condition = factor(rep(blockData$condition, each = length(kmGrid))), levels = conditions) %>%
  group_by(condition, time) %>%
  summarise(mean = mean(kmsc), se = sd(kmsc) / sqrt(length(kmsc)), min = mean - se, max = mean + se) %>% 
  ggplot(aes(time, mean, color = condition, fill = condition)) + 
  geom_ribbon(aes(ymin=min, ymax=max), colour=NA, alpha = 0.3)+
  geom_line(size = 1.5) + myTheme + scale_fill_manual(values = conditionColors) + 
  xlab("Elapsed time (s)") + ylab("Survival rate") + scale_color_manual(values = conditionColors)
ggsave("figures/expDataAnalysis/zTruc_kmsc_timecourse.png", width = 5, height = 4) 



# plot longtermR
data.frame(longtermR = as.vector(longtermR_), condition = rep(rep(c("LP", "HP"), each = nWindow), n), window = rep(1 : nWindow, n * nBlock)) %>%
  group_by(condition, window) %>% 
  summarise(muData = mean(longtermR, na.rm = T), seData = sd(longtermR, na.rm = T) / sqrt(sum(!is.na(longtermR))),
            min = muData - seData, max = muData + seData) %>%
  ggplot(aes(window, muData, color = condition, fill = condition)) + geom_line(size = 2) + 
  geom_ribbon(aes(ymin=min, ymax=max), colour=NA, alpha = 0.3) + 
  scale_color_manual(values = conditionColors) + scale_fill_manual(values = conditionColors)+
  myTheme + xlab("Time (s)") + ylab("E(R) / E(T)") +
    scale_x_continuous(labels = seq(1, 10, by = 2))
ggsave("figures/expDataAnalysis/zTruc_longtermR.png", width = 5, height = 4) 

data.frame(shorttermR= as.vector(shorttermR_), condition = rep(rep(c("LP", "HP"), each = nWindow), n), window = rep(1 : nWindow, n * nBlock)) %>%
  group_by(condition, window) %>% summarise(muData = mean(shorttermR, na.rm = T),
                                            seData = sd(shorttermR, na.rm = T) / sqrt(sum(!is.na(shorttermR))),
                                            min = muData - seData, max = muData + seData) %>%
  ggplot(aes(window, muData, color = condition, fill = condition)) + geom_line(size = 2) + 
  geom_ribbon(aes(ymin=min, ymax=max), colour=NA, alpha = 0.3) + 
  scale_color_manual(values = conditionColors) + scale_fill_manual(values = conditionColors)+
  myTheme + xlab("Time (s)")  + ylab("E(R / T)") + scale_x_continuous(labels = seq(1, 10, by = 2))
ggsave("figures/expDataAnalysis/zTruc_shorttermR.png", width = 5, height = 4) 

  
# plot LP AUC against HP AUC, to see adapation and correlation
data.frame(HPAUC = blockData[blockData$condition == "HP", "AUC"],
           LPAUC = blockData[blockData$condition == "LP", "AUC"]) %>% 
  ggplot(aes(LPAUC, HPAUC)) + geom_point(size = 2) +
  geom_abline(intercept = 0, slope = 1) +
  xlim(c(0, min(tMaxs)))+ ylim(c(0, min(tMaxs))) +
  xlab("LP AUC (s)") + ylab("HP AUC (s)") +
  myTheme
ggsave("figures/expDataAnalysis/zTruc_AUC_Cmp.png", width = 3.5, height = 3.5) 
                 