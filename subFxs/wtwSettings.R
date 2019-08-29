# in this script, we try to calculate the optimWaitingTimes and the optimRewardRates
# to be as close as to the normal analysis (like integrate the prob density)
# we assume rewards happen at the middle of the gap(therefore, the meanRewardDelay would be unbiased)
# yet in wtwSettingsEnd.R, to unfiy different algorithms, we assume rewards happen at the end of the gap
# however, results for LP still change with the stepDuration
# we do all the calculation by stepDuration = 0.5, and the optimWaitTime, you know is not that...

# we use this script to get stepDuration
# we don't use the reward rate here, it is close to the normal analysis, but not that good.
######## condition varibles #########
conditions = c("HP", "LP")
conditionColors = c("#008837", "#7b3294")

######## timing variables ########
tMaxs = c(16, 32) # trial durations
blockMins = 10 # block duration in mins
blockSecs = blockMins * 60 # block duration in secs
iti = 2 # iti duration in secs
tGrid = seq(0, blockSecs, 1)
kmGrid = seq(0, min(tMaxs), 0.1)

######### reward variable ########
tokenValue = 10 #value of the token
loseValue = 0
stepDuration = 1
########## supporting vairbales ########
junk = read.csv("data/scheduledWaitHP.csv")
seqHP = sort(as.double(unique(junk[[1]])))
junk = read.csv("data/scheduledWaitLP.csv")
seqLP = sort(as.double(unique(junk[[1]])))
nseq = length(seqLP)

optimWaitTimes = list(HP = 16, LP = seqLP[4])
HP = tokenValue / ((2 + 16) / 2 + iti)
LP = 0.5 * tokenValue / ((mean(seqLP[1:4]) + seqLP[4])/2 + iti)
optimRewardRates = c("HP" = HP, "LP" = LP) 
save("conditions", "conditionColors", "tMaxs", "blockMins", "blockSecs", "iti", "tGrid", 
     "tokenValue", "stepDuration", "optimRewardRates", 
     "optimWaitTimes", "loseValue", "kmGrid", file = "wtwSettings.RData")

# plot cdf
library('ggplot2')
source('subFxs/plotThemes.R')
library("tidyr"); library('dplyr')
dir.create('figures/exp')
data.frame(CDP = c(seq(0, 1, by = 1/8), 1, seq(0, 1, by = 1/8) ,1),
           index = c(0, seqHP, tMaxs[1], 0, seqLP, tMaxs[2]),
           cond =  rep(c('HP', 'LP'), c(length(seqHP) + 2, length(seqLP) + 2))) %>%
  ggplot(aes(index, CDP)) + geom_step(size = 2) + facet_grid(~cond) +
  ylim(c(0,1)) + 
  myTheme + xlab('Delay duration (s)') + ylab('CDF')
ggsave('figures/exp/cdp.png', width =6, height = 3)

# given seqHP we need smaller stepDuration
stepDuration = 0.1
seqHP = round(seqHP, 2)
seqLP = round(seqLP, 2)
# plot reward rates approximately
trialTicks = list(
  'HP' = round(seq(0, tMaxs[1], by = stepDuration), 1),
  'LP' = round(seq(0, tMaxs[2], by = stepDuration), 1)
)

trialGapValues = list(
  "HP" = round(seq(stepDuration, tMaxs[1], by = stepDuration), 2),
  "LP" = round(seq(stepDuration, tMaxs[2], by = stepDuration), 2)
)

trialGapIdxs = list(
  "HP" = 1 : length(trialGapValues$HP),
  "LP" = 1 : length(trialGapValues$LP)
)

HP = rep(0, length(trialGapIdxs$HP))
HP[sapply(1 :8, function(i) which.min(abs(trialGapValues$HP - seqHP[i])))] = 1 / 8
LP = rep(0, length(trialGapIdxs$LP))
LP[sapply(1 :8, function(i) which.min(abs(trialGapValues$LP - seqLP[i])))] = 1 / 8
rewardDelayPDF=list('HP'= HP, 'LP' = LP)
rewardDelayCDF = list('HP' = cumsum(rewardDelayPDF$HP),
                      'LP' = cumsum(rewardDelayPDF$LP))

# assume the rewards happen at the end of the gap
# E(t_reward | t_reward <= T) 
HP = cumsum((trialGapValues$HP - 0.5 * stepDuration) * rewardDelayPDF$HP) / cumsum(rewardDelayPDF$HP)
LP = cumsum((trialGapValues$LP - 0.5 * stepDuration) * rewardDelayPDF$LP) / cumsum(rewardDelayPDF$LP)
# no reward arrives before the first reward timing, so points before that turn to NAN
meanRewardDelay = list('HP' = HP, 'LP' = LP)


HP = tokenValue * rewardDelayCDF$HP /
  ((meanRewardDelay$HP * rewardDelayCDF$HP + trialGapValues$HP * (1 - rewardDelayCDF$HP)) + iti)
LP = tokenValue * rewardDelayCDF$LP /
  ((meanRewardDelay$LP * rewardDelayCDF$LP + trialGapValues$LP * (1 - rewardDelayCDF$LP)) + iti)
rewardRate = list('HP' = HP, 'LP' = LP)

policy = data.frame(cond = c("HP", "LP"), rewardRate = as.double(optimWaitTimes))
data.frame(rewardRate = c(0, rewardRate[[1]], 0, rewardRate[[2]]),
           time = c(trialTicks[[1]], trialTicks[[2]]),
           cond = rep(c("HP", "LP"), time = (tMaxs / stepDuration) + 1)) %>%
  ggplot(aes(time, rewardRate)) +
  geom_line(size = 3)  + myTheme + 
  ylab(expression(bold("Reward rate (cent s"^"-1"*")"))) + xlab("Waiting policy (s)")  +
  geom_vline(data = policy, aes(xintercept = rewardRate),
             linetype = "dashed", size = 1.5) + facet_grid(~cond)
ggsave("figures/exp/reward_rate.png", width = 6, height = 3)
