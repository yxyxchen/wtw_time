# in this script, we try to follow the wtwSettings analysis in other repos. 
# in this experiment, rewards arrive at discrete time points, Therefore, we actually can know 
# the exact optinal waiting policies and reward rates. 
# in this script, we just assume the temporal resolution is 0.1s. Therefore, our reward timings are rounded at 0.1s level
# It is not problematic in HP, since 2, 4, .. can be rounded exactly, however, it means rewards occur latter in LP, since we use ceiling.
# therefore, our reward rate function should be influenced in LP
# the distortion is regardless small, since the timestep is small here
# also, it is useful for some optimality analyses
######## condition varibles #########
conditions = c("HP", "LP")
conditionColors = c("#008837", "#7b3294")

######## timing variables ########
tMaxs = c(16, 32) # trial durations
blockMins = 10 # block duration in mins
blockSecs = blockMins * 60 # block duration in secs
iti = 2 # iti duration in secs
tGrid = seq(0, blockSecs, 0.1)

######### reward variable ########
tokenValue = 10 #value of the token
loseValue = 0
stepDuration = 0.01
########## supporting vairbales ########
# time ticks within a trial for timeEarnings or wtw analysis
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


########## additional  variables for optimal analysis ########
#  PDF of reward delays: p(t_reward = T)
junk = read.csv("data/scheduledWaitHP.csv")
seqHP = sort(as.double(unique(junk[[1]])))
junk = read.csv("data/scheduledWaitLP.csv")
seqLP = sort(as.double(unique(junk[[1]])))

lenSeq = length(seqHP)
HP = rep(0, length(trialGapValues$HP)); LP = rep(0, length(trialGapValues$LP))
HP[sapply(1 : lenSeq, function(i) min(which(trialGapValues$HP >= seqHP[[i]])))] = 1 / lenSeq
LP[sapply(1 : lenSeq, function(i) min(which(trialGapValues$LP >= seqLP[[i]])))] = 1 / lenSeq
rewardDelayPDF = list(
  "HP" = HP,
  "LP" = LP
)

# CDF of reward delays: p(t_reward <= T)
HP = cumsum(rewardDelayPDF$HP)
LP = cumsum(rewardDelayPDF$LP)
rewardDelayCDF = list(
  HP = HP,
  LP = LP
)

# assume the rewards happen at the middle of the gap
# E(t_reward | t_reward <= T) 
HP = sapply(1 : length(trialGapValues$HP), function(i){
  thisPDF = rewardDelayPDF$HP
  thisValues = trialGapValues$HP
  sum(thisPDF[1 : i] * thisValues[1 : i]) + (1 - sum(thisPDF[1 : i])) * i * stepDuration
})
LP = sapply(1 : length(trialGapValues$LP), function(i){
  thisPDF = rewardDelayPDF$LP
  thisValues = trialGapValues$LP
  sum(thisPDF[1 : i] * thisValues[1 : i]) + (1 - sum(thisPDF[1 : i])) * i * stepDuration
})
meanRewardDelay = list('HP' = HP, 'LP' = LP)

# rewardRate
HP = tokenValue * rewardDelayCDF$HP /
  ((meanRewardDelay$HP * rewardDelayCDF$HP + trialGapValues$HP * (1 - rewardDelayCDF$HP)) + iti)
LP = tokenValue * rewardDelayCDF$LP /
  ((meanRewardDelay$LP * rewardDelayCDF$LP + trialGapValues$LP * (1 - rewardDelayCDF$LP)) + iti)
rewardRate = list('HP' = HP, 'LP' = LP)

optimWaitTimes = list()
optimWaitTimes$HP = trialGapValues$HP[which.max(HP)]
optimWaitTimes$LP = trialGapValues$LP[which.max(LP)]

optimRewardRates = list()
optimRewardRates$HP = max(HP)
optimRewardRates$LP = max(LP)


# calculate the expected remaining time 
remainTime = rep(0, length = tMaxs[2] / stepDuration )
for(quitGap in 2 : (tMaxs[2] / stepDuration - 2)){
  select = (quitGap + 1) : length(rewardDelayCDF[[2]])
  remainTime[quitGap] = sum(rewardDelayPDF[[2]][select]  * (trialGapValues[[2]][select] - quitGap * stepDuration))
}
