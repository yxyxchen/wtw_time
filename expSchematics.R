# this script plots delay distributions and reward rates in two environments

# load experiment parameters 
load("expParas.RData")

seqHP = rewardDelays$HP
seqLP = rewardDelays$LP

## for display purposes, all variables on continous time
## are discretized into 0.1 secondition time bins
## which is lower than the time resoluation of reward delays in
## the experiment. 
bin = 0.1 # width of a time bin
time = list(
  HP = seq(bin, tMaxs[1], by = bin),
  LP = seq(bin, tMaxs[2], by = bin)
) 

## delays in both HP and LP are drawn repeatedly from 8 possible values,
## In HP, they are uniformly distributed from 0 to 16, 
## in LP, they are log-spaced from 0 to 32, therefore forming a heavy tail distribution

# delay PDFs
HP = rep(0, length(time$HP))
LP = rep(0, length(time$LP))
HP[sapply(1 :8, function(i) which.min(abs(time$HP - seqHP[i])))] = 1 / 8
LP[sapply(1 :8, function(i) which.min(abs(time$LP - seqLP[i])))] = 1 / 8
rewardDelayPDFs =list('HP'= HP, 'LP' = LP)

# delay CDFs
rewardDelayCDFs = list('HP' = cumsum(rewardDelayPDFs$HP),
                      'LP' = cumsum(rewardDelayPDFs$LP))

# average waiting durations given different policies
# Here we assume rewards occur at the middle of each time bin
HP = cumsum((time[['HP']] - 0.5 * bin) * rewardDelayPDFs$HP) / cumsum(rewardDelayPDFs$HP)
LP = cumsum((time[['LP']] - 0.5 * bin) * rewardDelayPDFs$LP) / cumsum(rewardDelayPDFs$LP)
meanRewardDelays = list('HP' = HP, 'LP' = LP)

# rewardRates given different policies
## noticably, since we lower the temporal resoluation, the reward rates 
## and the optimal behavior variables calculated here might be slightly
## different from those in expParas.R. These approximate values are only
## used for display purposes, and we use the exact values in all analysis scripts. 
HP = tokenValue * rewardDelayCDFs$HP /
  ((meanRewardDelays$HP * rewardDelayCDFs$HP + time[['HP']] * (1 - rewardDelayCDFs$HP)) + iti)
LP = tokenValue * rewardDelayCDFs$LP /
  ((meanRewardDelays$LP * rewardDelayCDFs$LP + time[['LP']] * (1 - rewardDelayCDFs$LP)) + iti)
rewardRates = list('HP' = HP, 'LP' = LP)

# optimal raward rates and optimal policies
optimWaitThresholds = list()
optimWaitThresholds$HP = time$HP[which.max(HP)]
optimWaitThresholds$LP = time$LP[which.max(LP)]
optimRewardRates = list()
optimRewardRates$HP = max(HP)
optimRewardRates$LP = max(LP)


# plot CDFs 
library('ggplot2')
source('subFxs/plotThemes.R')
library("tidyr"); library('dplyr')
dir.create('figures/expSchematics')
data.frame(CDF = c(0,c(rewardDelayCDFs$HP, rep(1, length(time$LP) - length(time$HP))), 0, rewardDelayCDFs$LP),
           time = c(0, time$LP, 0, time$LP),
           condition =  rep(c('HP', 'LP'), c(length(time$LP) + 1, length(time$LP) + 1))) %>%
  ggplot(aes(time, CDF)) + geom_line(size = 3, color = themeColor) + facet_grid(~condition) +
  ylim(c(0,1)) + scale_y_continuous(breaks = c(0,0.5,1)) + 
  scale_x_continuous(breaks = c(0, max(tMaxs)/ 2, max(tMaxs)),limits = c(0, max(tMaxs) * 1.1)) +
  myTheme + xlab('Delay duration (s)') + ylab('CDF') + ggtitle(expName) +
  annotate("text", x = 16, y = 0.4, label =  "+10¢", size = 6) + 
  theme(plot.title = element_text(hjust = 0.5, color = themeColor))
ggsave('figures/expSchematics/CDF.eps', width =4, height = 3)
ggsave('figures/expSchematics/CDF.png', width =4, height = 3)


# plot reward rates
optimData = data.frame(condition = c("HP", "LP"), waitThreshold = as.double(optimWaitThresholds))
data.frame(rewardRate = c(0, rewardRates[[1]], 0, rewardRates[[2]]),
           time = c(0, time[[1]], 0, time[[2]]),
           condition = rep(c("HP", "LP"), c(length(time$HP) + 1, length(time$LP) + 1))) %>%
  ggplot(aes(time, rewardRate)) +
  geom_line(size = 3, color = themeColor)  + myTheme + 
  scale_x_continuous(breaks = c(0, max(tMaxs)/ 2, max(tMaxs)),limits = c(0, max(tMaxs) * 1.1)) + 
  ylab(expression(bold("Reward rate (¢ s"^"-1"*")"))) + xlab("Waiting policy (s)")  + facet_grid(~condition) +
  scale_y_continuous(breaks = c(0, 0.4, 0.8, 1.2), limits = c(0, 1.3)) +
  ggtitle(expName) + theme(plot.title = element_text(hjust = 0.5, color = themeColor))
ggsave("figures/expSchematics/reward_rate.eps", width = 4, height = 3)
ggsave("figures/expSchematics/reward_rate.png", width = 4, height = 3)


