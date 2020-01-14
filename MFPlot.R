library('dplyr')
library("tidyr")
library("ggplot2")
library("ggpubr")
source("subFxs/plotThemes.R")
source('MFAnalysis.R')

# output dir
dir.create('figures/MFplot')

# load experiment parameters
load("expParas.RData")

# plot WTW timecourses in two environments
MFResults = MFAnalysis(isTrct = F)
sumStats = MFResults[['sumStats']]
timeWTW_ = MFResults[['timeWTW_']]
nSub = nrow(sumStats)
data.frame(wtw = unlist(timeWTW_),
           time = rep(tGrid, nSub * nBlock),
           condition = rep(rep(c("HP", "LP"), each = length(tGrid))), nSub) %>%
  group_by(condition, time) %>%
  dplyr::summarise(mu = mean(wtw, na.rm = T), se = sd(wtw, na.rm = T) / sqrt(sum(!is.na(wtw))),
                   min = mu- se, max = mu + se) %>%
  ggplot(aes(time, mu, linetype = condition)) +
  geom_ribbon(aes(ymin=min, ymax=max), fill = '#c7e9c0') +
  geom_line(color = themeColor, size = 1) +
  xlab("Elapsed time (s)") + ylab("WTW (s)") + 
  myTheme
ggsave("figures/MFPlot/wtw_timecourse.eps", width = 5, height = 4) 

# plot average WTWs in two environments
MFResults = MFAnalysis(isTrct = T)
sumStats = MFResults[['sumStats']]
wTest = wilcox.testc( sumStats[sumStats$condition == "HP", "muWTW"],
                     sumStats[sumStats$condition == "LP", "muWTW"],paired = T)
data.frame(muWTWHP = sumStats$muWTW[sumStats$condition == 'HP'],
           muWTWLP = sumStats$muWTW[sumStats$condition == 'LP']) %>%
  ggplot(aes(muWTWLP, muWTWHP)) +
  geom_point(color = themeColor, size = 5, shape = 21, fill = '#c7e9c0', stroke =1) +
  geom_abline(slope = 1, intercept = 0) + 
  annotate("text", x = 15, y = 3, label = sprintf('p < 0.001***', wTest$p.value)) +
  xlab("LP muAUC / (s)") + ylab("HP muAUC / (s)") + 
  myTheme + xlim(c(-1,17)) + ylim(c(-1,17)) 
ggsave("figures/MFPlot/muWTW_comparison.eps", width = 4, height = 3)

# plot an example participant 
AUC = sumStats$muWTW[1]
survCurve = MFResults$survCurve_[[1]]
data.frame(time = seq(0, tMaxs[1], by = 0.1), surveCurve = survCurve) %>% ggplot(aes(time, surveCurve)) +
  geom_step(size = 2) + myTheme + xlab("Elapsed time (s)") + ylab("Survival rate") +
  scale_x_continuous(limits = c(0, tMaxs[1]), breaks = c(0, tMaxs[1])) +
  annotate("text", x = 4, y = 0.6, label = sprintf("AUC = %.2f s", AUC), size = 5)
ggsave("figures/example_AUC.png", width = 4, height = 3)

data.frame(wtw = MFResults$timeWTW_[[1]], time = tGrid) %>% 
  ggplot(aes(time, wtw)) + geom_step(size = 1) + 
  xlab("Time (s)") + ylab("WTW (s)") + myTheme 
ggsave("figures/example_WTW.png", width = 4, height = 3)
