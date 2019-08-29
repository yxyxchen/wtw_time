conditions = c("HP", "LP")
tMaxs = c(16, 32) # max trial durations in secs
nBlock = 2
blockMin = 10 # block duration in mins
blockSec = blockMin * 60 # block duration in secs
iti = 2 # iti duration in secs
tokenValue = 10 # value of the matured token
# possible values of reward dlays
rewardDelays = list(HP = seq(2, 16, by = 2),
                    LP = c(0.2759766, 0.7589358, 1.6041142, 3.0831765,
                             5.6715355, 10.2011638, 18.1280133, 32)) 
# optimal waiting thresholds, unit = sec
optimWaitThresholds = list(HP = 16, LP = 3.0831765)
# optimal reward rates
optimRewardRates = c("HP" = 10 / 11, "LP" = 1.174574) 
# analyses parameters
tGrid = seq(0, blockSec, by = 2) # time grid for wtw time courses
kmGrid = seq(0, min(tMaxs), by = 0.1) # time grid for Kaplan-Meier survival curves
save("conditions" = conditions,
     "tMaxs" = tMaxs,
     "blockMin" = blockMin,
     "blockSec" = blockSec,
     "nBlock" = nBlock,
     "iti" = iti,
     "tokenValue" = tokenValue,
     "rewardDelays" = rewardDelays,
     "optimRewardRates" = optimRewardRates,
     "optimWaitThresholds" = optimWaitThresholds,
     'tGrid' = tGrid,
     'kmGrid' = kmGrid,
     file = "expParas.RData")