#load libraries
library('plyr'); library(dplyr); library(ggplot2);library('tidyr');
library("stringr")
library("loo")
library("coda") 
source('subFxs/modelFittingFxs.R') # for fitting each single participant
source('subFxs/loadFxs.R') # for load data
source("subFxs/helpFxs.R") # for getparaNames
load("wtwSettings.RData")
source("subFxs/analysisFxs.R")
decodes = c("RL1", "BL", "QL1", "QL2")
encodes = c("RL2", "QL2")

for(e in 1 : length(encodes)){
  encode = encodes[e]
  # load simData
  load(sprintf("genData/simulation/%s.RData", encode))
  ids = hdrData$ID
  nSub = length(ids)    
  for(d in 1 : length(decodes)){
    decode = decodes[d]
    # determine paraNames
    paraNames = getParaNames(decode)
    nPara = length(paraNames)
    if(paraNames == "wrong model name"){
      print(paraNames)
      break
    }
    # determine excID
    expPara = loadExpPara(paraNames,
                          sprintf("genData/simModelFitting/%s/%sdb", encode, decode))
    useID = getUseID(expPara, paraNames)
    excID = ids[!ids %in% useID]
    txt = sprintf("%s_%s:%s", encode, decode, length(useID))
    print(txt)
    readline("continue")
  }
}