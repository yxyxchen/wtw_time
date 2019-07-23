
# loadData.R
# varying-magnitude WTW

# each subject has:
#  header file (*hdr.txt)
#  main data file (ending both in .mat and .txt)
#  empty *keytimes.txt file, a vestige of the effort key-pressing version. 


loadAllData = function() {
  library("gtools")
  # define column names 
  colNames = c('blockNum', 'trialNum', 'trialStartTime', 'nKeyPresses', 'scheduledWait',
    'rewardTime', 'timeWaited', 'sellTime', 'trialEarnings','totalEarnings')
  
  # get fileNames
  dataDir = "data"
  fileNames = list.files(path= dataDir, pattern=('wtw-timing-fixed_[0-9]{3}_[0-9].txt'))
  fileNames = mixedsort(sort(fileNames))
  nFile = length(fileNames)
  if(any(duplicated(fileNames))){
    print("duplicated files!")
    break
  }else{
    sprintf("load %d files", nFile) 
  }
  
  # load data
  hdrData = matrix(ncol = 1, nrow = nFile)
  trialData = list()
  for(i in 1 : nFile){
    fileName = fileNames[i]
    id = substr(fileName, 18,20)
    junk = read.csv(sprintf("%s/%s", dataDir, fileName), col.names = colNames, header = F)
    junk$condition = ifelse(junk$blockNum == 1, "LP", "HP")
    hdrData[i,1] = id
    trialData[[id]] = junk
  }
  
  
 
  colnames(hdrData) = "ID"
  hdrData = as.data.frame(hdrData)
  outputData = list(hdrData=hdrData, trialData=trialData)
  return(outputData)
  
} 


loadExpPara = function(paraNames, dirName){
  # number of paraNames 
  nE = length(paraNames) + 1
  # number of files
  fileNames = list.files(path= dirName, pattern=("*_summary.txt"))
  library("gtools")
  fileNames = mixedsort(sort(fileNames))
  n = length(fileNames) 
  sprintf("load %d files", n)
  
  # initialize the outout variable 
  expPara = matrix(NA, n, nE * 4)
  idList = vector(length = n)
  # loop over files
  for(i in 1 : n){
    fileName = fileNames[[i]]
    address = sprintf("%s/%s", dirName, fileName)
    junk = read.csv(address, header = F)
    idIndexs = str_locate(fileName, "s[0-9]+")
    idList[i] = substr(fileName, idIndexs[1]+1, idIndexs[2])
    # delete the lp__ in the old version
    if(nrow(junk) > nE){
      junk = junk[1:nE,]
    }
    expPara[i, 1:nE] = junk[,1]
    expPara[i, (nE + 1) : (2 * nE)] = junk[,3]
    expPara[i, (2*nE + 1) : (3 * nE)] = junk[,9]
    expPara[i, (3 * nE + 1) : (4 * nE)] = junk[,10]
  }
  # transfer expPara to data.frame
  expPara = data.frame(expPara)
  junk = c(paraNames, "LL_all")
  colnames(expPara) = c(junk, paste0(junk, "SD"), paste0(junk, "Effe"), paste0(junk, "Rhat"))
  expPara$id = factor(idList, levels = levels(hdrData$ID)) # ensure the levels are consistent, usually not that problematic though
  return(expPara)
}

# I also need to load 2.5% and 97.5%
loadCVPara = function(paraNames, dirName, pattern){
  # number of paraNames 
  nE = length(paraNames) + 1
  # number of files
  fileNames = list.files(path= dirName, pattern= pattern)
  library("gtools")
  fileNames = mixedsort(sort(fileNames))
  n = length(fileNames) 
  sprintf("load %d files", n)
  
  # initialize the outout variable 
  expPara = matrix(NA, n, nE * 6)
  idList = vector(length = n)
  # loop over files
  for(i in 1 : n){
    fileName = fileNames[[i]]
    address = sprintf("%s/%s", dirName, fileName)
    junk = read.csv(address, header = F)
    sIndexs = str_locate(fileName, "s[0-9]+")
    s = substr(fileName, sIndexs[1]+1, sIndexs[2])
    fIndexs = str_locate(fileName, "f[0-9]+")
    f = substr(fileName, fIndexs[1]+1, fIndexs[2])
    idList[i] = sprintf("s%s_f%s", s, f)
    # delete the lp__ in the old version
    if(nrow(junk) > nE){
      junk = junk[1:nE,]
    }
    expPara[i, 1:nE] = junk[,1]
    expPara[i, (nE + 1) : (2 * nE)] = junk[,3]
    expPara[i, (2*nE + 1) : (3 * nE)] = junk[,9]
    expPara[i, (3 * nE + 1) : (4 * nE)] = junk[,10]
    expPara[i, (4*nE + 1) : (5 * nE)] = junk[,4]
    expPara[i, (5 * nE + 1) : (6 * nE)] = junk[,8]
  }
  # transfer expPara to data.frame
  expPara = data.frame(expPara)
  junk = c(paraNames, "LL_all")
  colnames(expPara) = c(junk, paste0(junk, "SD"), paste0(junk, "Effe"), paste0(junk, "Rhat"),
                        paste0(junk, "2.5"),paste0(junk, "97.5"))
  expPara$id = factor(idList, levels = idList)
  return(expPara)
}

