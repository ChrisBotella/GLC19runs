# Make run of lof300_85_4000_bgnPc6_v140419 for GLC19

originalVes = c('etp','chbio_12','chbio_1','chbio_5','alti','slope','awc_top','bs_top','clc')
# Load allowed_classes.txt
allowed = read.csv("P:/Partage GLC19/allowed_classes.txt",sep=";",header=F)
allowed = as.character(allowed[,1])

# Load lof model
saveDir = paste("C:/Users/",user,"/pCloud local/0_These/data/personalWorker/files/",sep="")
setwd(saveDir)
occ = readRDS('occ_lof300_85_4000')
lof.mod = readRDS('model_lof300_85_4000_bgnPc6_v140419')
glc19SpIdVeco = readRDS('glc19SpId_List_lof300')

coefficients = lof.mod$beta[,dim(lof.mod$beta)]

spIds = intersect(allowed,glc19SpIdVeco)

nOccPerSp = sapply(spIds,function(id) sum(occ$glc19SpId==id))
hist(nOccPerSp,breaks="fd")

# Load testSet 
testSet = read.csv("P:/Partage GLC19/testSet.csv",sep=";",header=T)


predMatrix = predict.lof.spLogRelativeIntensity(testSet,
                                                coefficients,
                                                Intensityformula = as.formula(" ~ etp + I(etp^2) + I(chbio_12-etp) + I((chbio_12-etp)^2) + chbio_1 + I(chbio_1^2) + chbio_5 + I(chbio_5^2) + alti + I(alti^2) + slope + I(slope^2) + awc_top + I(awc_top^2) + bs_top + I(bs_top^2) + spht + slope:I(chbio_12-etp)"),
                                                originalVes = originalVes, 
                                                taxasToPredict = spIds,
                                                glc19SpIdVec= glc19SpIdVeco,
                                                nOccPerSp = nOccPerSp,
                                                extractorScript = "C:/Users/user/pCloud local/0_These/Github/Run/all_functions.R",
                                                dataDir="C:/Users/user/pCloud local/0_These/data/")


# make run 
run = data.frame(glc19testOccId=NA,glc19SpId=NA,Rank=NA,Probability=NA)
run = run[-1,,drop=F]
for(i in 1:dim(testSet)[1]){
  tmp = data.frame(glc19testOccId=testSet$glc19TestOccId[i],glc19SpId=colnames(predMatrix),Rank=NA,Probability=exp(predMatrix[i,]))
  if(sum(is.na(tmp$Probability))==0){
    tmp$Probability = tmp$Probability / sum(tmp$Probability)
  }else{
    tmp$Probability = nOccPerSp / sum(nOccPerSp)
  }
  tmp = tmp[order(tmp$Probability,decreasing = T),]
  tmp = tmp[1:50,]
  tmp$Rank = 1:50
  run = rbind(run,tmp)
  if(i/1000==round(i/1000)){
    flush.console()
    cat('      \r       Process...',100*i/dim(testSet)[1],'%         \r ')
  }
}

setwd('C:/Users/user/pCloud local/0_These/Autres formes de travaux/19_04_30 Pascal GLC19')
write.table(run,"LOF_1.csv",sep=";",row.names = F, col.names=F)


