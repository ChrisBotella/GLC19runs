# Make run of MAXENT_300 for GLC19
library(raster)
library(rgdal)
library(Matrix)
user = "Christophe"
repoDir = paste("C:/Users/",user,"/pCloud local/0_These/Github/SamplingEffort/",sep="")
setwd(repoDir)
source('_functions.R')

originalVes = c('etp','chbio_12','chbio_1','chbio_5','alti','slope','awc_top','bs_top','clc',"spht")
Intensityformula= as.formula(" ~ etp + I(etp^2) + I(chbio_12-etp) + I((chbio_12-etp)^2) + chbio_1 + I(chbio_1^2) + chbio_5 + I(chbio_5^2) + alti + I(alti^2) + slope + I(slope^2) + awc_top + I(awc_top^2) + bs_top + I(bs_top^2) + spht + slope:I(chbio_12-etp)")

# Load allowed_classes.txt
allowed = read.csv("P:/Partage GLC19/allowed_classes.txt",sep=";",header=F)
allowed = as.character(allowed[,1])

# Load lof model
modelsDir = paste("C:/Users/",user,"/pCloud local/0_These/data/LOF data/19_02_19/",sep="")

saveDir = paste("C:/Users/",user,"/pCloud local/0_These/data/personalWorker/files/",sep="")
setwd(saveDir)
occ = readRDS('occ_lof300_85_4000')
glc19SpIdVeco = readRDS('glc19SpId_List_lof300')
# Load testSet 
testSet = read.csv("P:/Partage GLC19/testSet.csv",sep=";",header=T)
spIds = intersect(allowed,glc19SpIdVeco)

predMatrix = matrix(0,dim(testSet)[1],length(spIds))
nOccPerSp = sapply(spIds,function(id) sum(occ$glc19SpId==id))

Sp = "30248"
mod = readRDS(paste('maxent_',Sp,sep=""))
coefficients = mod$beta[,1]
coefficients[1:length(coefficients)] = 0

get.Mat = function(coos,
                   Intensityformula,
                   originalVes,
                   extractorScript = "C:/Users/Christophe/pCloud local/0_These/Github/Run/all_functions.R",
                   dataDir="C:/Users/Christophe/pCloud local/0_These/data/"){
  # coos : numeric matrix with 2 columns named Longitude and Latitude (coordinates in the WGS84 system)
  listVE = load_variables()
  avail = NULL
  for(i in 1:length(listVE)){avail =  c(avail,listVE[[i]])}
  avail = unique(avail)
  toGet = intersect(avail,originalVes)
  
  extra = setdiff(originalVes,avail)
  if(length(extra)!=1 || extra!="spht"){
    print('Error: Unknown environmental variables required in originalVes')
    return(NULL)
  }
  coos = get_variables(toGet,coos,dpath = dataDir,fix=F)
  if("clc"%in%toGet & length(extra)==1 & extra=="spht"){coos$spht = get_spht(coos$clc)
  }else{print('Error: spht required whereas clc (original source of it) not in originalVes')}
  
  Mat = sparse.model.matrix(Intensityformula,data = coos)
  Matos = matrix(NA,dim(coos)[1],dim(Mat)[2])
  Matos = data.frame(Matos)
  colnames(Matos) = colnames(Mat)
  Matos[complete.cases(coos),] = Mat
  Mat = Matos 
  Mat = as.matrix(Mat)
  rm(Matos)
  gc(reset=T)
  return(Mat)
}

Mat = get.Mat(testSet,
               Intensityformula,
               originalVes,
               extractorScript = "C:/Users/Christophe/pCloud local/0_These/Github/SamplingEffort/_functions.R",
               dataDir="C:/Users/Christophe/pCloud local/0_These/data/")


#####
# MAXENT lofModel 300
#####

i=1
for(Sp in spIds){
  setwd(modelsDir)
  mod = readRDS(paste('maxent_',Sp,sep=""))
  noNull = mod$betas
  coefficients[names(noNull)] = noNull
  
  pred = as.numeric(Mat %*% coefficients)
  predMatrix[,i] = log( exp(pred) * nOccPerSp[i] / sum(exp(pred),na.rm=T) )
  
  coefficients[1:length(coefficients)] = 0
  i=i+1
  if(i/10==round(i/10)){
    flush.console()
    cat('     \r      Process...',100*i/length(spIds),'%        \r           ')
  }
}
colnames(predMatrix) = spIds

if(F){
  tmp = as.vector(mSp)
  classes = cut(tmp,c(min(tmp,na.rm=T)-.1,quantile(tmp,c(.1,.4,.6,.9),na.rm=T),max(tmp,na.rm=T)+.1),dig.lab=2,na.rm=T)
  colo = colorRampPalette(c("darkorchid4","goldenrod"))(length(unique(classes)))
  lev= levels(classes)
  labels = paste(lev,paste('/ quantile',c('0 to .1','.1 to .4','.4 to .6','.6 to .9','.9 to 1')))
  p = ggplot() + geom_tile(data=grid,aes(x=x_l93,y=y_l93,fill=classes),size=1,alpha=.8)+scale_fill_manual(values=colo,name="log10 of species intensity",labels=labels)+theme_bw()
  p = p + xlab("Longitude in Lambert 93") + ylab("Latitude in Lambert 93")
  print(p)
}

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

setwd('C:/Users/Christophe/pCloud local/0_These/Autres formes de travaux/19_04_30 Pascal GLC19')
write.table(run,"Lot_Of_Lof_2.csv",sep=";",row.names = F, col.names=F)



setwd('C:/Users/Christophe/pCloud local/0_These/Autres formes de travaux/19_04_30 Pascal GLC19')
run=read.csv("Lot_Of_Lof_2.csv",sep=";",header = F)
run[,2] = as.numeric(run[,2])

#####
# TGB LofModel 300
#####

coefficients[1:length(coefficients)] = 0
i=1
for(Sp in spIds){
  setwd(modelsDir)
  fileName = paste('tgb_lofModel_',Sp,sep="")
  if(fileName%in%list.files()){
    mod = readRDS(fileName)
    noNull = mod$betas
    coefficients[names(noNull)] = noNull
    
    pred = as.numeric(Mat %*% coefficients)
    predMatrix[,i] = log( exp(pred) * nOccPerSp[i] / sum(exp(pred),na.rm=T) )
  }
  coefficients[1:length(coefficients)] = 0
  i=i+1
  if(i/10==round(i/10)){
    flush.console()
    cat('     \r      Process...',100*i/length(spIds),'%        \r           ')
  }
}
colnames(predMatrix) = spIds

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
run$glc19SpId = as.numeric(run$glc19SpId)

setwd('C:/Users/Christophe/pCloud local/0_These/Autres formes de travaux/19_04_30 Pascal GLC19')
write.table(run,"Lot_Of_Lof_3.csv",sep=";",row.names = F, col.names=F)
