# created on 18/02/2019
# Comparison of LOF, MAXENT, TGB (and TGBO) predictive performance on real data
library(glmnet)
library(maxnet)

user = "Christophe"
setwd(paste('C:/Users/',user,"/pCloud local/0_These/Github/R/_base/",sep=""))
source('functions.R')

originalVes = c('etp','chbio_12','chbio_1','chbio_5','alti','slope','awc_top','bs_top','clc')
predictor_formula = " ~ etp + I(etp^2) + I(chbio_12-etp) + I((chbio_12-etp)^2) + chbio_1 + I(chbio_1^2) + chbio_5 + I(chbio_5^2) + alti + I(alti^2) + slope + I(slope^2) + awc_top + I(awc_top^2) + bs_top + I(bs_top^2) + spht + slope:I(chbio_12-etp)"

saveDir = paste('C:/Users/',user,"/pCloud local/0_These/data/LOF data/19_02_19/",sep="")

MinOcc = 50

######
# Functions
######

# Juste la mise en forme de la table

prepare.for.maxnet=function(data,variables,factors){
  # select variables and convert factors
  data=data[,colnames(data)%in%variables]
  for(c in variables){
    if(c%in%factors){
      data[,colnames(data)==c]=factor(data[,colnames(data)==c])
    }else{
      data[,colnames(data)==c]=as.numeric(data[,colnames(data)==c])
    }
  }
  return(data)
}

# Change land cover to spht := simplified plant habitat type = urban / arable / grasses / forest / other
get_spht = function(clcVec){
  setwd(paste('C:/Users/',user,'/pCloud local/0_These/data/MTAP article/data/',sep=""))
  clc.variables = readRDS('clc.variables.RData')
  names = clc.variables[[1]]
  
  spht = rep(NA,length(clcVec))
  cd = clcVec %in% clc.variables[[2]][[which(names=='arable')]] | clcVec%in%c(12,13,15,16,20) 
  spht[cd] = "cultivated"
  cd = clcVec %in% clc.variables[[2]][[which(names=='pasture')]] | clcVec %in% clc.variables[[2]][[which(names=='nat_grass')]] | clcVec %in% clc.variables[[2]][[which(names=='moors')]] | clcVec %in% clc.variables[[2]][[which(names=='sclero')]]
  spht[cd] = "grasses"
  cd = clcVec %in% clc.variables[[2]][[which(names=='brl_for')]] | clcVec %in% clc.variables[[2]][[which(names=='coni_for')]] | clcVec %in% clc.variables[[2]][[which(names=='mixed_for')]] | clcVec %in% clc.variables[[2]][[which(names=='transi_wood')]]
  spht[cd] = "forest"
  cd = clcVec %in% clc.variables[[2]][[which(names=='arti')]] | clcVec %in% clc.variables[[2]][[which(names=='semi_arti')]] | clcVec%in%c(11)
  spht[cd] = "urban"
  spht[is.na(spht)] = "other"
  return(spht)
}


maxnet=function (p, data, f = maxnet.formula(p, data), regmult = 1, 
                 regfun = maxnet.default.regularization, ...) 
{
  mm <- model.matrix(f, data)
  reg <- regfun(p, mm) * regmult
  weights <- p + (1 - p) * 100
  glmnet::glmnet.control(pmin = 1e-08, fdev = 0)
  model <- glmnet::glmnet(x = mm, y = as.factor(p), family = "binomial", 
                          standardize = F, penalty.factor = reg, lambda = 10^(seq(4, 
                                                                                  0, length.out = 200)) * sum(reg)/length(reg) * sum(p)/sum(weights), 
                          weights = weights, ...)
  class(model) <- c("maxnet", class(model))
  if (length(model$beta) < 200) 
    stop("Error: glmnet failed to complete regularization path")
  bb <- model$beta[, 200]
  model$betas <- bb[bb != 0]
  model$alpha <- 0
  #rr <- predict.maxnet(model, data[p == 0, , drop = FALSE], 
  #                     type = "exponent", clamp = F)
  #raw <- rr/sum(rr)
  #model$entropy <- -sum(raw * log(raw))
  #model$alpha <- -log(sum(rr))
  model$penalty.factor <- reg
  model$featuremins <- apply(mm, 2, min)
  model$featuremaxs <- apply(mm, 2, max)
  vv <- (sapply(data, class) != "factor")
  model$varmin <- apply(data[, vv, drop = FALSE], 2, min)
  model$varmax <- apply(data[, vv, drop = FALSE], 2, max)
  means <- apply(data[p == 1, vv, drop = FALSE], 2, mean)
  majorities <- sapply(names(data)[!vv], function(n) which.max(table(data[p == 
                                                                            1, n, drop = FALSE])))
  names(majorities) <- names(data)[!vv]
  model$samplemeans <- c(means, majorities)
  model$levels <- lapply(data, levels)
  model
}

######
# Extract species occurrences data
######

setwd('P:/Partage GLC19/')
OCC= read.csv("PL_trusted.csv",sep=";",header=T)

# Get environmental variables 
OCC = get_variables(originalVes,OCC,dpath = paste('C:/Users/',user,"/pCloud local/0_These/data/",sep=""),fix=F)
OCC$spht = factor(get_spht(OCC$clc))
# remove NAs
tokeep = complete.cases(OCC[,colnames(OCC)%in%originalVes])
OCC = OCC[tokeep,,drop=F]

# remove species with less than MinOcc occurrences
length(unique(OCC$glc19SpId))
counts = table(OCC$glc19SpId)
spToKeep = names(counts)[counts>=MinOcc]
length(spToKeep)

occ = OCC[OCC$glc19SpId%in% spToKeep,]

print(length(spToKeep))

#occNoCity = occ[!occ$clc%in%c(1,6,10,11),,drop=F]

######
# Create LOF grid
# remove empty or scarce squares
######

squareSize = 4000

## Make raster of squares including all the metropolitan French territory  
# get france polygon 
setwd(paste('C:/Users/',user,"/pCloud local/0_These/Github/R/_base/",sep=""))
france = getData('GADM',country="FRA",level="0")
# Project to lambert 93 
frProj = spTransform(france,CRS("+proj=lcc +lat_1=49 +lat_2=44 +lat_0=46.5 +lon_0=3 +x_0=700000 +y_0=6600000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"))

plot(frProj)

frBuff= gBuffer(frProj,width=squareSize,joinStyle="ROUND") # create France polygon with 4km buffer area around borders
ext = extent(frBuff)
r = raster(ext= ext,resolution=c(squareSize,squareSize),crs = CRS("+proj=lcc +lat_1=49 +lat_2=44 +lat_0=46.5 +lon_0=3 +x_0=700000 +y_0=6600000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"))
#r[] = runif(ncell(r),0,1)
r_frBuf = rasterize(frBuff, r)  
# Values of the raster = Indices of squares  
r_frBuf[!is.na(r_frBuf[])] = 1:(sum(!is.na(getValues(r_frBuf))))
plot(r_frBuf)

# Number of cells
sum(!is.na(getValues(r_frBuf)))

## Remove squares with less than Min=5 occurrences 
Min= 5
# Count occurrences per cell
pts = SpatialPoints(occ[,c('Longitude','Latitude')],proj4string = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))
occProj = spTransform(pts,CRSobj = CRS("+proj=lcc +lat_1=49 +lat_2=44 +lat_0=46.5 +lon_0=3 +x_0=700000 +y_0=6600000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"))

occ = cbind(occ,as.data.frame(occProj@coords))
colnames(occ)[(dim(occ)[2]-1):dim(occ)[2]]= c('x_lamb93','y_lamb93')

occ$squareId = extract( r_frBuf, occ[,c('x_lamb93','y_lamb93')]  )

tab = table(occ$squareId)
tab = tab[order(tab,decreasing = T)] 
barplot(tab)
sum(tab>Min)

indicesToKeep = as.numeric(names(tab)[tab>Min])

r_lofCells = r_frBuf
r_lofCells[!r_lofCells[]%in%indicesToKeep] = NA
plot(r_lofCells)

sum(occ$squareId%in%indicesToKeep) # Number of occurrences for fitting LOF

######
# Prepare background points for MAXENT and LOF 
######

# We uniformly draw radom points inside the square extent of the raster and
# (i) MAXENT: Keep those falling inside the france raster (4 km buffer) for MAXENT, until at least 3 per cell
# (i) LOF: Only keep those falling inside LOF squares (until at least 5 per square)

nTmp = 30000
ext = extent(r_frBuf)
LOF_background = data.frame(x_lamb93 = NA , y_lamb93 = NA , q_hash = NA)
LOF_background = LOF_background[-1,,drop=F]
MAX_background = LOF_background
minNPerSquareLOF = 0
minNPerSquareMAX = 0 
while(minNPerSquareLOF<5 | minNPerSquareMAX<3){
  tmp = data.frame(x_lamb93 = runif(nTmp,ext[1],ext[2]) , y_lamb93 = runif(nTmp,ext[3],ext[4]) , q_hash = NA)
  
  if(minNPerSquareMAX<3){
    # Filter for MAXENT
    tmp$q_hash = extract(r_frBuf,tmp[,1:2])
    cdMAX = !is.na(tmp$q_hash)
    MAX_background = rbind(MAX_background,tmp[cdMAX,,drop=F])
    minNPerSquareMAX = min(table(MAX_background$q_hash))
  }
  
  if(minNPerSquareLOF<5){
    # Filter for LOF
    tmp$q_hash = extract(r_lofCells,tmp[,1:2])
    cdLOF = !is.na(tmp$q_hash)
    LOF_background = rbind(LOF_background,tmp[cdLOF,,drop=F])
    minNPerSquareLOF = min(table(LOF_background$q_hash))
  }
  
  flush.console()
  cat('     \r  Mini. n° points/LOF square:',minNPerSquareLOF,', tot. n° LOF pts:',dim(LOF_background)[1],'Mini. n° pts/MAX square:',minNPerSquareMAX,'tot. n° MAX pts:',dim(MAX_background)[1],'        \r')
}

## Get environmental variables
# MAXENT
pts = SpatialPoints(MAX_background[,1:2],proj4string = CRS("+proj=lcc +lat_1=49 +lat_2=44 +lat_0=46.5 +lon_0=3 +x_0=700000 +y_0=6600000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs") )
pts = spTransform(pts, CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))
MAX_background = cbind(MAX_background,data.frame(pts@coords))
colnames(MAX_background)[4:5] = c('Longitude','Latitude')
MAX_background = get_variables(originalVes,MAX_background,dpath = paste('C:/Users/',user,"/pCloud local/0_These/data/",sep=""),fix=F)
MAX_background$spht = factor(get_spht(MAX_background$clc))
nona = complete.cases(MAX_background)
MAX_background = MAX_background[nona,]
# LOF
pts = SpatialPoints(LOF_background[,1:2],proj4string = CRS("+proj=lcc +lat_1=49 +lat_2=44 +lat_0=46.5 +lon_0=3 +x_0=700000 +y_0=6600000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs") )
pts = spTransform(pts, CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))
LOF_background = cbind(LOF_background,data.frame(pts@coords))
colnames(LOF_background)[4:5] = c('Longitude','Latitude')
LOF_background = get_variables(originalVes,LOF_background,dpath = paste('C:/Users/',user,"/pCloud local/0_These/data/",sep=""),fix=F)
LOF_background$spht = factor(get_spht(LOF_background$clc))
nona = complete.cases(LOF_background)
LOF_background = LOF_background[nona,]

#####
# Prepare TGB background
#####

### from the 600 species
## Define the sites grid with the thinner raster (alti)
r = raster(paste('C:/Users/',user,'/pCloud local/0_These/data/0_mydata/alti/alti.tif',sep=""))
sum(is.na(getValues(r)))

r[] = 1:ncell(r)
rIds = extract(r,occ[,c('Longitude','Latitude')])
rIds = unique(rIds)
nas = setdiff(1:ncell(r),rIds)
r[nas] = NA
tgb = rasterToPoints(r)
tgb = tgb[,1:2]
colnames(tgb)=c('Longitude','Latitude')
tgb = as.data.frame(tgb)

tgb = get_variables(originalVes,tgb,dpath = paste('C:/Users/',user,"/pCloud local/0_These/data/",sep=""),fix=F)
tgb$spht = factor(get_spht(tgb$clc))
# remove NAs
tokeep = complete.cases(tgb[,colnames(tgb)%in%originalVes])
tgb = tgb[tokeep,,drop=F]

setwd(saveDir)
saveRDS(tgb,'bgTGB')

#####
# Fit TGB
#####

setwd(saveDir)
tgb = readRDS('bgTGB')

## Apply maxent with same variables as LOF
columns = c(originalVes,'spht')
for(sp in spToKeep){
  mod = NULL
  print(which(sp==spToKeep))
  occurrences = occ[occ$glc19SpId==sp,colnames(occ)%in%columns,drop=F]
  df = rbind(occurrences,tgb[,colnames(tgb)%in%columns])
  p = c(rep(1,dim(occurrences)[1]),rep(0,dim(tgb)[1]))
  df = prepare.for.maxnet(data = df,variables = columns,factors = c('spht'))
  
  deb = Sys.time()
  mod = try( maxnet(p,df,f = as.formula(paste("p",predictor_formula))) )
  print(Sys.time()-deb)
  if(all.equal(attr(mod,"class"),c("maxnet","lognet","glmnet"))){
    setwd(saveDir)
    saveRDS(mod,paste("tgb_lofModel_",sp,sep=""))
  }
  gc(reset=T)
}

## Apply maxent with its automatic variables transformations
# Change Pedologic variables to factors for MaxentModel
tgbPedoFact = tgb
tgbPedoFact$awc_top = factor(tgb$awc_top)
tgbPedoFact$bs_top = factor(tgb$bs_top)
MaxentVes= c(originalVes[originalVes!='clc'],'spht')
tgbPedoFact = tgb[,MaxentVes]

for(sp in spToKeep){
  mod = NULL
  print(which(sp==spToKeep))
  occurrences = occ[occ$glc19SpId==sp,colnames(occ)%in%columns,drop=F]
  # Change Pedologic variables to factors for MaxentModel
  occurrences$bs_top = factor(occurrences$bs_top,levels=levels(tgbPedoFact$bs_top))
  occurrences$awc_top = factor(occurrences$awc_top,levels=levels(tgbPedoFact$awc_top))
  occurrences = occurrences[,MaxentVes]
  
  df = rbind(occurrences,tgbPedoFact)
  p = c(rep(1,dim(occurrences)[1]),rep(0,dim(tgb)[1]))
  df = prepare.for.maxnet(data = df,variables = columns,factors = c('spht'))
  
  deb = Sys.time()
  mod = try( maxnet(p,df ))
  print(Sys.time()-deb)
  if(all.equal(attr(mod,"class"),c("maxnet","lognet","glmnet"))){
    setwd(saveDir)
    saveRDS(mod,paste("tgb_MaxModel_",sp,sep=""))
  }
  gc(reset=T)
} 


#####
# Fit MAXENT models
#####

columns = c(originalVes,'spht')
for(sp in spToKeep){
  print(which(sp==spToKeep))
  occurrences = occ[occ$glc19SpId==sp,colnames(occ)%in%columns,drop=F]
  df = rbind(occurrences,MAX_background[,colnames(MAX_background)%in%columns])
  p = c(rep(1,dim(occurrences)[1]),rep(0,dim(MAX_background)[1]))
  df = prepare.for.maxnet(data = df,variables = columns,factors = c('spht'))
  
  deb = Sys.time()
  mod = try( maxnet(p,df,f = as.formula(paste("p",predictor_formula))) )
  print(Sys.time()-deb)
  if(all.equal(attr(mod,"class"),c("maxnet","lognet","glmnet"))){
    setwd(saveDir)
    saveRDS(mod,paste("maxent_",sp,sep=""))
  }
}

## Apply maxent with its automatic variables transformations
# Change Pedologic variables to factors for MaxentModel
MAXbg= MAX_background
MAXbg$awc_top = factor(MAXbg$awc_top)
MAXbg$bs_top = factor(MAXbg$bs_top)
MaxentVes= c(originalVes[originalVes!='clc'],'spht')
MAXbg = MAXbg[,MaxentVes]
for(sp in spToKeep[341:length(spToKeep)]){
  print(which(sp==spToKeep))
  occurrences = occ[occ$glc19SpId==sp,colnames(occ)%in%columns,drop=F]
  # Change Pedologic variables to factors for MaxentModel
  occurrences$bs_top = factor(occurrences$bs_top,levels=levels(MAXbg$bs_top))
  occurrences$awc_top = factor(occurrences$awc_top,levels=levels(MAXbg$awc_top))
  occurrences = occurrences[,MaxentVes]
  
  df = rbind(occurrences,MAXbg)
  p = c(rep(1,dim(occurrences)[1]),rep(0,dim(MAXbg)[1]))
  df = prepare.for.maxnet(data = df,variables = columns,factors = c('spht'))
  
  deb = Sys.time()
  mod = try( maxnet(p,df ))
  print(Sys.time()-deb)
  if(length(attr(mod,"class"))==3 & all.equal(attr(mod,"class"),c("maxnet","lognet","glmnet"))){
    setwd(saveDir)
    saveRDS(mod,paste("maxent_MaxModel_",sp,sep=""))
  }
  gc(reset=T)
} 


#####
# Fit LOF models
#####

gc(reset=T)
occLOF = occ[occ$squareId%in%indicesToKeep,]
occLOF= occLOF[,c('glc19SpId',originalVes,'spht','x_lamb93','y_lamb93')]

### Format data for glmnet LOF 
occLOF$glc19SpId = factor(occLOF$glc19SpId)

LOF_background$glc19SpId = factor(NA,levels=levels(occLOF$glc19SpId))
LOF_background$q_hash = factor(LOF_background$q_hash)
LOF_background$spht = factor(LOF_background$spht)

occLOF$q_hash = extract(r_lofCells,occLOF[,c('x_lamb93','y_lamb93')])
occLOF$q_hash = factor(occLOF$q_hash,levels = levels(LOF_background$q_hash))
occLOF$spht = factor(occLOF$spht,levels=levels(LOF_background$spht))
occLOF = occLOF[complete.cases(occLOF),]

# Save LOF occ and BackgroundPts
setwd(saveDir)
saveRDS(occLOF,'occLOF')
saveRDS(LOF_background,'bgLOF')

### Prepare Design matrix
Bias_killer=100
Predictor = paste( " ~ glc19SpId*(",substr(predictor_formula,3,nchar(predictor_formula)),") + q_hash")
PTS = sparse.model.matrix(as.formula(Predictor),occLOF)

for(e in levels(occLOF$glc19SpId)){
  print(which(levels(occLOF$glc19SpId)==e))
  LOF_background$glc19SpId = factor(rep(e,dim(LOF_background)[1]),levels=levels(LOF_background$glc19SpId))
  
  tmpBG = sparse.model.matrix(as.formula(Predictor),LOF_background)
  if(which(levels(occLOF$glc19SpId)==e)==1){BG=tmpBG}else{BG=rbind(BG,tmpBG)}
  print(paste(as.numeric(object.size(BG))/(8*1024^2)," Mo"))
}
BG = rbind(PTS,BG)

# prepare weights and  pseudo_y output
w1 = rep(NA,dim(occLOF)[1])
y1 = rep(NA,dim(occLOF)[1])
for(e in levels(occLOF$glc19SpId)){
  cd = as.character(occLOF$glc19SpId)==e
  n = sum(cd)
  w1[cd] = 1/(Bias_killer*n)
  y1[cd] = Bias_killer*n
}
n_0 = dim(LOF_background)[1]
w0 = rep( (Bias_killer-1) / (n_0*Bias_killer) , dim(BG)[1]-dim(PTS)[1] ) 
y0 = rep( 0 , dim(BG)[1]-dim(PTS)[1] )
w = c(w1,w0)
y = c(y1,y0)

### Fit model
penaltyFactors = as.numeric(regexpr('q_hash',colnames(BG))<=0)
penaltyFactors[penaltyFactors==0] = .1

mod = glmnet(x=BG,y=y,family="poisson",weights = w,penalty.factor=penaltyFactors)
setwd(saveDir)
saveRDS(mod,'LOF_glmnet_model_sup350occ')

#####
# Select best Lambda value with CV
#####

setwd(saveDir)
occLOF = readRDS('occLOF')
LOF_background = readRDS('bgLOF')
BG = readRDS('DesSparseMatrix')

Bias_killer=100
# prepare weights and  pseudo_y output
w1 = rep(NA,dim(occLOF)[1])
y1 = rep(NA,dim(occLOF)[1])
for(e in levels(occLOF$glc19SpId)){
  cd = as.character(occLOF$glc19SpId)==e
  n = sum(cd)
  w1[cd] = 1/(Bias_killer*n)
  y1[cd] = Bias_killer*n
}
n_0 = dim(LOF_background)[1]
w0 = rep( (Bias_killer-1) / (n_0*Bias_killer) , dim(BG)[1]-dim(occLOF)[1] ) 
y0 = rep( 0 , dim(BG)[1]-dim(occLOF)[1] )
w = c(w1,w0)
y = c(y1,y0)

penaltyFactors = as.numeric(regexpr('q_hash',colnames(BG))<=0)
penaltyFactors[penaltyFactors==0] = .1

# Approx Max Bound of proba of 0 background points in at least one Cell because of CV partition
nfold= 15
print(((1-(1/nfold)^5)^(9000)*(1/nfold)^5)*9000)
CV = cv.glmnet(x=BG,y=y,weights=w,nfolds=nfold,penalty.factor=penaltyFactors,family="poisson")
setwd(saveDir)
saveRDS(CV,'LOF_CVglmnet_sup350occ')

#####
# Plot LOF observation effort 
#####

setwd(paste('C:/Users/',user,'/pCloud local/0_These/data/LOF data/19_02_19/',sep=""))
mod = readRDS('LOF_glmnet_model_sup350occ')
occLOF = readRDS('occLOF')
coefficients= mod$beta[,dim(mod$beta)[2],drop=F]


q_coef = coefficients[regexpr('q_hash',rownames(coefficients))>0,,drop=F]
q_coef = as.matrix(q_coef)
NumQ_hash = sapply(rownames(q_coef),function(n) as.numeric(substr(n,7,nchar(n))))

## Make raster of squares including all the metropolitan French territory  
setwd(paste('C:/Users/',user,"/pCloud local/0_These/Github/R/_base/",sep=""))
france = getData('GADM',country="FRA",level="0")
frProj = spTransform(france,CRS("+proj=lcc +lat_1=49 +lat_2=44 +lat_0=46.5 +lon_0=3 +x_0=700000 +y_0=6600000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"))
frBuff= gBuffer(frProj,width=squareSize,joinStyle="ROUND") # create France polygon with 4km buffer area around borders
ext = extent(frBuff)
r = raster(ext= ext,resolution=c(squareSize,squareSize),crs = CRS("+proj=lcc +lat_1=49 +lat_2=44 +lat_0=46.5 +lon_0=3 +x_0=700000 +y_0=6600000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"))
r_frBuf = rasterize(frBuff, r)  
r_frBuf[!is.na(r_frBuf[])] = 1:(sum(!is.na(getValues(r_frBuf))))
Min= 5
occLOF$squareId = extract( r_frBuf, occLOF[,c('x_lamb93','y_lamb93')]  )
tab = table(occLOF$squareId)
tab = tab[order(tab,decreasing = T)] 
indicesToKeep = as.numeric(names(tab)[tab>Min])
r_lofCells = r_frBuf
r_lofCells[!r_lofCells[]%in%indicesToKeep] = NA

# Observation effort
library(ggplot2)
df = as.data.frame(rasterToPoints(r_lofCells))
colnames(df)[3] = "q_hash"
df$coef = sapply(df$q_hash,function(hash){if(sum(NumQ_hash==hash)==1){q_coef[NumQ_hash==hash,1]}else{NA}})
ggplot(df,aes(x=x,y=y))+geom_tile(aes(fill=log10(exp(coef))))+scale_fill_gradient(low="darkblue",high="goldenrod3")+theme_bw()


df = data.frame(RastId= 1:ncell(r_lofCells), q_hash = getValues(r_lofCells),coef= NA)
df = df[!is.na(df$q_hash),]
df$coef = sapply(df$q_hash,function(hash){if(sum(NumQ_hash==hash)==1){q_coef[NumQ_hash==hash,1]}else{NA}})
r_effort = r_lofCells
r_effort[df$RastId] = log10(exp(df$coef))
plot(r_effort)



