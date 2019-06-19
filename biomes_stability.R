library(raster)
library(maptools)
library(rgdal)
library(randomForest)
library(rpart)
library(rpart.plot)	
library(irr)
library(diffeR)

setwd("./climate_data_1k/")
landmask <- raster("./001/landmask.asc")
bioclim_01 <- raster("./000/bioclim_01.asc")/landmask
bioclim_12 <- raster("./000/bioclim_12.asc")/landmask
bioclim_15 <- raster("./000/bioclim_15.asc")/landmask

present_vars_list <- list(bioclim_01,bioclim_12,bioclim_15);rm(bioclim_01,bioclim_12,bioclim_15)

WWFbiomes <- readOGR("./wwfbiomes/ecoregions_SA.shp")


#1 = Tropical moist forests
#2 = Dry broadleaf forests
#4 = Temperate forests
#7 = Cerrado (savanna)
#8 = Patagonian grasslands
#9 = Pantanal
#10 = Montane grasslands
#12 = Deserts
#13 = Xeric shrublands (Caatinga)
#14 = Mangrooves
#98 = Lake
#99 = Rock and Ice

WWFbiomes$BIOME[WWFbiomes$BIOME==2] <- 13
WWFbiomes$BIOME[WWFbiomes$BIOME==9] <- 7
WWFbiomes$BIOME[WWFbiomes$ECO_NAME %in% c("Sechura desert","Atacama desert")] <- 12
WWFbiomes$BIOME[WWFbiomes$ECO_NAME %in% c("Chilean matorral")] <- 13

biomes_num <- sort(unique(WWFbiomes$BIOME))[1:7]

biomes_rst <- list()

for(i in biomes_num){
  e <- extent(WWFbiomes[WWFbiomes$BIOME==i,])
  ebiome <- crop(present_vars_list[[1]],e,snap="out")
  crop  <-  setValues(ebiome, NA)
  myshp.r  <-  rasterize(WWFbiomes[WWFbiomes$BIOME==i,], crop,field=i)
  biomes_rst[[which(biomes_num==i)]] <- myshp.r
}



varsbybiomes <- list(list(),list(),list())
#slots <- cbind(c(sort(unique(WWFbiomes$BIOME))))

for(j in 1:length(present_vars_list)){
  for(i in biomes_num){
    e <- extent(WWFbiomes[WWFbiomes$BIOME==i,])
    ebiome <- crop(present_vars_list[[j]],e,snap="out")
    crop  <-  setValues(ebiome, NA)
    myshp.r  <-  rasterize(WWFbiomes[WWFbiomes$BIOME==i,], crop) 
    varsbybiomes[[j]][[which(biomes_num==i)]]  <-  mask(x=ebiome, mask=myshp.r)
  }
}


biome.df.list <- list()

for(i in biomes_num){
  biome.stack <- stack(varsbybiomes[[1]][[which(biomes_num==i)]],varsbybiomes[[2]][[which(biomes_num==i)]],varsbybiomes[[3]][[which(biomes_num==i)]])
  names(biome.stack) <- c("temperature","precipitation","cv_precipitation")
  biome.df.list[[which(biomes_num==i)]] <- na.omit(as.data.frame(biome.stack))
}


### Model calibration

biome.df.R <- lapply(biome.df.list,FUN=function(x) as.vector(na.omit(x[sample(1:round(nrow(x)),round(nrow(x)*0.1)),])))
names(biome.df.R) <- biomes_num

RF_data <- data.frame(do.call("rbind", biome.df.R),biomes=rep(as.numeric(names(biome.df.R)),t=lapply(biome.df.R,nrow)))

#RF_data$biomes[RF_data$biomes==13] <- 12
#RF_data$biomes[RF_data$biomes==4] <- 1
#RF_data$biomes[RF_data$biomes==7] <- 12
#RF_data$biomes[RF_data$biomes==10] <- 8

RF_data2 <- RF_data
RF_data2$biomes[RF_data$biomes==11] <- 13

mod_RF <- randomForest(as.factor(biomes)~temperature+precipitation+cv_precipitation,data=RF_data2,ntree=500)
mod_RF2 <- randomForest(as.factor(biomes)~temperature+precipitation+wet_quarter,data=RF_data,ntree=500)
mod_RF3 <- randomForest(as.factor(biomes)~temperature+precipitation+dry_quarter,data=RF_data,ntree=500)


#1 = Tropical moist forests
#4 = Temperate forests
#7 = Cerrado (savanna)
#8 = Patagonian grasslands
#10 = Montane grasslands
#12 = Deserts
#13 = Xeric shrublands (Caatinga)


partialPlot(mod_RF,pred.data=RF_data,x.var=precipitation,"1")

tree <- getTree(mod_RF,1, labelVar=T)

colnames(tree) <- sapply(colnames(tree),collapse)
rules <- getConds(tree)
print(rules)

plot(mod_RF$forest$treemap)

to.dendrogram(tree)


### Model prediction data

present_vars <- stack(present_vars_list);rm(present_vars_list)
names(present_vars) <- c("temperature","precipitation","cv_precipitation")#,"wet_quarter","dry_quarter")

SA_data <- na.omit(as.data.frame(present_vars))

pred_Rforest <- predict(mod_RF,SA_data)

SA_data_pred <- data.frame(SA_data,biomes_pred=pred_Rforest)

library(rpart)

tree.1 = rpart(as.factor(biomes_pred)~temperature+precipitation+cv_precipitation,data=SA_data_pred)
prp(tree.1)    			

new.tree.1  <-  prp(tree.1,snip=TRUE)$obj # interactively trim the tree
prp(new.tree.1)


biomes_pred <- present_vars[[3]]
biomes_pred[!is.na(biomes_pred[])] <- as.numeric(as.character(pred_Rforest))

plot(biomes_pred)

bio01_list <- lapply(list.files(rec=T,pattern=paste("bioclim_01",".asc",sep=""),full=T),raster)
bio12_list <- lapply(list.files(rec=T,pattern=paste("bioclim_12",".asc",sep=""),full=T),raster)
bio15_list <- lapply(list.files(rec=T,pattern=paste("bioclim_15",".asc",sep=""),full=T),raster)

vars_list <- list(bio01_list,bio12_list,bio15_list)
rm(bio01_list,bio12_list,bio15_list)

biomes.list <- list()

for(i in 1:length(vars_list[[1]])){
  past_vars <- stack(vars_list[[1]][[i]]/landmask,vars_list[[2]][[i]]/landmask,vars_list[[3]][[i]]/landmask)
  names(past_vars) <- c("temperature","precipitation","cv_precipitation")
  SA_past <- as.data.frame(past_vars)
  SA_past <- na.omit(SA_past)
  past_Rforest <- predict(mod_RF,SA_past)
  biomes_pred_past <- present_vars[[3]]
  biomes_pred_past[!is.na(biomes_pred[])] <- past_Rforest
  biomes.list[[i]] <- biomes_pred_past
}

plot(biomes.list[[1]])

stability <- landmask
stability[which(getValues(stability)==1)] <- 0

for(i in 1:c(length(biomes.list)-1)){
  stability[which(getValues(biomes.list[[i]])!=getValues(biomes.list[[i+1]]))] <- stability[which(getValues(biomes.list[[i]])!=getValues(biomes.list[[i+1]]))] + 1
}

plot(stability)
blue2red  <-  colorRampPalette(c("darkblue","darkblue","blue","yellow","orange","orange","red","red","darkred"),space = "rgb")

writeRaster(stability,file="stability.grd",overwrite=T)

### Model fit

biomes_pred <- present_vars[[3]]
biomes_pred[!is.na(biomes_pred[])] <- as.numeric(as.character(pred_Rforest))


biomes_SA_list <- list()

for(i in 1:7){
  biomes_SA_list[[i]]  <-  rasterize(WWFbiomes[WWFbiomes$BIOME==biomes_num[i],], landmask,field=biomes_num[i])
}

biomes_rst_full <- present_vars[[3]]
biomes_rst_full[!is.na(biomes_pred[])] <- NA
for(i in 1:7){
  biomes_rst_full[!is.na(biomes_SA_list[[i]])] <- biomes_SA_list[[i]][!is.na(biomes_SA_list[[i]])]
}

biomes_pred[is.na(biomes_rst_full)] <- NA
biomes_rst_full[is.na(biomes_pred)] <- NA

#biomes_rst_full[biomes_rst_full[]==11] <- 13

biomes_rst_full2 <- biomes_rst_full
biomes_rst_full2[is.na(biomes_acc)] <- NA

biomes_pred2 <- biomes_pred
biomes_pred2[is.na(biomes_acc)] <- NA

acc_data <- data.frame(biomes=as.factor(na.omit(biomes_rst_full2[])),pred=as.factor(na.omit(biomes_pred2[])), weights=na.omit(biomes_acc[]))

write.table(acc_data, file="acc_data.txt", row.names = F)

acc_data <- read.table("acc_data.txt", h=T)

library(irr)
kappa2(acc_data[,1:2], weight = acc_data[,3])

confusion_table <- table(data.frame(biomes=as.factor(na.omit(biomes_rst_full[])),pred=as.factor(na.omit(biomes_pred[]))))

#Error by biomes
round((rowSums(confusion_table)-diag(confusion_table))/rowSums(confusion_table),2)

#Overall error
1-sum(diag(confusion_table))/sum(confusion_table)


biomes_acc <- raster("C:/Users/gmazz/OneDrive/Documentos/Stability/biomes_acc_wwf.grd")


overallDiff(confusion_table)

overallQtyD(confusion_table)/overallDiff(confusion_table)
overallAllocD(confusion_table)/overallDiff(confusion_table)

overallQtyD(confusion_table)/sum(confusion_table)
overallAllocD(confusion_table)/sum(confusion_table)

overallExchangeD(confusion_table)/overallDiff(confusion_table)
overallShiftD(confusion_table)

sum(overallExchangeD(confusion_table),overallShiftD(confusion_table))


a <- confusion_table
diag(a) <- 0
sum(a)

round(confusion_table[,1]/sum(confusion_table[,c(1,3)]),2)
round(confusion_table[,3]/sum(confusion_table[,3]),3)
round(confusion_table[,7]/sum(confusion_table[,7]),2)

#######################################
RF_fit_table <- matrix(nr=nrow(SA_data),nc=100)

for(i in 1:100){
  biome.df.R_ <- lapply(biome.df.list,FUN=function(x) as.vector(na.omit(x[sample(1:round(nrow(x)),round(nrow(x)*0.1)),])))
  biomes_ <- list()
  for(j in 1:length(biome.df.R_)){
    biomes_[[j]] <- rep(biomes_num[j],nrow(biome.df.R_[[j]]) ) 
  }
  biomes_ <- do.call("c", biomes_)
  RF_data_ <- data.frame(do.call("rbind", biome.df.R_),as.factor(biomes_))
  mod_RF_ <- randomForest(as.factor(biomes_)~temperature+precipitation+cv_precipitation,data=RF_data_)
  pred_Rforest_ <- predict(mod_RF_,SA_data)
  RF_fit_table[,i] <- as.numeric(as.character(pred_Rforest_))
  cat(paste(i, "- "))
}

biomes_data <- biomes_rst_full[as.numeric(names(pred_Rforest_))]

RF_fit_table[biomes_data==13,2]

biomes_acc <- landmask
biomes_acc[] <- NA
biomes_acc[as.numeric(names(pred_Rforest_))] <- rowMeans(apply(RF_fit_table[,1:100],2,function(x) biomes_data==x))
writeRaster(biomes_acc, file="biomes_acc_wwf.grd")


biomes_acc[is.na(stability)] <- NA
stability[is.na(biomes_rst)] <- NA
biomes_acc[is.na(biomes_rst)] <- NA

plot(biomes_acc)
windows(14,7)
par(mfrow=c(1,2))
plot(biomes_rst,col=c("darkgreen","turquoise4","orange","yellow","orange4","red","indianred3"),main="Biomes",legend=F)
plot(stability,col=blue2red(21),main="Stability")

windows(14,7)
par(mfrow=c(1,2))
plot(biomes_rst,col=c("darkgreen","turquoise4","orange","yellow","orange4","red","indianred3"),main="Biomes",legend=F)
plot(stability,col=rich.colors(21),main="Stability")

#1 = Tropical moist forests
#2 = Temperate woodlands
#3 = Savannas
#4 = Patagonia grasslands
#5 = Montane grasslands
#6 = Deserts
#7 = Seasonal Dry Forests

plot(biomes_acc,col=rev(blue2red(25)),main="Accuracy")
plot(biomes_acc,col=rev(rich.colors(25)),main="Accuracy")


writeRaster(stability,"stability.asc")
writeRaster(biomes_rst,"biomes.asc")



##### Criar tabela com os thresholds climaticos

SA_data$biomes_pred <- pred_Rforest
threshold_table <- data.frame(matrix(nr=length(unique(SA_data$biomes_pred)),nc=9))

for(i in sort(unique(SA_data$biomes_pred))){
  qunt1 <- apply(SA_data[SA_data$biomes_pred==i,1:3],2,quantile,probs = c(0.975, 0.025))
  threshold_table[which(sort(unique(SA_data$biomes_pred))==i),1:7] <- c(i,diag(qunt1[c(1,2,1,2,1,2),c(1,1,2,2,3,3)])*c(0.1,0.1,1,1,0.1,0.1))
}  

threshold_table[,1]  <- c("Tropical moist forests",
                        "Temperate woodlands",
                        "Savannas",
                        "Patagonia grasslands",
                        "Montane grasslands",
                        "Deserts",
                        "Seasonal Dry Forests")

colnames(threshold_table) <- c("Biomes","Max Temp","Min Temp","Max Rain","Min Rain","Max Seasonality","Min Seasonality","N","Error rate")

threshold_table[,8] <- rowSums(mod_RF$confusion[,1:7])
threshold_table[,9] <- round(mod_RF$confusion[,8],3)


write.table(na.omit(threshold_table),"threshold_table.csv",row.names=F)


range_table <- data.frame(matrix(nr=11,nc=7))
for(i in 1:11){
  range1 <- apply(SA_data[SA_data$biomes_pred==i,1:3],2,range)
  range_table[i,] <- c(i,diag(range1[c(1,2,1,2,1,2),c(1,1,2,2,3,3)])*c(0.1,0.1,1,1,0.1,0.1))
}  

range_table[,1] <- threshold_table[,1]
colnames(range_table) <- colnames(threshold_table)
write.table(na.omit(range_table),"range_table.csv",row.names=F)

mean_table <- data.frame(matrix(nr=11,nc=4))
for(i in 1:11){
  mean_table[i,2:4] <- apply(SA_data[SA_data$biomes_pred==i,1:3],2,mean)
}  

mean_table[,1] <- threshold_table[,1]
colnames(mean_table) <- c("Biomes","Mean Temperature","Mean Rainfall","Mean Seasonality")
write.table(na.omit(mean_table),"mean_table.csv",row.names=F)

threshold_10_table <- data.frame(matrix(nr=11,nc=7))
for(i in 1:11){
  qunt1 <- apply(SA_data[SA_data$biomes_pred==i,1:3],2,quantile,probs = c(0.95, 0.05))
  threshold_10_table[i,] <- c(i,diag(qunt1[c(1,2,1,2,1,2),c(1,1,2,2,3,3)])*c(0.1,0.1,1,1,0.1,0.1))
}  

threshold_10_table[,1] <- threshold_table[,1]
colnames(threshold_10_table) <- colnames(threshold_table)
write.table(na.omit(threshold_10_table),"threshold_10_table.csv",row.names=F)


