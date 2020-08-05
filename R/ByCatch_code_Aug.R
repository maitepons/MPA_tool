############## Instructions #################################################
# Format of data entries:                                                   #
# There are 2 csv files that are needed,                                    #
# one called Data.csv and another one called Weights.csv                    #
# 1. data.csv: the first 5 columns are the same for all case studies:       # 
# Lat	Lon	Year Month Effort                                                 #
# The following columns may vary for each case study                        #
# put first the target species and then the bycatch by species              #
# you can include any number of target and bycatch species                  #
# or group of species as your convinience.                                  #
# IMPORTANT: each row in the data file corresponds to a quadrant (lat&long) #
# by year and month. There is only one record for each of these combinations# 
# The numbers in each column for target and bycatch species can be          #
# in numbers or biomass.                                                    #
# 2. Weights.csv: this csv has in the first row the names of the target     #
# species in the same order as in the Data.csv file and the second row      #
# is the weight for each target species. They should sum up to 1.           #
# The third row is the names of the Bycatch species                         #
# in the same order as in the Data.csv file. The fourth row is the weight   #
# for each bycatch species and they also should sum up to 1.                #
# Outputs:                                                                  #
# Running the code generates:                                               #
# 1. Plot with original proportions of catch by species to total catch      #
#    and proportions based in specified weights                             #
# 2. Correlation plot among Effort, target and bycatch species              #    
# 3. Maps with original data                                                #
# 4. running the function DoCalcs generates:                                #
# a) lines plot: changes in ByCatch, target catch, fishing efficiency and   # 
#   Effort for each % of area closure or number of months depending on      #
#   wether by.month=TRUE or FALSE. This is represented by the solid lines   # 
#   If by.month=F, the dashed lines represent the dinamic closures,         # 
#   which means closing a different area each year. In this case, each grey #
#   line represents a year (this can be turned off if tmp.lines=FALSE).     #
# b) barplot: showing relative changes in catch for each species            #
# c) maps: showing area closed in grey and effort, target and ByCatch CPUE  #
# d) RDS output with results to be use in future analysis                   #
# Outputs are saved in a folder called "Results"                            #
#############################################################################
rm(list=ls())
############################################
Fishery_name<-"Your_fishery_name" ## CHANGE!!! to a name that you can recognize
                                  # if you are working with multiple fisheries
####### Libraries needed
library(tidyverse)
library(reshape2)
library(grid)
library(gridExtra)
library(ggcorrplot)
library(dplyr)
###########################################################################
# Read data and weights csv files  ########################################
###########################################################################
# setwd("") # use thsi to set your working directory if needed#############
dir.create(path=paste(getwd(),"Results",sep="/"),showWarnings = FALSE)    # 
###########################################################################
D<-read.csv("Data.csv",header = T) 

Tweights<-read.csv("Weights.csv",nrows = 1)%>%
  select_if(~ !any(is.na(.))) 

BCweights<-read.csv("Weights.csv",nrows = 2,skip = 2)%>%
  select_if(~ !any(is.na(.))) 

NT<-ncol(Tweights) # number of target species 
NBC<-ncol(BCweights) # number of byCatch species

### Effort in thousands 
D$Effort<-D$Effort/1000

BCNames<-colnames(BCweights)
TNames<-colnames(Tweights)

col.BC <- which(colnames(D)%in%BCNames)
col.T <- which(colnames(D)%in%TNames)

####################################################
i=which(D$Lon<0)
D$Lon[i]<- D$Lon[i]+ 360 #  rescale to 360 degrees to have all positive values for longitude
D$Lat<- D$Lat+ 90  #  rescale to 180 degrees to have all positive values for latitude
####################################################
hr=.1   #a guess at the overall harvest rate for target species 
#######################################################
# sum(BCweights) # weights should sum up to 1
# sum(Tweights) # weights should sum up to 1
# multiply numbers of bycatch species by weights
mysums.B <- D %>% summarise(across(.cols=all_of(BCNames),sum))
mysums.T <- D %>% summarise(across(.cols=all_of(TNames),sum))

tmp_wb<-BCweights/mysums.B
tmp_wt<-Tweights/mysums.T

rel_wb<-tmp_wb/sum(tmp_wb)
# sum(rel_wb)
rel_wt<-tmp_wt/sum(tmp_wt)
# sum(rel_wt)

# plot weighted and unweighted proportions
D_nw<-D # unweighted

for (i in 1:nrow(D)){
  D[i,col.BC]<- D[i,col.BC]*rel_wb #asign bycatch weights
  D[i,col.T]<- D[i,col.T]*rel_wt   #asign target weights
} 

pdf(file=paste0("Results/","Proportions",Fishery_name,".pdf"),width = 8,height = 6)
par(mfrow=c(1,2),oma=c(5,1,1,1),mar=c(1,4,2,0))
plot(as.numeric(D_nw %>% summarise(across(.cols=all_of(TNames),sum)))/sum(D_nw[,TNames]),xaxt='n',
     col="blue",pch=19,xlab="",ylab="")
points(as.numeric(D %>% summarise(across(.cols=all_of(TNames),sum)))/sum(D[,TNames]),col="red",pch=19)
axis(side = 1, at = c(1:NT),labels = c(TNames),las=2)
title(main="Target",line = 0.5)
legend("topright",col=c("blue","red"),pch=19,legend = c("Unweighted","Weighted"))

plot(as.numeric(D_nw %>% summarise(across(.cols=all_of(BCNames),sum)))/sum(D_nw[,BCNames]),xaxt='n',
     col="blue",pch=19,xlab="",ylab="")
points(as.numeric(D %>% summarise(across(.cols=all_of(BCNames),sum)))/sum(D[,BCNames]),col="red",pch=19)
axis(side = 1, at = c(1:NBC),labels = c(BCNames),las=2)
title(main="Bycatch",line = 0.5)
dev.off()

D <- D %>% 
  mutate(TBC.w=rowSums(.[BCNames]),# total weighted By-Catch
         Target.w=rowSums(.[TNames])) %>% # total weighted By-Catch
  filter(Target.w>0)# eliminate sets with target catch =0 to be able to calculate the proportion of Bycatch/target

###### correlation analysis 

pdf(file=paste0("Results/","Corr_plot_",Fishery_name,".pdf"))
ggcorrplot(round(cor(D[,-c(1:4)]), 1), hc.order = FALSE, type = "lower",lab=TRUE,
           outline.col = "white")
dev.off()

# create a map with all years combined
D_sum<-aggregate(x=D[,5:ncol(D)],by=list(Lat=D$Lat,Lon=D$Lon),FUN=sum)
D_sum$Prop<-D_sum$TBC.w/D_sum$Target.w # % of ByCatch to target species

Tmp<-melt(D_sum,id.vars = c("Lat","Lon"))
head(Tmp)
summary(Tmp)
str(Tmp)

# Maps
pdf(paste0("Results/","Maps_BCW_N",Fishery_name,".pdf"),width=10, height = 9,onefile=F)
Tmp %>% group_by(variable) %>%
  do(gg = {ggplot(., aes(Lon, Lat, fill = value)) +
      geom_tile() + facet_grid(~variable) +
      scale_fill_gradient(low = "green", high = "red")+
      theme(legend.position = "top")}) %>%
  .$gg %>% arrangeGrob(grobs = ., nrow = 3) %>% grid.arrange()
dev.off()

###############################################
GlobalTC<-sum(D[,col.T])
GlobalEffort<-sum(D$Effort)
###############################################
# this routine searchers over spatial locations of closed areas 
# to find the pattern of closed areas that will have the max reduction in: 
# 1) "N": bycatch in numbers; 2) "rate": bycatch CPUE; 3) "prop": ByCatch/target   
# it will not be a square, but the areas closest to a centroid if mosaic=F
# if Mosaic =T, closes quadrants where bycatch is minimized, they don't have to be conected
# By.month=T, closes one month instead of an area
FindClose<-function(NClose,D,minimize.by="prop",mosaic=F,by.month=FALSE){
  if (by.month==FALSE){
    
    col.BC<-seq(from=ncol(D)- NBC -1 ,to=ncol(D)-2)
    
    Narea<-nrow(D) # Number of areas: this is just the number of rows in the data frame
    
    Bycatch<-D %>% dplyr::select(BCNames) # weighted numbers
    #now loop over over centroids
    BestByCatch=0 #the target bycatch to beat
    if (NClose>0) {
      if(mosaic==F){
        for (a in 1:Narea){
          lat<-D$Lat[a]
          lon<-D$Lon[a]
          TempClosed<-array(dim=Narea,FALSE)
          if (mosaic==F){
            d<-sqrt((D$Lat-lat)^2+(D$Lon-lon)^2) 
            ord<-order(d)
            TempClosed[ord[1:NClose]]<-TRUE }# closed areas for each cell centroid
          
          i<-which(TempClosed==TRUE)
          
          if(minimize.by=="prop"){
            
            Tbycatch<-sum(Bycatch[i,])/sum(D$Target.w[i]) #proportion
          } else if (minimize.by=="N"){
            Tbycatch<-sum(Bycatch[i,])
          } else {
            Tbycatch<-sum(Bycatch[i,])/sum(D$Effort[i]) #CPUE
          }
          if (Tbycatch>BestByCatch){ # it only replaces the Tbycatch if it is higher than the BestByCatch calculated for the previous area [a]
            BestByCatch<-Tbycatch;Closed<-TempClosed;BestLat<-lat;BestLon<-lon}
           
        } #end areas
        
      } else {
        # closed area mosaic
        TempClosed<-array(dim=Narea,FALSE)
        if(minimize.by=="prop"){
          j<-order(D$TBC.w/D$Target.w,decreasing = TRUE)[1:NClose] # #proportion
        } else if (minimize.by=="N"){
          j<-order(D$TBC.w,decreasing = TRUE)[1:NClose] # numbers
        } else {
          j<-order(D$TBC.w/D$Effort,decreasing = TRUE)[1:NClose] #CPUE
        }
        # closed area mosaic
        TempClosed[j]<-TRUE
        Closed<-TempClosed
      }
    }
    else{
      Closed<-NA
      for (a in 1:Narea){Closed[a]=FALSE}
    }
    M<-NULL
  }
  if (by.month == TRUE){
    
     TempClosed<-array(dim=nrow(D),FALSE)
    
    if(minimize.by=="N"){
      Table.month<-aggregate(x=D$TBC.w,by=list(Month=D$Month),FUN=sum) # bycatch in numbers
    } else if (minimize.by=="prop"){
      Table.month<-aggregate(x=D$TBC.w/D$Target.w,by=list(Month=D$Month),FUN=sum)  # bycatch in numbers
    } else {
      Table.month<-aggregate(x=D$TBC.w/D$Effort,by=list(Month=D$Month),FUN=sum) 
    }
    M<-Table.month$Month[which(Table.month$x==max(Table.month$x))]
   
    M<-M[1]  # sometimes it doesn't work if there is the same number of bycatch in different months when minimized.by="N"
    i<-which(D$Month==M)
    TempClosed[i]<-TRUE
    Closed<-TempClosed
    
  }
  return(list(Closed,M)) # matrix of TRUES and FALSES
}

###################################################################################
#FishToTC=T , fishing to reach the same total target catch
#FishToTC=F , total taget catch can change but effort remains the same, 
# effort inside closed area move to open areas proportional to the effort already in those open quadrants 
#FishEfficiency=T  CPUE of target species can change in open areas based in arbitrary exploitation rate
#FishEfficiency=F  CPUE of target species is the same in each open quadrant 
Calculate<-function(Closed,D,FishToTC,FishEfficiency,hr,by.month) { #Closed (matrix)comes from the previous function
  
  CPUE<-D$Target.w/D$Effort # CPUE target species
  
  T_CPUE<-D %>% 
    rowwise() %>% 
    transmute(across(TNames,function(x)x/Effort)) %>% 
    data.frame
  i<-which(Closed[[1]]==FALSE) # open area or open months
  
  
  TotalE<-sum(D$Effort)
  Narea<-length(Closed[[1]]) # Total number of areas (grids)
  q<-Narea*(-log(1-hr)/TotalE)  #this is the q for the target species 
  InitialCatch<-D$Target.w 
  InitialB<-InitialCatch/(1-exp(-q*D$Effort))  #the initial biomass in each cell for the target species  
  
  
  
  BeforeBC <- D %>% summarise(across(.cols=all_of(BCNames),sum))
  
  if( FishToTC==FALSE ) {  # total catch can change but not total effort, sum(NewEffort)==sum(D$Effort)
    OldEffort<-sum(D$Effort[i]) # old effort in open areas
    if(by.month==TRUE){
      DisplacedEffort<-GlobalEffort-OldEffort # effort in closed months
    }
    
    if(by.month==FALSE){
      DisplacedEffort<-sum(D$Effort[-i])  # total catch for the target species should be the same, effort can change 
    }
    
    NewEffort<-array(dim=length(D$Effort),0)
    NewEffort[i]<-DisplacedEffort*D$Effort[i]/sum(D$Effort[i])+D$Effort[i] 
    
  }   else { #FishToTC==TRUE
    
    # fishing to reach same TC 
    if(by.month==TRUE){
      TC<-GlobalTC
    }
    
    if(by.month==FALSE){
      TC<-sum(D$Target.w) # catch for each target species should be the same, but effort can change 
    }
    
    CatchFromOpen<- sum(D$Target.w[i])  # catch in open areas 
   
    CatchFromOpen_byspecies<- D[i,] %>% summarise(across(.cols=all_of(TNames),sum))
    
    HookIncrease<-TC/CatchFromOpen 
    NewEffort<-array(dim=length(D$Effort),0)
    NewEffort[i]<-HookIncrease*D$Effort[i] # New effort outside 
    OldEffort<-sum(D$Effort[i]) # old effort outside 
    
  }  # end of fishing to TC
  NewCatch<-NewEffort*CPUE
  NewCatch_byspecies<-NewEffort*T_CPUE
  Tcatch<-sum(NewEffort*CPUE)  # if CPUE (fishing Efficiency) doesn't change  sum(D$target[i])/sum(D$Effort[i]) = = TC/sum(NewEffort)
  #######################################################################################
  #at this point we have NewEffort (check, NewEffort when FishToTC= TRUE or FALSE is the same)
  if (FishEfficiency == TRUE){ # CPUE can change 
    Catch<- D$Target.w   # catch by area before closure
    
    Catch_byspecies<- D %>% dplyr::select(TNames)   # catch by area before closure by species
    Abu<-Catch/(1-exp(-q*D$Effort)) #Abundance by area
    Abu_byspecies<-Catch_byspecies/(1-exp(-q*D$Effort))
    NewCatch<-Abu*(1-exp(-q*NewEffort))  #new catch by area
    NewCatch_byspecies<-Abu_byspecies*(1-exp(-q*NewEffort))
    
    Tcatch<-sum(NewCatch)
    if  (FishToTC==TRUE) { #CPUE will change do an incremental calculation
      for (step in 1:20) {  #iterate
        if(by.month==TRUE){
          TC<-GlobalTC
        }
        if(by.month==FALSE){
          TC<-sum(D$Target.w) # total catch for the target species should be the same, effort can change 
        }
        ratio<-sum(NewCatch)/TC
        NewEffort<-NewEffort/ratio
        NewCatch<-Abu*(1-exp(-NewEffort*q))  #new catch by area
        NewCatch_byspecies<-Abu_byspecies*(1-exp(-NewEffort*q))  #new catch by area
        Tcatch<-sum(NewCatch)
      } 
    }
  }
  
  
  BcCPUE<-D %>% 
    rowwise() %>% 
    transmute(across(BCNames,function(x)x/Effort)) %>% 
    data.frame
  
  BcByArea<-BcCPUE*NewEffort 
  
  TotalBcByArea <- BcByArea %>% transmute(rowSums(.)) %>% pull
  
  ByCatch <- BcByArea %>% summarise(across(BCNames,sum))
  
 
  New_TCatch_byspecies <- NewCatch_byspecies %>% summarise(across(TNames,sum))
  
  TEffort<-sum(NewEffort)
  
  return(list(Tcatch=Tcatch,ByCatch=ByCatch,TEffort=TEffort,
              NewCatch=NewCatch,New_TCatch_byspecies=New_TCatch_byspecies,
              NewEffort=NewEffort,TotalBcByArea=TotalBcByArea))
}

################################################
# Do calculations and plots
# ClosedSeq=seq(0,.5,.1) (sequence of proportions of areas to closed) can changed and 
# Months_closed=seq(1,5,1) (sequence of number of months to close) can also be changed
# If maps are desired Maps=TRUE
# tmp.lines=TRUE, add grey lines for each year closed for bycatch species in case of months=F
DoCalcs<-function(D,FishToTC,FishEfficiency,hr,ClosedSeq=seq(0,.5,.1),Months_closed=seq(1,5,1),
                  minimize.by="prop",mosaic=FALSE,by.month=FALSE,Maps=TRUE,tmp.lines=TRUE){  # do the simulation
  if (by.month==FALSE){
    Closed_months<-NULL
    D<-aggregate(x=D[,5:ncol(D)],by=list(Lat=D$Lat,Lon=D$Lon,Year=D$Year),FUN=sum)
    
    BC<-list() # to store total bycatch by species, area closure and by year
    TC<-list() # to store total target by species, area closure and by year
    CPUESave<-list() # to store target cpue by closure and year
    TEffort<-list() 
    CatchSave<- list()
    ResArea<-list()
    
    Nyears<-length(unique(D$Year))
    
    for(y in 1:c(Nyears+1)){ # +1: in the last element of the list we will store all years combined
      if(y==Nyears+1){
        D2<- aggregate(x=D[,4:ncol(D)],by=list(Lat=D$Lat,Lon=D$Lon),FUN=sum)
      } else {
        D2<-D[D$Year==unique(D$Year)[y],] 
      }
      N<-nrow(D2) #number of areas
      
      BC[[y]]<-matrix(NA,nrow = length(ClosedSeq),ncol=length(col.BC))
      TC[[y]]<-matrix(NA,nrow = length(ClosedSeq),ncol=length(col.T))
      CPUESave[[y]]<- matrix(NA,ncol=length(ClosedSeq),nrow = 1)# for target species 
      TEffort[[y]]<- matrix(NA,ncol=length(ClosedSeq),nrow = 1) #for effort 
      CatchSave[[y]]<- matrix(NA,ncol=length(ClosedSeq),nrow = 1)
      ResArea[[y]]<-data.frame(Lat=NA,Lon=NA,NewEffort=NA,CPUE_Target=NA,
                               CPUE_ByCatch= NA,Closure=NA)
      # 
      for (i in 1:length(ClosedSeq) )   { #looping over the proportion of areas closed
        pClose<-ClosedSeq[i] # prop closed
        NClose<-round(pClose*N,0) # numbers of areas closed
        
        Closed<-FindClose(NClose,D=D2,minimize.by=minimize.by,mosaic = mosaic,by.month=FALSE)
        
        Res<-Calculate(Closed,D=D2,FishToTC,FishEfficiency,hr,by.month=FALSE)
        Res_a<-data.frame(Lat=D2$Lat,Lon=D2$Lon,NewEffort=Res$NewEffort,
                          CPUE_Target=Res$NewCatch/(Res$NewEffort*1E3),
                          CPUE_ByCatch= Res$TotalBcByArea/(Res$NewEffort),
                          Closure=ClosedSeq[i])
        BC[[y]][i,]<-as.matrix(Res$ByCatch) # weighted bycatch
        TC[[y]][i,]<-as.matrix(Res$New_TCatch_byspecies)
        CPUESave[[y]][i]<-Res$Tcatch/Res$TEffort
        TEffort[[y]][i]<-Res$TEffort
        CatchSave[[y]][i]<-Res$Tcatch
        ResArea[[y]]<-rbind(ResArea[[y]],Res_a)
      }
      ResArea[[y]]<-ResArea[[y]][-1,]
      
    }
    TBC<-apply(BC[[Nyears+1]],1,FUN=sum)
    BCScaled<-TBC/TBC[1] #relative to the first value (no closure)
    CPUEScaled<-CPUESave[[Nyears+1]]/CPUESave[[Nyears+1]][1]
    TEffortScaled<-TEffort[[Nyears+1]]/TEffort[[Nyears+1]][1]
    CatchScaled<-CatchSave[[Nyears+1]]/CatchSave[[Nyears+1]][1]
    
    # what if we sum up when areas were closed by year: 
    
    TBC_y<-0
    TEffort_y<-0
    CatchSave_y<-0
    TBC_tmp<-list()
    for(y in 1:Nyears){
      TBC_tmp[[y]]<-apply(BC[[y]],1,FUN=sum,na.rm=T)
      TBC_y<-TBC_y+TBC_tmp[[y]]
      TEffort_y<-TEffort_y+TEffort[[y]]
      CatchSave_y<-CatchSave_y+CatchSave[[y]]
      TBC_tmp[[y]]<-TBC_tmp[[y]]/TBC_tmp[[y]][1]
    }
    
    BCScaled_y<-TBC_y/TBC_y[1] #relative to the first value (no closure)
    CPUE_y<-CatchSave_y/TEffort_y
    CPUEScaled_y<-CPUE_y/CPUE_y[1]
    CatchScaled_y<-CatchSave_y/CatchSave_y[1]
    TEffortScaled_y <- TEffort_y/TEffort_y[1]
    
    ######### for barplot
    #################################################################################
    # data needed comes from BC, list where each element is a year, and the last one is the total (stationary closure)
    # rows are the closures and columns the species
    nPol<-length(BC[[Nyears+1]][,1])
    nSpec<-length(BC[[Nyears+1]][1,])+length(TC[[Nyears+1]][1,])
    BC_D<-matrix(0,nrow = nPol,ncol=nSpec) 
    for(i in 1:Nyears){
      Total<-cbind(TC[[i]],BC[[i]])
      BC_D<-BC_D  + Total 
    }
    
    BC_S<-cbind(TC[[Nyears+1]],BC[[Nyears+1]])
    BC_S_rel<-BC_S[-1,]
    BC_D_rel<-BC_D[-1,]
    for(i in 1:c(nPol-1)){
      BC_S_rel[i,]<-BC_S_rel[i,]/BC_S[1,]
      BC_D_rel[i,]<-BC_D_rel[i,]/BC_D[1,]
    }
    
    BC_D_rel<-as.data.frame(BC_D_rel); BC_S_rel<-as.data.frame(BC_S_rel)
    names(BC_S_rel)<-names(BC_D_rel)<-c(TNames,BCNames)
    BC_S_rel$Closure<-as.factor(ClosedSeq[-1]);BC_D_rel$Closure<-as.factor(ClosedSeq[-1])
    BC_S_rel$Area<-as.factor("Stationary");BC_D_rel$Area<-as.factor("Mobile")
    BC_tot<-rbind(BC_S_rel,BC_D_rel)
    BC_tot<-melt(BC_tot,value.name = "BCN")
    BC_tot$Relative_ByCatch<-BC_tot$BCN-1
    names(BC_tot)[3]<-"Species"
    tmp<-data.frame(Category="Target",Species=TNames)
    tmp<-rbind(tmp,data.frame(Category="ByCatch",Species=BCNames))
    BC_tot<-merge(BC_tot,tmp,by="Species",all = T)
    #########################################################################################
  } 
  if(by.month==TRUE){
    D2<-aggregate(x=D[,5:ncol(D)],by=list(Lat=D$Lat,Lon=D$Lon,Month=D$Month),FUN=sum)
    
    BC<-list() # to store total bycatch by species
    TC<-list() # to store total target by species, area closure and by year
    CPUESave<-list() # to store target cpue 
    TEffort<-list() # to store effort
    CatchSave<- list() # to store catch 
    #ResArea<-list()
    
    Closed<-list(array(dim=nrow(D2),FALSE),NULL)
    
    Res<-Calculate(Closed,D=D2,FishToTC,FishEfficiency,hr,by.month = TRUE) #when no closure
    BC[[1]]<-as.data.frame(Res$ByCatch) # weighted bycatch 
    TC[[1]]<-as.data.frame(Res$New_TCatch_byspecies) # weighted target species 
    CPUESave[[1]]<-Res$Tcatch/Res$TEffort
    TEffort[[1]]<-Res$TEffort
    CatchSave[[1]]<-Res$Tcatch
    Closed_months<-NULL
    # we can close 1 month to 5
    Months_closed<-Months_closed
    Nmonths<-length(Months_closed)
    
    for (i in 2:c(Nmonths+1)) { #looping over months
      
      Closed<-FindClose(NClose,D=D2,minimize.by=minimize.by,mosaic = mosaic,by.month=TRUE) #NClose does not matter for months
      Closed_months<-c(Closed_months,Closed[[2]])
      Res<-Calculate(Closed,D=D2,FishToTC,FishEfficiency,hr,by.month=TRUE)
      
      BC[[i]]<-as.data.frame(Res$ByCatch) # weighted bycatch 
      TC[[i]]<-as.data.frame(Res$New_TCatch_byspecies)
      CPUESave[[i]]<-Res$Tcatch/Res$TEffort
      TEffort[[i]]<-Res$TEffort
      CatchSave[[i]]<-Res$Tcatch
      
      #remove the month with the higest bycatch to calculate the next one
      D2<-D2[D2$Month!=Closed[[2]],]
    }
   
    TBC<-NULL
    CPUEScaled<-NULL
    TEffortScaled<-NULL
    CatchScaled<-NULL
    for(i in 1:length(BC))  {
      tmp<-c(apply(BC[[i]],1,FUN=sum))
      TBC<-c(TBC,tmp)
      tmp2<-CPUESave[[i]]/CPUESave[[1]]
      CPUEScaled<-c(CPUEScaled,tmp2)
      tmp3<-TEffort[[i]]/TEffort[[1]]
      TEffortScaled<-c(TEffortScaled,tmp3)
      tmp4<-CatchSave[[i]]/CatchSave[[1]]
      CatchScaled<-c(CatchScaled,tmp4)
    }
    BCScaled<-TBC/TBC[1] #relative to the first value (no closure)
    
    ######### for barplot
    #################################################################################
    # data needed comes from BC, list where each element is a year, and the last one is the total (stationary closure)
    # rows are the closures and columns the species
    nPol<-length(BC)
    nSpec<-ncol(BC[[1]])+ncol(TC[[1]])
    BC_S<-matrix(0,nrow = nPol,ncol=nSpec) 
    for(i in 1:c(Nmonths+1)){
      BC_S[i,]<- t(rbind(t(TC[[i]]),t(BC[[i]])))
    }
    rel<-BC_S[1,]
    BC_S_rel<-BC_S[-1,]
    
    for(i in 1:c(nPol-1)){
      BC_S_rel[i,]<-BC_S_rel[i,]/rel
    }
    
    BC_S_rel<-as.data.frame(BC_S_rel)
    colnames(BC_S_rel)<-c(TNames,BCNames)
    BC_S_rel$Closure<-as.factor(Months_closed)
    BC_S_rel$Month<-as.factor("Stationary")
    BC_tot<-melt(BC_S_rel,value.name = "BCN")
    BC_tot$Relative_ByCatch<-BC_tot$BCN-1 
    names(BC_tot)[3]<-"Species"
    tmp<-data.frame(Category="Target",Species=TNames)
    tmp<-rbind(tmp,data.frame(Category="ByCatch",Species=BCNames))
    BC_tot<-merge(BC_tot,tmp,by="Species",all = T)
    
    ### NULL, only outputs when year closures
    BCScaled_y=NULL;TEffortScaled_y=NULL;CatchScaled_y=NULL;CPUEScaled_y=NULL
  }
  
  ##############################################################
  par(mfrow=c(1,1),mar=(c(5, 4, 4, 2) + 0.1),oma=c(1,1,1,1))
  if (by.month==FALSE){
    plot(ClosedSeq,ClosedSeq,col="#BB000099",ylim=c(0,1.6),xlab="% of areas closed",type="n",
         ylab="Relative Amount",cex.main=0.8,
         main=paste("Mosaic=",mosaic,"Minimized by=",minimize.by,",", "Fish to TC =",FishToTC,",","Fishing Efficiency changes =",FishEfficiency))
    
    #   
    lines(ClosedSeq,BCScaled,col="#BB000099",lwd=3)
    lines(ClosedSeq,CatchScaled,lwd=3,col="#80008099",type="b",pch=16) # purple
    lines(ClosedSeq,CPUEScaled,lwd=3,col="#458B0099") #green
    lines(ClosedSeq,TEffortScaled,lwd=3,col="#0000FF99") #blue
    
    legend("bottomleft",pch=c(-1,16,-1,-1),bty = "n",
           legend=c("By-Catch","Target-Catch","Fishing Efficiency changes","Effort"),
           col=c("#BB000099","#80008099","#458B0099","#0000FF99"),lwd=2,cex=0.8)
    
    lines(ClosedSeq,BCScaled_y,col="#BB000099",lwd=3,lty=2)
    lines(ClosedSeq,CatchScaled_y,lwd=3,col="#80008099",type="b",pch=16,lty=2) # purple
    lines(ClosedSeq,CPUEScaled_y,lwd=3,col="#458B0099",lty=2) #green
    lines(ClosedSeq,TEffortScaled_y,lwd=3,col="#0000FF99",lty=2) #blue
    
    legend("topright",bty = "n",
           legend=c("Stationary closure","Mobile closure"),
           col="grey30",lty=c(1,2),cex=0.8,lwd=2)
    if(tmp.lines==TRUE){
      for(i in 1:Nyears){
        lines(ClosedSeq,TBC_tmp[[i]],lwd=1,col="#00000020",lty=1)
      }
    } 
  }
  
  if (by.month==TRUE){
    plot(c(0,Months_closed),c(0,Months_closed),col="#BB000099",ylim=c(0,1.6),xlab="ID Months closed",type="n",
         ylab="Relative Amount",cex.main=0.8,xaxt='n',
         main=paste("Minimized by=",minimize.by,",", "Fish to TC =",FishToTC,",","Fishing Efficiency changes =",FishEfficiency))
     axis(side = 1,at = Months_closed,labels = Closed_months)
       
     lines(c(0,Months_closed),BCScaled,col="#BB000099",lwd=3)
    lines(c(0,Months_closed),CatchScaled,lwd=3,col="#80008099",type="b",pch=16) # purple
    lines(c(0,Months_closed),CPUEScaled,lwd=3,col="#458B0099") #green
    lines(c(0,Months_closed),TEffortScaled,lwd=3,col="#0000FF99") #blue

    legend("bottomleft",pch=c(-1,16,-1,-1),bty = "n",
           legend=c("By-Catch","Target-Catch","Fishing Efficiency changes","Effort"),
           col=c("#BB000099","#80008099","#458B0099","#0000FF99"),lwd=2,cex=0.8)

    # mydat <- data.frame(Months_closed=c(0,Months_closed),Closed_months=c(0,Closed_months),CatchScaled,CPUEScaled,TEffortScaled) %>% 
    #   inner_join(data.frame(var=names(BCScaled),values=data.frame(BCScaled)) %>% 
    #   group_by(var) %>% 
    #   mutate(Months_closed=c(0,Months_closed)) %>% 
    #   spread(var,BCScaled) %>% 
    #   gather(var,values,-Months_closed))
    # 
    # plot(mydat$Months_closed,mydat$Months_closed,col="#BB000099",ylim=c(0,1.6),xlab="ID Months closed",type="n",
    #      ylab="Relative Amount",cex.main=0.8,xaxt='n',
    #      main=paste("Minimized by=",minimize.by,",", "Fish to TC =",FishToTC,",","Fishing Efficiency changes =",FishEfficiency))
    # axis(side = 1,at = Months_closed,labels = Closed_months)
    # #   
    # lines(c(0,Months_closed),BCScaled,col="#BB000099",lwd=3)
    # lines(c(0,Months_closed),CatchScaled,lwd=3,col="#80008099",type="b",pch=16) # purple
    # lines(c(0,Months_closed),CPUEScaled,lwd=3,col="#458B0099") #green
    # lines(c(0,Months_closed),TEffortScaled,lwd=3,col="#0000FF99") #blue
    # 
     legend("bottomleft",pch=c(-1,16,-1,-1),bty = "n",
            legend=c("By-Catch","Target-Catch","Fishing Efficiency changes","Effort"),
            col=c("#BB000099","#80008099","#458B0099","#0000FF99"),lwd=2,cex=0.8)
    
    
  }
  # barplot with ggplot
  if(by.month==FALSE){
    g<-ggplot(BC_tot, aes(x=Species, y=Relative_ByCatch, fill=Category, alpha=Area)) +
      geom_bar(stat="identity",position=position_dodge())+theme_minimal()+
      facet_wrap(~Closure,ncol=1)+ theme(legend.position="top",axis.text.x = element_text(angle = 90)) +
      scale_fill_manual(values = c("mediumorchid3", "turquoise3"))+
      scale_alpha_manual(values = c(0.5, 1))+ ylab("Relative Catch")
    print(g)
    
  }
  if(by.month==TRUE){
    g<-ggplot(BC_tot, aes(x=Species, y=Relative_ByCatch, fill=Category)) +
      geom_bar(stat="identity",position=position_dodge())+theme_minimal()+
      facet_wrap(~Closure,ncol=1)+ theme(legend.position="top",axis.text.x = element_text(angle = 90)) +
      scale_fill_manual(values = c("mediumorchid3", "turquoise3")) + ylab("Relative Catch")
    print(g)
  }
  
  # ############################################
  if (Maps==TRUE){
    N<-length(ResArea)
    for(y in N:1){
      ##################### Mapping
      ResArea[[y]][ResArea[[y]]$NewEffort==0,"NewEffort"]<-NA
      Tmp<-melt(ResArea[[y]],id.vars = c("Lat","Lon","Closure"))
      Tmp$Closure<-as.factor(Tmp$Closure)
      #
      g<-ggplot(Tmp[Tmp$variable=="NewEffort",], aes(Lon, Lat, fill = value)) +
        geom_tile() + facet_grid(Closure~variable) +
        scale_fill_gradient(low = "green", high = "red")+
        theme(legend.title = element_blank(),legend.position = "top",legend.text = element_text(size = 8),strip.text.y = element_blank())
      gg<-ggplot(Tmp[Tmp$variable=="CPUE_Target",], aes(Lon, Lat, fill = value)) +
        geom_tile() + facet_grid(Closure~variable) +
        scale_fill_gradient(low = "green", high = "red")+
        theme(legend.title = element_blank(),legend.position = "top",axis.ticks = element_blank(),axis.title.y = element_blank(),axis.text.y=element_blank(),legend.text = element_text(size = 8),strip.text.y = element_blank())
      ggg<-ggplot(Tmp[Tmp$variable=="CPUE_ByCatch",], aes(Lon, Lat, fill = value)) +
        geom_tile() + facet_grid(Closure~variable) +
        scale_fill_gradient(low = "green", high = "red")+
        theme(legend.title = element_blank(),legend.position = "top",axis.ticks = element_blank(),axis.title.y = element_blank(),axis.text.y=element_blank(),legend.text = element_text(size = 8))
      if (by.month==FALSE){
        grid.arrange(g,gg,ggg,ncol=3,top=paste("Year",unique(D$Year)[y],",","Mosaic=",mosaic,"Minimized by=",minimize.by,",", "Fish to TC =",FishToTC,",","Fishing Efficiency changes =",FishEfficiency))
      } else {
        grid.arrange(g,gg,ggg,ncol=3,top=paste("Month",unique(D$Month)[y],",","Mosaic=",mosaic,"Minimized by=",minimize.by,",", "Fish to TC =",FishToTC,",","Fishing Efficiency changes =",FishEfficiency))
      }
    }
  }
  saveRDS(object=list(Closed_Seq=ClosedSeq,Closed_months=Closed_months,# ResArea,BC=BC,CPUE_t=CPUESave,TEffort=TEffort,Catch=CatchSave,
                      BCScaled=BCScaled,TEffortScaled=TEffortScaled,CatchScaled=CatchScaled,CPUEScaled=CPUEScaled,
                      BCScaled_y=BCScaled_y,TEffortScaled_y=TEffortScaled_y,CatchScaled_y=CatchScaled_y,CPUEScaled_y=CPUEScaled_y,
                      Changes_bySpecies=BC_tot,ByArea=Tmp),
          file=paste0("Results/",Fishery_name,"_Mosaic=",mosaic,"_FishToTC=",FishToTC,"_FishEfficiency=",FishEfficiency,"_minimize.by=",minimize.by,
                      "_by.month=",by.month,".rds"))
} #end of function
############################################################################################################
# the following set of lines of codes run all posible combinations ##############
#################################################################################
pdf(file=paste0("Results/","Outs_prop_",Fishery_name,".pdf"))
DoCalcs(D,FishToTC=TRUE,FishEfficiency=FALSE,hr,minimize.by="prop",Maps=FALSE,by.month=TRUE,Months_closed=seq(1,5,1))
DoCalcs(D,FishToTC=TRUE,FishEfficiency=FALSE,hr,minimize.by="prop",Maps=FALSE,by.month=FALSE)
DoCalcs(D,FishToTC=TRUE,FishEfficiency=TRUE,hr,minimize.by="prop",Maps=FALSE,by.month=TRUE)
DoCalcs(D,FishToTC=TRUE,FishEfficiency=TRUE,hr,minimize.by="prop",Maps=FALSE,by.month=FALSE)
DoCalcs(D,FishToTC=FALSE,FishEfficiency=FALSE,hr,minimize.by="prop",Maps=FALSE,by.month=TRUE)
DoCalcs(D,FishToTC=FALSE,FishEfficiency=TRUE,hr,minimize.by="prop",Maps=FALSE,by.month=FALSE)
DoCalcs(D,FishToTC=FALSE,FishEfficiency=TRUE,hr,minimize.by="prop",Maps=FALSE,by.month=TRUE)
DoCalcs(D,FishToTC=FALSE,FishEfficiency=FALSE,hr,minimize.by="prop",Maps=TRUE,by.month=FALSE)
DoCalcs(D,FishToTC=TRUE,FishEfficiency=FALSE,hr,minimize.by="prop",mosaic=TRUE,Maps=FALSE)
DoCalcs(D,FishToTC=TRUE,FishEfficiency=TRUE,hr,minimize.by="prop",mosaic=TRUE,Maps=FALSE)
DoCalcs(D,FishToTC=FALSE,FishEfficiency=FALSE,hr,minimize.by="prop",mosaic=TRUE,Maps=FALSE)
DoCalcs(D,FishToTC=FALSE,FishEfficiency=TRUE,hr,minimize.by="prop",mosaic=TRUE,Maps=TRUE)
dev.off()
###############################################################################
pdf(file=paste0("Results/","Outs_N_",Fishery_name,".pdf"))
DoCalcs(D,FishToTC=TRUE,FishEfficiency=FALSE,hr,minimize.by="N",Maps=FALSE,by.month=TRUE)
DoCalcs(D,FishToTC=TRUE,FishEfficiency=FALSE,hr,minimize.by="N",Maps=FALSE,by.month=FALSE)
DoCalcs(D,FishToTC=TRUE,FishEfficiency=TRUE,hr,minimize.by="N",Maps=FALSE,by.month=TRUE)
DoCalcs(D,FishToTC=TRUE,FishEfficiency=TRUE,hr,minimize.by="N",Maps=FALSE,by.month=FALSE)
DoCalcs(D,FishToTC=FALSE,FishEfficiency=FALSE,hr,minimize.by="N",Maps=FALSE,by.month=TRUE)
DoCalcs(D,FishToTC=FALSE,FishEfficiency=TRUE,hr,minimize.by="N",Maps=FALSE,by.month=FALSE)
DoCalcs(D,FishToTC=FALSE,FishEfficiency=TRUE,hr,minimize.by="N",Maps=FALSE,by.month=TRUE)
DoCalcs(D,FishToTC=FALSE,FishEfficiency=FALSE,hr,minimize.by="N",Maps=TRUE,by.month=FALSE)
DoCalcs(D,FishToTC=TRUE,FishEfficiency=FALSE,hr,minimize.by="N",mosaic=TRUE,Maps=FALSE)
DoCalcs(D,FishToTC=TRUE,FishEfficiency=TRUE,hr,minimize.by="N",mosaic=TRUE,Maps=FALSE)
DoCalcs(D,FishToTC=FALSE,FishEfficiency=FALSE,hr,minimize.by="N",mosaic=TRUE,Maps=FALSE)
DoCalcs(D,FishToTC=FALSE,FishEfficiency=TRUE,hr,minimize.by="N",mosaic=TRUE,Maps=TRUE)
dev.off()
###############################################################################
pdf(file=paste0("Results/","Outs_BCrate_",Fishery_name,".pdf"))
DoCalcs(D,FishToTC=TRUE,FishEfficiency=FALSE,hr,minimize.by="rate",Maps=FALSE,by.month=TRUE)
DoCalcs(D,FishToTC=TRUE,FishEfficiency=FALSE,hr,minimize.by="rate",Maps=FALSE,by.month=FALSE)
DoCalcs(D,FishToTC=TRUE,FishEfficiency=TRUE,hr,minimize.by="rate",Maps=FALSE,by.month=TRUE)
DoCalcs(D,FishToTC=TRUE,FishEfficiency=TRUE,hr,minimize.by="rate",Maps=FALSE,by.month=FALSE)
DoCalcs(D,FishToTC=FALSE,FishEfficiency=FALSE,hr,minimize.by="rate",Maps=FALSE,by.month=TRUE)
DoCalcs(D,FishToTC=FALSE,FishEfficiency=TRUE,hr,minimize.by="rate",Maps=FALSE,by.month=FALSE)
DoCalcs(D,FishToTC=FALSE,FishEfficiency=TRUE,hr,minimize.by="rate",Maps=FALSE,by.month=TRUE)
DoCalcs(D,FishToTC=FALSE,FishEfficiency=FALSE,hr,minimize.by="rate",Maps=TRUE,by.month=FALSE)
DoCalcs(D,FishToTC=TRUE,FishEfficiency=FALSE,hr,minimize.by="rate",mosaic=TRUE,Maps=FALSE)
DoCalcs(D,FishToTC=TRUE,FishEfficiency=TRUE,hr,minimize.by="rate",mosaic=TRUE,Maps=FALSE)
DoCalcs(D,FishToTC=FALSE,FishEfficiency=FALSE,hr,minimize.by="rate",mosaic=TRUE,Maps=FALSE)
DoCalcs(D,FishToTC=FALSE,FishEfficiency=TRUE,hr,minimize.by="rate",mosaic=TRUE,Maps=TRUE)
dev.off()
