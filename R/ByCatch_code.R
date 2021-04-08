############## Instructions ####################################################
# Format of data entries:                                                      #
# There are 2 csv files that are needed:                                       #
# one called Data.csv and another one called Weights.csv                       #
# 1. data.csv: the first 5 columns are the same for all case studies:          # 
# Lat	Lon	Year Month Effort                                                    #
# The following columns may vary for each case study                           #
# put first the target species and then the bycatch by species                 #
# you can include any number of target and bycatch species                     #
# or group of species as your convenience.                                     #
# IMPORTANT: each row in the data file corresponds to a quadrant (lat&long)    #
# by year and month. There is only one record for each of these combinations   # 
# The numbers in each column for target and by-catch species can be            #
# in numbers or biomass but it has to be in the same units for each group.     #
# For example: all target species in weight and all by-catch species in numbers#
# 2. Weights.csv: this csv has in the first row the names of the target        #
# species in the same order as in the Data.csv file and the second row         #
# is the weight for each target species. They should sum up to 1.              #
# The third row is the names of the by-catch species                           #
# in the same order as in the Data.csv file. The fourth row is the weight      #
# for each by-catch species and they also should sum up to 1.                  #
# Outputs:                                                                     #
# Running the code generates:                                                  #
# 1. Plot with original proportions of catch by species to the total catch     #
#    in blue, and proportions based in specified weights in red                #
# 2. Correlation plot among Effort, target and by-catch species,               #
#    including totals (weighted and unweighted)                                #    
# 3. Maps with original data                                                   #
# 4. running the function DoCalcs generates:                                   #
# a) lines plot: changes in by-catch, target catch, fishing efficiency and     # 
#   Effort for each % of area closure or number of months depending on         #
#   whether by.month=TRUE or FALSE. This is represented by the solid lines     # 
#   If by.month=F, the dashed lines represent the dynamic closures,            # 
#   which means closing a different area each year. In this case, each gray    #
#   line represents a year (this can be turned off if tmp.lines=FALSE).        #
# b) barplot: showing relative changes in catch for each species               #
# c) maps: showing area closed in gray and effort, target and by-catch CPUE    #
# d) RDS output with results to be use for other plots and analysis of all     #
#    case studies combined. There is no raw data in the outputs                #
# All Outputs (plots and rds) are saved in a folder called "Results"           #
################################################################################
# Default names and locations
defaults = c(Fishery_name="Your_fishery_name", # Set your fishery name here! only change this entry (not the others)
             results_dir="Results", # default
             data_fn = "Data.csv", # default
             weights_fn = "Weights.csv") # default

# Set variables to default values if present in workspace
for (vn in names(defaults)) {
  if (!exists(vn)) # Only set to default if not already in workspace
    assign(vn,defaults[vn])
}

# Remove all but these parameter values
vns = ls()
vns = vns[!(vns %in% c("defaults",names(defaults)))]
rm(list=vns)

############################################
# if you don't have the following libraries needed to run the code do:
# install.packages("name of the library") 
# examples:
# install.packages(tidyverse)
# install.packages(reshape2)
# install.packages(grid)
# install.packages(gridExtra)
# install.packages(ggcorrplot)
# install.packages(dplyr)
# install.packages(maps)

####### Libraries needed
library(tidyverse)
library(reshape2)
library(grid)
library(gridExtra)
library(ggcorrplot)
library(dplyr) #make shure you have the last version of this package or you might have problems with the function accross not being recognized 
library(maps)

###########################################################################
# Read data and weights csv files  ########################################
###########################################################################
dir.create(results_dir,showWarnings = FALSE)    # create folder called "Results" in your directory
###########################################################################
# reading data
D <- read.csv(data_fn, header = T) 
D <- D[complete.cases(D),]
# reading weights
Tweights <- read.csv(weights_fn, nrows = 1) %>%
  select_if(~ !any(is.na(.))) 

BCweights <- read.csv(weights_fn, nrows = 2, skip = 2) %>%
  select_if(~ !any(is.na(.))) 

# sum(BCweights) # weights should sum up to 1
# sum(Tweights) # weights should sum up to 1

NT <- ncol(Tweights) # number of target species 
NBC <- ncol(BCweights) # number of by-catch species

### Effort in thousands (any unit)
D$Effort <- D$Effort/1000

BCNames <- colnames(BCweights) # name of the by-catch species
TNames <- colnames(Tweights)   # name of the target species

col.BC <- which(colnames(D)%in%BCNames)  # columns where the by-catch species are
col.T <- which(colnames(D)%in%TNames)    # columns where the target species are

####################################################
i = which(D$Lon<0)
D$Lon[i] <- D$Lon[i] + 360 #  re-scale to 360 degrees to have all positive values for longitude
D$Lat <- D$Lat + 90       #  re-scale to 180 degrees to have all positive values for latitude
####################################################
hr=.1   # a guess at the overall harvest rate for target species 
#######################################################
# sum all by-catch and target catch 
mysums.B <- D %>% summarise(across(.cols=all_of(BCNames),sum))
mysums.T <- D %>% summarise(across(.cols=all_of(TNames),sum))
# If you have problems with the function across not being recognized, please re-install the package 'dplyr' and re-start R. Thanks 

tmp_wb <- BCweights/mysums.B
tmp_wt <- Tweights/mysums.T
rel_wb <- tmp_wb/sum(tmp_wb)
rel_wt <- tmp_wt/sum(tmp_wt)

# sum(rel_wb)
# sum(rel_wt)

# plot weighted and un-weighted proportions
D_w <- D # keep D un-weighted

for (i in 1:nrow(D_w)){
  D_w[i,col.BC] <- D_w[i,col.BC]*rel_wb # assign by-catch weights
  D_w[i,col.T] <- D_w[i,col.T]*rel_wt   # assign target weights
} 

pdf(file=file.path(results_dir,paste0("Proportions",Fishery_name,".pdf")),width = 8,height = 6)
par(mfrow=c(1,2),oma=c(5,1,1,1),mar=c(1,4,2,0))
plot(as.numeric(D %>% summarise(across(.cols=all_of(TNames),sum)))/sum(D[,TNames]),xaxt='n',
     col="blue",pch=19,xlab="",ylab="Relative proportion")
points(as.numeric(D_w %>% summarise(across(.cols=all_of(TNames),sum)))/sum(D_w[,TNames]),col="red",pch=19)
axis(side = 1, at = c(1:NT),labels = c(TNames),las=2)
title(main="Target",line = 0.5)
legend("topright",col=c("blue","red"),pch=19,legend = c("Un-weighted","Weighted"))

plot(as.numeric(D %>% summarise(across(.cols=all_of(BCNames),sum)))/sum(D[,BCNames]),xaxt='n',
     col="blue",pch=19,xlab="",ylab="")
points(as.numeric(D_w %>% summarise(across(.cols=all_of(BCNames),sum)))/sum(D_w[,BCNames]),col="red",pch=19)
axis(side = 1, at = c(1:NBC),labels = c(BCNames),las=2)
title(main="Bycatch",line = 0.5)
dev.off()

D_w <- D_w %>% 
  mutate(TBC.w=rowSums(.[BCNames]), # total weighted By-Catch
         Target.w=rowSums(.[TNames])) %>% # total weighted Catch
  filter(Target.w>0)# eliminate sets with target catch = 0 to be able to calculate the proportion of By-catch/target
########################################################################
D <- D %>% 
  mutate(TBC=rowSums(.[BCNames]), # total weighted By-Catch
         Target=rowSums(.[TNames])) %>% # total weighted Catch
  filter(Target>0)
# add to D weighted columns for total by-catch and total target catch
D$TBC.w <- D_w$TBC.w
D$Target.w <- D_w$Target.w
###### correlation plot
pdf(file=file.path(results_dir,paste0("Corr_plot_",Fishery_name,".pdf")))
ggcorrplot(round(cor(D[,-c(1:4)]), 1), hc.order = FALSE, type = "lower",lab=TRUE,
           outline.col = "white")
dev.off()

# create a map with all years combined
D_sum <- aggregate(x=D[,5:ncol(D)],by=list(Lat=D$Lat,Lon=D$Lon),FUN=sum)
D_sum$Prop <- D_sum$TBC/D_sum$Target # proportion of By-Catch to target species
D_sum$Prop.w <- D_sum$TBC.w/D_sum$Target.w # proportion of By-Catch to target species weighted

Tmp<-melt(D_sum,id.vars = c("Lat","Lon"))
head(Tmp)
summary(Tmp)
str(Tmp)

# Maps
world <- map_data("world")
head(world)
i = which(world$long<0)
world$long[i] <- world$long[i]+ 360 #  re-scale to 360 degrees to have all positive values for longitude
world$lat <- world$lat+ 90  #  re-scale to 180 degrees to have all positive values for latitude

world.cut <- world[world$lat > min(D$Lat)-5 & world$lat < max(D$Lat)+5 ,]
world.cut <- world.cut[world.cut$long > min(D$Lon)-5 & world.cut$long < max(D$Lon)+5 ,]

# This creates a quick and dirty world map - playing around with the themes, aesthetics, and device dimensions is recommended!
worldmap <- ggplot(world.cut, aes(x=long, y=lat)) +
  geom_polygon(aes(group=group)) +
  theme(panel.background = element_rect(fill="skyblue", colour="skyblue"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  coord_equal()

pdf(file.path(results_dir,paste0("Maps_data",Fishery_name,".pdf")),width=10, height = 9,onefile=F)
Tmp %>% group_by(variable) %>%
  do(gg = {
    worldmap +
      geom_tile(., mapping=aes(Lon, Lat, fill = value)) +
      facet_grid(~variable) +
      scale_fill_gradient(low = "white", high = "red")+
      theme(panel.background = element_rect(fill="skyblue", colour="skyblue"),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),legend.position = "top")
  }) %>%
  .$gg %>% arrangeGrob(grobs = ., nrow = 3) %>% grid.arrange()
dev.off()

###############################################
GlobalTC <- sum(D[,col.T])  # total Catch for target species 
GlobalEffort <- sum(D$Effort) # total Effort 
###############################################
# the next function searchers over spatial locations of closed areas 
# to find the pattern of closed areas that will have the max reduction in: 
# 1) "N": by-catch in numbers; 2) "rate": by-catch CPUE; 3) "prop": By-Catch/target   
# it will not be a square, but the areas closest to a centroid if mosaic=F
# if Mosaic = T, quadrants where by-catch is minimized can be closed independently, they don't have to be connected
# if By.month = T, one to five months are closed instead of an area

FindClose <- function(NClose, D, weights=TRUE, minimize.by="prop", mosaic=F, by.month=FALSE){
  if (weights == TRUE){
    Bycatch <- D$TBC.w
    Target <- D$Target.w
  } else {
    Bycatch <- D$TBC
    Target <- D$Target
  }
  if (by.month==FALSE){
    
    Narea <- nrow(D) # Number of areas: this is just the number of rows in the data frame
    
    # now loop over centroids
    BestByCatch = 0 # the target by-catch to beat
    if (NClose>0) {
      if(mosaic==F){
        for (a in 1:Narea){
          lat <- D$Lat[a]
          lon <- D$Lon[a]
          TempClosed<-array(dim=Narea,FALSE)
          if (mosaic==F){
            d <- sqrt((D$Lat-lat)^2+(D$Lon-lon)^2) 
            ord <- order(d)
            TempClosed[ord[1:NClose]] <- TRUE } # closed areas for each cell centroid
          
          i<-which(TempClosed==TRUE)
          
          if(minimize.by=="prop"){
            Tbycatch <- sum(Bycatch[i])/sum(Target[i]) # proportion
            
          } else if (minimize.by=="N"){
            Tbycatch <- sum(Bycatch[i])  # Numbers
            
          } else {
            Tbycatch <- sum(Bycatch[i])/sum(D$Effort[i]) # CPUE
          }
          
          if (Tbycatch > BestByCatch){ # it only replaces the Tbycatch if it is higher than the BestByCatch calculated for the previous area [a]
            BestByCatch <- Tbycatch; Closed <- TempClosed; BestLat <- lat; BestLon <- lon}
          
        } # end areas
        
      } else {
        # closed area mosaic
        TempClosed <- array(dim=Narea,FALSE)
        if(minimize.by == "prop"){
          j <- order(Bycatch/Target, decreasing = TRUE)[1:NClose] # proportion
        } else if (minimize.by == "N"){
          j <- order(Bycatch, decreasing = TRUE)[1:NClose] # numbers
        } else {
          j <- order(Bycatch/D$Effort,decreasing = TRUE)[1:NClose] # CPUE
        }
        # closed area mosaic
        TempClosed[j] <- TRUE
        Closed <- TempClosed
      }
    }
    else{
      Closed <- NA
      for (a in 1:Narea){Closed[a] = FALSE}
    }
    M <- NULL
  }
  if (by.month == TRUE){
    
    TempClosed <- array(dim=nrow(D), FALSE)
    
    if(minimize.by == "N"){
      Table.month <- aggregate(x=Bycatch, by=list(Month=D$Month), FUN=sum) # by-catch in numbers
    } else if (minimize.by == "prop"){
      Table.month <- aggregate(x=Bycatch/Target, by=list(Month=D$Month), FUN=sum)  # prop
    } else {
      Table.month <- aggregate(x=Bycatch/D$Effort, by=list(Month=D$Month), FUN=sum) # CPUE
    }
    M <- Table.month$Month[which(Table.month$x==max(Table.month$x))]
    
    M <- M[1]  # sometimes it doesn't work if there is the same number of by-catch in different months when minimized.by = "N"
    i <- which(D$Month == M)
    TempClosed[i] <- TRUE
    Closed <- TempClosed
    
  }
  return(list(Closed, M)) # matrix of TRUES and FALSES
}
#Closed<-FindClose(NClose=50,weights=TRUE, D,minimize.by="N",mosaic=F,by.month=FALSE);sum(Closed[[1]])
###################################################################################
# FishToTC = T , fish to reach the same total target catch, effort can change
# FishToTC = F , total target catch can change but effort remains the same, 
# effort inside the closed area moves to open areas proportional to the effort already in those open quadrants 
# FishEfficiency = T  CPUE of the target species can change in open areas based in arbitrary exploitation rate
# FishEfficiency = F  CPUE of the target species is the same in each open quadrant, CPUe doesn't change 
Calculate <- function(Closed,D,FishToTC,FishEfficiency, hr, by.month) { #Closed (matrix)comes from the previous function
  
  CPUE <- D$Target/D$Effort # CPUE target species
  
  T_CPUE <- D %>%   # CPUE for each target species
    rowwise() %>% 
    transmute(across(all_of(TNames),function(x)x/Effort)) %>% 
    data.frame
  i <- which(Closed[[1]]==FALSE) # open area or open months
  
  
  TotalE <- sum(D$Effort)
  Narea <- length(Closed[[1]]) # Total number of areas (grids)
  q <- Narea*(-log(1-hr)/TotalE)  # this is the q (catchability) for the target species 
  InitialCatch <- D$Target 
  InitialB <- InitialCatch/(1-exp(-q*D$Effort))  #the initial biomass in each cell for the target species  
  
  BeforeBC <- D %>% summarise(across(.cols=all_of(BCNames),sum)) # before closure
  
  if( FishToTC == FALSE ) {  # total catch can change but not total effort, sum(NewEffort)==sum(D$Effort)
    OldEffort <- sum(D$Effort[i]) # old effort in open areas
    if(by.month == TRUE){
      DisplacedEffort <- GlobalEffort-OldEffort # effort in closed months
    }
    
    if(by.month == FALSE){
      DisplacedEffort <- sum(D$Effort[-i])  # total catch for the target species should be the same, effort can change 
    }
    
    NewEffort <- array(dim=length(D$Effort),0)
    NewEffort[i] <- DisplacedEffort*D$Effort[i]/sum(D$Effort[i])+D$Effort[i] 
    
  }   else { #FishToTC == TRUE
    
    # fishing to reach same TC 
    if(by.month == TRUE){
      TC <- GlobalTC
    }
    
    if(by.month == FALSE){
      TC <- sum(D$Target) # catch for each target species should be the same, but effort can change 
    }
    
    CatchFromOpen <- sum(D$Target[i])  # catch in open areas 
    
    CatchFromOpen_byspecies <- D[i,] %>% summarise(across(.cols=all_of(TNames),sum))
    
    HookIncrease <- TC/CatchFromOpen 
    NewEffort <- array(dim=length(D$Effort),0)
    NewEffort[i] <- HookIncrease*D$Effort[i] # New effort outside 
    OldEffort <- sum(D$Effort[i]) # old effort outside 
    
  }  # end of fishing to TC
  NewCatch <- NewEffort*CPUE
  NewCatch_byspecies <- NewEffort*T_CPUE
  Tcatch <- sum(NewEffort*CPUE) # if CPUE (fishing Efficiency) doesn't change  sum(D$target[i])/sum(D$Effort[i]) = = TC/sum(NewEffort)
  #######################################################################################
  #at this point we have NewEffort (check, NewEffort when FishToTC= TRUE or FALSE is the same)
  if (FishEfficiency == TRUE){ # CPUE can change 
    Catch <- D$Target   # catch by area before closure
    
    Catch_byspecies<- D %>% dplyr::select(TNames)   # catch by area before closure by species
    Abu <- Catch/(1-exp(-q*D$Effort)) #Abundance by area
    Abu_byspecies <- Catch_byspecies/(1-exp(-q*D$Effort))
    NewCatch <- Abu*(1-exp(-q*NewEffort))  #new catch by area
    NewCatch_byspecies <- Abu_byspecies*(1-exp(-q*NewEffort))
    
    Tcatch<-sum(NewCatch)
    if  (FishToTC == TRUE) { #CPUE will change do an incremental calculation
      for (step in 1:20) {  #iterate
        if(by.month == TRUE){
          TC <- GlobalTC
        }
        if(by.month == FALSE){
          TC <- sum(D$Target) # total catch for the target species should be the same, effort can change 
        }
        ratio <- sum(NewCatch)/TC
        NewEffort <- NewEffort/ratio
        NewCatch <- Abu*(1-exp(-NewEffort*q))  #new catch by area
        NewCatch_byspecies <- Abu_byspecies*(1-exp(-NewEffort*q))  #new catch by area
        Tcatch <- sum(NewCatch)
      } 
    }
  }
  
  
  BcCPUE <-D %>% 
    rowwise() %>% 
    transmute(across(BCNames, function(x)x/Effort)) %>% 
    data.frame
  
  BcByArea <- BcCPUE*NewEffort 
  
  TotalBcByArea <- BcByArea %>% transmute(rowSums(.)) %>% pull
  
  ByCatch <- BcByArea %>% summarise(across(BCNames,sum))
  
  
  New_TCatch_byspecies <- NewCatch_byspecies %>% summarise(across(TNames,sum))
  
  TEffort <- sum(NewEffort)
  
  return(list(Tcatch=Tcatch, ByCatch=ByCatch, TEffort=TEffort,
              NewCatch=NewCatch, New_TCatch_byspecies=New_TCatch_byspecies,
              NewEffort=NewEffort, TotalBcByArea=TotalBcByArea))
}
# Calculate(Closed=Closed,D,FishToTC=TRUE,FishEfficiency=FALSE, hr=0.1, by.month=FALSE)
################################################
# Do calculations and plots
# ClosedSeq = seq(0,.5,.1) (sequence of proportions of areas to close) can change
# Months_closed = seq(1,5,1) (sequence of number of months to close) can also be changed
# If maps are desired Maps = TRUE
# tmp.lines = TRUE, add gray lines for each year closed for by-catch species in case of months=F
DoCalcs <- function(D, FishToTC, FishEfficiency, hr, weights=TRUE, ClosedSeq=seq(0,.5,.1), Months_closed=seq(1,5,1),
                  minimize.by="prop", mosaic=FALSE, by.month=FALSE, Maps=TRUE, tmp.lines=TRUE){  
  Corr <- cor(D[,-c(1:4)])
  if (by.month == FALSE){
    Closed_months<-NULL
    D <- aggregate(x=D[,5:ncol(D)],by=list(Lat=D$Lat,Lon=D$Lon,Year=D$Year),FUN=sum)
    
    BC<-list() # to store total by-catch by species, area closure and by year
    TC<-list() # to store total target by species, area closure and by year
    CPUESave<-list() # to store target cpue by closure and year
    TEffort<-list() 
    CatchSave<- list()
    ResArea<-list()
    
    Nyears<-length(unique(D$Year))
    
    for(y in 1:c(Nyears+1)){ # +1: in the last element of the list we will store all years combined
      if(y == Nyears+1){
        D2 <- aggregate(x=D[,4:ncol(D)],by=list(Lat=D$Lat,Lon=D$Lon),FUN=sum)
      } else {
        D2 <- D[D$Year == unique(D$Year)[y],] 
      }
      N <- nrow(D2) #number of areas
      
      BC[[y]] <- matrix(NA, nrow = length(ClosedSeq), ncol=length(col.BC))
      TC[[y]] <- matrix(NA, nrow = length(ClosedSeq), ncol=length(col.T))
      CPUESave[[y]] <- matrix(NA, ncol=length(ClosedSeq), nrow = 1) # for target species 
      TEffort[[y]] <- matrix(NA, ncol=length(ClosedSeq), nrow = 1) # for effort 
      CatchSave[[y]] <- matrix(NA, ncol=length(ClosedSeq), nrow = 1)
      ResArea[[y]] <- data.frame(Lat=NA, Lon=NA, NewEffort=NA, CPUE_Target=NA,
                               CPUE_ByCatch= NA, Closure=NA)
      # 
      for (i in 1:length(ClosedSeq) )   { # looping over the proportion of areas closed
        pClose <- ClosedSeq[i] # prop closed
        NClose <- round(pClose*N,0) # numbers of areas closed
        
        Closed <- FindClose(NClose, D=D2, weights = weights, minimize.by=minimize.by, mosaic = mosaic,by.month=FALSE)
        
        Res <- Calculate(Closed, D=D2, FishToTC, FishEfficiency, hr, by.month=FALSE)
        Res_a <- data.frame(Lat=D2$Lat, Lon=D2$Lon, NewEffort=Res$NewEffort,
                          CPUE_Target=Res$NewCatch/(Res$NewEffort*1E3),
                          CPUE_ByCatch= Res$TotalBcByArea/(Res$NewEffort),
                          Closure=ClosedSeq[i])
        BC[[y]][i,] <- as.matrix(Res$ByCatch) 
        TC[[y]][i,] <- as.matrix(Res$New_TCatch_byspecies)
        CPUESave[[y]][i] <- Res$Tcatch/Res$TEffort
        TEffort[[y]][i] <- Res$TEffort
        CatchSave[[y]][i] <- Res$Tcatch
        ResArea[[y]] <- rbind(ResArea[[y]],Res_a)
      }
      ResArea[[y]] <- ResArea[[y]][-1,]
      
    }
    
    if (weights == TRUE){
      BC_w <- list()
      TC_w <- list()
      
      for(i in 1:c(Nyears+1)){
        BC_w[[i]] <- sweep(BC[[i]], MARGIN=2, as.numeric(rel_wb), `*`) # weighted by-Catch
        TC_w[[i]] <- sweep(TC[[i]], MARGIN=2, as.numeric(rel_wt), `*`) # weighted Target Catch
      }
      
    } else {
      BC_w <- BC
      TC_w <- TC
    }
    
    TBC <- apply(BC_w[[Nyears+1]],1,FUN=sum)
    BCScaled <- TBC/TBC[1] # relative to the first value (no closure)
    CPUEScaled <- CPUESave[[Nyears+1]]/CPUESave[[Nyears+1]][1]
    TEffortScaled <- TEffort[[Nyears+1]]/TEffort[[Nyears+1]][1]
    CatchScaled <- CatchSave[[Nyears+1]]/CatchSave[[Nyears+1]][1]
    
    # what if we sum up when areas were closed by year: 
    
    TBC_y<-0
    TEffort_y<-0
    CatchSave_y<-0
    TBC_tmp<-list()
    for(y in 1:Nyears){
      TBC_tmp[[y]]<-apply(BC_w[[y]],1,FUN=sum,na.rm=T)
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
    nPol<-length(BC_w[[Nyears+1]][,1])
    nSpec<-length(BC_w[[Nyears+1]][1,])+length(TC_w[[Nyears+1]][1,])
    BC_D<-matrix(0,nrow = nPol,ncol=nSpec) 
    for(i in 1:Nyears){
      Total<-cbind(TC_w[[i]],BC_w[[i]])
      BC_D<-BC_D  + Total 
    }
    
    BC_S<-cbind(TC_w[[Nyears+1]],BC_w[[Nyears+1]])
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
      
      Closed<-FindClose(NClose,D=D2,minimize.by=minimize.by, weights=weights,mosaic = mosaic,by.month=TRUE) #NClose does not matter for months
      Closed_months<-c(Closed_months,Closed[[2]])
      Res<-Calculate(Closed,D=D2,FishToTC,FishEfficiency,hr,by.month=TRUE)
      
      BC[[i]]<-as.data.frame(Res$ByCatch) 
      TC[[i]]<-as.data.frame(Res$New_TCatch_byspecies)
      CPUESave[[i]]<-Res$Tcatch/Res$TEffort
      TEffort[[i]]<-Res$TEffort
      CatchSave[[i]]<-Res$Tcatch
      
      #remove the month with the highest bycatch to calculate the next one
      D2<-D2[D2$Month!=Closed[[2]],]
    }
    
    if (weights == TRUE){
      BC_w <- list()
      TC_w <- list()
      
      for(i in 1:c(Nmonths+1)){
        BC_w[[i]] <- sweep(BC[[i]], MARGIN=2, as.numeric(rel_wb), `*`) # weighted by-Catch
        TC_w[[i]] <- sweep(TC[[i]], MARGIN=2, as.numeric(rel_wt), `*`) # weighted Target Catch
      }
      
    } else {
      BC_w <- BC
      TC_w <- TC
    }
    
    TBC<-NULL
    CPUEScaled<-NULL
    TEffortScaled<-NULL
    CatchScaled<-NULL
    for(i in 1:length(BC_w))  {
      tmp<-c(apply(BC_w[[i]],1,FUN=sum))
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
    nSpec<-ncol(BC_w[[1]])+ncol(TC_w[[1]])
    BC_S<-matrix(0,nrow = nPol,ncol=nSpec) 
    for(i in 1:c(Nmonths+1)){
      BC_S[i,]<- t(rbind(t(TC_w[[i]]),t(BC_w[[i]])))
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
      g<-worldmap+
        geom_tile(data=Tmp[Tmp$variable=="NewEffort",], mapping=aes(Lon, Lat, fill = value)) +
        facet_grid(Closure~variable) +
        scale_fill_gradient(low = "white", high = "red")+
        theme(legend.title = element_blank(),legend.position = "top",legend.text = element_text(size = 8),strip.text.y = element_blank())
      gg<-worldmap+
        geom_tile(data=Tmp[Tmp$variable=="CPUE_Target",], mapping=aes(Lon, Lat, fill = value)) +
        facet_grid(Closure~variable) +
        scale_fill_gradient(low = "white", high = "red")+
        theme(legend.title = element_blank(),legend.position = "top",axis.ticks = element_blank(),axis.title.y = element_blank(),axis.text.y=element_blank(),legend.text = element_text(size = 8),strip.text.y = element_blank())
      ggg<-worldmap+
        geom_tile(data=Tmp[Tmp$variable=="CPUE_ByCatch",], mapping=aes(Lon, Lat, fill = value)) +
        facet_grid(Closure~variable) +
        scale_fill_gradient(low = "white", high = "red")+
        theme(legend.title = element_blank(),legend.position = "top",axis.ticks = element_blank(),axis.title.y = element_blank(),axis.text.y=element_blank(),legend.text = element_text(size = 8))
      if (by.month==FALSE){
        grid.arrange(g,gg,ggg,ncol=3,top=paste("Year",unique(D$Year)[y],",","Mosaic=",mosaic,"Minimized by=",minimize.by,",", "Fish to TC =",FishToTC,",","Fishing Efficiency changes =",FishEfficiency))
      } else {
        grid.arrange(g,gg,ggg,ncol=3,top=paste("Month",unique(D$Month)[y],",","Mosaic=",mosaic,"Minimized by=",minimize.by,",", "Fish to TC =",FishToTC,",","Fishing Efficiency changes =",FishEfficiency))
      }
    }
  }
  saveRDS(object=list(Closed_Seq=ClosedSeq,Closed_months=Closed_months, BC=BC, BC_w=BC_w, TC=TC, TC_w=TC_w,
                      BCScaled=BCScaled,TEffortScaled=TEffortScaled,CatchScaled=CatchScaled,CPUEScaled=CPUEScaled,
                      BCScaled_y=BCScaled_y,TEffortScaled_y=TEffortScaled_y,CatchScaled_y=CatchScaled_y,CPUEScaled_y=CPUEScaled_y,
                      Changes_bySpecies=BC_tot, Corr=Corr),
          file=file.path(results_dir,paste0(Fishery_name,"_Weights=",weights,"_Mosaic=",mosaic,"_FishToTC=",FishToTC,"_FishEfficiency=",FishEfficiency,"_minimize.by=",minimize.by,
                      "_by.month=",by.month,".rds")))
} # end of function
##################################################################################################################
# the following set of lines run all possible scenarios, and with and without weights for sensitivity analysis
##################################################################################################################
# with weights identified by the "W"
pdf(file=file.path(results_dir,paste0("Outs_prop_W_",Fishery_name,".pdf")))
DoCalcs(D,FishToTC=TRUE,FishEfficiency=FALSE,weights=TRUE,hr,minimize.by="prop",Maps=FALSE,by.month=TRUE,Months_closed=seq(1,5,1))
DoCalcs(D,FishToTC=TRUE,FishEfficiency=FALSE,weights=TRUE,hr,minimize.by="prop",Maps=FALSE,by.month=FALSE)
DoCalcs(D,FishToTC=TRUE,FishEfficiency=TRUE,weights=TRUE,hr,minimize.by="prop",Maps=FALSE,by.month=TRUE)
DoCalcs(D,FishToTC=TRUE,FishEfficiency=TRUE,weights=TRUE,hr,minimize.by="prop",Maps=FALSE,by.month=FALSE)
DoCalcs(D,FishToTC=FALSE,FishEfficiency=FALSE,weights=TRUE,hr,minimize.by="prop",Maps=FALSE,by.month=TRUE)
DoCalcs(D,FishToTC=FALSE,FishEfficiency=TRUE,weights=TRUE,hr,minimize.by="prop",Maps=FALSE,by.month=FALSE)
DoCalcs(D,FishToTC=FALSE,FishEfficiency=TRUE,weights=TRUE,hr,minimize.by="prop",Maps=FALSE,by.month=TRUE)
DoCalcs(D,FishToTC=FALSE,FishEfficiency=FALSE,weights=TRUE,hr,minimize.by="prop",Maps=TRUE,by.month=FALSE)
DoCalcs(D,FishToTC=TRUE,FishEfficiency=FALSE,weights=TRUE,hr,minimize.by="prop",mosaic=TRUE,Maps=FALSE)
DoCalcs(D,FishToTC=TRUE,FishEfficiency=TRUE,weights=TRUE,hr,minimize.by="prop",mosaic=TRUE,Maps=FALSE)
DoCalcs(D,FishToTC=FALSE,FishEfficiency=FALSE,weights=TRUE,hr,minimize.by="prop",mosaic=TRUE,Maps=FALSE)
DoCalcs(D,FishToTC=FALSE,FishEfficiency=TRUE,weights=TRUE,hr,minimize.by="prop",mosaic=TRUE,Maps=TRUE)
dev.off()
###############################################################################
pdf(file=file.path(results_dir,paste0("Outs_N_W_",Fishery_name,".pdf")))
DoCalcs(D,FishToTC=TRUE,FishEfficiency=FALSE,weights=TRUE,hr,minimize.by="N",Maps=FALSE,by.month=TRUE)
DoCalcs(D,FishToTC=TRUE,FishEfficiency=FALSE,weights=TRUE,hr,minimize.by="N",Maps=FALSE,by.month=FALSE)
DoCalcs(D,FishToTC=TRUE,FishEfficiency=TRUE,weights=TRUE,hr,minimize.by="N",Maps=FALSE,by.month=TRUE)
DoCalcs(D,FishToTC=TRUE,FishEfficiency=TRUE,weights=TRUE,hr,minimize.by="N",Maps=FALSE,by.month=FALSE)
DoCalcs(D,FishToTC=FALSE,FishEfficiency=FALSE,weights=TRUE,hr,minimize.by="N",Maps=FALSE,by.month=TRUE)
DoCalcs(D,FishToTC=FALSE,FishEfficiency=TRUE,weights=TRUE,hr,minimize.by="N",Maps=FALSE,by.month=FALSE)
DoCalcs(D,FishToTC=FALSE,FishEfficiency=TRUE,weights=TRUE,hr,minimize.by="N",Maps=FALSE,by.month=TRUE)
DoCalcs(D,FishToTC=FALSE,FishEfficiency=FALSE,weights=TRUE,hr,minimize.by="N",Maps=TRUE,by.month=FALSE)
DoCalcs(D,FishToTC=TRUE,FishEfficiency=FALSE,weights=TRUE,hr,minimize.by="N",mosaic=TRUE,Maps=FALSE)
DoCalcs(D,FishToTC=TRUE,FishEfficiency=TRUE,weights=TRUE,hr,minimize.by="N",mosaic=TRUE,Maps=FALSE)
DoCalcs(D,FishToTC=FALSE,FishEfficiency=FALSE,weights=TRUE,hr,minimize.by="N",mosaic=TRUE,Maps=FALSE)
DoCalcs(D,FishToTC=FALSE,FishEfficiency=TRUE,weights=TRUE,hr,minimize.by="N",mosaic=TRUE,Maps=TRUE)
dev.off()
###############################################################################
pdf(file=file.path(results_dir,paste0("Outs_BCrate_W_",Fishery_name,".pdf")))
DoCalcs(D,FishToTC=TRUE,FishEfficiency=FALSE,weights=TRUE,hr,minimize.by="rate",Maps=FALSE,by.month=TRUE)
DoCalcs(D,FishToTC=TRUE,FishEfficiency=FALSE,weights=TRUE,hr,minimize.by="rate",Maps=FALSE,by.month=FALSE)
DoCalcs(D,FishToTC=TRUE,FishEfficiency=TRUE,weights=TRUE,hr,minimize.by="rate",Maps=FALSE,by.month=TRUE)
DoCalcs(D,FishToTC=TRUE,FishEfficiency=TRUE,weights=TRUE,hr,minimize.by="rate",Maps=FALSE,by.month=FALSE)
DoCalcs(D,FishToTC=FALSE,FishEfficiency=FALSE,weights=TRUE,hr,minimize.by="rate",Maps=FALSE,by.month=TRUE)
DoCalcs(D,FishToTC=FALSE,FishEfficiency=TRUE,weights=TRUE,hr,minimize.by="rate",Maps=FALSE,by.month=FALSE)
DoCalcs(D,FishToTC=FALSE,FishEfficiency=TRUE,weights=TRUE,hr,minimize.by="rate",Maps=FALSE,by.month=TRUE)
DoCalcs(D,FishToTC=FALSE,FishEfficiency=FALSE,weights=TRUE,hr,minimize.by="rate",Maps=TRUE,by.month=FALSE)
DoCalcs(D,FishToTC=TRUE,FishEfficiency=FALSE,weights=TRUE,hr,minimize.by="rate",mosaic=TRUE,Maps=FALSE)
DoCalcs(D,FishToTC=TRUE,FishEfficiency=TRUE,weights=TRUE,hr,minimize.by="rate",mosaic=TRUE,Maps=FALSE)
DoCalcs(D,FishToTC=FALSE,FishEfficiency=FALSE,weights=TRUE,hr,minimize.by="rate",mosaic=TRUE,Maps=FALSE)
DoCalcs(D,FishToTC=FALSE,FishEfficiency=TRUE,weights=TRUE,hr,minimize.by="rate",mosaic=TRUE,Maps=TRUE)
dev.off()

# Un-weighted
pdf(file=file.path(results_dir,paste0("Outs_prop_",Fishery_name,".pdf")))
DoCalcs(D,FishToTC=TRUE,FishEfficiency=FALSE,weights=FALSE,hr,minimize.by="prop",Maps=FALSE,by.month=TRUE,Months_closed=seq(1,5,1))
DoCalcs(D,FishToTC=TRUE,FishEfficiency=FALSE,weights=FALSE,hr,minimize.by="prop",Maps=FALSE,by.month=FALSE)
DoCalcs(D,FishToTC=TRUE,FishEfficiency=TRUE,weights=FALSE,hr,minimize.by="prop",Maps=FALSE,by.month=TRUE)
DoCalcs(D,FishToTC=TRUE,FishEfficiency=TRUE,weights=FALSE,hr,minimize.by="prop",Maps=FALSE,by.month=FALSE)
DoCalcs(D,FishToTC=FALSE,FishEfficiency=FALSE,weights=FALSE,hr,minimize.by="prop",Maps=FALSE,by.month=TRUE)
DoCalcs(D,FishToTC=FALSE,FishEfficiency=TRUE,weights=FALSE,hr,minimize.by="prop",Maps=FALSE,by.month=FALSE)
DoCalcs(D,FishToTC=FALSE,FishEfficiency=TRUE,weights=FALSE,hr,minimize.by="prop",Maps=FALSE,by.month=TRUE)
DoCalcs(D,FishToTC=FALSE,FishEfficiency=FALSE,weights=FALSE,hr,minimize.by="prop",Maps=TRUE,by.month=FALSE)
DoCalcs(D,FishToTC=TRUE,FishEfficiency=FALSE,weights=FALSE,hr,minimize.by="prop",mosaic=TRUE,Maps=FALSE)
DoCalcs(D,FishToTC=TRUE,FishEfficiency=TRUE,weights=FALSE,hr,minimize.by="prop",mosaic=TRUE,Maps=FALSE)
DoCalcs(D,FishToTC=FALSE,FishEfficiency=FALSE,weights=FALSE,hr,minimize.by="prop",mosaic=TRUE,Maps=FALSE)
DoCalcs(D,FishToTC=FALSE,FishEfficiency=TRUE,weights=FALSE,hr,minimize.by="prop",mosaic=TRUE,Maps=TRUE)
dev.off()
###############################################################################
pdf(file=file.path(results_dir,paste0("Outs_N_",Fishery_name,".pdf")))
DoCalcs(D,FishToTC=TRUE,FishEfficiency=FALSE,weights=FALSE,hr,minimize.by="N",Maps=FALSE,by.month=TRUE)
DoCalcs(D,FishToTC=TRUE,FishEfficiency=FALSE,weights=FALSE,hr,minimize.by="N",Maps=FALSE,by.month=FALSE)
DoCalcs(D,FishToTC=TRUE,FishEfficiency=TRUE,weights=FALSE,hr,minimize.by="N",Maps=FALSE,by.month=TRUE)
DoCalcs(D,FishToTC=TRUE,FishEfficiency=TRUE,weights=FALSE,hr,minimize.by="N",Maps=FALSE,by.month=FALSE)
DoCalcs(D,FishToTC=FALSE,FishEfficiency=FALSE,weights=FALSE,hr,minimize.by="N",Maps=FALSE,by.month=TRUE)
DoCalcs(D,FishToTC=FALSE,FishEfficiency=TRUE,weights=FALSE,hr,minimize.by="N",Maps=FALSE,by.month=FALSE)
DoCalcs(D,FishToTC=FALSE,FishEfficiency=TRUE,weights=FALSE,hr,minimize.by="N",Maps=FALSE,by.month=TRUE)
DoCalcs(D,FishToTC=FALSE,FishEfficiency=FALSE,weights=FALSE,hr,minimize.by="N",Maps=TRUE,by.month=FALSE)
DoCalcs(D,FishToTC=TRUE,FishEfficiency=FALSE,weights=FALSE,hr,minimize.by="N",mosaic=TRUE,Maps=FALSE)
DoCalcs(D,FishToTC=TRUE,FishEfficiency=TRUE,weights=FALSE,hr,minimize.by="N",mosaic=TRUE,Maps=FALSE)
DoCalcs(D,FishToTC=FALSE,FishEfficiency=FALSE,weights=FALSE,hr,minimize.by="N",mosaic=TRUE,Maps=FALSE)
DoCalcs(D,FishToTC=FALSE,FishEfficiency=TRUE,weights=FALSE,hr,minimize.by="N",mosaic=TRUE,Maps=TRUE)
dev.off()
###############################################################################
pdf(file=file.path(results_dir,paste0("Outs_BCrate_",Fishery_name,".pdf")))
DoCalcs(D,FishToTC=TRUE,FishEfficiency=FALSE,weights=FALSE,hr,minimize.by="rate",Maps=FALSE,by.month=TRUE)
DoCalcs(D,FishToTC=TRUE,FishEfficiency=FALSE,weights=FALSE,hr,minimize.by="rate",Maps=FALSE,by.month=FALSE)
DoCalcs(D,FishToTC=TRUE,FishEfficiency=TRUE,weights=FALSE,hr,minimize.by="rate",Maps=FALSE,by.month=TRUE)
DoCalcs(D,FishToTC=TRUE,FishEfficiency=TRUE,weights=FALSE,hr,minimize.by="rate",Maps=FALSE,by.month=FALSE)
DoCalcs(D,FishToTC=FALSE,FishEfficiency=FALSE,weights=FALSE,hr,minimize.by="rate",Maps=FALSE,by.month=TRUE)
DoCalcs(D,FishToTC=FALSE,FishEfficiency=TRUE,weights=FALSE,hr,minimize.by="rate",Maps=FALSE,by.month=FALSE)
DoCalcs(D,FishToTC=FALSE,FishEfficiency=TRUE,weights=FALSE,hr,minimize.by="rate",Maps=FALSE,by.month=TRUE)
DoCalcs(D,FishToTC=FALSE,FishEfficiency=FALSE,weights=FALSE,hr,minimize.by="rate",Maps=TRUE,by.month=FALSE)
DoCalcs(D,FishToTC=TRUE,FishEfficiency=FALSE,weights=FALSE,hr,minimize.by="rate",mosaic=TRUE,Maps=FALSE)
DoCalcs(D,FishToTC=TRUE,FishEfficiency=TRUE,weights=FALSE,hr,minimize.by="rate",mosaic=TRUE,Maps=FALSE)
DoCalcs(D,FishToTC=FALSE,FishEfficiency=FALSE,weights=FALSE,hr,minimize.by="rate",mosaic=TRUE,Maps=FALSE)
DoCalcs(D,FishToTC=FALSE,FishEfficiency=TRUE,weights=FALSE,hr,minimize.by="rate",mosaic=TRUE,Maps=TRUE)
dev.off()