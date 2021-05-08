####### Some libraries needed
library(tidyverse)
library(reshape2)
library(grid)
library(gridExtra)
library(ggcorrplot)
library(shiny)
library(maps)
library(shinycssloaders)

###############################################
# this routine searches over spatial locations of closed areas 
# to find the pattern of closed areas that will have the max reduction in: 
# 1) "N": bycatch in numbers; 2) "rate": bycatch CPUE; 3) "prop": ByCatch/target   
# it will not be a square, but the areas closest to a centroid
FindClose<-function(NClose,D,Tweights,BCweights,minimize.by="prop",mosaic=F,by.month=FALSE){
  # NT<-length(Tweights) # number of target species 
  # NBC<-length(BCweights) # number of byCatch species
  # 
  # col.BC<-seq(from=ncol(D2)- NBC-3,to=ncol(D2)-4) # columns for by-catch species 
  # col.T<-seq(from=ncol(D2)- NBC-4-NT,to=1+
  # BCNames<-names(D[col.BC])
  # TNames<-names(D[col.T])
  #BCNames<-colnames(BCweights)
  Bycatch <- D$TBC.w
  Target <- D$Target.w
  
  if (by.month==FALSE){
    
    #col.BC<-seq(from=ncol(D)- NBC -1 ,to=ncol(D)-4)
    
    Narea<-nrow(D) # Number of areas: this is just the number of rows in the data frame
    
    #Bycatch<-D2 %>% dplyr::select(all_of(BCNames)) # weighted numbers
    #now loop over over centroids
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
  return(list(Closed,M)) # matrix of TRUES and FALSES
}
###################################################################################
Calculate<-function(Closed,D,Tweights,BCweights,GlobalTC=GlobalTC,GlobalEffort=GlobalEffort,FishToTC,FishEfficiency,hr=0.1,by.month) { #Closed (matrix)comes from the previous function
  NT<-length(Tweights) # number of target species 
  NBC<-length(BCweights) # number of byCatch species
  
  BCNames<-names(BCweights)
  TNames<-names(Tweights)
  
  CPUE<-D$Target/D$Effort # CPUE target species
  #T_CPUE<-D[,Tcols]/D$Effort
  
  T_CPUE<-D %>%
    rowwise() %>%
    transmute(across(all_of(TNames),function(x)x/Effort)) %>%
    data.frame
  i<-which(Closed[[1]]==FALSE) # open area or open months
  
  # TotAbu=sum(D$target*D$Effort)/hr
  TotalE<-sum(D$Effort)
  Narea<-length(Closed[[1]]) # Total number of areas (grids)
  q<-Narea*(-log(1-hr)/TotalE)  #this is the q for the target species 
  InitialCatch<-D$Target
  InitialB<-InitialCatch/(1-exp(-q*D$Effort))  #the initial biomass in each cell for the target species  
  
  
  #set Bycatch to what it would be in the open areas with existing effort
  
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
    #Displaced effort allocated proportional to the effort outside the closed area + the effort already existing in that cell
    
  }   else { #FishToTC==TRUE
    
    # fishing to reach same TC 
    if(by.month==TRUE){
      TC<-GlobalTC
    }
    
    if(by.month==FALSE){
      TC<-sum(D$Target) # catch for each target species should be the same, but effort can change 
    }
    
    CatchFromOpen<- sum(D$Target[i])  # catch in open areas 
    
    CatchFromOpen_byspecies<- D[i,] %>% summarise(across(.cols=all_of(TNames),sum))
    #CatchFromOpen_byspecies<- colSums(D[i,Tcols])
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
    Catch<- D$Target  # catch by area before closure
    #Catch_byspecies<- D[,Tcols] 
    Catch_byspecies<- D %>% dplyr::select(all_of(TNames))   # catch by area before closure by species
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
          TC<-sum(D$Target) # total catch for the target species should be the same, effort can change 
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
    transmute(across(all_of(BCNames),function(x)x/Effort)) %>% 
    data.frame
  
  BcByArea<-BcCPUE*NewEffort 
  
  TotalBcByArea <- BcByArea %>% transmute(rowSums(.)) %>% pull
  
  ByCatch <- BcByArea %>% summarise(across(all_of(BCNames),sum))
  
  New_TCatch_byspecies <- NewCatch_byspecies %>% summarise(across(all_of(TNames),sum))
  
  TEffort<-sum(NewEffort)
  
  return(list(Tcatch=Tcatch,ByCatch=ByCatch,TEffort=TEffort,
              NewCatch=NewCatch,New_TCatch_byspecies=New_TCatch_byspecies,
              NewEffort=NewEffort,TotalBcByArea=TotalBcByArea))
}
################################################
DoCalcs<-function(D,Tweights,BCweights,GlobalTC=GlobalTC,GlobalEffort=GlobalEffort,FishToTC,FishEfficiency,hr=0.1,ClosedSeq=seq(0,.5,.1),Months_closed=seq(1,5,1),
                  minimize.by="prop",mosaic=FALSE,by.month=FALSE,Maps=TRUE,rel_wb, rel_wt){  # do the simulation
  NT<-length(Tweights) # number of target species 
  NBC<-length(BCweights) # number of byCatch species
  
  # col.BC<-seq(from=ncol(D)- NBC-3,to=ncol(D)-4) # columns for by-catch species 
  # col.T<-seq(from=5+1,to=5+NT)
  
  BCNames<-names(BCweights)
  TNames<-names(Tweights)
  
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
      
      BC[[y]]<-matrix(NA,nrow = length(ClosedSeq),ncol=length(BCweights))
      TC[[y]]<-matrix(NA,nrow = length(ClosedSeq),ncol=length(Tweights))
      CPUESave[[y]]<- matrix(NA,ncol=length(ClosedSeq),nrow = 1)# for target species 
      TEffort[[y]]<- matrix(NA,ncol=length(ClosedSeq),nrow = 1) #for effort 
      CatchSave[[y]]<- matrix(NA,ncol=length(ClosedSeq),nrow = 1)
      ResArea[[y]]<-data.frame(Lat=NA,Lon=NA,NewEffort=NA,ByCatch=NA,CPUE_Target=NA,
                               CPUE_ByCatch= NA,Proportion = NA, Closure=NA)
      # 
      for (i in 1:length(ClosedSeq) )   { #looping over the proportion of areas closed
        pClose<-ClosedSeq[i] # prop closed
        NClose<-round(pClose*N,0) # numbers of areas closed
        
        Closed<-FindClose(NClose,D=D2,minimize.by=minimize.by,mosaic = mosaic,by.month=FALSE,Tweights=Tweights,BCweights=BCweights)
        
        Res<-Calculate(Closed,D=D2,FishToTC=FishToTC,FishEfficiency=FishEfficiency,hr=hr,by.month=FALSE,
                       Tweights=Tweights,BCweights=BCweights)
        Res_a<-data.frame(Lat = D2$Lat,Lon = D2$Lon, NewEffort = Res$NewEffort,
                          ByCatch = (Res$TotalBcByArea/(Res$NewEffort))*(Res$NewEffort),
                          CPUE_Target=Res$NewCatch/(Res$NewEffort),
                          CPUE_ByCatch = Res$TotalBcByArea/(Res$NewEffort),
                          Proportion = Res$TotalBcByArea/Res$NewCatch,
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
    
    # if I want to plot the weighted bycatch I have to multiply again for the tmp.weights becasue in "Calculate" we use raw numbers
    for(i in 1:c(Nyears+1)){
      BC[[i]] <- sweep(BC[[i]], MARGIN=2, as.numeric(rel_wb), `*`) # weighted by-Catch
      TC[[i]] <- sweep(TC[[i]], MARGIN=2, as.numeric(rel_wt), `*`) # weighted Target Catch
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
    # data needed comes from BC, list where each element is a year, and the last one is the total (Static closure)
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
    BC_S_rel$Area<-as.factor("Static");BC_D_rel$Area<-as.factor("Mobile")
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
    Nyears=NULL
    TBC_tmp=NULL
    ResArea=NULL
    D2<-aggregate(x=D[,5:ncol(D)],by=list(Lat=D$Lat,Lon=D$Lon,Month=D$Month),FUN=sum)
    
    BC<-list() # to store total bycatch by species
    TC<-list() # to store total target by species, area closure and by year
    CPUESave<-list() # to store target cpue 
    TEffort<-list() # to store effort
    CatchSave<- list() # to store catch 
    #ResArea<-list()
    
    Closed<-list(array(dim=nrow(D2),FALSE),NULL)
    
    Res<-Calculate(Closed,D=D2,GlobalTC=GlobalTC,GlobalEffort=GlobalEffort,
                   FishToTC=FishToTC,FishEfficiency=FishEfficiency,hr=hr,by.month = TRUE,
                   Tweights=Tweights,BCweights=BCweights) #when no closure
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
      
      Closed<-FindClose(NClose,D=D2,minimize.by=minimize.by,mosaic = mosaic,by.month=TRUE,Tweights=Tweights,BCweights=BCweights) #NClose does not matter for months
      Closed_months<-c(Closed_months,Closed[[2]])
      Res<-Calculate(Closed,D=D2,GlobalTC=GlobalTC,GlobalEffort=GlobalEffort,
                     FishToTC=FishToTC,FishEfficiency=FishEfficiency,
                     hr=hr,by.month=TRUE,Tweights=Tweights,BCweights=BCweights)
      
      BC[[i]]<-as.data.frame(Res$ByCatch) # weighted bycatch 
      TC[[i]]<-as.data.frame(Res$New_TCatch_byspecies)
      CPUESave[[i]]<-Res$Tcatch/Res$TEffort
      TEffort[[i]]<-Res$TEffort
      CatchSave[[i]]<-Res$Tcatch
      
      #remove the month with the higest bycatch to calculate the next one
      D2<-D2[D2$Month!=Closed[[2]],]
    }
    
    for(i in 1:c(Nmonths+1)){
      BC[[i]] <- sweep(BC[[i]], MARGIN=2, as.numeric(rel_wb), `*`) # weighted by-Catch
      TC[[i]] <- sweep(TC[[i]], MARGIN=2, as.numeric(rel_wt), `*`) # weighted Target Catch
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
    # data needed comes from BC, list where each element is a year, and the last one is the total (Static closure)
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
    BC_S_rel$Month<-as.factor("Static")
    BC_tot<-melt(BC_S_rel,value.name = "BCN")
    BC_tot$Relative_ByCatch<-BC_tot$BCN-1 
    names(BC_tot)[3]<-"Species"
    tmp<-data.frame(Category="Target",Species=TNames)
    tmp<-rbind(tmp,data.frame(Category="ByCatch",Species=BCNames))
    BC_tot<-merge(BC_tot,tmp,by="Species",all = T)
    
    ### NULL, only outputs when year closures
    BCScaled_y=NULL;TEffortScaled_y=NULL;CatchScaled_y=NULL;CPUEScaled_y=NULL
  }
  
  return(list(Nyears=Nyears,ClosedSeq = ClosedSeq, Closed_months=Closed_months,Months_closed=Months_closed,
              BC=BC,CPUE_t=CPUESave,TEffort=TEffort,Catch=CatchSave, TBC_tmp=TBC_tmp,ResArea = ResArea,
              BCScaled=BCScaled,TEffortScaled=TEffortScaled,CatchScaled=CatchScaled,CPUEScaled=CPUEScaled,BC_tot=BC_tot,
              BCScaled_y=BCScaled_y,TEffortScaled_y=TEffortScaled_y,CatchScaled_y=CatchScaled_y,CPUEScaled_y=CPUEScaled_y))
}
############################################
# Define IU for the shiny application
ui <- fluidPage(
  h1("Spatial-temporal closures as a tool to reduce bycatch interactions (under development)"),
  h4("Instructions"),
  h5("You'll need 2 csv files to run this application. 
        One should be called 'Data.csv' and the other one 'Weights.csv'"),
  h5("After uploading the csv files, wait until the plots are generated. This could take a few 
     seconds or minutes based on how fast is your computer and how much data you are uploading. If you decide to send the model outputs to
     mpons@uw.edu, please do not change the default name when downloading the RDS file. 
     You can also send to the same email address any questions and suggestions to improve the app. Thanks!"),
  
  
  sidebarLayout(
    sidebarPanel(
      h6("This is only used to name the output RDS file. Example: US swordfish"),
      textInput("Fishery_name", "Name of the Fishery", ""),
      # Input: Select a file ----
      h6("The first 5 columns are the same for all case studies:
      'Lat',	'Lon',	'Year', 'Month',	and 'Effort'. The following columns may vary for each case study. 
       First input the target and then the bycatch by species. You can include any number of target and bycatch species or group of species. 
      IMPORTANT: each row in the data file corresponds to a quadrant (lat&long) by year and month combination.
      There is only one record for each of these combinations. 
      The values in each column for the target and bycatch species can be in numbers or biomass."),
      
      fileInput("file1", "Choose csv Data File",
                multiple = TRUE,
                accept = c("text/csv",
                           "text/comma-separated-values,text/plain",
                           ".csv")),
      
      radioButtons("sep1", "Separator",
                   choices = c(Comma = ",",
                               Semicolon = ";"),
                   selected = ","),
      h6("This csv needs the first row to be the names of the target species in the same order as in the data file. The second row
         should be the weight for each target species. These should sum up to 1.
         The third row should be the names of the Bycatch species in the same order as those in the data file. 
         The fourth row should be the weights for those species. These also should sum up to 1."),
      # Input: Select a file ----
      fileInput("file2", "Choose Weights File",
                multiple = TRUE,
                accept = c("text/csv",
                           "text/comma-separated-values,text/plain",
                           ".csv")),
      radioButtons("sep2", "Separator",
                   choices = c(Comma = ",",
                               Semicolon = ";"),
                   selected = ","),
      
      br(),
      h6("The areas will be closed by minimizing the following:"),
      radioButtons(inputId = "minimize.by", label="Minimize by:",
                   choices=c("ByCatch/Target", "ByCatch rate" ,"Absolute bycatch")),
      
      br(),
      h6("Fishing to Total Catch means that the effort outside the closed areas will increase in order to get the same Total Catch as it was before the closure.
         Total Catch before and after the closure will be the same.
         Fishing to Total Effort means that the total effort inside the closed area will be reallocated outside and the total catch might change. 
         The total effort before and after the closure is the same.
         In all cases we assume that when areas are closed, effort moves to open areas proportional to the amount of effort that was in the open areas before closure."),
      radioButtons(inputId = "FishToTC", label= "Fish to:", 
                   choices=c("Total Catch" , "Total Effort")),
      
      br(),
      h6("If Fishing efficiency does not change, it means that the CPUE for the target species in each quadrant 
         outside the area remains the same no matters the increase in effort. If fishing efficiency can change, 
         we assume that as fishing pressure increases in open areas, due to displaced effort, 
         the CPUE of the target species will change due to a decline in biomass"),
      radioButtons(inputId = "FishEfficiency" , label = "Fishing efficiency", 
                   choices=c("Does not change", "Changes")) 
    ),
    
    mainPanel(
      
      tabsetPanel(type = "tabs",

                  tabPanel("Main Results", actionButton("plot2", "Run", icon("running")) %>% withSpinner(color="#0dc5c1"),
                           plotOutput("plot2"),
                            downloadButton(outputId='downloadLPlot',label='download.LPlot'),
                           downloadButton("SaveRDS", "Save Output as a RDS file")),
                  tabPanel("Barplot_Centroid", plotOutput("plot3"),
                            downloadButton(outputId='downloadBPlotC',label='download.BPlot.C')),
                  tabPanel("Barplot_Mosaic", plotOutput("plot8"),
                           downloadButton(outputId='downloadBPlotM',label='download.BPlot.M')),
                  tabPanel("Barplot_Temporal", plotOutput("plot9"),
                           downloadButton(outputId='downloadBPlotT',label='download.BPlot.T')),
                  tabPanel("Static centroid closures", plotOutput("plot4"),
                           downloadButton(outputId='downloadCMap',label='download.CMap')),
                  tabPanel("Static mosaic closures", plotOutput("plot7"),
                           downloadButton(outputId='downloadMMap',label='download.MMap')),
                  tabPanel("Correlation", plotOutput("plot5"),
                           downloadButton(outputId='downloadCorPlot',label='download.CorPlot')),
                  tabPanel("Weights", plotOutput("plot6"),
                           downloadButton(outputId='downloadWPlot',label='download.WPlot')),
                  tabPanel("Map raw_data", plotOutput("plot1"),
                           downloadButton(outputId='downloadMap',label='download.Map'))
                  
      ),
      
    )
  )
)

server <- function(input, output, session) {
  dataframe<-reactive({
    req(input$file1,input$file2)
    
    if (is.null(input$file1))
      return(NULL)                
    D<-read.csv(input$file1$datapath,
                header = TRUE,
                sep =input$sep1)
    D<-D[complete.cases(D),]
    if (is.null(input$file2))
      return(NULL)
    
    Tweights<-read.csv(input$file2$datapath,sep=input$sep2,nrows = 1)%>%
      select_if(~ !any(is.na(.))) 
    
    BCweights<-read.csv(input$file2$datapath,sep=input$sep2,nrows = 2,skip = 2)%>%
      select_if(~ !any(is.na(.))) 
    
    NT<-ncol(Tweights) # number of target species 
    NBC<-ncol(BCweights) # number of byCatch species
    
    ### Effort in thousands 
    D$Effort<-D$Effort/1000
    
    BCNames<-colnames(BCweights)
    TNames<-colnames(Tweights)
    
    col.BC <- which(colnames(D)%in%all_of(BCNames))
    col.T <- which(colnames(D)%in%all_of(TNames))
    
    #Grid_res<-input$Grid_res
    ####################################################
    i=which(D$Lon<0)
    D$Lon[i]<- D$Lon[i]+ 360 #  rescale to 360 degrees to have all positive values for longitude
    D$Lat<- D$Lat+ 90  #  rescale to 180 degrees to have all positive values for latitude
    
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
    D_w<-D # keep D unweighted
    
    for (i in 1:nrow(D_w)){
      D_w[i,col.BC]<- D_w[i,col.BC]*rel_wb #assign bycatch weights
      D_w[i,col.T]<- D_w[i,col.T]*rel_wt   #assign target weights
    } 
    
    D_w <- D_w %>% 
      mutate(TBC.w=rowSums(.[all_of(BCNames)]),# total weighted By-Catch
             Target.w=rowSums(.[all_of(TNames)])) %>% # total weighted By-Catch
      filter(Target.w>0)# eliminate sets with target catch =0 to be able to calculate the proportion of Bycatch/target
    
    D <- D %>% 
      mutate(TBC=rowSums(.[all_of(BCNames)]),# total weighted By-Catch
             Target=rowSums(.[all_of(TNames)])) %>% # total weighted By-Catch
      filter(Target>0)
    # add to D weighted columns for total bycatch and total target
    D$TBC.w<-D_w$TBC.w
    D$Target.w<-D_w$Target.w
    
    # create a map with all years combined
    D_sum<-aggregate(x=D[,5:ncol(D)],by=list(Lat=D$Lat,Lon=D$Lon),FUN=sum)
    D_sum$Prop<-D_sum$TBC.w/D_sum$Target.w # % of ByCatch to target species
    
    Tmp<-melt(D_sum,id.vars = c("Lat","Lon"))
    # for mapping
    world <- ggplot2::map_data('world')
   
    i=which(world$long<0)
    world$long[i]<- world$long[i]+ 360 #  rescale to 360 degrees to have all positive values for longitude
    world$lat<- world$lat+ 90  #  rescale to 180 degrees to have all positive values for latitude
    
    world.cut<-world[world$lat > min(D$Lat)-5 & world$lat < max(D$Lat)+5 ,]
    world.cut<-world.cut[world.cut$long > min(D$Lon)-5 & world.cut$long < max(D$Lon)+5 ,]
    
    #This creates a quick and dirty world map - playing around with the themes, aesthetics, and device dimensions is recommended!
    worldmap <- ggplot(world.cut, aes(x=long, y=lat)) +
      geom_polygon(aes(group=group)) +
      theme(panel.background = element_rect(fill="skyblue", colour="skyblue"),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank()) +
      coord_equal()
    ##################### 
    GlobalTC<-sum(D[,col.T])
    GlobalEffort<-sum(D$Effort)
    
    return(list(Tmp=Tmp,D=D,D_w=D_w,col.BC=col.BC,col.T=col.T,NBC=NBC,NT=NT,BCNames=BCNames,TNames=TNames,
                Tweights=Tweights,BCweights=BCweights,worldmap=worldmap,rel_wb=rel_wb,rel_wt=rel_wt,
                GlobalTC=GlobalTC,GlobalEffort=GlobalEffort))
    
  })
  myplot1<-function(){
    
    dataframe()$Tmp %>% group_by(variable) %>%
      do(gg = {
        dataframe()$worldmap +
          geom_tile(., mapping=aes(Lon, Lat, fill = value)) +
          facet_grid(~variable) +
          scale_fill_gradient(low = "white", high = "red")+

        theme(legend.title = element_blank(),legend.position = "top",
              legend.text = element_text(size = 10),strip.text.y = element_blank(),
              plot.margin = unit(c(0, 0, 0, 0), "lines"),
              text = element_text(size=15))
      }) %>%
      .$gg %>% arrangeGrob(grobs = ., nrow = 3) %>% grid.arrange()

  }
  output$plot1<-renderPlot({
    myplot1()
  })
  output$downloadMap <- downloadHandler(  
    filename = function() {
      # specify the filename
      paste('Map_data', 'png',sep=".")
    },
    content = function(file) {
      
      png(file,width = 1500, height = 1000)
      # create the plot
        myplot1()
      # close the device
      dev.off()
    }
  )
  
  ### Run functions
  Out<-reactive({
    
    if(input$minimize.by=="ByCatch/Target"){
      minimize.by="prop"
    }
    if(input$minimize.by=="Absolute bycatch"){
      minimize.by="N"
    }
    if(input$minimize.by=="ByCatch rate"){
      minimize.by="rate"
    }

    if(input$FishToTC=="Total Catch"){
      FishToTC=TRUE
    }
    if(input$FishToTC=="Total Effort"){
      FishToTC=FALSE
    }
    if(input$FishEfficiency=="Changes"){
      FishEfficiency=TRUE
    }
    if(input$FishEfficiency=="Does not change"){
      FishEfficiency=FALSE
    }
    OutID <- paste0("_minimize.by=",minimize.by, "_FishToTC=",FishToTC,
                 "_FishEfficiency=",FishEfficiency)
    
    Res1<-DoCalcs(D=dataframe()$D,GlobalTC=dataframe()$GlobalTC,GlobalEffort=dataframe()$GlobalEffort,
                 BCweights=dataframe()$BCweights,Tweights=dataframe()$Tweights,
                 ClosedSeq=seq(0,.5,.1),Months_closed=seq(1,5,1),
                 mosaic = FALSE, rel_wb=dataframe()$rel_wb, rel_wt=dataframe()$rel_wt, 
                 FishToTC=FishToTC,FishEfficiency=FishEfficiency,hr=0.1,
                 minimize.by=minimize.by,Maps=FALSE,by.month=FALSE)
    Res2<-DoCalcs(D=dataframe()$D,GlobalTC=dataframe()$GlobalTC,GlobalEffort=dataframe()$GlobalEffort,
                  BCweights=dataframe()$BCweights,Tweights=dataframe()$Tweights,
                  ClosedSeq=seq(0,.5,.1),Months_closed=seq(1,5,1),
                  mosaic = TRUE,rel_wb=dataframe()$rel_wb, rel_wt=dataframe()$rel_wt,
                  FishToTC=FishToTC,FishEfficiency=FishEfficiency,hr=0.1,
                  minimize.by=minimize.by,Maps=FALSE,by.month=FALSE)
    Res3<-DoCalcs(D=dataframe()$D,GlobalTC=dataframe()$GlobalTC,GlobalEffort=dataframe()$GlobalEffort,
                  BCweights=dataframe()$BCweights,Tweights=dataframe()$Tweights,
                  ClosedSeq=seq(0,.5,.1),Months_closed=seq(1,5,1),
                  mosaic = TRUE,rel_wb=dataframe()$rel_wb, rel_wt=dataframe()$rel_wt,
                  FishToTC=FishToTC,FishEfficiency=FishEfficiency,hr=0.1,
                  minimize.by=minimize.by,Maps=FALSE,by.month=TRUE)

    return(list(Centroid=Res1, Mosaic=Res2, Temporal=Res3, OutID=OutID))
  })

   Out_plot<-eventReactive(input$plot2, {
    df <- NULL
    df <- rbind(df, data.frame(Type.closure = rep("Static Centroid",6), #FishToTC= rep(FishToTC, 6),
                               #FishEfficiency= rep(FishEfficiency, 6),
                               Closure = Out()[[1]]$ClosedSeq, Effort = t(Out()[[1]]$TEffortScaled), ByCatch= Out()[[1]]$BCScaled,
                               Catch= t(Out()[[1]]$CatchScaled)))
    df<- rbind(df, data.frame(Type.closure = rep("Dynamic Centroid",6), #FishToTC= rep(FishToTC, 6),
                              #FishEfficiency= rep(FishEfficiency, 6),
                              Closure = Out()[[1]]$ClosedSeq, Effort = t(Out()[[1]]$TEffortScaled_y), ByCatch= Out()[[1]]$BCScaled_y,
                              Catch= t(Out()[[1]]$CatchScaled_y)))
    df <- rbind(df, data.frame(Type.closure = rep("Static Mosaic",6), #FishToTC= rep(FishToTC, 6),
                               #FishEfficiency= rep(FishEfficiency, 6),
                               Closure = Out()[[2]]$ClosedSeq, Effort = t(Out()[[2]]$TEffortScaled), ByCatch= Out()[[2]]$BCScaled,
                               Catch= t(Out()[[2]]$CatchScaled)))
    df<- rbind(df, data.frame(Type.closure = rep("Dynamic Mosaic",6), #FishToTC= rep(FishToTC, 6),
                              #FishEfficiency= rep(FishEfficiency, 6),
                              Closure = Out()[[2]]$ClosedSeq, Effort = t(Out()[[2]]$TEffortScaled_y), ByCatch= Out()[[2]]$BCScaled_y,
                              Catch= t(Out()[[2]]$CatchScaled_y)))
    df<- rbind(df, data.frame(Type.closure = rep("Temporal",6), #FishToTC= rep(FishToTC, 6),
                              #FishEfficiency= rep(FishEfficiency, 6),
                              Closure = Out()[[3]]$ClosedSeq, Effort = Out()[[3]]$TEffortScaled, ByCatch= Out()[[3]]$BCScaled,
                              Catch= Out()[[3]]$CatchScaled))
    df2 <- gather(df, "ByCatch", "Effort", "Catch", key= "quant", value= "value")
    
    df2$Type.closure <- as.factor(df2$Type.closure)
    df2$Type.closure <- factor(df2$Type.closure, levels = c("Static Centroid", "Dynamic Centroid",
                                                            "Static Mosaic", "Dynamic Mosaic" , "Temporal" ))
    return(df2)
   })
  # 

  myplot2<-function(){
    
        if(!is.null(Out_plot()))
        par(mfrow=c(1,1),mar=(c(5, 4, 4, 2) + 0.1),oma=c(1,1,1,1))
      cbPalette <- c("#CC79A7","#009E73","#56B4E9","#E69F00","#999999")
    
    p <- ggplot(Out_plot(), aes(x = Closure, y = value,  color=Type.closure))
    p <- p + geom_hline(yintercept=1, linetype="dashed",
                        color = "grey40", size=0.4)
    p <- p + geom_line(alpha=0.8, size=1) +
      scale_color_manual(values = cbPalette ) +
      scale_fill_manual(values = cbPalette ) +
      theme_bw() +
      theme(legend.position = "top", legend.title = element_blank(), panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(), legend.text=element_text(size=9), strip.text.y = element_text(angle=0),
            axis.text.x = element_text(angle = 45, hjust=1, vjust=1)) +
      ylim(0,1.5)
    p <- p +  facet_grid ( ~  quant,
                          labeller = label_wrap_gen(width = 17,multi_line = TRUE))
    p + labs(y="Relative change") + guides(colour=guide_legend(nrow=2,byrow=TRUE))+
      scale_x_continuous(sec.axis = sec_axis(~ . *10, name = "Number of months closed"),
                         name="Proportion of area closed") + theme(axis.text.x.top = element_text(angle = 0, hjust = 0.5))
    if(input$minimize.by=="ByCatch/Target"){
      minimize.by="ByCatch/Target"
    }
    if(input$minimize.by=="Absolute bycatch"){
      minimize.by="Absolute ByCatch"
    }
    if(input$minimize.by=="ByCatch rate"){
      minimize.by="Bycatch rate"
    }
    if(input$FishToTC=="Total Catch"){
      FishToTC="Total Catch"
    }
    if(input$FishToTC=="Total Effort"){
      FishToTC="Total Effort"
    }
    if(input$FishEfficiency=="Changes"){
      FishEfficiency="Changes"
    }
    if(input$FishEfficiency=="Does not change"){
      FishEfficiency="Does not change"
    }
    Title = paste("Minimized by:",minimize.by, ";", "Fishing to reach same:", FishToTC, ";", 
                  "Fishing efficiency:", FishEfficiency, sep=" ")
    
    p<-p + ggtitle(Title)
    print(p)
  }
  output$plot2 <- renderPlot({
    if(!is.null(Out_plot()))
      myplot2()
  })


  output$downloadLPlot <- downloadHandler(
    filename = function() {
      # specify the filename
      paste('L.plot', 'png',sep=".")
    },
    content = function(file) {
      # open the device
      ggsave(file, device = "png", width=10, height=6)
      # create the plot
      myplot2()
      # close the device
      dev.off()
    }
  )
  
   myplot3<-function(){
    
      g<-ggplot(Out()[[1]]$BC_tot, aes(x=Species, y=Relative_ByCatch, fill=Category, alpha=Area)) +
        geom_bar(stat="identity",position=position_dodge())+theme_minimal()+
        facet_wrap(~Closure,ncol=1)+ theme(legend.position="top",axis.text.x = element_text(angle = 90)) +
        scale_fill_manual(values = c("mediumorchid3", "turquoise3"))+
        scale_alpha_manual(values = c(0.5, 1))+ ylab("Relative Catch")

    print(g)
  }
  output$plot3<-renderPlot({
    myplot3()
  })
  output$downloadBPlotC <- downloadHandler(
    filename = function() {
      paste('Barplot_Centroid', '.png',sep="")
    },
    content = function(file) {
      #ggsave(file, device = "png", width=8, height=7)

      png(file,width = 500, height = 500, pointsize = 22)
      myplot3()
      dev.off()
    }
  )
  myplot8<-function(){

    g<-ggplot(Out()$Mosaic$BC_tot, aes(x=Species, y=Relative_ByCatch, fill=Category, alpha=Area)) +
      geom_bar(stat="identity",position=position_dodge())+theme_minimal()+
      facet_wrap(~Closure,ncol=1)+ theme(legend.position="top",axis.text.x = element_text(angle = 90)) +
      scale_fill_manual(values = c("mediumorchid3", "turquoise3"))+
      scale_alpha_manual(values = c(0.5, 1))+ ylab("Relative Catch")

      print(g)
  }
  output$plot8<-renderPlot({
    myplot8()
  })
  output$downloadBPlotM <- downloadHandler(
    filename = function() {
      paste('Barplot_Mosaic', '.png',sep="")
    },
    content = function(file) {
      #ggsave(file, device = "png", width=8, height=7)
      png(file,width = 500, height = 500, pointsize = 22)
      myplot8()
      dev.off()
    }
  )

  myplot9<-function(){

      g<-ggplot(Out()$Temporal$BC_tot, aes(x=Species, y=Relative_ByCatch, fill=Category)) +
        geom_bar(stat="identity",position=position_dodge())+theme_minimal()+
        facet_wrap(~Closure,ncol=1)+ theme(legend.position="top",axis.text.x = element_text(angle = 90)) +
        scale_fill_manual(values = c("mediumorchid3", "turquoise3")) + ylab("Relative Catch")

    print(g)
  }
  output$plot9<-renderPlot({
    myplot9()
  })
  output$downloadBPlotT <- downloadHandler(
    filename = function() {
      paste('Barplot_Temporal', '.png',sep="")
    },
    content = function(file) {
      #ggsave(file, device = "png", width=8, height=7)
      png(file,width = 500, height = 500, pointsize = 22)
      myplot9()
      dev.off()
    }
  )
  ###################### Mapping
  Tmp<-reactive({
    N<-length(Out()[[1]]$ResArea)
    ResArea<-Out()[[1]]$ResArea[[N]]
    ResArea[ResArea$NewEffort==0,"NewEffort"]<-NA
    tmp<-melt(ResArea,id.vars = c("Lat","Lon","Closure"))
    tmp$Closure<-as.factor(tmp$Closure)
    #tmp$ByCatch
    return(tmp)
  })
  myplot4<-function(){
    g<-dataframe()$worldmap +
      geom_tile(data=Tmp()[Tmp()$variable=="NewEffort",], mapping=aes(Lon, Lat, fill = value)) +
      facet_grid(Closure~variable) +
      scale_fill_gradient(low = "white", high = "red")+
      theme(legend.title = element_blank(),legend.position = "top",
            legend.text = element_text(size = 10),strip.text.y = element_blank(),
            plot.margin = unit(c(0, 0, 0, 0), "lines"),
            text = element_text(size=15))
    gg<-dataframe()$worldmap+
      geom_tile(data=Tmp()[Tmp()$variable=="ByCatch",], mapping=aes(Lon, Lat, fill = value)) +
      facet_grid(Closure~variable) +
      scale_fill_gradient(low = "white", high = "red")+
      theme(legend.title = element_blank(),legend.position = "top",
            axis.ticks = element_blank(),axis.title.y = element_blank(),
            axis.text.y=element_blank(),legend.text = element_text(size = 10),
            strip.text.y = element_blank(), 
            plot.margin = unit(c(0, 0, 0, 0), "lines"),
            text = element_text(size=15))
    ggg<-dataframe()$worldmap+
      geom_tile(data=Tmp()[Tmp()$variable=="CPUE_ByCatch",], mapping=aes(Lon, Lat, fill = value)) +
      facet_grid(Closure~variable) +
      scale_fill_gradient(low = "white", high = "red")+
      theme(legend.title = element_blank(),legend.position = "top",
            axis.ticks = element_blank(),axis.title.y = element_blank(),
            axis.text.y=element_blank(),legend.text = element_text(size = 10),
            plot.margin = unit(c(0, 0, 0, 0), "lines"),
            text = element_text(size=15))
    gggg<-dataframe()$worldmap+
      geom_tile(data=Tmp()[Tmp()$variable=="Proportion",], mapping=aes(Lon, Lat, fill = value)) +
      facet_grid(Closure~variable) +
      scale_fill_gradient(low = "white", high = "red")+
      theme(legend.title = element_blank(),legend.position = "top",
            axis.ticks = element_blank(),axis.title.y = element_blank(),
            axis.text.y=element_blank(),legend.text = element_text(size = 10),
            plot.margin = unit(c(0, 0, 0, 0), "lines"),
            text = element_text(size=15))
    G<-grid.arrange(g,gg,ggg,gggg,ncol=4)
  }
  output$plot4<-renderPlot({
    #
    myplot4()
  })
  output$downloadCMap <- downloadHandler(
    filename = function() {
      paste('OutMapsC', '.png',sep="")
    },
    content = function(file) {
      #ggsave(file, device = "png", width=11, height=8.5)
      png(file, width = 1000, height = 1000, pointsize = 12)
      myplot4()
      dev.off()
    }
  )
  #######################################################
  Tmp2<-reactive({
    N<-length(Out()[[2]]$ResArea)
    ResArea<-Out()[[2]]$ResArea[[N]]
    ResArea[ResArea$NewEffort==0,"NewEffort"]<-NA
    tmp<-melt(ResArea,id.vars = c("Lat","Lon","Closure"))
    tmp$Closure<-as.factor(tmp$Closure)
    #tmp$ByCatch
    return(tmp)
  })
  myplot7<-function(){
    g<-dataframe()$worldmap +
      geom_tile(data=Tmp2()[Tmp2()$variable=="NewEffort",], mapping=aes(Lon, Lat, fill = value)) +
      facet_grid(Closure~variable) +
      scale_fill_gradient(low = "white", high = "red")+
      theme(legend.title = element_blank(),legend.position = "top",
            legend.text = element_text(size = 10),strip.text.y = element_blank(),
            plot.margin = unit(c(0, 0, 0, 0), "lines"),
            text = element_text(size=15))
    gg<-dataframe()$worldmap+
      geom_tile(data=Tmp2()[Tmp2()$variable=="ByCatch",], mapping=aes(Lon, Lat, fill = value)) +
      facet_grid(Closure~variable) +
      scale_fill_gradient(low = "white", high = "red")+
      theme(legend.title = element_blank(),legend.position = "top",
            axis.ticks = element_blank(),axis.title.y = element_blank(),
            axis.text.y=element_blank(),legend.text = element_text(size = 10),
            strip.text.y = element_blank(), 
            plot.margin = unit(c(0, 0, 0, 0), "lines"),
            text = element_text(size=15))
    ggg<-dataframe()$worldmap+
      geom_tile(data=Tmp2()[Tmp2()$variable=="CPUE_ByCatch",], mapping=aes(Lon, Lat, fill = value)) +
      facet_grid(Closure~variable) +
      scale_fill_gradient(low = "white", high = "red")+
      theme(legend.title = element_blank(),legend.position = "top",
            axis.ticks = element_blank(),axis.title.y = element_blank(),
            axis.text.y=element_blank(),legend.text = element_text(size = 10),
            plot.margin = unit(c(0, 0, 0, 0), "lines"),
            text = element_text(size=15))
    gggg<-dataframe()$worldmap+
      geom_tile(data=Tmp2()[Tmp2()$variable=="Proportion",], mapping=aes(Lon, Lat, fill = value)) +
      facet_grid(Closure~variable) +
      scale_fill_gradient(low = "white", high = "red")+
      theme(legend.title = element_blank(),legend.position = "top",
            axis.ticks = element_blank(),axis.title.y = element_blank(),
            axis.text.y=element_blank(),legend.text = element_text(size = 10),
            plot.margin = unit(c(0, 0, 0, 0), "lines"),
            text = element_text(size=15))
    G<-grid.arrange(g,gg,ggg,gggg,ncol=4)
  }
  output$plot7<-renderPlot({
    #
    myplot7()
  })
  output$downloadMMap <- downloadHandler(
    filename = function() {
      paste('OutMapsM', '.png',sep="")
    },
    content = function(file) {
      #ggsave(file, device = "png", width=11, height=8.5)
      png(file,width = 1000, height = 1000, pointsize = 12)
      myplot7()
      dev.off()
    }
  )
  #######################################################
  myplot5<-function(){
    ggcorrplot(round(cor(dataframe()$D[,-c(1:4)]), 1), hc.order = FALSE, type = "lower",lab=TRUE,
               outline.col = "white")

  }
  output$plot5<-renderPlot({
    myplot5()
  })
  output$downloadCorPlot <- downloadHandler(
    filename = function() {
      paste('Corrplot', '.png',sep="")
    },
    content = function(file) {
      ggsave(file, device = "png", width=11, height=8.5)

      #png(file,width = 1000, height = 1000, pointsize = 22)
      myplot5()
      #dev.off()
    }
  )
  #######################################################
  myplot6<-function(){
    par(mfrow=c(1,2),oma=c(5,1,1,1),mar=c(1,4,2,0))
    plot(as.numeric(dataframe()$D %>% summarise(across(.cols=all_of(dataframe()$TNames),sum)))/
           sum(dataframe()$D[,dataframe()$TNames]),xaxt='n',
         col="blue",pch=19,xlab="",ylab="")
    points(as.numeric(dataframe()$D_w %>% summarise(across(.cols=all_of(dataframe()$TNames),sum)))/
             sum(dataframe()$D_w[,dataframe()$TNames]),col="red",pch=19)
    axis(side = 1, at = c(1:dataframe()$NT),labels = c(dataframe()$TNames),las=2)
    title(main="Target",line = 0.5)
    legend("topright",col=c("blue","red"),pch=19,legend = c("Unweighted","Weighted"))

    plot(as.numeric(dataframe()$D %>% summarise(across(.cols=all_of(dataframe()$BCNames),sum)))/
           sum(dataframe()$D[,dataframe()$BCNames]),xaxt='n',
         col="blue",pch=19,xlab="",ylab="")
    points(as.numeric(dataframe()$D_w %>% summarise(across(.cols=all_of(dataframe()$BCNames),sum)))/
             sum(dataframe()$D_w[,dataframe()$BCNames]),col="red",pch=19)
    axis(side = 1, at = c(1:dataframe()$NBC),labels = c(dataframe()$BCNames),las=2)
    title(main="Bycatch",line = 0.5)

  }
  output$plot6<-renderPlot({
    myplot6()
  })
  output$downloadWPlot <- downloadHandler(
    filename = function() {
      paste('Weigths_plot', '.png',sep="")
    },
    content = function(file) {
      png(file,width = 1500, height = 1000, pointsize = 22)
      myplot6()
      dev.off()
    }
  )
  # Save state
  output$SaveRDS <- downloadHandler(
    filename = function() {

      paste0(input$Fishery_name,Out()$OutID,".rds")
    },
    content = function(file) {
      data_out <- isolate(Out())

      saveRDS(data_out, file)
    }
  )

}

shinyApp(ui=ui, server=server)