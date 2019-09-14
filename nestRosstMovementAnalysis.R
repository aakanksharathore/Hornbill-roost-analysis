##Code to analyse hornbill movement patterns from nest to roost
##Written by AA on 14 Sept 2019

##Import movement data from csv file (raw data file)
fname <- file.choose()   #choose rawDat.csv in RoostAnalysis folder
dat = read.csv(fname, header=TRUE)
attach(dat)

#Impot roost information from csv file
fname <- file.choose()   #choose Roost.csv in RoostAnalysis folder
datR = read.csv(fname,header=TRUE)

##Remove garbage data
dat<-dat[dat$lon != 0 & dat$lat != 0,]


##Distribution of movement distances and angles
datN=dat[dat$code=="bill",]    ##Set hornbill name here

#------------------------------------------------------------------------------#

  ##Plotting distributions for the pooled data

hist(datN$dist.moved, xlab="Step length",ylab="Frequency",col="grey",main="") #confirm with Rohit the unit of this distance
       
      ##Calculate angles for each step
cnt=1
angD=numeric()
for(i in 2:nrow(datN)){
  
  if(datN$date.time.1[i] == datN$date.time.1[i-1]){
    angD[cnt] = atan((datN$lon[i]-datN$lon[i-1])/(datN$lat[i]-datN$lat[i-1]))*180/pi
    cnt=cnt+1
  }
  
}

hist(angD, xlab="Step length",ylab="Frequency",col="grey",main="") #confirm with Rohit the unit of this distance

#---------------------------------------------------------------------------------------------#

###################Path Directedness############################################
###Use sliding window approach
##Analysis for each day separately and save the time of day

dats=unique(datN$date.time.1)
drtness = numeric()
dtime=vector()
winS =    1    #Sliding window steps
winL =    10    #Sliding window length
cnt=1
for(i in 1:length(dats)){
  temp=datN[datN$date.time.1==dats[i],]
  if(nrow(temp)==0){
    next
  }
  j=1
  while(1){
    if(j>=nrow(temp)){
      break
    }
  
    #convert displacement to the unit in which distance was calculated
  drtness[cnt] = ( sqrt(temp$lon[min(nrow(temp),(j+winL))]-temp$lon[j])^2 + (temp$lat[min(nrow(temp),(j+winL))]-temp$lat[j])^2)/(sum(temp$dist.moved[j:min(nrow(temp),(j+winL))]))
  dtime[cnt]= as.character(temp$rounded[j])
  cnt=cnt+1
  j=j+winS
  }
}
hist(drtness)
drct=na.omit(as.data.frame(cbind(dtime,as.numeric(as.character(drtness)))))
##Directness as a function of the time of day
plot(drct$dtime,drct$drtness)

library(magrittr)
library(ggpubr)

ggerrorplot(drct, x = "dtime", y = "drtness", 
            desc_stat = "mean_sd", color = "black",
            add = "jitter", add.params = list(color = "darkgray")
)

ggerrorplot(drct, x = "dtime", y = "drtness", 
            desc_stat = "mean_sd", color = "black",
            add = "violin", add.params = list(color = "darkgray")
)


#-------------------------------------------------------------------------------------------------------#

###Heading of hornbill wrt to nest-roost vector at all time-points
