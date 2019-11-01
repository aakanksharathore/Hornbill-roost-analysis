##Code to analyse hornbill movement patterns from nest to roost
##Written by AA on 14 Sept 2019

library(geosphere)
library(magrittr)
library(ggplot2)
library(ggpubr)
##Import movement data from csv file (raw data file)
fname <- file.choose()   #choose rawDat.csv in RoostAnalysis folder
dat = read.csv(fname, header=TRUE)

#Impot roost information from csv file
fname <- file.choose()   #choose Roost.csv in RoostAnalysis folder
datR = read.csv(fname,header=TRUE)

#File - mean daily distance moved meanDist.csv
fname <- file.choose()
meanDD = read.csv(fname, header=TRUE)

##Remove garbage data
dat<-dat[dat$lon != 0 & dat$lat != 0,]


##Distribution of movement distances and angles
datN=dat[dat$code=="gabbar",]    ##Set hornbill name here

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
##Distance from roost at each time-point
dats=unique(datN$date.time.1)
dtime=vector()
ddate=vector()
datRi = datR[datR$code=="gabbar",]   

distR=numeric()

#towards and away movement
toaway=vector()
cnt=1
for(i in 2:nrow(datN)){
  
  roost=datRi[which(as.character(datRi$date.time.1) == datN$date.time.1[i]),c("lon","lat")]
  if(nrow(roost)==0){
    next
  }
  if(datN$date.time.1[i] == datN$date.time.1[i-1]){
    distR[cnt] = distHaversine(roost,c(datN$lon[i],datN$lat[i]))
    dtime[cnt]= as.character(datN$rounded[i])
    ddate[cnt]= as.character(datN$date.time.1[i])
    
    if(cnt>1){
    if(distR[cnt]>distR[cnt-1]){
      toaway[cnt]= "Away"
    }else{
      toaway[cnt]="Towards"
    }
    }
    cnt=cnt+1
  }
  
}

#Plot distance from roost
dist2R = na.omit(as.data.frame(cbind(as.numeric(distR),dtime,ddate)))
dtimeO=c("3:30","3:45","4:0","4:15","4:30","4:45","5:0","5:15","5:30","5:45","6:0","6:15","6:30","6:45","7:0","7:15","7:30","7:45","8:0","8:15","8:30","8:45","9:0","9:15","9:30","9:45","10:0","10:15","10:30","10:45","11:0","11:15","11:30","11:45","12:0","12:15","12:30","12:45","13:0","13:15","13:30","13:45","14:0","14:15","14:30","14:45","15:0","15:15","15:30","15:45","16:0","16:15","16:30","16:45","17:0","17:15","17:30","17:45","18:0","18:15","18:30","18:45","19:0")
ggerrorplot(dist2R, x = "dtime", y = "distR", 
            desc_stat = "mean_sd", color = "black",
            add = "jitter", add.params = list(color = "darkgray"),order=dtimeO,
            xlab="Time of the day",ylab="Distance in meters",title="Gabbar"
)+ rotate_x_text(90)+
geom_point(data=meanDD,mapping=aes(y=gabbar,x=Row.Labels),shape=23, fill="blue",size=3)+scale_x_discrete(limits=dtimeO)

# Plot distance from roost - daily separately
ggplot(dist2R, aes(x=dtime, y=distR, 
                 group=ddate,
                 col=factor(ddate))) +
  geom_line()+scale_x_discrete(limits=dtimeO)+ rotate_x_text(90)+
  ggtitle("Gabbar")

#Plot movement towards and away from nest
dist2away = na.omit(as.data.frame(cbind(toaway,dtime,ddate)))
dtimeO=c("3:30","3:45","4:0","4:15","4:30","4:45","5:0","5:15","5:30","5:45","6:0","6:15","6:30","6:45","7:0","7:15","7:30","7:45","8:0","8:15","8:30","8:45","9:0","9:15","9:30","9:45","10:0","10:15","10:30","10:45","11:0","11:15","11:30","11:45","12:0","12:15","12:30","12:45","13:0","13:15","13:30","13:45","14:0","14:15","14:30","14:45","15:0","15:15","15:30","15:45","16:0","16:15","16:30","16:45","17:0","17:15","17:30","17:45","18:0","18:15","18:30","18:45","19:0")

counts <- table(dist2away$toaway, dist2away$dtime)
# barplot(counts, main="Gabbar",
#         xlab="Time of the day", ylab="Distance in meters",col=c("blue","gray"),
#         legend = c("Towards","Away"), beside=TRUE)

counts1=as.data.frame(counts)
ggbarplot(counts1, x = "Var2", y = "Freq", fill = "Var1",order=dtimeO,palette=c("blue","gray"))+ rotate_x_text(90)+ggtitle("Gabbar")

##Separately for each day

###################Path Directedness############################################
###Use sliding window approach
##Analysis for each day separately and save the time of day


drtness = numeric()
dtime=vector()
ddate=vector()
winS =    1    #Sliding window steps
winL =    4    #Sliding window length - one hour
cnt=1

dats <- unique(datN$date.time.1)
for(i in 1:length(dats)){
  temp=datN[datN$date.time.1==dats[i],]
  if(nrow(temp)==0){
    next
  }
  j=1
  while(1){
    if(j>=(nrow(temp)-winL)){
      break
    }
  
    #convert displacement to the unit in which distance was calculated
  drtness[cnt] = ( sqrt(temp$lon[min(nrow(temp),(j+winL))]-temp$lon[j])^2 + (temp$lat[min(nrow(temp),(j+winL))]-temp$lat[j])^2)/(sum(temp$dist.moved[j:min(nrow(temp),(j+winL))]))
  dtime[cnt]= as.character(temp$rounded[j])
  ddate[cnt]= as.character(temp$date.time.1[j])
  cnt=cnt+1
  j=j+winS
  }
}
hist(drtness)
drct=na.omit(as.data.frame(cbind(dtime,ddate,as.numeric(as.character(drtness)))))
##Directness as a function of the time of day

ggerrorplot(drct, x = "dtime", y = "drtness", 
            desc_stat = "mean_sd", color = "black",xlab="Time of the day",ylab="Path directedness",
            add = "violin", add.params = list(color = "darkgray"),order=dtimeO
)+ rotate_x_text(90)+ggtitle("Gabbar")

##Not a very clear pattern for pooled data but we look at each day separately, e might see some trends
##Directedness for each day

ggplot(drct, aes(x=dtime, y=drtness, 
                   group=ddate,xlab="Time of the day",ylab="Path directedness",
                   col=factor(ddate))) +
  geom_line()+scale_x_discrete(limits=dtimeO)+ rotate_x_text(90)+
  ggtitle("Gabbar")

