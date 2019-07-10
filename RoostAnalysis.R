##Code for analysing see dispersal at roosts by hornbill
##Written on 10 July 2019 by AA

library(lubridate)
library(sp)
library(geosphere)
library(plotrix)

##Import movement data from csv file (raw data file)
fname <- file.choose()   #choose rawdata.csv in RoostAnalysis folder
dat = read.csv(fname, header=TRUE)
attach(dat)

#Impot roost information from csv file
fname <- file.choose()   #choose Roost.csv in RoostAnalysis folder
datR = read.csv(fname,header=TRUE)

##Remove garbage data
dat<-dat[dat$lon != 0 & dat$lat != 0,]

##Import gut passage data
gpt<-read.csv("regur2.csv",header=TRUE)

#Activity budget file - we will use it for sampling during bootstrapping
act=read.csv("proportions_foraging.csv",header=TRUE)


##Output data table
final = data.frame(matrix(ncol = 4, nrow = 0))
colnames(final)=c("type","lat","long","dist")
##Specify the individual name and its nest location
#Change these values for which individual results are required
Hname="gabbar"
nest_lat=26.936533
nest_long=92.967067  

##Only consider data from days when more than 30 locations are available (to avoid underestimationof dispersal)
#and the dates for each Hornbill which are available in Roost.csv
idate = datR[datR$code==Hname,"date.time.1"]
datF = dat[dat$date.time.1 %in% idate,]


##Vectors to save nest and roost values
datF$nest = numeric(length=nrow(datF))
datF$roost =numeric(length=nrow(datF))
k = 1

#Find out nest locations and remove them
for(i in 1:nrow(datF))
{
  lat = datF[i,"lat"]
  long = datF[i,"lon"]
  dist_nest = distHaversine(c(long,lat), c(nest_long,nest_lat))
  if(dist_nest < 20)
    datF$nest[i] = 1
  #dist_roost = distHaversine(c(long,lat), c(roost_long,roost_lat))
  #if(dist_roost < 50)
  #  data$roost[i] = 1
}
data1 <- datF[datF$nest != 1,]

if(nrow(data1) == 0)
{
  message <- "All points lie within nest area"
  message
}

##Create a vector to represent each row of data frame as many times as activity level in that hour
times=numeric(length=nrow(data1))
for(i in 1:nrow(data1)){
  times[i]=round(act$Breeding[max(which((strptime(act$Time.of.the.day,format="%H:%M:%S") == strptime(data1$rounded[i],format="%H:%M"))))])

}
times[is.na(times)] = 0
samp_vec = rep(1:nrow(data1),times=times)

##calculation for every time bin
#final$nest = 0
#final$roost = 0
##Loop over all gpts
for(i in 1:length(gpt$time))
{
  #Select first bin from GPT
  tim = gpt$bin[i]
  j = 1
  while(j <= 10)
  {
    #Sampling according to activity levels
    ran <- sample(samp_vec,1)
    ##Select source location and time
    latx = data1[ran,"lat"]
    lonx = data1[ran,"lon"]  
    day = data1[ran,"date.time.1"] 
    reftime = strptime(data1$date.time[ran],format = "%d/%m/%y %H:%M")
    reftime = format(reftime,"%d/%m/%y %H:%M")
    if(reftime>=format(strptime("17:00","%H:%M"),"%H:%M"))
      next
    t = 0
    temp = which(datF$no == data1$no[ran])
    if(temp <= nrow(datF))
      while(datF[temp,"date.time.1"]  == day)
      {
        temp = temp+1
        if(temp >= nrow(datF))
          break
        
        dif = as.duration(strptime(datF$date.time[temp],format = "%d/%m/%y %H:%M") - strptime(reftime,format="%d/%m/%y %H:%M")) 
        dif = as.numeric(dif/60)
        if(!is.na(dif) & length(dif) != 0 & dif <= tim & dif > (tim-15))
        {
          t = temp
          break
        }
        
        
      }
    if(t == 0)                        ##take last value if destination was not found within given interval
      t = temp-1
    
    laty = datF[t,"lat"]
    lony = datF[t,"lon"]
    #dist = distHaversine(c(lony,laty), c(lonx,latx))
    
    dist = distHaversine(c(lony,laty), c(lonx,latx))
    
    ## Fins the distance from that day's roost and if it is smaller than 50 count in as on roost dropped seed
    roost_lat= datR[match(day,datR$date.time.1),"lat"]
    roost_lon= datR[match(day,datR$date.time.1),"lon"]
    
    dist_roost = distHaversine(c(lony,laty), c(roost_lon,roost_lat))
    
    if(dist_roost < 50){
      final[k,"type"] = "Roost"
    }else{
      final[k,"type"] = "Other"
    }
    
    final[k,"lat"] = latx
    final[k,"long"] = lonx
    final[k,"dist"] = dist

    k = k + 1
    j=j+1
  }
}
write.csv(final,"Roost_Seed_Dispersal_Gabbar.csv")

tiff(file = paste("Roost_Seed_Dispersal_",Hname,".tiff"),width = 1420, height = 840,pointsize = 18)
barplot(table(final$type),beside = TRUE,xlab = "Category (Inside roost or other locations)",ylab = "Probability of seed arrival", main=Hname)
dev.off()
