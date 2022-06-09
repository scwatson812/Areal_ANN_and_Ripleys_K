##########################
##### Load Libraries #####
##########################

library(rgdal)
library(rgeos)
library(spatialEco)
library(spatstat)
library(maptools)

#########################
#####  User Inputs  #####
#########################

###Specify the file path to the working directory.  The working directory should contain a subfolder containing the shapefiles holding the data
path = 'C:/Users/'
path = 'C:/Users/scwatson/Stella Self Dropbox/Stella Watson/Synced Files/Research/Spatial Clustering/Spatial Clustering Comparison/Transactions in GIS Submission/Simulations/Null Distributions/'

###Specify the data.folder.  This should be the name of the subfolder containing the shapefile which contains the data
###This shapefile must contain a variable called 'Observations' with a value 1 if a unit is an observed location and 0 otherwise.
data.folder = 'cb_2020_us_county_20m'

###Specify the number of Monte Carlo Simulations
nsim = 999

###Specify the level of the hypothesis test
alpha = 0.25

###Specify the direction of the test, direction must be one of 'lower', 'upper', 'two.sided'
direction = 'two.sided'

###Specify the minimum radius at which to evaluate Ripley's K Function
r.min = 0.1

###Specify the maximum radius at which to evaluate Ripley's K Function
r.max = 1

###Specify the number of radii at which to evaluate Ripley's K Function
n.r = 5

################################
##### PreProcess the Data  #####
################################

###Set the working directory
setwd(path)

###Read in the shapefile containing the data
###data.folder should be the name of the folder containing the shapefile
data.folder = 'cb_2020_us_county_20m'
shp <-readOGR(data.folder)

###Pre compute some objects
n.units = dim(shp)[1]
n.selected.points = sum(shp$Observations)
shp1 <-unionSpatialPolygons(shp,IDs = rep(1,n.units))

##############################################################
##### Compute Ripley's K Function for the Observed Data  #####
##############################################################

radii = c(0,seq(r.min, r.max,length.out = n.r))
centroids <- gCentroid(shp[shp$Observations ==1,],byid=TRUE)
cord <- as.data.frame(coordinates(centroids))
centroids.ppp <- ppp(cord$x,cord$y,window = as.owin(shp1))
RK.data <- Kest(centroids.ppp,r=radii, correction="iso")$iso 

################################################
##### Perform the Monte Carlo Simulations  #####
################################################

RK.mat <- matrix(NA, nsim,(n.r+1))
for (i in 1:nsim){
  ind <- sample(x = seq(1,n.units),n.selected.points,replace = FALSE)
  shp$selected<- 0
  shp$selected[ind] = 1
  centroids <- gCentroid(shp[shp$selected ==1,],byid=TRUE)
  cord <- as.data.frame(coordinates(centroids))
  centroids.ppp <- ppp(cord$x,cord$y,window = as.owin(shp1))
  RK <- Kest(centroids.ppp,r=radii, correction="iso")$iso 
  RK.mat[i,] <- RK
}

########################################
##### Perform the Hypothesis Test  #####
########################################

if(direction == 'lower'){
  cut.off = apply(RK.mat,2,quantile,prob = alpha)
  for(i in 2:(n.r+1)){
    if(RK.data[i] < cut.off[i]){ 
      print(paste0("Significant Dispersion Detected at radius ",radii[i]))
    }else{
      print(paste0("Fail to Reject the Null Hypothesis at radius ",radii[i]))
    }
  }
  plot(radii,cut.off,type = 'o',ylim = range(c(cut.off,RK.data)), ylab = "K(r)")
  points(radii,RK.data,col = 'blue',type = 'o')
  legend('topleft',legend = c("Upper Bound of Rejection Rejection","Observed Data"),col = c('black','blue'),lty = 1)
}

if(direction == 'upper'){
  cut.off = apply(RK.mat,2,quantile,prob = 1-alpha)
  for(i in 2:(n.r+1)){
    if(RK.data[i] > cut.off[i]){ 
      print(paste0("Significant Clustering Detected at radius ",radii[i]))
    }else{
      print(paste0("Fail to Reject the Null Hypothesis at radius ",radii[i]))
    }
  }
  plot(radii,cut.off,type = 'o',ylim = range(c(cut.off,RK.data)), ylab = "K(r)")
  points(radii,RK.data,col = 'blue',type = 'o')
  legend('topleft',legend = c("Lower Bound of Rejection Rejection","Observed Data"),col = c('black','blue'),lty = 1)
}


if(direction == 'two.sided'){
  cut.off.lower = apply(RK.mat,2,quantile,prob = alpha/2)
  cut.off.upper = apply(RK.mat,2,quantile,prob = 1-alpha/2)
  for(i in 2:(n.r+1)){
    if(RK.data[i] < cut.off.lower[i] |RK.data[i] > cut.off.upper[i] ){ 
      print(paste0("Significant Clustering or Dispersion Detected at radius ",radii[i]))
    }else{
      print(paste0("Fail to Reject the Null Hypothesis at radius ",radii[i]))
    }
  }
  plot(radii,cut.off.lower,type = 'o',ylim = range(c(cut.off.lower,cut.off.upper,RK.data)), ylab = "K(r)")
  points(radii,cut.off.upper,type = 'o')
  points(radii,RK.data,col = 'blue',type = 'o')
  legend('topleft',legend = c("Bounds of Rejection Rejection","Observed Data"),col = c('black','blue'),lty = 1)
  
}


