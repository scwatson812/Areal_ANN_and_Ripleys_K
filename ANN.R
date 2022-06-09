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

###Specify the data.folder.  This should be the name of the subfolder containing the shapefile which contains the data
###This shapefile must contain a variable called 'Observations' with a value 1 if a unit is an observed location and 0 otherwise.
data.folder = 'cb_2020_us_county_20m'

###Specify the number of Monte Carlo Simulations
nsim = 999

###Specify the level of the hypothesis test
alpha = 0.05

###Specify the direction of the test, direction must be one of 'lower', 'upper', 'two.sided'
direction = 'two.sided'


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

########################################################
##### Compute the ANN Ratio for the Observed Data  #####
########################################################

centroids <- gCentroid(shp[shp$Observations ==1,],byid=TRUE)
cord <- as.data.frame(coordinates(centroids))
centroids.ppp <- ppp(cord$x,cord$y,window = as.owin(shp1))
ANN.z.data <- nni(centroids ,win="extent")$z.score 

################################################
##### Perform the Monte Carlo Simulations  #####
################################################

ann.mat <- matrix(NA, nsim)
for (i in 1:nsim){
  ind <- sample(x = seq(1,n.units),n.selected.points,replace = FALSE)
  shp$selected<- 0
  shp$selected[ind] = 1
  centroids <- gCentroid(shp[shp$selected ==1,],byid=TRUE)
  cord <- as.data.frame(coordinates(centroids))
  centroids.ppp <- ppp(cord$x,cord$y,window = as.owin(shp1))
  ANN.z <- nni(centroids ,win="extent")$z.score 
  ann.mat[i] <- ANN.z
}

########################################
##### Perform the Hypothesis Test  #####
########################################

if(direction == 'lower'){
  cut.off = quantile(ann.mat,prob = alpha)
  if(ANN.z.data < cut.off){ 
    print(paste0("Significant Clustering Detected"))
  }else{
    print("Fail to Reject the Null Hypothesis")
  }
}

if(direction == 'upper'){
  cut.off = quantile(ann.mat,prob = 1-alpha)
  if(ANN.z.data > cut.off){ 
    print(paste0("Significant Dispersion Detected"))
  }else{
    print("Fail to Reject the Null Hypothesis")
  }
}

if(direction == 'two.sided'){
  cut.off.lower = quantile(ann.mat,prob = alpha/2)
  cut.off.upper = quantile(ann.mat,prob = 1-alpha/2)
  if( ANN.z.data < cut.off.lower | ANN.z.data > cut.off.upper ){ 
    print(paste0("Significant Clustering or Dispersion Detected"))
  }else{
    print("Fail to Reject the Null Hypothesis")
  }
}
