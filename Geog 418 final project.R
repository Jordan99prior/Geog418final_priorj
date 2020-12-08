#Libraries
library(spgwr)
library(spatstat)
library(tmap)
library(gstat)
library(sf)
library(raster)
library(rgdal)
library(e1071)
library(spdep)

library(gtable)
library(gridExtra)
library(grid)
library(ggplot2)
#Set working directory
dir <- "C:/Users/jorda/OneDrive/Documents/UVIC year 4/Goeg 418/Final project"
setwd(dir)

#Reading in elevation dataset
elev <- readOGR("ElevSample.shp") #Read in data
elev <- spTransform(elev, CRS("+init=epsg:26910"))

#Reading in VRI data
VRI <- readOGR("WatershedVRI.shp") #Read in shapefile
VRI <- spTransform(VRI, CRS("+init=epsg:26910"))
head(VRI@data)

#chooses which variables to keep
vriCleanCols <- c("FID_VEG_CO", "POLYGON_ID", 
                  "SITE_INDEX","FOLIAGE_BI",
                  "Shape_Leng","Shape_Area")

# clips data to selected columns
vriClean <- VRI[,vriCleanCols] 

#Create choropleth map of height
map_Bio <- tm_shape(vriClean) +
              tm_polygons(col = "FOLIAGE_BI",
              title = "Foliage Biomass (tonnes/ha)",
              style = "jenks",
              palette = "-RdBu", n = 6) +
  tm_legend(legend.position = c("LEFT", "BOTTOM"))+
  tm_shape(elev)+ tm_dots(size = 0.1)+
  tm_compass(north = 0, cardinal.directions = c("N","E","S","W"), position = c("LEFT", "BOTTOM"))+
  tm_scale_bar(position = c("LEFT", "BOTTOM"))+
  tm_layout(title = "Map of Foliage Biomass in Greater Victoria Water Supply Area", title.position = c("LEFT", "TOP"))

map_Bio
plot(elev)

#Create a grid called grd to use in your interpolation
# Create an empty grid where n is the total number of cells
grd <- as.data.frame(spsample(elev, "regular", n=50000))
names(grd)       <- c("X", "Y")
coordinates(grd) <- c("X", "Y")
gridded(grd)     <- TRUE  # Create SpatialPixel object
fullgrid(grd)    <- TRUE  # Create SpatialGrid object
proj4string(grd) <- proj4string(elev)
##################################################

#Descriptive stats:
#Mean
meanFOLIAGE <- mean(vriClean$FOLIAGE_BI, na.rm = TRUE) #Use na.rm = TRUE to ignore NA values in calculation
meanelev <- mean(elev$grid_code, na.rm = TRUE)

#Standard Deviation
sdFOLIAGE <- sd(vriClean$FOLIAGE_BI, na.rm = TRUE) #Calculate the SD, ignoring NA values
sdelev <- sd(elev$grid_code, na.rm = TRUE)

#Mode
modeFOLIAGE <- as.numeric(names(sort(table(vriClean$FOLIAGE_BI), decreasing = TRUE))[1]) #make frequency table of fire size variable and sort it in desending order and extract the first row (Most Frequent)
modeelev <- as.numeric(names(sort(table(elev$grid_code), decreasing = TRUE))[1])

#Median
medFOLIAGE <- median(vriClean$FOLIAGE_BI, na.rm = TRUE)
medelev <- median(elev$grid_code, na.rm = TRUE)

#Skewness
skewFOLIAGE <- skewness(vriClean$FOLIAGE_BI, na.rm = TRUE)[1]
skewelev <- skewness(elev$grid_code, na.rm = TRUE)

#Kurtosis
kurtFOLIAGE <- kurtosis(vriClean$FOLIAGE_BI, na.rm = TRUE)[1]
kurtelev <- kurtosis(elev$grid_code, na.rm = TRUE)
#CoV
CoVFOLIAGE <- (sdFOLIAGE / meanFOLIAGE) * 100
CoVelev <- (sdelev / meanelev) * 100

#Normal distribution test
qqnorm(vriClean$FOLIAGE_BI, pch = 1, frame = FALSE)
qqline(vriClean$FOLIAGE_BI, col = "steelblue", lwd = 2)
qqnorm(elev$grid_code, pch = 1, frame = FALSE)
qqline(elev$grid_code, col = "steelblue", lwd = 2)


#####
#Create a table of descriptive stats 

samples = c("Foilage Biomass", "Elevation") #Create an object for the labels
means = c(meanFOLIAGE, meanelev) #Create an object for the means
sd = c(sdFOLIAGE, sdelev)
median = c(medFOLIAGE, medelev) #Create an object for the medians
mode <- c(modeFOLIAGE, modeelev) #Create an object for the modes
skewness <- c(skewFOLIAGE, skewelev) #Create an object for the skewness
kurtosis <- c(kurtFOLIAGE, kurtelev)#Create an object for the kurtosis
CoV <- c(CoVFOLIAGE, CoVelev) #Create an object for the CoV

##Check table values for sigfigs?
means <- round(means, 3)
sd <- round(sd, 3)
median <- round(median, 3)
mode<- round(mode, 3)
skewness <- round(skewness, 3)
kurtosis <- round(kurtosis, 3)
CoV <- round(CoV, 3)

data.for.table1 = data.frame(samples, means, sd, median, mode)
data.for.table2 = data.frame(samples, skewness, kurtosis, CoV)


#Make table 1
table1 <- tableGrob(data.for.table1, rows = c("","")) #make a table "Graphical Object" (GrOb) 
t1Caption <- textGrob("Table 1: Table shows means, sd, median, and mode 
for Foilage Biomass(tonnes/ha) and Elevation (m)
in Greater Victoria Water Supply Area", gp = gpar(fontsize = 09))
padding <- unit(5, "mm")

table1 <- gtable_add_rows(table1, 
                          heights = grobHeight(t1Caption) + padding, 
                          pos = 0)

table1 <- gtable_add_grob(table1,
                          t1Caption, t = 1, l = 2, r = ncol(data.for.table1) + 1)


table2 <- tableGrob(data.for.table2, rows = c("",""))
t2Caption <- textGrob("Table 2: Table shows skewness, kurtosis, CoV, and normality of 
Foilage Biomass(tonnes/ha)and Elevation(m)
in Greater Victoria Water Supply Area", gp = gpar(fontsize = 09))
padding <- unit(5, "mm")

table2 <- gtable_add_rows(table2, 
                          heights = grobHeight(t2Caption) + padding, 
                          pos = 0)

table2 <- gtable_add_grob(table2,
                          t2Caption, t = 1, l = 2, r = ncol(data.for.table2) + 1)

grid.arrange(table1, newpage = TRUE)

grid.arrange(table2, newpage = TRUE)

##################################################
#obj 1: Global morans/ local morans- forest variable

#define our neighbours for the moran's I using defult or Queens
vri.nb <- poly2nb(vriClean)
vri.net <- nb2lines(vri.nb, coords=coordinates(vriClean))
crs(vri.net) <- crs(vriClean)

#map weight matrix using Queens
tm_shape(vriClean) + tm_borders(col='lightgrey') + 
  tm_shape(vri.net) + tm_lines(col='red') +
  tm_layout(title = "Map of queen's neighbours for Greater Victoria Water Supply Area", title.position = c("LEFT", "BOTTOM"))

#creates weight matrix from neighbours
vri.lw <- nb2listw(vri.nb, zero.policy = TRUE, style = "W")
print.listw(vri.lw, zero.policy = TRUE)
#creates global moran's I
mi <- moran.test(vriClean$FOLIAGE_BI, vri.lw, zero.policy = TRUE) #uses weight matrix and your chosen variable
mi

#creates the range of GLobal Moran's i assuming random or dispersed
moran.range <- function(lw) {
  wmat <- listw2mat(lw)
  return(range(eigen((wmat + t(wmat))/2)$values))
}
moran.range(vri.lw) 

#calculations/values derived from GLobal Moran's i
mI <- mi$estimate[[1]] #moran's I global value
eI <- mi$estimate[[2]] #Expected Moran's I
var <- mi$estimate[[3]] #variance
z <- (mI - eI)/sqrt(var) # Z-value equation

mI <- signif(mI,3)
eI <- signif(eI,3)
var <- signif(var,3)
z <- signif(z,3)

samples2 <- c("Elevation")
data.for.table3 <-data.frame(samples2, mI, eI, var, z)

table3 <- tableGrob(data.for.table3, rows = c(""))
t2Caption <- textGrob("Table 3: Table shows Moran's global value, expected moran's i 
Variance and Z-score of Elevation(m)
in Greater Victoria Water Supply Area", gp = gpar(fontsize = 09))
padding <- unit(5, "mm")

table3 <- gtable_add_rows(table3, 
                          heights = grobHeight(t2Caption) + padding, 
                          pos = 0)

table3 <- gtable_add_grob(table3,
                          t2Caption, t = 1, l = 2, r = ncol(data.for.table3) + 1)
grid.arrange(table3, newpage = TRUE)
#create local Moran's I test
lisa.test <- localmoran(vriClean$FOLIAGE_BI, vri.lw, zero.policy = TRUE)

vriClean$Ii <- lisa.test[,1] #local moran's I
vriClean$E.Ii<- lisa.test[,2] #expected moran's I
vriClean$Var.Ii<- lisa.test[,3] #variance
vriClean$Z.Ii<- lisa.test[,4] #z-value
vriClean$P<- lisa.test[,5] # p-value

map_LISA <- tm_shape(vriClean) + 
  tm_polygons(col = "Z.Ii", #use Z-value to help discern whether regions is significantly correalted 
              title = "Local Moran's I",
              style = "fixed", breaks = c(-Inf, -1.96, 0, 1.96, +Inf), #uses these values because they are reperesentative of significance 
              palette = "-Spectral", n = 5) + #Palette with N=5 due to above breaks
  tm_layout(title = "Map of Local Moran's I for Foliage Biomass in Greater Victoria Water Supply Area", title.position = c("LEFT", "TOP"))

map_LISA
########################
# creates a plot showing values relative to thier spatially lagged variable
moran.plot(vriClean$FOLIAGE_BI, vri.lw, zero.policy=TRUE, spChk=NULL, labels=NULL, xlab="Foliage Biomas", 
           ylab="Spatially Lagged Variable", quiet=NULL)
##################################################
#Obj 2: Spatial interpolation -IDW- elevation data
##Spatial Interpolation with IDW

#IDW Interpolation
P.idw <- gstat::idw(grid_code ~ 1, elev, newdata=grd, idp= 9)
r       <- raster(P.idw)
r.m     <- mask(r, vriClean)

tm_shape(r.m) + 
  tm_raster(n=10,palette = "-Spectral",
            title="Prediced elevation \n(in m)") + 
  tm_shape(elev) + tm_dots(size=0.2) +
  tm_legend(legend.outside=TRUE)+
  tm_layout(title = "Map of Greater Victoria Water \nSupply Area and elevation \nusing IDW", title.position = c("LEFT", "TOP"))

#################################################
# Leave-one-out validation routine
IDW.out <- vector(length = length(elev))
for (i in 1:length(elev)) {
  IDW.out[i] <- idw(grid_code ~ 1, elev[-i,], elev[i,], idp=9)$var1.pred
}

# Plot the differences
OP <- par(pty="s", mar=c(4,3,0,0))
plot(IDW.out ~ elev$grid_code, asp=1, xlab="Observed", ylab="Predicted", pch=16,
     col=rgb(0,0,0,0.5)) 
abline(lm(IDW.out ~ elev$grid_code), col="red", lw=2,lty=2)
abline(0,1)
title(main ="\n\nPlot of Validation routine for Greater Victoria \nWater Supply Area using IDW")
par(OP)
sqrt( sum((IDW.out - elev$grid_code)^2) / length(elev$grid_code))

#################################################
# Implementation of a jackknife technique to estimate a confidence interval at each unsampled point.
# Create the interpolated surface
img <- gstat::idw(grid_code~1, elev, newdata=grd, idp=9)
n   <- length(elev)
Zi  <- matrix(nrow = length(img$var1.pred), ncol = n)

# Remove a point then interpolate (do this n times for each point)
st <- stack()
for (i in 1:n){
  Z1 <- gstat::idw(grid_code~1, elev[-i,], newdata=grd, idp=9)
  st <- addLayer(st,raster(Z1,layer=1))
  # Calculated pseudo-value Z at j
  Zi[,i] <- n * img$var1.pred - (n-1) * Z1$var1.pred
}

# Jackknife estimator of parameter Z at location j
Zj <- as.matrix(apply(Zi, 1, sum, na.rm=T) / n )

# Compute (Zi* - Zj)^2
c1 <- apply(Zi,2,'-',Zj)            # Compute the difference
c1 <- apply(c1^2, 1, sum, na.rm=T ) # Sum the square of the difference

# Compute the confidence interval
CI <- sqrt( 1/(n*(n-1)) * c1)

# Create (CI / interpolated value) raster
img.sig   <- img
img.sig$v <- CI /img$var1.pred 

# Clip the confidence raster to Southern California
r2 <- raster(img.sig, layer="v")
r.m2 <- mask(r2, vriClean)

# Plot the mape
tm_shape(r.m2) + tm_raster(n=7, title="95% confidence interval \n(in m)") +
  tm_legend(legend.outside=TRUE)+
  tm_layout(title = "Map of 95% confidence interval\nGreater Victoria Water \nSupply Area and elevation\n using IDW", title.position = c("LEFT", "TOP"))

##################################################
#obj 3: Point pattern- point sampling
kma <- elev
kma$x <- coordinates(kma)[,1]
kma$y <- coordinates(kma)[,2]

#check for and remove duplicated points
#first, finds zero distance among points to see if there are any duplicates
zd <- zerodist(kma)
zd

#if there are duplicates, remove them
kma <- remove.duplicates(kma)

#create an "extent" object which can be used to create the observation window for spatstat
kma.ext <- as.matrix(extent(kma)) 

#observation window
window <- as.owin(list(xrange = kma.ext[1,], yrange = kma.ext[2,])) #use extent of city boundary
#create ppp oject from spatstat
kma.ppp <- ppp(x = kma$x, y = kma$y, window = window)
nearestNeighbour <- nndist(kma.ppp)
##Convert the nearestNeighbor object into a dataframe.
nearestNeighbour=as.data.frame(as.numeric(nearestNeighbour))
##Change the column name to "Distance"
colnames(nearestNeighbour) = "Distance"


##Calculate the nearest neighbor statistic to test for a random spatial distribution.
#mean nearest neighbour
Foliage_sum = NROW(kma)
nnd_sum = sum(nearestNeighbour)
nnd = nnd_sum/Foliage_sum

#mean nearest neighbour for random spatial distribution

studyArea <- (kma.ext[1,1] - kma.ext[1,2]) * (kma.ext[2,1] - kma.ext[2,2])
#sqmeters
pointDensity <-(Foliage_sum/studyArea)

r.nnd = 1/(2*sqrt(pointDensity))

d.nnd = 1.07453/sqrt(pointDensity)

R = nnd/d.nnd

SE.NND <- 0.26136/(sqrt(Foliage_sum*(pointDensity)))

z = (nnd-r.nnd)/SE.NND

nnd <- signif(nnd, digits = 3)
r.nnd <- signif(r.nnd, digits = 3)  
d.nnd <- signif(d.nnd, digits = 3)  
R <- signif(R, digits = 3) 
SE.NND <- signif(SE.NND, digits = 3)  
z <- signif(z, digits = 3)  
data.for.table4 <-data.frame(nnd, r.nnd, d.nnd, R, SE.NND, z)
samples3 = c("Nearest neighbour") 

data.for.table4 <-data.frame(samples3, nnd, r.nnd, d.nnd, R, z)

table4 <- tableGrob(data.for.table4, rows = c("")) #make a table "Graphical Object" (GrOb) 
t2Caption <- textGrob("Table 4: Table shows NND, random and distributed NND, R and 
z-score in Greater Victoria Water Supply Area", gp = gpar(fontsize = 09))
padding <- unit(5, "mm")

table4 <- gtable_add_rows(table4, 
                          heights = grobHeight(t2Caption) + padding, 
                          pos = 0)

table4 <- gtable_add_grob(table4,
                          t2Caption, t = 1, l = 2, r = ncol(data.for.table4) + 1)
grid.arrange(table4, newpage = TRUE)
##################################################
#Combining VRI data
#These steps will help you combine the outputs 
#from your spatial interpolation with your income data.
#Convert your interpolation into a raster and map it:
sufaceMap <- tm_shape(r.m) + 
  tm_raster(n=5,palette = "viridis",
            title="Elev (m)") +
  tm_shape(elev) + tm_dots(size=0.2)
sufaceMap

#If you have too many cells, 
#you can reduce the number by aggregating values
#agg <- aggregate(yourRasterFromKriging, fact=??, fun=mean)

#Extract average elev for each polygon
vriClean$elev <- extract(r.m, vriClean, fun = mean)[,1]
##################################################

######Linear Regression##########
#Let's say your dataset with both Elev and Height are stored in a dataset called VRI.
#Plot Height and Elev from the VRI dataset you created
plot(vriClean$FOLIAGE_BI ~ vriClean$elev)

#Notice that there are a lot of 0's in this dataset. If you decide to remove them, use the following line:
VRI.no0 <-  vriClean[which(vriClean$FOLIAGE_BI > 0), ]
VRI.no0 <-  vriClean[which(vriClean$elev > 0), ]

#Now plot the data again
plot(VRI.no0$FOLIAGE_BI ~ VRI.no0$elev)

#Perform a linear regression on the two variables. You should decide which one is dependent.
lm.model <- lm(VRI.no0$FOLIAGE_BI ~ VRI.no0$elev)

#Add the regression model to the plot you created
plot(VRI.no0$FOLIAGE_BI ~ VRI.no0$elev)
abline(lm.model, col = "red")
title(main ="\n\nPlot of Foliage biomas versus elevation for the lienar regresion model")

#Get the summary of the results
summary(lm.model)

#add the fitted values to your spatialpolygon dataframe
VRI.no0$predictlm <- lm.model$fitted.values

#You want to determine if the model residuals are spatially clustered. 
#add the residuals to your spatialpolygon dataframe
VRI.no0$residuals <- residuals.lm(lm.model)

#Observe the result to make sure it looks correct
head(VRI.no0@data)

#Now, create choropleth map of residuals
map_resid <- tm_shape(VRI.no0) +
  tm_polygons(col = "residuals",
              title = "Foliage Biomass Residuals",
              style = "jenks",
              palette = "viridis", n = 6)+
  tm_layout(title = "Map of Linear regression residuals for GVWSA", title.position = c("LEFT", "BOTTOM"))

map_resid

#global moran's I- independance of residuals (check assumptions)
#obj 1: Global morans/ local morans- forest variable

#define our neighbours for the moran's I using defult or Queens
vri.nb1 <- poly2nb(VRI.no0)
vri.net1 <- nb2lines(vri.nb, coords=coordinates(VRI.no0))
crs(vri.net1) <- crs(VRI.no0)

#map weight matrix using Queens
tm_shape(VRI.no0) + tm_borders(col='lightgrey') + 
  tm_shape(vri.net1) + tm_lines(col='red') +
  tm_layout(title = "Map of queen's neighbours for Greater Victoria Water Supply Area of residuals", title.position = c("LEFT", "BOTTOM"))

#creates weight matrix from neighbours
vri.lw1 <- nb2listw(vri.nb, zero.policy = TRUE, style = "W")
print.listw(vri.lw, zero.policy = TRUE)
#creates global moran's I
mi1 <- moran.test(VRI.no0$residuals, vri.lw, zero.policy = TRUE) #uses weight matrix and your chosen variable
mi1

#creates the range of GLobal Moran's i assuming random or dispersed
moran.range1 <- function(lw) {
  wmat <- listw2mat(lw)
  return(range(eigen((wmat + t(wmat))/2)$values))
}
moran.range(vri.lw1) 

#calculations/values derived from GLobal Moran's i
mI1 <- mi1$estimate[[1]] #moran's I global value
eI1 <- mi1$estimate[[2]] #Expected Moran's I
var1 <- mi1$estimate[[3]] #variance
z1 <- (mI1 - eI1)/sqrt(var) # Z-value equation

mI1 <- signif(mI1,3)
eI1 <- signif(eI1,3)
var1 <- signif(var1,3)
z1 <- signif(z1,3)

samples2 <- c("Elevation")
data.for.table5 <-data.frame(samples2, mI1, eI1, var1, z1)

table5 <- tableGrob(data.for.table5, rows = c(""))
t2Caption <- textGrob("Table 5: Table shows Moran's global value, expected moran's i 
Variance and z-value of Elevation(m) residuals
in Greater Victoria Water Supply Area", gp = gpar(fontsize = 09))
padding <- unit(5, "mm")

table5 <- gtable_add_rows(table5, 
                          heights = grobHeight(t2Caption) + padding, 
                          pos = 0)

table5 <- gtable_add_grob(table5,
                          t2Caption, t = 1, l = 2, r = ncol(data.for.table5) + 1)
grid.arrange(table5, newpage = TRUE)
##################################################
####Geographically Weighted Regression
#Let's say you are continuing with 
#your data from the regression analysis. 
#The first thing you need to do is to add the 
#polygon coordinates to the spatialpolygondataframe.
#You can obtain the coordinates using the 
#"coordinates" function from the sp library
VRI.no0.coords <- sp::coordinates(VRI.no0)
#Observe the result:
head(VRI.no0.coords)
#Now add the coordinates back to the spatialpolygondataframe
VRI.no0$X <- VRI.no0.coords[,1]
VRI.no0$Y <- VRI.no0.coords[,2]

###Determine the bandwidth for GWR: this will take a while
GWRbandwidth <- gwr.sel(VRI.no0$FOLIAGE_BI ~ VRI.no0$elev, 
                        data=VRI.no0, coords=cbind(VRI.no0$X,VRI.no0$Y),adapt=T) 

###Perform GWR on the two variables with the bandwidth determined above
###This will take a looooooong while
gwr.model = gwr(VRI.no0$FOLIAGE_BI ~ VRI.no0$elev, 
                data=VRI.no0, coords=cbind(VRI.no0$X,VRI.no0$Y), 
                adapt=GWRbandwidth, hatmatrix=TRUE, se.fit=TRUE) 

#Print the results of the model
gwr.model

#Look at the results in detail
results<-as.data.frame(gwr.model$SDF)
head(results)

#Now for the magic. Let's add our local r-square values to the map
VRI.no0$localr <- results$localR2

#Create choropleth map of r-square values
map_r2 <- tm_shape(VRI.no0) +
  tm_polygons(col = "localr",
              title = "R2 values",
              style ="fixed", breaks = c(-Inf, 0.170, 0.253, 0.365, 0.620),
              palette = "viridis", n = 6)+
  tm_layout(title = "Map of r-squared values in Greater Victoria Water Supply Area", title.position = c("LEFT", "TOP"))

map_r2

#Time for more magic. Let's map the coefficients
VRI.no0$coeff <- results$VRI.no0.elev
#Create choropleth map of the coefficients
map_coef <- tm_shape(VRI.no0) +
  tm_polygons(col = "elev",
              title = "Coefficients",
              style = "jenks",
              palette = "viridis", n = 6)+
  tm_layout(title = "Map Coefficients in Greater Victoria Water Supply Area", title.position = c("LEFT", "TOP"))

map_coef
