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
library(gridExtra)
library(grid)

#?rgdal::set_thin_PROJ6_warnings()

#Set working directory
dir <- "C:/Users/rnawa/Desktop/GEOG 418/Final/Census"
setwd(dir)

#Reading in particulate matter dataset
#Read in PM2.5 data:
pm2.5 <- readOGR(dsn =dir, layer = "Pm25Sample") 
pm2.5 <- spTransform(pm2.5, CRS("+init=epsg:26910"))

#Reading in dissemination tract and income data
#Read in census income data:
income <- read.csv("Income.csv", header = T, sep = ",")  
#Select only ID and Income columns:
colnames(income) <- c("DAUID", "Income")

#Read in dissemination tract shapefile:
census.tracts <- readOGR(dsn = dir, layer = "BC_DA")
#Merge income and dissemination data:
income.tracts <- merge(census.tracts,income, by = "DAUID") 
#Determine the number of columns in the dataframe:
nrow(income.tracts)



#Remove NA values:
#have to change this >>>>>
income.tracts <- income.tracts[!is.na(income.tracts$Income),]


#Reproject the data:
income.tracts <- spTransform(income.tracts, CRS("+init=epsg:26910"))

#Study Area
tm_shape(income.tracts)+tm_polygons()+ tm_layout(title = "Study Area", title.size = 1)+
  tm_compass(position=c("right","top"), size = 0.5)+
  tm_scale_bar(position = c("left","bottom"), size = 0.3)


tmap_mode("plot")
#Create choropleth map of income:
map_Income <- tm_shape(income.tracts) +
  tm_polygons(col = "Income",
              title = "Median Income",
              style = "jenks",
              palette = "inferno", n = 6) +
  tm_legend(legend.outside=TRUE)
map_Income



pm2.5 <-  pm2.5[which(pm2.5$PM25 > 0), ]
pm2.5
#This is similar to what we did in Lab4

#Create a grid called grd to use in your interpolation
# Create an empty grid where n is the total number of cells
grd <- as.data.frame(spsample(pm2.5, "regular", n=5000))
names(grd)       <- c("X", "Y")
coordinates(grd) <- c("X", "Y")
# Create SpatialPixel object:
gridded(grd)     <- TRUE  
# Create SpatialGrid object:
fullgrid(grd)    <- TRUE  
#Reproject the grid:
proj4string(grd) <- proj4string(income.tracts)





#Step 1: (Analysis #2)
#Have to do moran's I (Global and Local) only on income for objective 1 
#and result of regression for objective 2.
#To do local still (LISA TEST)
#Spatial Seg of Income (Global/Local moran's i) (lab3)







#Stepp 2:
#Creates an interpolated surface for income and Pm2.5 data using IDW
#Spatial Interp. (lab4)
#Don't have to use all methods just use best one you like
#Justify for it
#uses grid created before

P.idw <- gstat::idw(PM25~ 1, pm2.5, newdata=grd, idp=3)
r       <- raster(P.idw)
idwrast     <- mask(r, income.tracts)

tm_shape(idwrast) + 
 tm_raster(n=6,palette = "-RdBu",
            title="Predicted PM2.5 \n(in ppm)") + 
  tm_shape(income.tracts) + tm_dots(size=0.05) +
  tm_legend(legend.outside=TRUE)




#################################################
# Leave-one-out validation routine
IDW.out <- vector(length = length(income.tracts))
for (i in 1:length(income.tracts)) {
  IDW.out[i] <- idw(Income ~ 1, income.tracts[-i,], income.tracts[i,], idp=2)$var1.pred
}

# Plot the differences
OP <- par(pty="s", mar=c(4,3,0,0))
plot(IDW.out ~ income.tracts$Income, asp=1, xlab="Observed", ylab="Predicted", pch=16,
     col=rgb(0,0,0,0.5))
abline(lm(IDW.out ~ income.tracts$Income), col="red", lw=2,lty=2)
abline(0,1)
par(OP)
sqrt( sum((IDW.out - income.tracts$Income)^2) / length(income.tracts))



#These steps will help you combine the outputs 
#from your spatial interpolation with your income data.
# Convert your interpolation into a raster and map it:
sufaceMap <- tm_shape(idwrast) + 
  tm_raster(n=5,palette = "viridis",
            title="PM 2.5 \n(in ppm)") +
  tm_shape(income.tracts) + tm_dots(size=0.2)
sufaceMap
#If you have too many cells, 
#you can reduce the number by aggregating values
#agg <- aggregate(yourRasterFromKriging, fact=??, fun=mean)

#Extract average pm2.5 for each polygon
income.tracts$Pm2.5 <- round(extract(idwrast, income.tracts, fun = mean)[,1], 5)




#Step 3
######Linear Regression Analysis##########
#Let's say your dataset with both PM2.5 and Income 
#are stored in a dataset called income.tracts.
#Plot income and PM2.5 from the income.tracts dataset you created
plot(income.tracts$Income~income.tracts$Pm2.5)

#Notice that there are a lot of 0's in this dataset. If you decide to remove them, use the following line:
income.tracts.no0 <-  income.tracts[which(income.tracts$Pm2.5 > 0), ]

#Now plot the data again
plot(income.tracts.no0$Income~income.tracts.no0$Pm2.5)

#Perform a linear regression on the two variables. You should decide which one is dependent.
lm.model <- lm(income.tracts.no0$Income~income.tracts.no0$Pm2.5)
#Add the regression model to the plot you created
plot(income.tracts.no0$Income~income.tracts.no0$Pm2.5)
abline(lm.model, col = "red")
#Get the summary of the results
summary(lm.model)

#add the fitted values to your spatialpolygon dataframe
income.tracts.no0$predictlm <- lm.model$fitted.values

#You want to determine if the model residuals are spatially clustered. 
#add the residuals to your spatialpolygon dataframe
income.tracts.no0$residuals <- residuals.lm(lm.model)

#Observe the result to make sure it looks correct
head(income.tracts.no0)
income.tracts.no0 <- income.tracts.no0[!is.na(income.tracts.no0$residuals),]


#Now, create choropleth map of residuals
map_resid <- tm_shape(income.tracts.no0) +
  tm_polygons(col = "residuals",
              title = "Residuals Choropleth",
              style = "jenks",
              palette = "PiYG", n = 6) + tm_legend(legend.outside=TRUE) 

map_resid




#Step 4 Use the residuals from above and do global moran's i
residuals.nb <- poly2nb(income.tracts.no0, queen = FALSE)
residuals.lw <- nb2listw(residuals.nb, zero.policy = TRUE, style = "W")

miIncome <- moran.test(income.tracts.no0$residuals, residuals.lw, zero.policy = TRUE)
miIncome

moran.range <- function(lw) {
  wmat <- listw2mat(lw)
  return(range(eigen((wmat + t(wmat))/2)$values))
}   
moran.range(residuals.lw)


mI <- miIncome$estimate[[1]]
eI <- miIncome$estimate[[2]]
var <- miIncome$estimate[[3]]

z <- (mI - eI)/sqrt(var)
z




#Steps 2-4 are for objective 2

#Step 5 do GWR (Takes a while)
####Geographically Weighted Regression
#Let's say you are continuing with 
#your data from the regression analysis. 
#The first thing you need to do is to add the 
#polygon coordinates to the spatialpolygondataframe.
#You can obtain the coordinates using the 
#"coordinates" function from the sp library
income.tracts.no0.coords <- sp::coordinates(income.tracts.no0)
#Observe the result:
head(income.tracts.no0.coords)
#Now add the coordinates back to the spatialpolygondataframe
income.tracts.no0$X <- income.tracts.no0.coords[,1]
income.tracts.no0$Y <- income.tracts.no0.coords[,2]

###Determine the bandwidth for GWR: this will take a while
GWRbandwidth <- gwr.sel(income.tracts.no0$Income~income.tracts.no0$Pm2.5, 
                        data=income.tracts.no0, coords=cbind(income.tracts.no0$X,income.tracts.no0$Y),adapt=T) 

###Perform GWR on the two variables with the bandwidth determined above
###This will take a looooooong while
gwr.model = gwr(income.tracts.no0$Income~income.tracts.no0$Pm2.5, 
                data=income.tracts.no0, coords=cbind(income.tracts.no0$X,income.tracts.no0$Y), 
                adapt=GWRbandwidth, hatmatrix=TRUE, se.fit=TRUE) 

#Print the results of the model
gwr.model

#Look at the results in detail
results<-as.data.frame(gwr.model$SDF)
head(results)

#Now for the magic. Let's add our local r-square values to the map
income.tracts.no0$localr <- results$localR2

#Create choropleth map of r-square values
map_r2 <- tm_shape(income.tracts.no0) +
  tm_polygons(col = "localr", midpoint = NA,
              title = "R2 values",
              style = "jenks",
              palette = "PiYG", n = 5) + tm_legend(legend.outside=TRUE)
map_r2

#Time for more magic. Let's map the coefficients
income.tracts.no0$coeff <- results$income.tracts.no0.Pm2.5
#Create choropleth map of the coefficients
map_coef <- tm_shape(income.tracts.no0) +
  tm_polygons(col = "coeff",
              title = "Coefficients",
              style = "jenks",
              palette = "-RdBu", n = 6)+tm_legend(legend.outside=TRUE)
map_coef

library("gtable")
#Descriptive Stats:
samples = c("Income", "Pm2.5")
#Mean:
meanI = mean(income.tracts.no0$Income)
meanP = mean(income.tracts.no0$Pm2.5)
means = c(meanI, meanP)
#SD
sdI = sd(income.tracts.no0$Income)
sdP = sd(income.tracts.no0$Pm2.5)
sds = c(sdI, sdP)
#Mode
modeI = as.numeric(names(sort(table(income.tracts.no0$Income), decreasing = TRUE))[1])
modeP = as.numeric(names(sort(table(income.tracts.no0$Pm2.5), decreasing = TRUE))[1])
modes = c(modeI, modeP)
#Median
medI = median(income.tracts.no0$Income)
medP = median(income.tracts.no0$Pm2.5)
medians = c(medI, medP)
#Skewness
skewI = skewness(income.tracts.no0$Income)
skewP = skewness(income.tracts.no0$Pm2.5)
skewness = c(skewI, skewP)
#Kurtosis
kurtI = kurtosis(income.tracts.no0$Income)
kurtP = kurtosis(income.tracts.no0$Pm2.5)
kurtosis = c(kurtI, kurtP)
#CoV
CoVI = (sdI/meanI) * 100
CoVP = (sdP/meanP) * 100
CoV = c(CoVI, CoVP)
#Normal Distribution Test
normI = shapiro.test(income.tracts.no0$Income)$p.value
normP = shapiro.test(income.tracts.no0$Pm2.5)$p.value
normality = c(normI, normP)
normality
#Have to round all values still
means = round(means, 3)
medians = round(medians, 3)
modes = round(modes, 3)
sds = round(sds, 3)
skewness = round(skewness, 3)
kurtosis = round(kurtosis, 3)
CoV = round(CoV, 3)
normality = signif(normality, digits = 3)

data.for.table1 = data.frame(samples, means, medians, modes,sds, skewness, kurtosis, CoV, normality)
data.for.table2 = data.frame(samples, skewness, kurtosis, CoV, normality)

#Make table 1
table1 <- tableGrob(data.for.table1, rows = c("","")) #make a table "Graphical Object" (GrOb) 
t1Caption <- textGrob("Table 1: Descriptive Data for Pm2.5 and Income", gp = gpar(fontsize = 09))
padding <- unit(5, "mm")

table1 <- gtable_add_rows(table1, 
                          heights = grobHeight(t1Caption) + padding, 
                          pos = 0)

table1 <- gtable_add_grob(table1,
                          t1Caption, t = 1, l = 2, r = ncol(data.for.table1) + 1)


grid.arrange(table1, newpage = TRUE)


#Printing a table (You can use the same setup for printing other types of objects (see ?png))
png("Output_Table1.png") #Create an object to print the table to
grid.arrange(table1, newpage = TRUE)
dev.off() #Print table


#After Description perform analysis on income (Analysis #2)
income.nb <- poly2nb(income.tracts.no0)
income.lw <- nb2listw(income.nb, zero.policy = TRUE, style = "W")

miIncome <- moran.test(income.tracts.no0$Income, income.lw, zero.policy = TRUE)
miIncome

moran.range <- function(lw) {
  wmat <- listw2mat(lw)
  return(range(eigen((wmat + t(wmat))/2)$values))
}   
moran.range(income.lw)


mI <- miIncome$estimate[[1]]
eI <- miIncome$estimate[[2]]
var <- miIncome$estimate[[3]]

z <- (mI - eI)/sqrt(var)
z

income.net <- nb2lines(income.nb, coords=coordinates(income.tracts.no0))
crs(income.net) <- crs(income.tracts.no0)

tm_shape(income.tracts.no0) + tm_borders(col='lightgrey') + 
  tm_shape(income.net) + tm_lines(col='red')









lisa.test <- localmoran(income.tracts.no0$Income, income.lw, zero.policy = TRUE)

income.tracts.no0$Ii <- lisa.test[,1]
income.tracts.no0$E.Ii<- lisa.test[,2]
income.tracts.no0$Var.Ii<- lisa.test[,3]
income.tracts.no0$Z.Ii<- lisa.test[,4]
income.tracts.no0$P<- lisa.test[,5]
########################
tmaptools::palette_explorer()
map_LISA <- tm_shape(income.tracts.no0) + 
  tm_polygons(col = "Z.Ii", midpoint = NA,
              title = "Z Score", 
              style = "fixed", breaks = c(-Inf, -1.96, 1.96, Inf), 
              palette = "RdBu") + tm_legend(legend.outside=TRUE) 


map_LISA



#Step 6 perform point pattern analysis

kma <- pm2.5
kma$x <- coordinates(kma)[,1] 
kma$y <- coordinates(kma)[,2]
kma <- remove.duplicates(kma)
kma.ext <- as.matrix(extent(kma)) 
window <- as.owin(list(xrange = kma.ext[1,], yrange = kma.ext[2,]))
kma.ppp <- ppp(x = kma$x, y = kma$y, window = window)
sarea = (kma.ext[3]-kma.ext[1])*(kma.ext[4]-kma.ext[2]) 
points = nrow(kma)
points
quads <- 4

qcount <- quadratcount(kma.ppp, nx = quads, ny = quads)
qcount
#Run these two together 
plot(kma.ppp, pch = "+", cex = 0.5)
plot(qcount, add = T, col = "red")

qcount.df <- as.data.frame(qcount)

##Second, count the number of quadrats with a distinct number of points.
qcount.df <- plyr::count(qcount.df,'Freq')
##Change the column names so that x=number of points and f=frequency of quadrats with x point.
colnames(qcount.df) <- c("x","f")

sum.f.x2 <- sum((qcount.df["x"]^2)*(qcount.df["f"]))
M <- sum(qcount.df["f"])
N <- sum(qcount.df["x"]* qcount.df["f"])
sum.fx.2 <- (N)^2
VAR <- (sum.f.x2-(sum.fx.2/M))/(M-1)
MEAN <- N/M
VMR <- VAR/MEAN

chi.square = VMR*(M-1)

p = 1 - pchisq(chi.square, (M - 1))

XS = qcount.df["x"]
FS = qcount.df["f"]
Variables <- c("fx^2", "(fx)^2", "Mean", "VMR", "p value")
VAR
VMR
MEAN
N
M
P

Results <- c(sum.f.x2, sum.fx.2, MEAN, VMR, p)
Results
