################################################################################
#Calculate AGB of tree individuals (different stage) using BIOMASS package
################################################################################
# Load necessary libraries
library(sf)
library(terra)

## Linking to GEOS 3.10.1, GDAL 3.4.0, PROJ 8.2.0; sf_use_s2() is TRUE
##Load data
#First load the LiDAR canopy heigh model (CHM):
CHMlidar = rast("E:/ALS_processing2017/02_processed/CHM_comb/CHM2017.tif")
#CHMlidar = rast("E:/ALS_processing2017/02_processed/CHM_comb/CHM_pitfree2017.tif")
CHMlidar[CHMlidar<0]=0
plot(CHMlidar)

#Second load the field AGB estimates
AGBplots_0.25ha=st_read("Output/polygons_res0.25ha.shp")
AGBplots_1ha=st_read("Output/polygons_res1ha.shp")



##Build a LiDAR and field AGB dataset
#First extract LiDAR data for each subplots:
meanTCH = extract(CHMlidar, AGBplots_0.25ha, mean, na.rm=T)
H50 = extract(CHMlidar, AGBplots_0.25ha, fun = function(x) quantile(x, probs = 0.5, na.rm = T))

AGBplots_0.25ha$meanTCH = meanTCH$CHM_pitfree
AGBplots_0.25ha$H50 = H50$CHM_pitfree

meanTCH = extract(CHMlidar, AGBplots_1ha, mean, na.rm=T)
H50 = extract(CHMlidar, AGBplots_1ha, fun = function(x) quantile(x, probs = 0.5, na.rm = T))

AGBplots_1ha$meanTCH = meanTCH$CHM_pitfree
AGBplots_1ha$H50 = H50$CHM_pitfree

##Choose data 
AGBplots = AGBplots_0.25ha
AGBplots = AGBplots_1ha

# Rename a specific column in AGBplots
colnames(AGBplots)[colnames(AGBplots) == "AGB_s__"] <- "AGB_ha"

#Check relationship between LiDAR and field AGB estimates
plot(AGB_ha~meanTCH, data=AGBplots,log="xy", xlab="Mean TCH (m)",
     ylab="Field AGB (Mg/ha)", pch=20)

## Build a LiDAR-AGB model
mod=lm(log(AGB_ha)~log(meanTCH), data=AGBplots)
summary(mod)

rse=summary(mod)$sigma

AGBpredicted=exp(0.5*rse^2+predict(mod))

plot(AGBplots$AGB_ha,AGBpredicted,pch=20,xlab="observed AGB", ylab="predicted AGB",
     cex.lab=1.3, xlim=range(c(AGBplots$AGB_ha,AGBpredicted)),
     ylim=range(c(AGBplots$AGB_ha,AGBpredicted)),asp=1)
abline(a=0,b=1)

RMSE= sqrt(sum((AGBpredicted-AGBplots$AGB_ha)^2)/length(AGBpredicted))
RMSErel=100*RMSE/mean(AGBplots$AGB_ha)

##Apply the model to the whole CHM
FunPredAGB=function(x) exp(mod$coefficients[1]+0.5*rse^2+mod$coefficients[2]*log(x))

Rast_TCHm=aggregate(CHMlidar,50,mean,na.rm=T)

AGBmap = lapp(Rast_TCHm, fun = FunPredAGB)

plot(AGBmap,plg=list(title='Predicted AGB\n(Mg/ha)', title.cex=0.9))

writeRaster(AGBmap,"Output/AGBmap_2017.tif",overwrite=TRUE)

sum(AGBmap$totCHM)
