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

################################################################################
# Define the function
CHM_metrics <- function(AGBplots, CHMlidar) {
  colnames(AGBplots)[colnames(AGBplots) == "AGB_s__"] <- "AGB_ha"
  
  ##Build a LiDAR and field AGB dataset
  #First extract LiDAR data for each subplots:
  meanTCH <- terra::extract(CHMlidar, AGBplots, fun = mean, na.rm = TRUE)
  H50 <- terra::extract(CHMlidar, AGBplots, 
                        fun = function(x) quantile(x, probs = 0.5, na.rm = TRUE))
  
  # Assign the extracted values back to the data
  AGBplots$meanTCH <- meanTCH[, 2]  # Second column has the extracted values
  AGBplots$H50 <- H50[, 2]
  
  return(AGBplots)  # Return the updated data
}
#Use function
AGBplots_0.25ha <- CHM_metrics(AGBplots_0.25ha, CHMlidar)
AGBplots_1ha <- CHM_metrics(AGBplots_1ha, CHMlidar)

#Check relationship between LiDAR and field AGB estimates
png("Output/AGBfield_LiDARmetrics.png", width = 4000, height = 3000, res = 400)
par(mfrow = c(2, 2))

plot(AGB_ha ~ meanTCH, data = AGBplots_0.25ha, log = "xy", 
     xlab = "Mean TCH (m)", ylab = "Field AGB (Mg/ha)", pch = 20, 
     col = as.factor(Stage), main = "0.25 ha")
legend("bottomright", legend = levels(as.factor(AGBplots_0.25ha$Stage)), 
       col = 1:length(unique(AGBplots_0.25ha$Stage)), pch = 20, title = "Stage")

plot(AGB_ha ~ meanTCH, data = AGBplots_1ha, log = "xy", 
     xlab = "Mean TCH (m)", ylab = "Field AGB (Mg/ha)", pch = 20, 
     col = as.factor(Stage), main = "1 ha")
legend("bottomright", legend = levels(as.factor(AGBplots_1ha$Stage)), 
       col = 1:length(unique(AGBplots_1ha$Stage)), pch = 20, title = "Stage")

plot(AGB_ha ~ H50, data = AGBplots_0.25ha, log = "xy", 
     xlab = "H50 (m)", ylab = "Field AGB (Mg/ha)", pch = 20, 
     col = as.factor(Stage), main = "0.25 ha")
legend("bottomright", legend = levels(as.factor(AGBplots_0.25ha$Stage)), 
       col = 1:length(unique(AGBplots_0.25ha$Stage)), pch = 20, title = "Stage")

plot(AGB_ha ~ H50, data = AGBplots_1ha, log = "xy", 
     xlab = "H50 (m)", ylab = "Field AGB (Mg/ha)", pch = 20, 
     col = as.factor(Stage), main = "1 ha")
legend("bottomright", legend = levels(as.factor(AGBplots_1ha$Stage)), 
       col = 1:length(unique(AGBplots_1ha$Stage)), pch = 20, title = "Stage")

dev.off()

Predict_AGB <- function(AGBplots, metric) {
  
## Build a LiDAR-AGB model
mod=lm(log(AGB_ha)~log(metric), data=AGBplots)
summary(mod)

rse=summary(mod)$sigma

AGBpredicted=exp(0.5*rse^2+predict(mod))

RMSE= sqrt(sum((AGBpredicted-AGBplots$AGB_ha)^2)/length(AGBpredicted))
RMSErel=100*RMSE/mean(AGBplots$AGB_ha)

# Print RMSE and relative RMSE
print(paste("RSE:", rse))
print(paste("RMSE:", RMSE))
print(paste("Relative RMSE:", RMSErel))

return(AGBpredicted)

}

#Use function
AGBplots_0.25ha$AGBpredicted <- Predict_AGB(AGBplots_0.25ha, AGBplots_0.25ha$meanTCH)
AGBplots_1ha$AGBpredicted <- Predict_AGB(AGBplots_1ha, AGBplots_1ha$meanTCH)

plot(AGBplots_0.25ha$AGB_ha,AGBplots_0.25ha$AGBpredicted,pch=20,xlab="observed AGB", ylab="predicted AGB",
     cex.lab=1.3, xlim=range(c(AGBplots_0.25ha$AGB_ha,AGBplots_0.25ha$AGBpredicted)),
     ylim=range(c(AGBplots_0.25ha$AGB_ha,AGBplots_0.25ha$AGBpredicted)),asp=1)
abline(a=0,b=1)

plot(AGBplots_1ha$AGB_ha,AGBplots_1ha$AGBpredicted,pch=20,xlab="observed AGB", ylab="predicted AGB",
     cex.lab=1.3, xlim=range(c(AGBplots_1ha$AGB_ha,AGBplots_1ha$AGBpredicted)),
     ylim=range(c(AGBplots_1ha$AGB_ha,AGBplots_1ha$AGBpredicted)),asp=1)
abline(a=0,b=1)

################################################################################
## Build a LiDAR-AGB model
mod=lm(log(AGB_ha)~log(meanTCH), data=AGBplots_0.25ha)
summary(mod)

rse=summary(mod)$sigma

##Apply the model to the whole CHM
FunPredAGB=function(x) exp(mod$coefficients[1]+0.5*rse^2+mod$coefficients[2]*log(x))

Rast_TCHm=aggregate(CHMlidar,50,mean,na.rm=T)

AGBmap = lapp(Rast_TCHm, fun = FunPredAGB)

plot(AGBmap,plg=list(title='Predicted AGB\n(Mg/ha)', title.cex=0.9))

writeRaster(AGBmap,"Output/AGBmap_2017.tif",overwrite=TRUE)
################################################################################