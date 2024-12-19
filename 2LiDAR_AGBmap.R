################################################################################
#Calculate AGB of tree individuals (different stage) using BIOMASS package
################################################################################
# Load necessary libraries
library(sf)
library(terra)

## Linking to GEOS 3.10.1, GDAL 3.4.0, PROJ 8.2.0; sf_use_s2() is TRUE
##Load data
#First load the LiDAR canopy heigh model (CHM):
CHMlidar2017 = rast("E:/ALS_processing2017/02_processed/CHM_comb/CHM2017.tif")
#CHMlidar2017 = rast("E:/ALS_processing2017/02_processed/CHM_comb/CHM_pitfree2017.tif")
CHMlidar2017[CHMlidar2017<0]=0
plot(CHMlidar2017)

CHMlidar2024 = rast("E:/ALS_processing2024/02_processed/CHM_comb/totCHM.tif")
#CHMlidar2024 = rast("E:/ALS_processing2024/02_processed/CHM_comb/CHM_pitfree.tif")
CHMlidar2024[CHMlidar2024<0]=0
plot(CHMlidar2024)

#par(mfrow = c(1, 2))

#Crop CHM
CHMlidar2024 = crop(CHMlidar2024,CHMlidar2017)
CHMlidar2024 = mask(CHMlidar2024,CHMlidar2017)
plot(CHMlidar2024)

#Second load the field AGB estimates
AGBplots_0.25ha=st_read("Output/AGBplots_0.25ha.shp")
AGBplots_1ha=st_read("Output/AGBplots_1ha.shp")

################################################################################
# Define the function for extract CHM metrics
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
AGBplots_0.25ha <- CHM_metrics(AGBplots_0.25ha, CHMlidar2017)
AGBplots_1ha <- CHM_metrics(AGBplots_1ha, CHMlidar2017)
################################################################################
#Check relationship between LiDAR and field AGB estimates
#png("Output/AGBfield_LiDARmetrics_lm.png", width = 4000, height = 3000, res = 400)
par(mfrow = c(2, 2))

plot(AGB_ha ~ meanTCH, data = AGBplots_0.25ha, #log = "xy", 
     xlab = "Mean TCH (m)", ylab = "Field AGB (Mg/ha)", pch = 20, 
     col = as.factor(Stage), main = "0.25 ha")
legend("bottomright", legend = levels(as.factor(AGBplots_0.25ha$Stage)), 
       col = 1:length(unique(AGBplots_0.25ha$Stage)), pch = 20, title = "Stage")

plot(AGB_ha ~ meanTCH, data = AGBplots_1ha, #log = "xy", 
     xlab = "Mean TCH (m)", ylab = "Field AGB (Mg/ha)", pch = 20, 
     col = as.factor(Stage), main = "1 ha")
legend("bottomright", legend = levels(as.factor(AGBplots_1ha$Stage)), 
       col = 1:length(unique(AGBplots_1ha$Stage)), pch = 20, title = "Stage")

plot(AGB_ha ~ H50, data = AGBplots_0.25ha, #log = "xy", 
     xlab = "H50 (m)", ylab = "Field AGB (Mg/ha)", pch = 20, 
     col = as.factor(Stage), main = "0.25 ha")
legend("bottomright", legend = levels(as.factor(AGBplots_0.25ha$Stage)), 
       col = 1:length(unique(AGBplots_0.25ha$Stage)), pch = 20, title = "Stage")

plot(AGB_ha ~ H50, data = AGBplots_1ha, #log = "xy", 
     xlab = "H50 (m)", ylab = "Field AGB (Mg/ha)", pch = 20, 
     col = as.factor(Stage), main = "1 ha")
legend("bottomright", legend = levels(as.factor(AGBplots_1ha$Stage)), 
       col = 1:length(unique(AGBplots_1ha$Stage)), pch = 20, title = "Stage")

#dev.off()

################################################################################
# Function to predict AGB
Predict_AGB <- function(AGBplots, metric) {
  
  ## Build a LiDAR-AGB model
  mod=lm(log(AGB_ha)~log(metric), data=AGBplots)
  summary(mod)

  rse=summary(mod)$sigma

  AGBpredicted=exp(0.5*rse^2+predict(mod))

  RMSE= sqrt(sum((AGBpredicted-AGBplots$AGB_ha)^2)/length(AGBpredicted))
  RMSErel=100*RMSE/mean(AGBplots$AGB_ha)
  
  AGBplots$AGBpredicted <- AGBpredicted
  
  # Print RMSE and relative RMSE
  print(paste("res:", rse))
  print(paste("RMSE:", RMSE))
  print(paste("Relative RMSE:", RMSErel))

  return(list(AGBplots = AGBplots, rse = rse, RMSE = RMSE, RMSErel = RMSErel))
}

#Use function
AGBpred_meanTCH_0.25ha = Predict_AGB(AGBplots_0.25ha, AGBplots_0.25ha$meanTCH)
AGBpred_meanTCH_1ha = Predict_AGB(AGBplots_1ha, AGBplots_1ha$meanTCH)
AGBpred_H50_0.25ha = Predict_AGB(AGBplots_0.25ha, AGBplots_0.25ha$H50)
AGBpred_H50_1ha = Predict_AGB(AGBplots_1ha, AGBplots_1ha$H50)

# Function to plot AGB predicted
AGB_plot <- function(AGBplots_dat, title) {
  
  data = AGBplots_dat$AGBplots

  # Add text to the plot
  text_label <- paste(
    "RSE:", round(AGBplots_dat$rse, 2), "\n",
    "RMSE:", round(AGBplots_dat$RMSE, 2), "\n",
    "Relative RMSE:", round(AGBplots_dat$RMSErel, 2), "%")
  
  plot(data$AGB_ha,data$AGBpredicted,pch=20,xlab="observed AGB", ylab="predicted AGB",
       cex.lab=1.3, xlim=range(c(data$AGB_ha,data$AGBpredicted)),
       ylim=range(c(data$AGB_ha,data$AGBpredicted)),asp=1, 
       main = title)
  text(x = min(data$AGB_ha, na.rm = TRUE), y = max(data$AGBpredicted, na.rm = TRUE), 
       labels = text_label, pos = 1, # Position text to the right of the coordinates
       cex = 0.9 # Text size
  )
  abline(a=0,b=1)
}

#png("Output/AGBfield_AGBpred.png", width = 4000, height = 3000, res = 400)
par(mfrow = c(2, 2))
##Use function
AGB_plot(AGBpred_meanTCH_0.25ha, "0.25 ha (mean TCH)")
AGB_plot(AGBpred_meanTCH_1ha, "1 ha (mean TCH)")
AGB_plot(AGBpred_H50_0.25ha, "0.25 ha (H50)")
AGB_plot(AGBpred_H50_1ha, "1 ha (H50)")
#dev.off()

################################################################################
## Build a LiDAR-AGB model
mod=lm(log(AGB_ha)~log(meanTCH), data=AGBplots_0.25ha)
#mod=lm(log(AGB_ha)~log(meanTCH), data=AGBplots_1ha)
summary(mod)

rse=summary(mod)$sigma

##Apply the model to the whole CHM
FunPredAGB=function(x) exp(mod$coefficients[1]+0.5*rse^2+mod$coefficients[2]*log(x))

Rast_TCHm=aggregate(CHMlidar2017,50,mean,na.rm=T)
Rast_TCHm_2024=aggregate(CHMlidar2024,50,mean,na.rm=T)

AGBmap = lapp(Rast_TCHm, fun = FunPredAGB)
AGBmap_2024 = lapp(Rast_TCHm_2024, fun = FunPredAGB)

summary(AGBmap$totCHM)
summary(AGBmap_2024$totCHM)

par(mfrow = c(1, 2))
plot(AGBmap,plg=list(title='Predicted AGB\n(Mg/ha)', title.cex=0.9))
plot(AGBmap_2024,plg=list(title='Predicted AGB\n(Mg/ha)', title.cex=0.9))

#AGBmap_resampled <- resample(AGBmap, AGBmap_2024, method = "bilinear")
#plot(AGBmap_resampled)

AGBmapdif = AGBmap_2024-AGBmap
plot(AGBmapdif)

hist(AGBmapdif$totCHM)
summary(AGBmapdif$totCHM)

writeRaster(AGBmap,"Output/res0.25ha/AGBmap2017_0.25ha.tif",overwrite=TRUE)
writeRaster(AGBmap_2024,"Output/res0.25ha/AGBmap2024_0.25ha.tif",overwrite=TRUE)
writeRaster(AGBmapdif,"Output/res0.25ha/AGBmap_dif_0.25ha.tif",overwrite=TRUE)

writeRaster(AGBmap,"Output/res1ha/AGBmap2017_1ha.tif",overwrite=TRUE)
writeRaster(AGBmap_2024,"Output/res1ha/AGBmap2024_1ha.tif",overwrite=TRUE)
writeRaster(AGBmapdif,"Output/res1ha/AGBmap_dif_1ha.tif",overwrite=TRUE)
################################################################################
AGBmap = rast("Output/res0.25ha/AGBmap2017_0.25ha.tif")

Mean <- global(AGBmap, fun = "mean", na.rm = TRUE)[1, 1]

png("Output/AGB_predit_2017.png", width = 3000, height = 1500, res = 300)
par(mfrow=c(1,2))
plot(AGBmap,plg=list(title='Predicted AGB\n(Mg/ha)', title.cex=0.9))
hist_data <- hist(AGBmap$totCHM, breaks = 20,# Number of bins
  main = "Histogram of AGB predict in 2017", 
  xlab = "Predicted AGB (Mg/ha)")
abline(v = Mean, col = "red", lwd = 2, lty = 2) 
text(x = Mean+400, y = max(hist_data$counts) + 20,  # Position the text above the tallest bar
     labels = paste("Mean:", round(Mean, 2)), col = "red", cex = 0.9)
dev.off()


AGBmap_2024 = rast("Output/res0.25ha/AGBmap2024_0.25ha.tif")

Mean <- global(AGBmap_2024, fun = "mean", na.rm = TRUE)[1, 1]

png("Output/AGB_predit_2024.png", width = 3000, height = 1500, res = 300)
par(mfrow=c(1,2))
plot(AGBmap_2024,plg=list(title='Predicted AGB\n(Mg/ha)', title.cex=0.9))
hist_data <- hist(AGBmap_2024$totCHM, breaks = 20,# Number of bins
                  main = "Histogram of AGB predict in 2024", 
                  xlab = "Predicted AGB (Mg/ha)")
abline(v = Mean, col = "red", lwd = 2, lty = 2) 
text(x = Mean+400, y = max(hist_data$counts) + 20,  # Position the text above the tallest bar
     labels = paste("Mean:", round(Mean, 2)), col = "red", cex = 0.9)
dev.off()


AGBmapdif = rast("Output/res0.25ha/AGBmap_dif_0.25ha.tif")

Mean <- global(AGBmapdif, fun = "mean", na.rm = TRUE)[1, 1]

png("Output/AGB_predit_dif.png", width = 3000, height = 1500, res = 300)
par(mfrow=c(1,2))
plot(AGBmapdif,plg=list(title='Predicted AGB\n(Mg/ha)', title.cex=0.9))
hist_data <- hist(AGBmapdif$totCHM, breaks = 20,# Number of bins
                  main = "Histogram of AGB predict difference", 
                  xlab = "Predicted AGB (Mg/ha)")
abline(v = Mean, col = "red", lwd = 2, lty = 2) 
text(x = Mean+80, y = max(hist_data$counts) + 20,  # Position the text above the tallest bar
     labels = paste("Mean:", round(Mean, 2)), col = "red", cex = 0.9)
dev.off()

################################################################################