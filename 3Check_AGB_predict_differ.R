
# Load necessary libraries
library(sf)
library(terra)

CHMlidar2017 = rast("E:/ALS_processing2017/02_processed/CHM_comb/CHM2017.tif")
CHMlidar2017[CHMlidar2017<0]=0
plot(CHMlidar2017)

CHMlidar2024 = rast("E:/ALS_processing2024/02_processed/CHM_comb/totCHM.tif")
CHMlidar2024[CHMlidar2024<0]=0
plot(CHMlidar2024)

par(mfrow=c(1,2))
#Crop CHM
CHMlidar2024 = crop(CHMlidar2024,CHMlidar2017)
CHMlidar2024 = mask(CHMlidar2024,CHMlidar2017)
plot(CHMlidar2024)

#Second load the field AGB estimates
AGBplots_0.25ha=st_read("Output/AGBplots_0.25ha.shp")
AGBplots_1ha=st_read("Output/AGBplots_1ha.shp")

AGBplots_0.25ha = AGBplots_1ha
colnames(AGBplots_0.25ha)[colnames(AGBplots_0.25ha) == "AGB_s__"] <- "AGB_ha"

##Build a LiDAR and field AGB dataset
#First extract LiDAR data for each subplots:
meanTCH_2017 <- terra::extract(CHMlidar2017, AGBplots_0.25ha, fun = mean, na.rm = TRUE)
meanTCH_2024 <- terra::extract(CHMlidar2024, AGBplots_0.25ha, fun = mean, na.rm = TRUE)

AGBplots_0.25ha$meanTCH_2017 = meanTCH_2017$totCHM
AGBplots_0.25ha$meanTCH_2024 = meanTCH_2024$totCHM
################################################################################
mod=lm(log(AGB_ha)~log(meanTCH_2017), data=AGBplots_0.25ha)
summary(mod)

rse=summary(mod)$sigma

#AGBpredicted=exp(0.5*rse^2+predict(mod))
#AGBplots_0.25ha$AGBpredicted_2017 = AGBpredicted

#Raster from polygon
template_raster <- rast(ext(AGBplots_0.25ha), resolution = 1) 
plot(template_raster)

subplots_0.25ha <- rasterize(AGBplots_0.25ha, template_raster, field = "sbplt_d")
subplots_0.25ha$meanTCH_2017 <- rasterize(AGBplots_0.25ha, template_raster, field = "meanTCH_2017")
subplots_0.25ha$meanTCH_2024 <- rasterize(AGBplots_0.25ha, template_raster, field = "meanTCH_2024")

plot(subplots_0.25ha)

FunPredAGB=function(x) exp(mod$coefficients[1]+0.5*rse^2+mod$coefficients[2]*log(x))

subplots_0.25ha$AGBpredicted2017 = lapp(subplots_0.25ha$meanTCH_2017, fun = FunPredAGB)
subplots_0.25ha$AGBpredicted2024 = lapp(subplots_0.25ha$meanTCH_2024, fun = FunPredAGB)
plot(subplots_0.25ha)

plot(subplots_0.25ha$AGBpredicted2017,plg=list(title='Predicted AGB\n(Mg/ha)', title.cex=0.9))
plot(subplots_0.25ha$AGBpredicted2024,plg=list(title='Predicted AGB\n(Mg/ha)', title.cex=0.9))
################################################################################

subplots_0.25ha$AGBpredict_dif = subplots_0.25ha$AGBpredicted2024 - subplots_0.25ha$AGBpredicted2017
plot(subplots_0.25ha$AGBpredict_dif,plg=list(title='Predicted AGB \n difference (Mg/ha)', title.cex=0.9))

writeRaster(subplots_0.25ha$AGBpredict_dif,"Output/AGBdif_0.25ha.tif",overwrite=TRUE)

df <- as.data.frame(subplots_0.25ha, xy = TRUE)

df_unique <- df[!duplicated(df$sbplt_d), ]
summary(df_unique$AGBpredict_dif)

#hist(df_unique$AGBpredict_dif,main = "Predicted AGB difference (Mg/ha)",nclass=20)
#abline(v=quantile(df_unique$AGBpredict_dif, probs=c(0.25,0.5,0.75)),col="red")
#plot(df_unique$AGBpredict_dif,text= df_unique$sbplt_d)
#abline(h=quantile(df_unique$AGBpredict_dif, probs=c(0.25,0.5,0.75)),col="red")

################################################################################
ggplot(df_unique,aes(x= "", y = AGBpredict_dif)) +
  geom_boxplot() + 
  geom_jitter(position = position_jitter(seed = 1)) + 
  ggrepel::geom_text_repel(aes(label = sbplt_d),
            position = position_jitter(seed = 1), size = 2) +
  coord_flip()


Add_Q <- function(data, value) {
  if (value < summary(data$AGBpredict_dif)["1st Qu."]) {
    return("Q1")
  } else if (value >= summary(data$AGBpredict_dif)["1st Qu."] & 
             value < summary(data$AGBpredict_dif)["Median"]) {
    return("Q2")
  } else if (value >= summary(data$AGBpredict_dif)["Median"] & 
             value < summary(data$AGBpredict_dif)["3rd Qu."]) {
    return("Q3")
  } else if (value >= summary(data$AGBpredict_dif)["3rd Qu."]) {
    return("Q4")
  } else {
    return(NA)
  }
}

# Apply the function using sapply (applies the function to each row's AGBpredict_dif)
df_unique$Quartile <- sapply(df_unique$AGBpredict_dif, Add_Q, data = df_unique)

write.csv(df_unique, "Output/AGB_dif1ha_Quartile.csv")

#library(tidyverse)
#df_unique_Q = df_unique %>% group_by(Quartile) %>% sample_n(4)
################################################################################

df_unique = read.csv("Output/AGB_dif1ha_Quartile.csv")

library(ggplot2)

png("Output/AGBpreditdif_plot1ha.png", width = 3000, height = 1500, res = 300)
ggplot(df_unique,aes(x= "", y = AGBpredict_dif)) +
  geom_boxplot() + 
  geom_jitter(position = position_jitter(seed = 1)) + 
  ggrepel::geom_text_repel(aes(label = sbplt_d),
                           position = position_jitter(seed = 1), size = 2) +
  coord_flip()

dev.off()
