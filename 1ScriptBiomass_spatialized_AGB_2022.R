################################################################################
#Calculate AGB of tree individuals (different stage) using BIOMASS package
################################################################################
# Load necessary libraries
#remotes::install_github("umr-amap/BIOMASS" , ref="subplot_summary")
#library(devtools)
# Read data
Census_allplot_sub=read.csv("Data/All_Census_KhaoYai_clean_5cm_2022.csv",na.strings = c("NA","#VALUE!"))

str(Census_allplot_sub)

# New column with stage of each plot and change MST --> OGS stage)
Census_allplot_sub$Stage = substr(Census_allplot_sub$Plot, start = 1, stop = 3)
Census_allplot_sub$Stage[Census_allplot_sub$Plot=="MST"] <- "OGS"

summary(Census_allplot_sub$DBH_Cen2022)

# Check typo in taxonomy
genspcorr = correctTaxo(Census_allplot_sub$GenusSpecies)#,useCache=TRUE)
filt=genspcorr$nameModified=="TRUE"
unique(cbind(Census_allplot_sub$GenusSpecies[filt],genspcorr[filt,]))
Census_allplot_sub$GenusCorr=genspcorr$genusCorrected
Census_allplot_sub$SpeciesCorr=genspcorr$speciesCorrected

# Get wood density from the global wood density database
#Census_allplot_sub$Plot="MST30ha"
genspcorr = getWoodDensity(Census_allplot_sub$GenusCorr,Census_allplot_sub$SpeciesCorr,stand = Census_allplot_sub$Plot)
Census_allplot_sub$meanWD=genspcorr$meanWD
Census_allplot_sub$sdWD=genspcorr$sdWD

# Assign tree height
# Define the function to calculate Hpred in different stage
Calculate_Hpred <- function(Stage, dbh) {
  if (Stage == "OGS") {
    return(exp(0.221 + 0.5 * 0.324^2 + 1.115 * log(dbh) - 0.075 * log(dbh)^2))
  } else if (Stage == "SES") {
    return(exp(0.919 + 0.5 * 0.286^2 + 0.727 * log(dbh) - 0.0189 * log(dbh)^2))
  } else if (Stage == "SIS") {
    return(exp(1.034 + 0.5 * 0.493^2 + 0.564 * log(dbh) - 0.0323 * log(dbh)^2))
  } else {
    return(NA)
  }
}
#Use function 
Census_allplot_sub$Hpred <- mapply(Calculate_Hpred, Census_allplot_sub$Stage, Census_allplot_sub$DBH_Cen2022)
summary(Census_allplot_sub$Hpred)


# Estimate AGB 
Census_allplot_sub$AGB=computeAGB(D = Census_allplot_sub$DBH_Cen2022, 
                                  WD = Census_allplot_sub$meanWD, 
                                  H = Census_allplot_sub$Hpred)

# Save data in order to load these data in a new script :
save(Census_allplot_sub , file = "Output/Field2022/AGB_calculate_2022.RData")
################################################################################
### Spatialized AGB 
################################################################################
# Load necessary libraries
library(ggplot2)

# Load data
load("Output/AGB_calculate_2022.RData")
coord <- read.csv(file = "Data/Corner_allplots.csv")

str(Census_allplot_sub)

# Split the data by plots in order to apply the check_plot_coord() function
split_coord <- split(coord , f = coord$Plot)

#par(mfrow=c(2,2))
res_check <- lapply(split_coord , function(dat) { #dat = split_coord[[1]] for an example
  print(unique(dat$Plot))
  check_plot <- check_plot_coord(
    proj_coord = dat[c("X_UTM","Y_UTM")],
    rel_coord = dat[c("X_loc","Y_loc")], 
    trust_GPS_corners = FALSE, 
    tree_df = Census_allplot_sub[Census_allplot_sub$Plot==unique(dat$Plot),], 
    tree_coords = c("PX","PY")
  )
})

# There's a warning message telling us that some trees are not inside the plot -> it concerns the plot MST !
################################################################################
# Adding projected coordinates to the tree data.frame
# for one plot, we got the projected coordinates in check_plot$tree_proj_coord.
# So for all plots, we have to loop over res_check, extract $tree_proj_coord, and then rbind the list
tree_projected_coord <- do.call(rbind,lapply(res_check , function(x) x$tree_proj_coord ))
#Census_allplot_sub[c("X_UTM","Y_UTM")] <- tree_projected_coord

tree_projected_coord$Plot <- rownames(tree_projected_coord)
tree_projected_coord$Plot = substr(tree_projected_coord$Plot, start = 1, stop = 4)
tree_projected_coord$Plot[tree_projected_coord$Plot=="MST."]<- "MST"

plot_names = unique(tree_projected_coord$Plot)

all_trees <- data.frame()

# Loop through each plot name

for (i in 1:length(plot_names)) {
  
  pos <- tree_projected_coord[tree_projected_coord$Plot == plot_names[i], ]
  tree <- Census_allplot_sub[Census_allplot_sub$Plot == plot_names[i], ]
  
  tree[c("X_UTM", "Y_UTM")] <- pos
  
  all_trees <- rbind(all_trees, tree)
}
Census_allplot_sub=all_trees
################################################################################
##### Dividing plots into subplots #####

# You can apply the function to all the plot but the grid dimensions will be the same for all the plots :
DividePlot = function (coord,grid_size,dat){
  
  plot_division <- divide_plot(
    rel_coord = coord[c("X_loc","Y_loc")],
    proj_coord = coord[c("X_UTM","Y_UTM")], 
    grid_size = grid_size, 
    corner_plot_ID = coord$Plot,
    grid_tol = 0.5, # =0.1 by defaut, so we increase it to 0.5 : if more than 50% of the plot area is not include in the grid area, it returns an error
    tree_df = dat,
    tree_coords = c("PX","PY"),
    tree_plot_ID = dat$Plot,
  )
}
#Use function for 50x50 m (all plots)
plot_division_0.25ha = DividePlot(coord, grid_size = c(50,50), Census_allplot_sub)

# The first 16 warning messages tell us that the dimensions of the grid does not fit perfectly the dimensions of the plots, but its OK
# The 17th warning message tell us that some trees are not assigned to a subplot :
# It happens (for example) with the 80mx60m plots when a tree as an X or a Y coordinate > 50

# The function returns a list containing 2 elements : 
corner_subplot_0.25ha <- plot_division_0.25ha$sub_corner_coord
tree_subplot_0.25ha <- plot_division_0.25ha$tree_df

##Cut SIS4 because not 50*50 m
corner_subplot_0.25ha <- corner_subplot_0.25ha[!corner_subplot_0.25ha$subplot_id %in% c("SIS4_0_0","SIS4_0_1","SIS4_1_0","SIS4_1_1"),]
tree_subplot_0.25ha <- tree_subplot_0.25ha[!tree_subplot_0.25ha$subplot_id %in% c("SIS4_0_0","SIS4_0_1","SIS4_1_0","SIS4_1_1"),]

##### Summarising tree metrics #####
#source("~/BIOMASS/R/subplot_summary.R")

subplots_0.25ha <- subplot_summary(subplots = plot_division_0.25ha, value = "AGB")

# By defaut, the summary is the sum :
AGB0.25ha = subplots_0.25ha$polygon
AGB0.25ha$Stage = substr(AGB0.25ha$subplot_id, start = 1, stop = 3)
AGB0.25ha$Stage[AGB0.25ha$Stage=="MST"] <- "OGS"

##Cut SIS4 because not 50*50 m
AGB0.25ha <- AGB0.25ha[!AGB0.25ha$subplot_id %in% c("SIS4_0_0","SIS4_0_1","SIS4_1_0","SIS4_1_1"),]
##################################################################################
##Use function for SIS4 plot (9 subplots = 0.25 ha)

plot_division_SIS4 = DividePlot(coord[coord$Plot=="SIS4",], 
                                grid_size = c(46.66, 46.66),
                                Census_allplot_sub[Census_allplot_sub$Plot=="SIS4",])

##### Summarising tree metrics #####
subplots_SIS4 <- subplot_summary(subplots = plot_division_SIS4, value = "AGB")

## The function returns a list containing 2 elements : 
corner_subplot_0.25ha = rbind(corner_subplot_0.25ha, plot_division_SIS4$sub_corner_coord)
tree_subplot_0.25ha <- rbind(tree_subplot_0.25ha, plot_division_SIS4$tree_df)

## By defaut, the summary is the sum :
AGB0.25ha_SIS4 = subplots_SIS4$polygon
AGB0.25ha_SIS4$Stage = substr(AGB0.25ha_SIS4$subplot_id, start = 1, stop = 3)

AGB0.25ha = rbind(AGB0.25ha, AGB0.25ha_SIS4)
##################################################################################

##Use function for all plot (1 ha)
# Function to find max of X_loc and Y_loc for each plot
get_max_values <- function(split_coord) {
  max_values <- lapply(split_coord, function(df) {
    X_loc <- max(df$X_loc)
    Y_loc <- max(df$Y_loc)
    # Limit X_loc and Y_loc to a maximum value of 100
    X_loc <- pmin(X_loc, 100)
    Y_loc <- pmin(Y_loc, 100)
    return(c(X_loc = X_loc, Y_loc = Y_loc))
  }
  )
}
plot_names <- as.data.frame(get_max_values(split_coord))
plot_names <- t(plot_names)
plot_names <- as.data.frame(plot_names)
plot_names$Plot <- rownames(plot_names)


corner_subplot_list <- list()
tree_subplot_list <- list()
subplots_list <- list()

# Loop over each plot name
for (i in 1:length(plot_names$Plot)) {
  # Filter the data for the current plot
  coord1ha <- coord[coord$Plot == plot_names$Plot[i], ]
  Census_allplot_sub_1ha <- Census_allplot_sub[Census_allplot_sub$Plot == plot_names$Plot[i], ]
  
  plot_division <- DividePlot(coord1ha, 
                              grid_size=c(plot_names$X_loc[i],plot_names$Y_loc[i]), 
                              Census_allplot_sub_1ha)
  
  corner_subplot_list[[i]] <- plot_division$sub_corner_coord
  tree_subplot_list[[i]] <- plot_division$tree_df
  
  subplots_list[[i]] <- subplot_summary(subplots = plot_division, value = "AGB")
} 

# Combine results into data frames
corner_subplot_1ha <- do.call(rbind, corner_subplot_list)
tree_subplot_1ha <- do.call(rbind, tree_subplot_list)
# By defaut, the summary is the sum :
AGB1ha <- do.call(rbind, lapply(subplots_list, function(x) x$polygon))
AGB1ha$Stage = substr(AGB1ha$subplot_id, start = 1, stop = 3)
AGB1ha$Stage[AGB1ha$Stage=="MST"] <- "OGS"

##View AGB
subplots_list[[10]]
subplots_list[[12]]
################################################################################
## Save polygon as shape file for analyse CHM metrics 
library(sf)

st_crs(AGB0.25ha) <- 32647
st_crs(AGB1ha) <- 32647

st_write(AGB0.25ha,"Output/Field2022/AGBplots_0.25ha.shp", delete_layer = TRUE)
st_write(AGB1ha,"Output/Field2022/AGBplots_1ha.shp", delete_layer = TRUE)

library(terra)
writeRaster(AGB0.25ha$AGB_summary_per_ha,"Output/Field2022/AGBdif_0.25ha.tif",overwrite=TRUE)
writeRaster(AGB1ha$AGB_summary_per_ha,"Output/Field2022/AGBdif_0.25ha.tif",overwrite=TRUE)
################################################################################
################################################################################

AGBplots_0.25ha_2022=st_read("Output/Field2022/AGBplots_0.25ha.shp")
AGBplots_1ha_2022=st_read("Output/Field2022/AGBplots_1ha.shp")

AGBplots_0.25ha_2017=st_read("Output/AGBplots_0.25ha.shp")
AGBplots_1ha_2017=st_read("Output/AGBplots_1ha.shp")

#Raster from polygon
template_raster_0.25ha <- rast(ext(AGBplots_0.25ha_2017), resolution = 1) 
plot(template_raster_0.25ha)

subplots_0.25ha <- rasterize(AGBplots_0.25ha_2017, template_raster_0.25ha, field = "sbplt_d")
subplots_0.25ha$AGB_2017 <- rasterize(AGBplots_0.25ha_2017, template_raster_0.25ha, field = "AGB_s__")
subplots_0.25ha$AGB_2022 <- rasterize(AGBplots_0.25ha_2022, template_raster_0.25ha, field = "AGB_s__")

subplots_0.25ha

template_raster_1ha <- rast(ext(AGBplots_1ha_2017), resolution = 1) 
plot(template_raster_1ha)

subplots_1ha <- rasterize(AGBplots_1ha_2017, template_raster_1ha, field = "sbplt_d")
subplots_1ha$AGB_2017 <- rasterize(AGBplots_1ha_2017, template_raster_1ha, field = "AGB_s__")





subplots_0.25ha

AGBplots_0.25ha = AGBplots_0.25ha_2017 %>% 
  select(sbplt_d,AGB_s__) %>% 
  st_drop_geometry() %>%
  left_join(AGBplots_0.25ha_2022, by="sbplt_d")%>%
  mutate(AGB_dif=AGB_s__.y-AGB_s__.x)

AGBplots_1ha = AGBplots_1ha_2017 %>% 
  select(sbplt_d,AGB_s__) %>% 
  st_drop_geometry() %>%
  left_join(AGBplots_1ha_2022, by="sbplt_d")%>%
  mutate(AGB_dif=AGB_s__.y-AGB_s__.x)

#Raster from polygon
template_raster_0.25ha <- rast(ext(AGBplots_0.25ha_2022), resolution = 1) 
plot(template_raster_0.25ha)


#Raster from polygon
template_raster <- rast(ext(AGBplots_0.25ha), resolution = 1) 
plot(template_raster)

subplots_0.25ha <- rasterize(AGBplots_0.25ha, template_raster, field = "sbplt_d")
subplots_0.25ha$meanTCH_2017 <- rasterize(AGBplots_0.25ha, template_raster, field = "meanTCH_2017")
subplots_0.25ha$meanTCH_2024 <- rasterize(AGBplots_0.25ha, template_raster, field = "meanTCH_2024")

plot(subplots_0.25ha)
library(terra)
library(sf)

# Use the original sf object with geometry
template_raster <- rast(ext(AGBplots_0.25ha), resolution = 1)

# Check the template raster
print(template_raster)


writeRaster(subplots_0.25ha$AGBpredict_dif,"Output/AGBdif_0.25ha.tif",overwrite=TRUE)

