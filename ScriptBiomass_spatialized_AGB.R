################################################################################
#Calculate AGB of tree individuals (different stage) using BIOMASS package
################################################################################
# Load necessary libraries
#remotes::install_github("umr-amap/BIOMASS" , ref="subplot_summary")
#library(devtools)
library(BIOMASS)

# Set the working directory to the folder containing the data file
#setwd("D:/KhaoYai_Biomass/")

# Read data
Census_allplot <- read.csv("Data/All_Census2017-2022_KhaoYai.csv")
str(Census_allplot)

# New column with stage of each plot and change MST --> OGS stage)
Census_allplot$Stage = substr(Census_allplot$Plot, start = 1, stop = 3)
Census_allplot$Plot[Census_allplot$Plot=="MST30ha"] <- "MST"
Census_allplot$Stage[Census_allplot$Plot=="MST"] <- "OGS"

# Keep only individuals > 5 cm dbh and remove NA 
Census_allplot_sub <- Census_allplot[Census_allplot$DBH_Cen2017 >= 5 & !is.na(Census_allplot$DBH_Cen2017),]
summary(Census_allplot_sub$DBH_Cen2017)

# Check typo in taxonomy
genspcorr = correctTaxo(Census_allplot_sub$GenusSpecies,useCache=TRUE)
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
Census_allplot_sub$Hpred <- mapply(Calculate_Hpred, Census_allplot_sub$Stage, Census_allplot_sub$DBH_Cen2017)
summary(Census_allplot_sub$Hpred)


# Estimate AGB 
Census_allplot_sub$AGB=computeAGB(D = Census_allplot_sub$DBH_Cen2017, 
                                  WD = Census_allplot_sub$meanWD, 
                                  H = Census_allplot_sub$Hpred)

#AGB in each quadrat (20x20 m = 0.04 ha)
#Res_Quadrat=summaryByPlot(Census_allplot_sub$AGB, plot=Census_allplot_sub$Quadrat , drawPlot = T)
#Res_Quadrat$AGBperha=Res_Quadrat$AGB/0.04

#### AGBmonteCarlo ####
resultMC <- AGBmonteCarlo(
  D = Census_allplot_sub$DBH_Cen2017, WD = Census_allplot_sub$meanWD, 
  errWD = Census_allplot_sub$sdWD, H=Census_allplot_sub$Hpred,errH=0
)
str(resultMC)

# Save data in order to load these data in a new script :
save(Census_allplot_sub , file = "Output/AGB_calculate.RData")

################################################################################
### Spatialized AGB 
################################################################################
# Load necessary libraries
library(ggplot2)

# Load data
load("Output/AGB_calculate.RData")
coord <- read.csv(file = "Data/Corner_allplots.csv")

# Excluding SIS4 and SIS5 because there are no trees associated in Census_allplot_sub
coord <- coord[ ! coord$Plot %in% c("SIS4","SIS5"),]

# Split the data by plots in order to apply the check_plot_coord() function
split_coord <- split(coord , f = coord$Plot)

##### Ckecking plot coordinates #####
# OGS1 example 
OGS1_coord <- split_coord[["OGS1"]]
range(OGS1_coord$X_loc)
range(OGS1_coord$Y_loc)

OGS1_trees <- Census_allplot_sub[Census_allplot_sub$Plot=="OGS1",]
range(OGS1_trees$PX2)
range(OGS1_trees$PY2)
### !!! There has been an inversion of X and Y axis (either on the corner coordinates, or on the tree data frame)

# Let's suppose that the inversion has been made on the trees (but not for MST plot) :
Census_allplot_sub[Census_allplot_sub$Plot!="MST", c("PX2","PY2")] <- Census_allplot_sub[Census_allplot_sub$Plot!="MST", c("PY2","PX2")]

res_check <- lapply(split_coord , function(dat) { #dat = split_coord[[1]] for an example
  print(unique(dat$Plot))
  check_plot <- check_plot_coord(
    proj_coord = dat[c("X_UTM","Y_UTM")],
    rel_coord = dat[c("X_loc","Y_loc")], 
    trust_GPS_corners = FALSE, 
    tree_df = Census_allplot_sub[Census_allplot_sub$Plot==unique(dat$Plot),], 
    tree_coords = c("PX2","PY2")
  )
})

# There's a warning message telling us that some trees are not inside the plot -> it concerns the plot MST !

# Adding projected coordinates to the tree data.frame
# for one plot, we got the projected coordinates in check_plot$tree_proj_coord.
# So for all plots, we have to loop over res_check, extract $tree_proj_coord, and then rbind the list
tree_projected_coord <- do.call(rbind,lapply(res_check , function(x) x$tree_proj_coord ))
Census_allplot_sub[c("X_UTM","Y_UTM")] <- tree_projected_coord

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
    tree_coords = c("PX2","PY2"),
    tree_plot_ID = dat$Plot,
    )
}
#Use function for 50x50 m (all plots)
plot_division_0.25ha = DividePlot(coord, grid_size = 50, Census_allplot_sub)

# The first 16 warning messages tell us that the dimensions of the grid does not fit perfectly the dimensions of the plots, but its OK
# The 17th warning message tell us that some trees are not assigned to a subplot :
# It happens (for example) with the 80mx60m plots when a tree as an X or a Y coordinate > 50

# The function returns a list containing 2 elements : 
corner_subplot_0.25ha <- plot_division_0.25ha$sub_corner_coord
tree_subplot_0.25ha <- plot_division_0.25ha$tree_df
##################################################################################
##Use function for SIS4 plot (9 subplots = 0.25 ha)
#plot_division_SIS4 = DividePlot(coord[coord$Plot=="SIS4",], 
#                                grid_size = 46.66, 
#                                Census_allplot_sub[Census_allplot_sub$Plot=="SIS4",])
## The function returns a list containing 2 elements : 
#corner_subplot_SIS4 <- plot_division_SIS4$sub_corner_coord
#tree_subplot_SIS4 <- plot_division_SIS4$tree_df

## The function returns a list containing 2 elements : 
#corner_subplot_0.25ha = rbind(corner_subplot_0.25ha, corner_subplot_SIS4)
#tree_subplot_0.25ha <- rbind(tree_subplot_0.25ha, tree_subplot_SIS4)
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
} 
  
# Combine results into data frames
corner_subplot_df <- do.call(rbind, corner_subplot_list)
tree_subplot_df <- do.call(rbind, tree_subplot_list)

################################################################################
##### Summarising tree metrics #####

#source("~/BIOMASS/R/subplot_summary.R")

subplots <- subplot_summary(subplots = plot_division_0.25ha, value = "AGB")

# By defaut, the summary is the sum :
AGB0.25ha = subplots$tree_summary

# But you can supply any function you want (even your own function) : 
# subplots <- subplot_summary(subplots = plot_division ,
#                             value = "AGB",
#                             fun = median , probs = 0.9)

# Adding tree representation (we have to loop on plot's name) : 
treeplot = lapply(names(subplots$plot_design) , function(plot_name) { #plot_name = "OGS1" #for an example
  plot_representation <- subplots$plot_design[[plot_name]]
  plot_representation + 
    geom_point(data = Census_allplot_sub[Census_allplot_sub$Plot==plot_name,] ,
               mapping = aes(x = X_UTM , y = Y_UTM) , 
               shape = 3) # you can represent what you want, as in the vignette
})

#############################################################################################################
# If you want to center the subplots for the 80mx60m plots : 
plot_division2 <- divide_plot(
  rel_coord = coord[c("X_loc","Y_loc")],
  proj_coord = coord[c("X_UTM","Y_UTM")], 
  grid_size = 50, 
  corner_plot_ID = coord$Plot,
  grid_tol = 0.5,
  tree_df = Census_allplot_sub,
  tree_coords = c("PX2","PY2"),
  tree_plot_ID = Census_allplot_sub$Plot,
  centred_grid = TRUE # add this parameter = TRUE 
)

subplots2 <- subplot_summary(subplots2 = plot_division2, value = "AGB")

# Adding tree representation (we have to loop on plot's name) : 
treeplot2 = lapply(names(subplots2$plot_design) , function(plot_name) { #plot_name = "OGS1" #for an example
  plot_representation <- subplots2$plot_design[[plot_name]]
  plot_representation + 
    geom_point(data = Census_allplot_sub[Census_allplot_sub$Plot==plot_name,] ,
               mapping = aes(x = X_UTM , y = Y_UTM) , 
               shape = 3) # you can represent what you want, as in the vignette
})
#############################################################################################################