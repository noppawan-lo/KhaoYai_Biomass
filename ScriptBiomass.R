################################################################################
#Calculate AGB of tree individuals (different stage) using BIOMASS package
################################################################################
# Load necessary libraries
#install_github("umr-amap/BIOMASS")
library(devtools)
library(BIOMASS)

# Set the working directory to the folder containing the data file
setwd("D:/KhaoYai_Biomass/")

# Read data
Census_allplot <- read.csv("Data/All_Census2017_2022_KhaoYai_SC.CSV")
str(Census_allplot)

# New column with stage of each plot and change MST --> OGS stage)
Census_allplot$Stage = substr(Census_allplot$Plot, start = 1, stop = 3)
Census_allplot$Stage[Census_allplot$Plot=="MST30ha"] <- "OGS"

# Keep only individuals > 5 cm dbh and remove NA 
Census_allplot_sub <- Census_allplot[Census_allplot$DBH_C2017 >= 5 & !is.na(Census_allplot$DBH_C2017),]
summary(Census_allplot_sub$DBH_C2017)

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
Census_allplot_sub$Hpred <- mapply(Calculate_Hpred, Census_allplot_sub$Stage, Census_allplot_sub$DBH_C2017)
summary(Census_allplot_sub$Hpred)


# Estimate AGB 
Census_allplot_sub$AGB=computeAGB(D = Census_allplot_sub$DBH_C2017, 
                                  WD = Census_allplot_sub$meanWD, 
                                  H = Census_allplot_sub$Hpred)

#AGB in each quadrat (20x20 m = 0.04 ha)
Res=summaryByPlot(Census_allplot_sub$AGB,plot=Census_allplot_sub$Quadrat)
Res$AGBperha=Res$AGB/0.04


#### AGBmonteCarlo ####
resultMC <- AGBmonteCarlo(
  D = Census_allplot_sub$DBH_C2017, WD = Census_allplot_sub$meanWD, 
  errWD = Census_allplot_sub$sdWD, H=Census_allplot_sub$Hpred,errH=0
)
str(resultMC)

################################################################################