###Prepare datasheet for recensus###
##OGS stage at 10 cm and secondary forest at 5 cm

library(BIOMASS)
library(tidyverse)

# Read data
coord <- read.csv(file = "Data/Corner_allplots.csv")
Census_allplot <- read.csv("Data/All_Census_KhaoYai.csv")
str(Census_allplot)

#Change chr to num
Census_allplot$PX <- as.numeric(Census_allplot$PX)
Census_allplot$PY <- as.numeric(Census_allplot$PY)

#DBH difference between 2017 and 2022
Census_allplot$DBH_dif = Census_allplot$DBH_Cen2022-Census_allplot$DBH_Cen2017
hist(Census_allplot$DBH_dif,  breaks = 100)
summary(Census_allplot$DBH_dif)

#Select plot =<  1 ha (recensus all plot)
Tree_second_1ha <- Census_allplot[Census_allplot$Plot != "MST" & Census_allplot$Plot != "SIS4", ]
Tree_second_1ha$subplot_id <- Tree_second_1ha$Plot

#Function for divine plot
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

#Use function for plot arae > 1 ha
plot_division_MST_SIS4 = DividePlot(coord[coord$Plot %in% c("MST", "SIS4"),], 
                                grid_size = c(100, 100),
                               Census_allplot[Census_allplot$Plot %in% c("MST", "SIS4"),])

#Datafame of each subplot
Tree_MST_SIS4_1ha = plot_division_MST_SIS4$tree_df
Tree_MST_SIS4_1ha$plot_id = NULL
Tree_MST_SIS4_1ha = Tree_MST_SIS4_1ha[!is.na(Tree_MST_SIS4_1ha$subplot_id),]

#Combine data
Tree_1ha = rbind(Tree_MST_SIS4_1ha,Tree_second_1ha)
################################################################################

Plotnames = unique(Tree_1ha$subplot_id)

library(writexl)
# Loop through all unique subplot IDs
for (plot_id in Plotnames) {
  data <- Tree_1ha %>% filter(subplot_id == plot_id) %>%
    mutate(DBH_2024 = NA, Status = NA, Note = NA)
  
  data5 <- data[!data$DBH_Cen2022<5,]
  data5 = data5[!is.na(data5$subplot_id),]
  
  data10 <- data[!data$DBH_Cen2022<10,]
  data10 = data10[!is.na(data10$subplot_id),]

  # Write each file
  write_xlsx(data, paste0("Recensus_Sheets/1 cm/", plot_id, ".xlsx"))
  write_xlsx(data5, paste0("Recensus_Sheets/5 cm/", plot_id, ".xlsx"))
  write_xlsx(data10, paste0("Recensus_Sheets/10 cm/", plot_id, ".xlsx"))
}
################################################################################
