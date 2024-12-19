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

Tree_1ha_MST = Tree_1ha[Tree_1ha$Plot=="MST",]
Tree_1ha_other = Tree_1ha[Tree_1ha$Plot!="MST",]

Tree_1ha_MST$Tree_Tag <- Tree_1ha_MST$Tag
Tree_1ha_other$Tree_Tag <- sapply(strsplit(Tree_1ha_other$Tag, "_"), function(x) paste(x[-1], collapse = "_"))

Tree_1ha = rbind(Tree_1ha_MST,Tree_1ha_other)
################################################################################

Plotnames = data.frame(subplot_id = unique(Tree_1ha$subplot_id))
Plotnames$Stage = substr(Plotnames$subplot_id, start = 1, stop = 3)
Plotnames$Stage[Plotnames$Stage=="MST"] <- "OGS"

Plotnames_OGS_SES = Plotnames[Plotnames$Stage!="SIS",]
Plotnames_SIS = Plotnames[Plotnames$Stage=="SIS",]

library(tidyverse)
library(writexl)
library(openxlsx)

# Loop through all unique subplot IDs
for (plot_id in Plotnames_OGS_SES$subplot_id) {
  
  data <- Tree_1ha %>% filter(subplot_id == plot_id) %>%
    mutate(DBH_2024 = NA, Status = NA, Note = NA)
  
  dataNA <- data %>% filter(is.na(DBH_Cen2022))
  
  data10 <- data[!data$DBH_Cen2022<9,]
  data10 = data10[!is.na(data10$subplot_id),]
  
  data10 = rbind(data10,dataNA)
  
  wb <- createWorkbook()
  
  # Add sheets and write data to them
  addWorksheet(wb, "9cm")
  writeData(wb, "9cm", data10)
  
  addWorksheet(wb, "1cm")
  writeData(wb, "1cm", data)

  # Save the workbook
  saveWorkbook(wb, file = paste0("Recensus_Sheets/", plot_id, ".xlsx"), overwrite = TRUE)
}

# Loop through all unique subplot IDs
for (plot_id in Plotnames_SIS$subplot_id) {
  data <- Tree_1ha %>% filter(subplot_id == plot_id) %>%
    mutate(DBH_2024 = NA, Status = NA, Note = NA)
  
  dataNA <- data %>% filter(is.na(DBH_Cen2022))
  
  data5 <- data[!data$DBH_Cen2022<7,]
  data5 = data5[!is.na(data5$subplot_id),]
  
  wb <- createWorkbook()
  
  # Add sheets and write data to them
  addWorksheet(wb, "7cm")
  writeData(wb, "7cm", data5)
  
  addWorksheet(wb, "1cm")
  writeData(wb, "1cm", data)
  
  # Save the workbook
  saveWorkbook(wb, file = paste0("Recensus_Sheets/", plot_id, ".xlsx"), overwrite = TRUE)
}
################################################################################