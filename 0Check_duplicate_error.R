
# Load field data
# Read data
Census_allplot=read.csv("Data/All_Census_KhaoYai_clean.csv",na.strings = c("NA","#VALUE!"))
str(Census_allplot)

# Keep only individuals > 5 cm dbh and remove NA 
Census_allplot_sub <- Census_allplot[Census_allplot$DBH_Cen2017 >= 5 & !is.na(Census_allplot$DBH_Cen2017),]
summary(Census_allplot_sub$DBH_Cen2017)

#Census_allplot_sub <- Census_allplot[Census_allplot$DBH_Cen2022 >= 5 & !is.na(Census_allplot$DBH_Cen2022),]
#summary(Census_allplot_sub$DBH_Cen2022)

Census_allplot_sub$tag_stem = paste(Census_allplot_sub$Tag, Census_allplot_sub$StemTag, sep="_")
Census_allplot_sub$tag_stem_Q = paste(Census_allplot_sub$tag_stem, Census_allplot_sub$Quadrat, sep="_")

# Find rows that are completely duplicated
duplicates_same <- Census_allplot_sub[duplicated(Census_allplot_sub), ]
all_duplicates_same <- Census_allplot_sub[Census_allplot_sub$tag_stem_Q %in% duplicates_same$tag_stem_Q, ]
#remove duplicate
Census_allplot_sub <- Census_allplot_sub[!duplicated(Census_allplot_sub), ]

# Find rows that are individual duplicated
duplicates <- Census_allplot_sub[duplicated(Census_allplot_sub$tag_stem_Q), ]

duplicates_notsame = Census_allplot_sub[Census_allplot_sub$tag_stem_Q %in% Census_allplot_sub$tag_stem_Q[duplicated(Census_allplot_sub$tag_stem_Q)], ]

Census_allplot_rm_dup <- Census_allplot_sub %>% 
  filter(!tag_stem_Q %in% duplicates$tag_stem_Q)

ordered_Tag <- duplicates_notsame %>%
  arrange(Tag)

write.csv(ordered_Tag,"Data/duplicates_check_5cm.csv",row.names = F)
#write.csv(ordered_Tag,"Data/duplicates_check_5cm_2022.csv",row.names = F)
################################################################################
#manual remove duplicate
ordered_Tag_rm_dup = read.csv("Data/duplicates_check_5cm_rm_dup.csv")
ordered_Tag_rm_dup = read.csv("Data/duplicates_check_5cm_2022_rm_dup.csv")

Census_allplot_clean = rbind(Census_allplot_rm_dup,ordered_Tag_rm_dup)

write.csv(Census_allplot_clean,"Data/All_Census_KhaoYai_clean_5cm_2022.csv",row.names = F)
