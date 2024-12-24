# Load field data
# Read data
Census_allplot=read.csv("Data/All_Census_KhaoYai_clean.csv",na.strings = c("NA","#VALUE!"))
str(Census_allplot)

# Keep only individuals > 5 cm dbh and remove NA 
Census_allplot_sub <- Census_allplot[Census_allplot$DBH_Cen2017 >= 1 & !is.na(Census_allplot$DBH_Cen2017),]
summary(Census_allplot_sub$DBH_Cen2017)

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
Census_allplot_sub <- Census_allplot_sub[!duplicated(Census_allplot_sub$tag_stem_Q), ]

A <- Census_allplot_sub[!duplicated(Census_allplot_sub$Tag), ]

unique(A$StemTag)

