library(dplyr)


expression <- readRDS("data/expression/merged.dgecounts.rds")
BL6 <- read.csv("data/allelic/merged.BL6_fractional_UMIs.csv")
CAST <- read.csv(file = "data/allelic/merged.CAST_fractional_UMIs.txt", sep = "\t")

ex <- as.matrix(expression$umicount$inex$all)


df_bl6 <- as.data.frame(BL6)
df_CAST <- as.data.frame(CAST)
df_ex <- as.data.frame(ex)

#sort allelic columns to match expression df
excol <- colnames(ex)
df_ex <- df_ex[,which(excol %in% colnames(df_bl6))]
reordered_BL6 <- df_bl6[,match(colnames(df_ex), colnames(df_bl6))]
df_ex <- df_ex[which(rownames(df_ex) %in% rownames(df_bl6)),]


reordered_CAST <- df_CAST[,match(colnames(df_ex), colnames(df_CAST))]


#fix_allelic_umi <- function(ex, allelic_umi){
#  if((ex > 0  | is.na(ex)) & (allelic_umi == 0 | is.na(allelic_umi))){
#    allelic_umi <- NA
#  }
#  else if(ex == 0){
#    allelic_umi <- 0
#  }
#  return(allelic_umi)
#}


fix_allelic_umi <- function(ex, allelic_umi_1, allelic_umi_2 ){
  allelic_umi <- allelic_umi_1
  #catch (0,NA) condition
  v <- c(allelic_umi_1, allelic_umi_2)
  v[which(is.na(v))] <- 0
  if((sum(v) == 0)){
    if((ex > 0  | is.na(ex))){
      allelic_umi <- NA
    }
  } else if(ex == 0){
    allelic_umi <- 0
  }
    return(allelic_umi)
}

fixed.reordered_BL6 <- reordered_BL6
fixed.reordered_CAST <- reordered_CAST

for(col in colnames(reordered_BL6)){
  fixed.reordered_BL6[,col] <- mapply(fix_allelic_umi, df_ex[,col], reordered_BL6[,col], reordered_CAST[,col])
}
write.csv(fixed.reordered_BL6, file = "data/allelic/merged.fixed.BL6_fractional_UMIs.csv")


for(col in colnames(reordered_CAST)){
  fixed.reordered_CAST[,col] <- mapply(fix_allelic_umi, df_ex[,col], reordered_CAST[,col], reordered_BL6[,col])
}
write.csv(fixed.reordered_CAST, file = "data/allelic/merged.fixed.CAST_fractional_UMIs.csv")



## Creating log files
diff_df <- is.na(reordered_CAST) != is.na(fixed.reordered_CAST)
CAST_original_difference <- as.character(reordered_CAST[diff_df])
CAST_changed_difference <- as.character(fixed.reordered_CAST[diff_df])
diff_df_indices <- which(is.na(reordered_CAST) != is.na(fixed.reordered_CAST), arr.ind = TRUE)
CAST.log <- paste(names(diff_df_indices[,1]), "at column index", diff_df_indices[,2], ":", CAST_original_difference, "->", CAST_changed_difference)
header <- paste("## Total number of changes:", length(diff_df_indices[,1]))
cat("## CAST", header, CAST.log, file = "data/allelic/CAST_changed_log.txt", sep = "\n")

diff_df <- is.na(reordered_BL6) != is.na(fixed.reordered_BL6)
original_difference <- as.character(reordered_BL6[diff_df])
changed_difference <- as.character(fixed.reordered_BL6[diff_df])
diff_df_indices <- which(is.na(reordered_BL6) != is.na(fixed.reordered_BL6), arr.ind = TRUE)
BL6.log <- paste(names(diff_df_indices[,1]), "at column index", diff_df_indices[,2], ":", original_difference, "->", changed_difference)
header <- paste("## Total number of changes:", length(diff_df_indices[,1]))
cat("## C57_BL6", header, BL6.log, file = "data/allelic/BL6_changed_log.txt", sep = "\n")

