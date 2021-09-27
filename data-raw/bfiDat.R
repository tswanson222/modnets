# BFI data

library(psychTools)
library(usethis)

bfiDat <- na.omit(psychTools::bfi)[, 1:26]

# Reverse-scoring variables
rv <- c("A1", "C4", "C5", "E1", "E2", "O2", "O5")
bfiDat[, rv] <- (bfiDat[, rv] * (-1)) + (max(bfiDat) + 1)

# Create composite scores and re-code gender
bfiDat <- structure(data.frame(
  sapply(strsplit('ACENO', '')[[1]], function(z){
    rowMeans(bfiDat[, grep(z, colnames(bfiDat))])
  }),
  gender = bfiDat$gender - 1
), row.names = 1:nrow(bfiDat))

write.csv(bfiDat, 'data-raw/bfiData.csv')
usethis::use_data(bfiDat, overwrite = TRUE)
