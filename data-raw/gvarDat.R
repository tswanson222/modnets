# Simulated GGM data from simNet

library(modnets)

set.seed(666)
gvarDat <- simNet(N = 5000, p = 5, m = TRUE, m2 = .15, m1 = .5,
                  modType = 'partial', lags = 1)

atts <- gvarDat[setdiff(names(gvarDat), 'data')]
gvarDat <- gvarDat$data
attributes(gvarDat)[names(atts)] <- atts

write.csv(gvarDat, 'data-raw/gvarDat.csv', row.names = FALSE)
usethis::use_data(gvarDat, overwrite = TRUE)
