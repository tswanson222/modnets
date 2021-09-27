# Simulated GGM data from simNet

library(modnets)

set.seed(666)
ggmDat <- simNet(N = 5000, p = 5, m = TRUE, m2 = .3, m1 = .5,
                 modType = 'partial', nCores = 16)

atts <- ggmDat[setdiff(names(ggmDat), 'data')]
ggmDat <- ggmDat$data
attributes(ggmDat)[names(atts)] <- atts

write.csv(ggmDat, 'data-raw/ggmDat.csv', row.names = FALSE)
usethis::use_data(ggmDat, overwrite = TRUE)
