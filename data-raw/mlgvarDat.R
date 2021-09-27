# Simulated multi-level network data

library(modnets)

set.seed(666)
mlgvarDat <- mlGVARsim(nTime = 500, nPerson = 100, nNode = 5, m = TRUE,
                       m2 = .2, m1 = .5, modType = 'partial')

atts <- mlgvarDat[setdiff(names(mlgvarDat), 'data')]
mlgvarDat <- mlgvarDat$data
attributes(mlgvarDat)[names(atts)] <- atts

write.csv(mlgvarDat, 'data-raw/mlgvarDat.csv', row.names = FALSE)
usethis::use_data(mlgvarDat, overwrite = TRUE)
