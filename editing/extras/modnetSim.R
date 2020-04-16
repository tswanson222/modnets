tri <- function(x, d = 0, tt = 0){
  stopifnot(all(c(d, tt) %in% 0:1))
  dtri <- switch(d + 1, lower.tri(x), upper.tri(x))
  if(tt){x <- t(x)}
  return(x[dtri])
}


if(FALSE){
  x1 <- expand.grid(0:1, 0:1)
  x2 <- x1[rep(1:4, each = 4), ]
  x2 <- lapply(split(x2, rep(1:4, each = 4)), as.matrix)
  x <- list(poop, shit)
  names(x) <- paste0('p', 1:2)
  unname(sapply(x2, function(z){
    oo <- sapply(1:4, function(i){
      t1 <- tri(x$p1, z[i, 1], z[i, 2])
      t2 <- tri(x$p2, x1[i, 1], x1[i, 2])
      return(cor0(t1, t2))
    })
    return(which.max(abs(oo)))
  }))
}


################################################################################
################################################################################
out <- simNet(1000, 10, T)
dat <- out$data

fit0 <- fitNetwork(dat, 11)
v1 <- varSelect(dat, 11, 'AIC')
v2 <- varSelect(dat, 11, 'EBIC')
v3 <- varSelect(dat, 1:11, 'EBIC')
fit1 <- fitNetwork(dat, 11, v1)
fit2 <- fitNetwork(dat, 11, v2)
fit3 <- fitNetwork(dat, 1:11, v3)
fit4 <- mgm2(dat[, -11])
fit5 <- mgm2(dat, 11)
fit6 <- mgm2(dat, 1:11)
fits <- nlist('fit', T)
names(fits) <- c(paste0('fit', c(0, 'AIC', 'EBIC', 'ALL')), paste0('mgm', 0:2))


fits1 <- lapply(fits, function(z) net(z)[1:10, 1:10])
fits2 <- lapply(fits, function(z) net(z, T)[1:10, 1:10])
fits3 <- lapply(fits, function(z) netInts(z, r = 11))
fits4 <- lapply(fits, function(z) netInts(z, threshold = T, r = 11))
fits5 <- lapply(fits, function(z) netInts(z, avg = T, r = 11))
fits6 <- lapply(fits, function(z) netInts(z, threshold = T, avg = T, r = 11))
fits7 <- nlist('fits[1-6]', T)


fits8 <- setNames(lapply(seq_along(fits7), function(z){
  fun <- rep(c(net, netInts), c(2, 4))[[z]]
  out1 <- sapply(fits7[[z]], matrixDist, fun(out))
  out2 <- do.call(rbind, lapply(fits7[[z]][sapply(fits7[[z]], length) > 0], performance, fun(out)))
  return(list(cor = out1, perform = out2))
}), c('net', 'netThresh', 'ints', 'intsThresh', 'intsAVG', 'intsAVGthresh'))


res1 <- do.call(rbind, lapply(fits8, '[[', 1))
res1 <- structure(cbind.data.frame(matrix = rownames(res1), res1), row.names = 1:nrow(res1))
res1 <- data.frame(matrix = factor(rep(res1[, 'matrix'], ncol(res1) - 1)), 
                   model = factor(rep(colnames(res1)[-1], each = nrow(res1))),
                   value = unname(unlist(c(res1[, 2:ncol(res1)]))))

require(ggplot2)
p1 <- ggplot(na.omit(res1), aes(x = matrix, y = value, colour = model)) + 
  geom_line() + geom_point() + theme_bw()



res2 <- do.call(rbind, lapply(fits8, '[[', 2))
res2 <- cbind.data.frame(do.call(rbind, strsplit(rownames(res2), '[.]')), res2)
colnames(res2)[1:2] <- c('matrix', 'model')
res2$matrix <- factor(res2$matrix)
res2$model <- factor(res2$model)
rownames(res2) <- 1:nrow(res2)





