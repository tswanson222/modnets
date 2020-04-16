### Initial look at types of moderation
mt <- rep(c('none', 'full', 'partial', 'full2', 'partial2'), each = 100)
fits <- pbapply::pblapply(seq_along(mt), function(i){
  simNet(100, 5, T, .1, onlyNets = TRUE, m1 = .5, modType = mt[i])
}, cl = 4)



fits_out <- data.frame(do.call(rbind, lapply(fits, function(z){
  mtz <- attr(z, 'modType')
  `%b1%` <- ifelse(grepl('full', mtz), function(x, y){x == y}, function(x, y){x != y})
  `%m1%` <- ifelse(mtz == 'full2', `%b1%`, function(x, y){x != y})
  FUN1 <- switch(2 - (mtz == 'none'), function(x){sum(x)/2}, match.fun('all'))
  FUN2 <- match.fun(ifelse(grepl('2', mtz), 'all', 'sum'))
  z1 <- FUN1(z$b1[z$b2 != 0] %b1% 0)
  z2 <- FUN2(z$m1[apply(z$b2, 2, function(k) any(k != 0))] %m1% 0)
  z3 <- c(b1 = z1, m1 = z2, type = mtz)
  return(z3)
})))


# Multi-plotting
multiplot <- function(grobs, y = NULL, box = FALSE, legend = NULL, 
                      nrow = 2, group = TRUE, getLegend = FALSE, exclude = NULL){
  if(is(grobs, c('data.frame', 'matrix'))){
    stopifnot(!is.null(y))
    setup <- function(dat, y, box = FALSE, group = TRUE, legend = NULL){
      y <- toupper(match.arg(tolower(y), c('ecr', 'saam', 'other')))
      if(y == 'OTHER'){y <- 'mood|SE$|stress'}
      y <- colnames(dat)[grep(y, colnames(dat))]
      yt <- setNames(switch(
        length(y) - 1, c('ECR Anxious', 'ECR Avoidant'), 
        c('SAAM Anxiety', 'SAAM Avoidance', 'SAAM Secure'), 
        c('Positive Mood', 'Negative Mood', 'Self-esteem', 'Stress')), y)
      if(!is.null(exclude)){
        yt <- yt[setdiff(names(yt), exclude)]
        y <- setdiff(y, exclude)
      }
      p <- lapply(y, function(z){
        pdat(dat, y = z, ylab = FALSE, title = unname(yt[z]), group = group,
             legend = ifelse(group & !identical(legend, FALSE), 'bottom', FALSE), 
             box = box, plot = FALSE)
      })
      g_legend <- function(a.gplot){
        tmp <- ggplot_gtable(ggplot_build(a.gplot))
        leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
        legend <- tmp$grobs[[leg]]
        return(legend)
      }
      out <- list(p = p)
      if(group & !identical(legend, FALSE)){
        out$p <- lapply(p, function(z) z + theme(legend.position = 'none'))
        out$legend <- g_legend(p[[1]])
      }
      return(out)
    }
    leg <- switch(2 - is.null(legend), NULL, FALSE)
    grobs <- setup(grobs, y, box, group, leg)
    if(getLegend){return(grobs$legend)}
  }
  if(!is.null(names(grobs))){
    if('legend' %in% names(grobs)){legend <- grobs$legend}
    grobs <- grobs[[1]]
  }
  if(is.null(legend)){
    gridExtra::grid.arrange(gridExtra::arrangeGrob(grobs = grobs, nrow = nrow))
  } else {
    gridExtra::grid.arrange(gridExtra::arrangeGrob(grobs = grobs, nrow = nrow), 
                            legend, nrow = 2, heights = c(10, 1))
  }
}


