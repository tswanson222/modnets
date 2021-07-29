##### show which variables were selected in each model
#' Shows which variables were selected for each node of the network
#'
#' @param object fitNetwork object
#' @param threshold numeric or logical
#' @param mod not sure
#'
#' @return A table
#' @export
#'
#' @examples
#' 1 + 1
selected <- function(object, threshold = FALSE, mod = 1){
  ints <- NULL
  if(threshold != FALSE & !is.numeric(threshold)){threshold <- .05}
  if(isTRUE(attr(object, "mlGVAR"))){object <- object[[mod + 1]]}
  if(isTRUE(attr(object, "SURnet"))){
    if('SURnet' %in% names(object)){object <- object$SURnet}
    mods0 <- object$temporal$coefs[[ifelse(threshold == FALSE, 1, 2)]]
    mods <- lapply(seq_len(nrow(mods0)), function(z){
      if(threshold == FALSE){
        z <- colnames(mods0)[mods0[z, ] != 0]
      } else {
        z <- colnames(mods0)[mods0[z, ] <= threshold]
      }
      z <- replace(z, length(z) == 0, "")
      z <- replace(z, is.null(z), "")
      return(z[z != "(Intercept)"])
    })
    names(mods) <- rownames(mods0)
  } else {
    if(threshold != FALSE & "fitobj" %in% names(object)){
      mods <- lapply(object$fitobj, function(z){
        z <- summary(z)$coefficients[, 4]
        z <- ifelse(z <= threshold, 1, 0)
        z <- names(z)[which(z == 1)]
        z <- replace(z, length(z) == 0, "")
        z <- replace(z, is.null(z), "")
        return(z[z != "(Intercept)"])
      })
    } else {
      mods <- lapply(object$mods, function(z){
        z <- rownames(z$model)[-1]
        return(replace(z, length(z) == 0, ""))
      })
    }
  }
  if(any(grepl(":", mods))){
    ints <- lapply(mods, function(z){
      z <- z[grepl(":", z)]
      return(replace(z, length(z) == 0, ""))
    })
    ix <- max(sapply(ints, length))
    ints <- do.call(cbind.data.frame, lapply(ints, function(z){
      if(length(z) < ix){z <- c(z, rep("", ix - length(z)))}
      return(z)
    }))
    mods <- lapply(mods, function(z){
      z <- z[!grepl(":", z)]
      return(replace(z, length(z) == 0, ""))
    })
  }
  mx <- max(sapply(mods, length))
  mods <- do.call(cbind.data.frame, lapply(mods, function(z){
    if(length(z) < mx){z <- c(z, rep("", mx - length(z)))}
    return(z)
  }))
  out <- mods
  if(!is.null(ints)){out <- list(mods = mods, ints = ints)}
  return(out)
}
