#' Convert SpatRaster to Pixel Image
#' @description
#' Converts a SpatRaster to a pixel image
#' @param X a SpatRaster
#' @param ... additional arguments (unused)
#' @returns A pixel image
#' @rdname as.im.SpatRaster
#'
#' @importFrom spatstat.geom as.im
#' @export
#' @method as.im SpatRaster
as.im.SpatRaster <- function(X, ...) {
  X <- X[[1]]
  g <- as.list(X, geom=TRUE)

  isfact <- is.factor(X)
  if (isfact) {
    v <- matrix(as.data.frame(X)[, 1], nrow=g$nrows, ncol=g$ncols, byrow=TRUE)
  } else {
    v <- as.matrix(X, wide=TRUE)
  }
  vtype <- if(isfact) "factor" else typeof(v)
  if(vtype == "double") vtype <- "real"
  tv <- v[g$nrows:1, ]
  if(isfact) tv <- factor(tv, levels=levels(X))
  out <- list(
    v = tv,
    dim = c(g$nrows, g$ncols),
    xrange = c(g$xmin, g$xmax),
    yrange = c(g$ymin, g$ymax),
    xstep = g$xres[1],
    ystep = g$yres[1],
    xcol = g$xmin + (1:g$ncols) * g$xres[1] + 0.5 * g$xres,
    yrow = g$ymax - (g$nrows:1) * g$yres[1] + 0.5 * g$yres,
    type = vtype,
    units  = list(singular=g$units, plural=g$units, multiplier=1)
  )
  attr(out$units, "class") <- "unitname"
  attr(out, "class") <- "im"
  out
}
