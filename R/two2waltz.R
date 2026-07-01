#' Simulate associated virtual species (Linear)
#' @description
#' Simulates two virtual species whose log-intensity (occurrence probability)
#' responds *linearly* to an environmental predictor, in contrast to
#' [two2tango()] which uses a Gaussian niche curve. As in [two2tango()], the
#' two species are made statistically associated by adding correlated random
#' noise to their log-intensity surfaces before sampling point locations from
#' an inhomogeneous Poisson process, allowing you to generate paired species
#' distributions with a known, tunable degree of spatial co-occurrence.
#'
#' @param alpha1 Numeric. Intercept of the log-linear response for species 1
#'   (log-intensity when the predictor equals 0).
#' @param alpha2 Numeric. Intercept of the log-linear response for species 2.
#' @param beta1 Numeric. Slope of the log-linear response for species 1.
#'   Positive values mean species 1's intensity increases with the predictor;
#'   negative values mean it decreases; 0 means no environmental effect.
#' @param beta2 Numeric. Slope of the log-linear response for species 2.
#' @param var.e Numeric. Variance of the random noise (`e`) added to each
#'   species' log-intensity surface. Shared by both species so the noise
#'   covariance matrix stays symmetric.
#' @param cov.e Numeric. Covariance between the noise terms for species 1 and
#'   species 2. Controls the strength and sign of spatial association between
#'   the two species independent of any similarity in their environmental
#'   responses: positive values push the species toward co-occurring more
#'   than expected from their environmental responses alone, negative values
#'   push them apart.
#' @param predictor A `SpatRaster` (from \pkg{terra}) or a pixel image
#'   (`im` object from \pkg{spatstat.geom}) giving the environmental predictor
#'   surface. If a `SpatRaster` is supplied it is converted internally via
#'   [as.im.SpatRaster()].
#'
#' @returns A list of length 3:
#' \describe{
#'   \item{sp1}{An `sf` `POINT` object with the simulated occurrence locations
#'   of species 1.}
#'   \item{sp2}{An `sf` `POINT` object with the simulated occurrence locations
#'   of species 2.}
#'   \item{cor.pp}{An `htest` object (the result of [stats::cor.test()])
#'   giving the correlation between the two species' log-intensity surfaces at
#'   the level of the underlying Poisson process — i.e. the "true" spatial
#'   association driving the simulation, as opposed to any correlation you
#'   might later measure from the sampled points themselves. Access the
#'   correlation coefficient with `cor.pp$estimate`.}
#' }
#'
#' @seealso [two2waltz()] for a Gaussian (rather than linear) environmental response.
#'
#' @examples
#' set.seed(134)
#' size=10
#' env <- spatstat.geom::im(matrix(0, size, size),
#'           xrange=c(0,1), yrange=c(0,1))
#' xy <- expand.grid(x = env$xcol, y = env$yrow)
#' temp.dummy <- gstat::gstat(formula=z~1, locations=~x+y, dummy=TRUE, beta=0,
#'                     model=gstat::vgm(psill=0.1, range=50, model='Exp'), nmax=20)
#' temp <- predict(temp.dummy, newdata=xy, nsim=1, debug.level=0)
#' temp$sim1 <- scale(temp$sim1)
#'
#' example <- two2waltz(alpha1=1, alpha2=1,
#'                      beta1=0, beta2=0,
#'                      var.e = 1, cov.e=-0.9,
#'                      predictor = temp)
#'
#' sp1 <- example[[1]]
#' sp2 <- example[[2]]
#' cor.sp1.sp2 <- example[[3]]
#' env <- spatstat.geom::as.im(temp)
#'
#' plot(env,
#'      main= paste0('cov=',
#'                   round(cor.sp1.sp2$estimate,2)))
#' plot(sp1, col='white', add=TRUE)
#' plot(sp2, col='black', add=TRUE)
#'
#' @export
#'
two2waltz <- function(alpha1=1, alpha2=2,
                      beta1=0, beta2=0,
                      var.e=1, cov.e=0,
                      predictor) {

  crs_predictor <- terra::crs(terra::rast(predictor))

  # convert raster to pixel image
  if(!spatstat.geom::is.im(predictor)){
    predictor <- spatstat.geom::as.im(predictor)
  }

  # number of grid cells in the rasters
  # N.cell <- terra::ncell(predictor)
  N.cell <- prod(dim(predictor))

  # generate the correlated e-plurals
  sigma.e <- matrix(c(var.e, cov.e, cov.e, var.e), ncol=2)
  e <- mvtnorm::rmvnorm(n = N.cell, mean = c(0, 0), sigma = sigma.e)

  # True point pattern
  # point process intensity lambda as a function of environment
  lambda.sp1 <- exp(alpha1 + beta1*predictor + e[,1])
  lambda.sp2 <- exp(alpha2 + beta2*predictor + e[,2])

  # the correlation at the level of the Poisson process
  cor.pp <- stats::cor.test(log(lambda.sp1[]), log(lambda.sp2[]))

  # area of each pixel
  # pixel.area <- 1/(dim(predictor)[1]*dim(predictor)[2])
  pixel.area <- predictor$xstep * predictor$ystep

  # sample points using inhomogeneous poisson point process
  points.sp1 <- spatstat.random::rpoispp(lambda.sp1/pixel.area)
  points.sp2 <- spatstat.random::rpoispp(lambda.sp2/pixel.area)

  # generate counts
  abund.sp1 <- stats::rpois(n = N.cell, lambda.sp1[])
  abund.sp2 <- stats::rpois(n = N.cell, lambda.sp2[])

  # if(return=='coords'){
  #   sp1 <- sf::st_coordinates(points.sp1) %>%
  #     dplyr::as_tibble() %>%
  #     dplyr::rename(X1=X, Y1=Y)
  #   sp2 <- sf::st_coordinates(points.sp2) %>%
  #     dplyr::as_tibble() %>%
  #     dplyr::rename(X2=X, Y2=Y)
  # } else {
    sp1 <- sf::st_as_sf(points.sp1) %>%
      dplyr::filter(.data$label == 'point') %>%
      dplyr::select(-label) %>%
      sf::st_set_crs(crs_predictor)
    sp2 <- sf::st_as_sf(points.sp2) %>%
      dplyr::filter(.data$label == 'point') %>%
      dplyr::select(-label) %>%
      sf::st_set_crs(crs_predictor)
  # }

  # can also return the value of the predictor for each data point
  # env1 <- terra::extract(rast(predictor), vect(PTS.sp1)) %>% select(env1=2)
  # env2 <- terra::extract(rast(predictor), vect(PTS.sp2)) %>% select(env2=2)
  #return(list(sp1, sp2, env1, env2))

  return(list(sp1, sp2, cor.pp))
}


