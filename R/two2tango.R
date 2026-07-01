#' Simulate associated virtual species (Gaussian)
#' @description
#' Simulates two virtual species whose intensity (occurrence probability)
#' responds to an environmental predictor following a Gaussian (bell-shaped)
#' niche curve, each species peaking at its own environmental optimum. The two
#' species' spatial patterns are made statistically associated by adding
#' correlated random noise to their log-intensity surfaces before sampling
#' point locations from an inhomogeneous Poisson process. This lets you
#' generate paired species distributions with a known, tunable degree of
#' spatial co-occurrence, useful for testing species-association or
#' co-occurrence methods against a known ground truth.
#'
#'
#' @param peak1 Numeric. Maximum intensity (height of the Gaussian curve) for
#'   species 1, reached when the predictor equals `mu1`.
#' @param peak2 Numeric. Maximum intensity (height of the Gaussian curve) for
#'   species 2, reached when the predictor equals `mu2`.
#' @param mu1 Numeric. Environmental optimum for species 1 — the predictor
#'   value at which species 1's intensity is highest.
#' @param mu2 Numeric. Environmental optimum for species 2 — the predictor
#'   value at which species 2's intensity is highest.
#' @param sigma1 Numeric. Niche breadth (standard deviation of the Gaussian
#'   response curve) for species 1. Larger values produce a flatter, wider
#'   response to the environment; smaller values produce a narrower, more
#'   specialized response.
#' @param sigma2 Numeric. Niche breadth (standard deviation of the Gaussian
#'   response curve) for species 2.
#' @param var.e Numeric. Variance of the random noise (`e`) added to each
#'   species' log-intensity surface. Shared by both species so the noise
#'   covariance matrix stays symmetric.
#' @param cov.e Numeric. Covariance between the noise terms for species 1 and
#'   species 2. Controls the strength and sign of spatial association between
#'   the two species independent of any similarity in their environmental
#'   responses: positive values push the species toward co-occurring more
#'   than expected from their niches alone, negative values push them apart.
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
#'@seealso [two2waltz()] for a linear (rather than Gaussian) environmental response.
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
#' example <- two2tango(peak1=1, peak2=1,
#'                      mu1=-1, sigma1=0.5,
#'                      mu2=1, sigma2=0.5,
#'                      cov=0, predictor = temp)
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
two2tango <- function(peak1=1, peak2=1,
                      mu1=0.5, mu2=0.5,
                      sigma1=0.5, sigma2=0.5,
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
  lambda.sp1 <- peak1*exp(-0.5*((predictor-mu1)^2)/(sigma1^2) + e[,1])
  lambda.sp2 <- peak2*exp(-0.5*((predictor-mu2)^2)/(sigma2^2) + e[,2])

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


