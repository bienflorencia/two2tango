two2tango <- function(peak1, peak2,
                      mu1, mu2, sigma1, sigma2,
                      var=1, cov,
                      predictor, return='sf') {
  
  crs_predictor <- terra::crs(rast(predictor))
  
  # convert raster to pixel image
  if(!is.im(predictor)){
    predictor <- as.im.SpatRaster(predictor)
  }
  
  # number of grid cells in the rasters
  N.cell <- ncell(predictor)
  
  # generate the correlated e-plurals
  sigma <- matrix(c(var, cov, cov, var), ncol=2)
  e <- mvtnorm::rmvnorm(n = N.cell, mean = c(0, 0), sigma = sigma)
  
  # True point pattern
  # point process intensity lambda as a function of environment
  lambda.sp1 <- peak1*exp(-0.5*((predictor-mu1)^2)/(sigma1^2) + e[,1])
  lambda.sp2 <- peak2*exp(-0.5*((predictor-mu2)^2)/(sigma2^2) + e[,2])
  
  # sample points using inhomogeneous poisson point process
  PTS.sp1 <- rpoispp(lambda.sp1)
  PTS.sp2 <- rpoispp(lambda.sp2)
  
  if(return=='coords'){
    sp1 <- st_coordinates(PTS.sp1) %>% as_tibble() %>% rename(X1=X, Y1=Y)
    sp2 <- st_coordinates(PTS.sp2) %>% as_tibble() %>% rename(X2=X, Y2=Y)
  } else {
    sp1 <- st_as_sf(PTS.sp1) %>%
      filter(label=='point') %>% select(-label) %>%
      st_set_crs(crs_predictor)
    sp2 <- st_as_sf(PTS.sp2) %>%
      filter(label=='point') %>% select(-label) %>%
      st_set_crs(crs_predictor)
  }
  
  # can also return the value of the predictor for each data point
  # env1 <- terra::extract(rast(predictor), vect(PTS.sp1)) %>% select(env1=2)
  # env2 <- terra::extract(rast(predictor), vect(PTS.sp2)) %>% select(env2=2)
  #return(list(sp1, sp2, env1, env2))
  
  return(list(sp1, sp2))
}