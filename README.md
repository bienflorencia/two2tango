
# two2tango

<!-- badges: start -->
[![R-CMD-check](https://github.com/yourusername/two2tango/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/yourusername/two2tango/actions/workflows/R-CMD-check.yaml)
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
<!-- badges: end -->

`two2tango` simulates pairs of virtual species with a known, adjustable degree
of spatial association, for testing and benchmarking species co-occurrence
and association methods against a known ground truth.

Given an environmental predictor surface, `two2tango` generates two species'
intensity surfaces — one following a **Gaussian** (niche-shaped) response to
the environment, the other a **linear** response — and injects correlated
random noise into their log-intensities before sampling occurrence points
from an inhomogeneous Poisson process. Because the true correlation between
species is set directly as a simulation parameter, you always know the
"ground truth" association to validate a method against.

## Installation

You can install the development version of `two2tango` with:

``` r
# install.packages("pak")
pak::pak("bienflorencia/two2tango")
```

or with `devtools`:

``` r
# install.packages("devtools")
devtools::install_github("bienflorencia/two2tango")
```

## Functions

| Function | Description |
|---|---|
| `two2tango()` | Simulates two virtual species with a **Gaussian** (bell-shaped) response to an environmental predictor. |
| `two2waltz()` | Simulates two virtual species with a **linear** response to an environmental predictor. |

## Example

This basic example simulates two Gaussian-response species over a randomly
generated environmental surface, and plots their occurrence points against
the predictor:

``` r
library(two2tango)

set.seed(134)
size=10
env <- spatstat.geom::im(matrix(0, size, size),
          xrange=c(0,1), yrange=c(0,1))
xy <- expand.grid(x = env$xcol, y = env$yrow)
temp.dummy <- gstat::gstat(formula=z~1, locations=~x+y, dummy=TRUE, beta=0,
                    model=gstat::vgm(psill=0.1, range=50, model='Exp'), nmax=20)
temp <- predict(temp.dummy, newdata=xy, nsim=1, debug.level=0)
temp$sim1 <- scale(temp$sim1)

example <- two2tango(peak1=1, peak2=1,
                     mu1=-1, sigma1=0.5,
                     mu2=1, sigma2=0.5,
                     cov=0, predictor = temp)

sp1 <- example[[1]]
sp2 <- example[[2]]
cor.sp1.sp2 <- example[[3]]
env <- spatstat.geom::as.im(temp)

plot(env, main= paste0('cov=', round(cor.sp1.sp2$estimate,2)))
plot(sp1, col='white', add=TRUE)
plot(sp2, col='black', add=TRUE)
```

`two2waltz()` works the same way, but with `alpha1`/`alpha2` (intercepts)
and `beta1`/`beta2` (slopes) in place of `peak`/`mu`/`sigma`, since its
species respond linearly rather than following a Gaussian niche curve.

For a more detailed example — including simulating over a real
environmental raster — see the package vignette:

``` r
vignette("intro_to_two2tango", package = "two2tango")
```

## Dependencies

`two2tango` builds on:

- [terra](https://rspatial.github.io/terra/) and
  [spatstat](https://spatstat.org/) for raster/pixel-image handling and
  spatial point process simulation
- [sf](https://r-spatial.github.io/sf/) for point geometry output
- [mvtnorm](https://cran.r-project.org/package=mvtnorm) for correlated noise
  generation
- [dplyr](https://dplyr.tidyverse.org/) for internal data wrangling

## Citation

If you use `two2tango` in a publication, please cite it — see
`citation("two2tango")` once installed, or the `CITATION.cff` file in this
repository.

## Contributing

Bug reports, feature requests, and pull requests are welcome — please open an
[issue](https://github.com/bienflorencia/two2tango/issues) on GitHub.

## License

See [LICENSE.md](LICENSE.md) for details.

## Author 
Florencia Grattarola <a dir="ltr" href="http://orcid.org/0000-0001-8282-5732" target="_blank"><img class="is-rounded" src="https://upload.wikimedia.org/wikipedia/commons/0/06/ORCID_iD.svg" width="15"></a> 

## Acknowledgements
Thanks to Petr Keil and Alejandra Zarzo-Arias for early stages discussions about simulating virtual species.  
