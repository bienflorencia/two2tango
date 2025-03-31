# Two2tango

This repository contains the `R` function `two2tango()`, and other auxiliary ones, that will allow you to create two virtual species that have a certain response to the environment and degree of association between them. For instructions on running the function, including example data, please see the vignette tutorial.  

The model assumes an underlying Poisson point pattern, and a Gaussian response:

$\lambda_1 = peak_1 \times exp^{(-1/2 \times \dfrac{(temp - \mu_1)^2}{\sigma_1^2} + e_1)}$  
$\lambda_2 = peak_2 \times exp^{(-1/2 \times \dfrac{(temp - \mu_2)^2}{\sigma_2^2} + e_2)}$  

$e_{ij} \sim \textsf{MVN}(0, \Sigma)$  

Where, $\Sigma = \begin{bmatrix} var_{1,1} & cov_{1,2} \\ cov_{2,1} & var_{2,2} \end{bmatrix}$  

The inverse of the covariance matrix is called the precision matrix, denoted by $\tau = {\Sigma}^{-1}$.

## Author 
Florencia Grattarola <a dir="ltr" href="http://orcid.org/0000-0001-8282-5732" target="_blank"><img class="is-rounded" src="https://upload.wikimedia.org/wikipedia/commons/0/06/ORCID_iD.svg" width="15"></a>. 

## Acknowledgements
Thanks to Petr Keil and Alejandra Zarzo-Arias.  