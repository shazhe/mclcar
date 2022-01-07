# mclcar

## Aim
The ```mclcar``` provide two implementations of the Monte Carlo likelihood methods for estimating parameters in conditional auto-regressive (CAR) models for spatially correlated data on areal supports. The Monte Carlo likelihood methods are based on [Geyer and Thompson (1992)](http://www.jstor.org/stable/2345852) and two maximisation procedures are developed in [Sha (2016)](https://ora.ox.ac.uk/objects/uuid:6cc56943-2b4d-4931-895a-f3ab67e48e3a/). The first procedures iteratively maximise the Monte Carlo likelihood with a variance constraints on the step size and the second procedure use a response surface desgin method (rsm). 

The package is available on [CRAN](https://cran.r-project.org/) at https://CRAN.R-project.org/package=mclcar

This github repo contains latest update to the package and also maintains the previous version. 

## Installation
* Latest github version: 0.2.0
* Latest CRAN version: 0.2.0

### 1 From github
The latest github version can be installed by 
```R
library(devtools)
install_github("shazhe/mclcar/mclcar_0.2-0", build_vignettes=F, dependencies=T)
```

### 2 From package .tar.gz
Alternatively you can download the ```mclcar_x.x-x.tar.gz``` and install from the source. You can also download the ```pkg_skeleton``` folder and build the package and install.

### 3 From CRAN
Stable version can be installed directly from CRAN by calling ```install.packages("mclcar")``` in R.

## Features
1. Simulate data from the CAR models.
2. Estimate the CAR models using the Monte Carlo likelihood methods with either the iterative procedure or the rsm procedure.
3. Plot the likelihood surface evolution when using the rsm procedure.
4. The CAR models can be fiited by the package currently include the followings:
   * a direct CAR model -- linear model with CAR errors, on a torus or a map with given neighbourhood matrix,
   * a Binomial/Poisson GLM with latent CAR variables,
   * a hierachical CAR models with individual and district level CAR correlations.

