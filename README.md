R package for "Improved doubly robust inference with nonprobability survey samples using finite mixture models: Application to Health Monitoring SMS Survey data"

Implement functions to make an inference of a continuous/binary outcome for nonprobability sample.

**Installation**
```
source(main_function.R)
source(simu_data.R)
```

**Example for application based on simuated data**  

Estimate the mean and 95%CI of a continuous outcome with simulated data using the NHDR method.   

1.Generate the simulated data using `simucode()`
```
library(NHDR)
data <- simu_data(N=20000,n_A=500,n_B=1000,R2=0.5,n_LC=3,icc=0.5,n_Z=2,n_P=1,Ytype='gaussian',dist=1)
data_SA <- data$nonprobability_data
data_SB <- data$reference_data
```
2.Estimate the mean and the 95%CI, with the upper limit of latent classes set as 5
```
NHDR <- NPinfer_ALC(data_NP=data_SA,data_P=data_SB,
                    CIV=c('Z1','Z2','Z3','Z4','Z5'),
                    type_CIV=c('binomial','gaussian','gaussian','gaussian','gaussian'),
                    DV_NP='Y',DV_type='gaussian',
                    sample_weight='sample_weight',maxclassnum=5)
mean <- NHDR$point_kwdr
CI <- NHDR$CI_kwdr
```
**Example for simulation study**  

1.Load the `simucode`
```
source(simulation.R)
```

2.Set the simulation parameters and run the simulation. For example, take the following as the illustrative parameter setup and the code shown below:
```
simucode(
  B = 3000,        # number of replications
  N = 100000,      # finite population size
  n_A = 500,        # size of the nonprobability sample
  n_B = 2000,       # size of the reference sample
  R2 = 0.5,          # correlation parameter between key variables
  n_LC = 2,           # number of latent classes
  icc = 0.5,          # intraclass correlation coefficient
  n_Z = 1,            # number of covariates with sign differences (example value)
  n_P = 1,            # number of subpopulations with opposing covariate effects
  Ytype = 'gaussian',  # outcome distribution
  dist = 1,             # how Z_i distributions vary across subpopulations
  maxclassnum = 4       # upper limit for latent classes
)

```
