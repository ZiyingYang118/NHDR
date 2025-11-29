R package for "Improved doubly robust inference with nonprobability survey samples using finite mixture models: Application to Health Monitoring SMS Survey data"

Implement functions to make an inference of a continuous/binary outcome for nonprobability sample.

**Installation**
```
if (!require("devtools")) install.packages("devtools")
devtools::install_github("ZiyingYang118/NHDR")
```

**Example**  

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
