Hello guys, this is a code respiritory for the NHDR method proposed by our research team. We propose a novel inference method for nonprobability samples that accounts for population heterogeneity when adjusting for selection bias. 
Our approach first identifies latent heterogeneous subpopulations using finite mixture models and incorporate these groups as a high-level variable within the mixed-effects models. A doubly robust estimator is then constructed to improve reliability and accuracy of inferences.
To use our proposed NHDR approach, please first run the main_function.R for preparation.
The main function of NHDR is NPinfer_ALC.
Suppose there is a nonprobability sample called dnp, and a reference sample called dref6_1, the common coavriates are job, region, sex, age, race, edu, health, and marry, with the corresponding variable types are multinomial, multinomial, bionimal, gaussian, multinomial, gaussian, and multinomial, respectively.
And the target outcome is Y, you can obtain an adjusted estimate of Y by our NHDR method by the following codes:
 NHDR <- NPinfer_ALC(data_NP=dnp,data_P=dref6_1,
          CIV=c('job','region','sex','age','race','edu','health','marry'),
          type_CIV=c('multinomial','multinomial','binomial','gaussian','multinomial','multinomial','gaussian','multinomial'),
          DV_NP='Y',DV_type='binomial',
          sample_weight='sample_weight',maxclassnum=5)
          
