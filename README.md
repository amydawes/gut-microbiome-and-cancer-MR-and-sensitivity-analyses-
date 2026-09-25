########### Readme for the analytical pipeline of running MR to assess the causal role of microbial traits on cancer outcomes##########

#This read me contains an example MR pipeline, for assessing the causal role of the microbiome on cancer,
#using prostate cancer as the cancer outcome of interest.

#The pipeline should be followed sequentially from scripts 1-6, where:

#script 1: Extract the variants required from the exposure and the outcome MR
#script 2: Perform main MR analysis across a primary outcome, and replication outcomes
#script 3: perform a manual exploration of pleiotropy
#scripts 4: scripts a-d should be followed sequentially to perform colocalisation analysis
#script 5: A lenient threshold MR analysis, where the P-value threshold of P<1x10-5 is used and traits where:
#a) directional consistency is identified and b) all traits at P<1x10-5 are tested as a sensitivity analysis
#script 6: Reverse MR is used to test for reverse causality of the cancer on the microbiome
#script 7: to assess the extent for which heterogeneity due to differences in genetic ancestry or population stratification may have an effect #on 
#multi-ancestry meta-analysis used as the exposure, we repeat main MR analyses for the MiBioGen microbiome traits, restricting only to #individuals and cohorts where all
#individuals are of European ancestry, and we compare the directionality of the results with the main analysis. 

