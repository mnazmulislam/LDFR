
title : Description of the required functions for fitting "Longitudinal dynamic functional regression" (LDFR)


1. data generator.R : 
mufun() : defines true mean function

phi_fun() : defines true eigen functions

onesubj() : generates data for single subject

sparsefun() : creates sparse data set

data_func() : generates data for all subjects
 
2. data fit.R       
Yestim() : fits the model for training data and test data


3. LDFR.R :
LDFR() : applies LDFR model and predicting the full trajectory
