

############################################################################################################

  title : Description of the required functions for fitting "Longitudinal dynamic functional regression" (LDFR)
   date : 2017/11

############################################################################################################


1. data generator.R : 

     [ i ]       mufun() : defines true mean function
     [ ii ]    phi_fun() : defines true eigen functions
     [ iii ]   onesubj() : generates data for single subject
     [ iv ]  sparsefun() : creates sparse data set
     [ v ]   data_func() : generates data for all subjects
 

2. data fit.R :  
     
     [ i ]      Yestim() : fits the model for training data and test data



3. LDFR.R :

     [ i ]        LDFR() : applies LDFR model and predicting the full trajectory
