#!/bin/bash

# Script 1: ChileanTest_Regressor_KL_Spline.R
Rscript experiment/EmpiricalTests/ChileanTest_Regressor_KL_Spline.R
Rscript experiment/EmpiricalTests/Japan_Regressor_KL_Spline.R


# Script 2: ChileanTest_Regressor_KL_AR1.R
Rscript experiment/EmpiricalTests/ChileanTest_Regressor_KL_AR1.R
Rscript experiment/EmpiricalTests/JapanTest_Regressor_KL_AR1.R

# Script 3: ChileanTest_Regressor_KL_AR1_Spline.R
Rscript experiment/EmpiricalTests/ChileanTest_Regressor_KL_AR1_Spline.R
Rscript experiment/EmpiricalTests/JapanTest_Regressor_KL_AR1_Spline.R


# Script 4: ChileanTest_Regressor_K.R
Rscript experiment/EmpiricalTests/ChileanTest_Regressor_K.R
Rscript experiment/EmpiricalTests/JapanTest_Regressor_K.R

# Script 5: ChileanTest_Regressor_K_AR1.R
Rscript experiment/EmpiricalTests/ChileanTest_Regressor_K_AR1.R
Rscript experiment/EmpiricalTests/JapanTest_Regressor_K_AR1.R

# Script 6: ChileanTest_AR1.R
Rscript experiment/EmpiricalTests/ChileanTest_AR1.R
Rscript experiment/EmpiricalTests/JapanTest_AR1.R


# Script 7: Plain 
Rscript experiment/EmpiricalTests/ChileanTest.R
Rscript experiment/EmpiricalTests/JapanTest.R


# Script 8: K and L
Rscript experiment/EmpiricalTests/ChileanTest_Regressor_KL.R
Rscript experiment/EmpiricalTests/JapanTest_Regressor_KL.R