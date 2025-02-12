#!/bin/bash

# python -u empirical_test.py plain_mixture > output_empirical_plain_mixture.log 2>&1 &
# python -u empirical_test.py k_mixture > output_empirical_k_mixture.log 2>&1 &
# python -u empirical_test.py kl_mixture > output_empirical_kl_mixture.log 2>&1 &

# python -u empirical_test.py kl_mixture_spline > output_empirical_kl_mixture_spline.log 2>&1 &

# python -u empirical_test.py ar1_plain_mixture > output_empirical_ar1_plain_mixture.log 2>&1 &

# python -u empirical_test.py ar1_k_mixture > output_empirical_ar1_k_mixture.log 2>&1 

# python -u empirical_test.py ar1_kl_mixture > output_empirical_ar1_kl_mixture.log 2>&1 

# python -u empirical_test.py ar1_kl_mixture_spline > output_empirical_ar1_kl_mixture_spline.log 2>&1 

# -----
# python -u empirical_test.py plain > output_empirical_plain.log 2>&1 
# python -u empirical_test.py k > output_empirical_k.log 2>&1 & 
# python -u empirical_test.py kl > output_empirical_kl.log 2>&1 &
# python -u empirical_test.py kl_spline > output_empirical_kl_spline.log 2>&1 &

# -----
# python -u empirical_test.py ar1_k > output_empirical_ar1_k.log 2>&1 &
# python -u empirical_test.py ar1_kl > output_empirical_ar1_kl.log 2>&1 &
# python -u empirical_test.py ar1_kl_spline > output_empirical_ar1_kl_spline.log 2>&1 &


# Comparison case 1: need to make sure mixture version gives smaller number
# ----
# python -u empirical_test.py ar1_plain_mixture > output_empirical_ar1_plain_mixture.log 2>&1 & 
# python -u empirical_test.py ar1_plain > output_empirical_ar1_plain.log 2>&1 &
# # ----
# [*] Comparison case 2: ???Still don't know why
# python -u empirical_test.py ar1_k_mixture > output_empirical_ar1_k_mixture.log 2>&1 & 
# python -u empirical_test.py ar1_k > output_empirical_ar1_k.log 2>&1 &
# # ----
# # Comparison case 3:
# python -u empirical_test.py ar1_plain > output_empirical_ar1_plain.log 2>&1 &
# python -u empirical_test.py ar1_k > output_empirical_ar1_k.log 2>&1 &

# # ----
# # Comparison case 4:
# python -u empirical_test.py plain_mixture > output_empirical_plain_mixture.log 2>&1 & 
# python -u empirical_test.py k_mixture > output_empirical_k_mixture.log 2>&1 & 


python -u empirical_test.py plain_mixture_3 5 > output_empirical_plain_mixture_3.log 2>&1 &
python -u empirical_test.py k_mixture_3 5> output_empirical_k_mixture_3.log 2>&1 &

python -u empirical_test.py plain 3 > output_empirical_plain.log 2>&1 &
python -u empirical_test.py plain_mixture 3 > output_empirical_plain_mixture.log 2>&1 &

python -u empirical_test.py k 3 > output_empirical_k.log 2>&1 &
python -u empirical_test.py k_mixture 3 > output_empirical_k_mixture.log 2>&1 &

python -u empirical_test.py ar1_plain 3 > output_empirical_ar1_plain.log 2>&1 &
python -u empirical_test.py ar1_plain_mixture 3 > output_empirical_ar1_plain_mixture.log 2>&1 & 


python -u empirical_test.py ar1_k 3 > output_empirical_ar1_k.log 2>&1 &
python -u empirical_test.py ar1_k_mixture 3 > output_empirical_ar1_k_mixture.log 2>&1 & 

