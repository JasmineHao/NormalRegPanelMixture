#!/bin/bash

# python -u empirical_test.py plain_mixture > output_empirical_plain_mixture.log 2>&1 

# python -u empirical_test.py k_mixture > output_empirical_k_mixture.log 2>&1 

# python -u empirical_test.py kl_mixture > output_empirical_kl_mixture.log 2>&1 &

# python -u empirical_test.py kl_mixture_spline > output_empirical_kl_mixture_spline.log 2>&1 &

# python -u empirical_test.py ar1_plain_mixture > output_empirical_ar1_plain_mixture.log 2>&1 

# python -u empirical_test.py ar1_k_mixture > output_empirical_ar1_k_mixture.log 2>&1 

# python -u empirical_test.py ar1_kl_mixture > output_empirical_ar1_kl_mixture.log 2>&1 

# python -u empirical_test.py ar1_kl_mixture_spline > output_empirical_ar1_kl_mixture_spline.log 2>&1 

python -u empirical_test.py plain > output_empirical_plain.log 2>&1 
python -u empirical_test.py k > output_empirical_k.log 2>&1 & 
python -u empirical_test.py kl > output_empirical_kl.log 2>&1 &
python -u empirical_test.py kl_spline > output_empirical_kl_spline.log 2>&1 &

python -u empirical_test.py ar1_plain > output_empirical_ar1_plain.log 2>&1 &
python -u empirical_test.py ar1_k > output_empirical_ar1_k.log 2>&1 &
python -u empirical_test.py ar1_kl > output_empirical_ar1_kl.log 2>&1 &
python -u empirical_test.py ar1_kl_spline > output_empirical_ar1_kl_spline.log 2>&1 &


