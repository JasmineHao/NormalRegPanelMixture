#!/bin/bash

python -u empirical_test.py plain_mixture > output_empirical_plain_mixture.log 2>&1 &

python -u empirical_test.py k_mixture > output_empirical_k_mixture.log 2>&1 &

python -u empirical_test.py kl_mixture > output_empirical_kl_mixture.log 2>&1 &

python -u empirical_test.py kl_mixture_spline > output_empirical_kl_mixture_spline.log 2>&1 &

python -u empirical_test.py ar1_plain_mixture > output_empirical_plain_mixture.log 2>&1 &

python -u empirical_test.py ar1_k_mixture > output_empirical_k_mixture.log 2>&1 &

python -u empirical_test.py ar1_kl_mixture > output_empirical_kl_mixture.log 2>&1 &

python -u empirical_test.py ar1_kl_mixture_spline > output_empirical_kl_mixture_spline.log 2>&1 &