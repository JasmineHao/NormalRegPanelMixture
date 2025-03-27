#!/bin/bash


python -u empirical_test.py plain 3 no_ciiu mY_share > output_empirical_plain.log 2>&1 &
python -u empirical_test.py plain_mixture 3 no_ciiu mY_share > output_empirical_plain_mixture.log 2>&1 &

python -u empirical_test.py plain_mixture_3 3 no_ciiu mY_share > output_empirical_plain_mixture_3.log 2>&1 &

python -u empirical_test.py plain_mixture_4 3 no_ciiu mY_share > output_empirical_plain_mixture_4.log 2>&1 &


python -u empirical_test.py k 3 ciiu mY_share > output_empirical_k_ciiu.log 2>&1 &
python -u empirical_test.py k_mixture 3 ciiu mY_share > output_empirical_k_mixture_ciiu.log 2>&1 &
python -u empirical_test.py k_mixture_3 3 ciiu mY_share > output_empirical_k_mixture_3_ciiu.log 2>&1 &
python -u empirical_test.py k_mixture_4 3 ciiu mY_share > output_empirical_k_mixture_4_ciiu.log 2>&1 &

# AR1 
python -u empirical_test.py ar1_plain 3 no_ciiu mY_share > output_empirical_ar1_plain.log 2>&1 &
python -u empirical_test.py ar1_plain_mixture 3 no_ciiu mY_share > output_empirical_ar1_plain_mixture.log 2>&1 & 

python -u empirical_test.py ar1_plain_mixture_3 3 no_ciiu mY_share > output_empirical_ar1_plain_mixture_3.log 2>&1 & 


python -u empirical_test.py ar1_k 3 ciiu mY_share > output_empirical_ar1_k.log 2>&1 &
python -u empirical_test.py ar1_k_mixture 3 ciiu mY_share > output_empirical_ar1_k_mixture.log 2>&1 & 
python -u empirical_test.py ar1_k_mixture_3 3 ciiu mY_share > output_empirical_ar1_k_mixture_3.log 2>&1 & 


# python -u empirical_test.py k 3 no_ciiu mY_share > output_empirical_k.log 2>&1 &
# python -u empirical_test.py k_mixture 3 no_ciiu mY_share > output_empirical_k_mixture.log 2>&1 &

# python -u empirical_test.py k_mixture_3 3 no_ciiu mY_share > output_empirical_k_mixture_3.log 2>&1 &
# python -u empirical_test.py k_mixture_4 3 no_ciiu mY_share > output_empirical_k_mixture_4.log 2>&1 &


# python -u empirical_test.py ar1_k 3 no_ciiu mY_share > output_empirical_ar1_k.log 2>&1 &
# python -u empirical_test.py ar1_k_mixture 3 no_ciiu mY_share > output_empirical_ar1_k_mixture.log 2>&1 & 
# python -u empirical_test.py ar1_k_mixture_3 3 no_ciiu mY_share > output_empirical_ar1_k_mixture_3.log 2>&1 & 

# --
# include ciiu
# python -u empirical_test.py plain 3 ciiu mY_share > output_empirical_plain_ciiu.log 2>&1 &
# python -u empirical_test.py plain_mixture 3 ciiu mY_share > output_empirical_plain_mixture_ciiu.log 2>&1 &
# python -u empirical_test.py plain_mixture_3 3 ciiu mY_share > output_empirical_plain_mixture_3_ciiu.log 2>&1 &
# python -u empirical_test.py plain_mixture_4 3 ciiu mY_share > output_empirical_plain_mixture_4_ciiu.log 2>&1 &



# *******************************

# kmshare

# python -u empirical_test.py kmshare 3 no_ciiu mY_share > output_empirical_kmshare.log 2>&1 &
# python -u empirical_test.py kmshare_mixture 3 no_ciiu mY_share > output_empirical_kmshare_mixture.log 2>&1 &
# python -u empirical_test.py kmshare_mixture_3 3 no_ciiu mY_share > output_empirical_kmshare_mixture_3.log 2>&1 &
# python -u empirical_test.py kmshare_mixture_4 3 no_ciiu mY_share > output_empirical_kmshare_mixture_4.log 2>&1 &



# python -u empirical_test.py kmshare 3 ciiu mY_share > output_empirical_kmshare.log 2>&1 &
# python -u empirical_test.py kmshare_mixture 3 ciiu mY_share > output_empirical_kmshare_mixture.log 2>&1 &
# python -u empirical_test.py kmshare_mixture_3 3 ciiu mY_share > output_empirical_kmshare_mixture_3.log 2>&1 &
# python -u empirical_test.py kmshare_mixture_4 3 ciiu mY_share > output_empirical_kmshare_mixture_4.log 2>&1 &

# -------
python -u empirical_test.py nonpar 3 no_ciiu mY_share