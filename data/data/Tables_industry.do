***********************
* Descriptive Statistics
***********************


****** (1) Industry 3220 ********

* Creat Summary Statistics
clear
set memory 500m
set matsize 800
 
use "C:\Documents and Settings\hkasahar\My Documents\work\project\exportimport\data\DATA_KL_EXIM5.dta"
drop if year<90
replace L = 0 if L ==.
sort padron
by padron: egen min_M = min(M) if L>0
drop if min_M == 0
by padron: egen min_Y = min(Y) if L>0
drop if min_Y == 0
by padron: egen sum_L = sum(L)
drop if sum_L == 0

by padron: egen mciiu = mode(ciiu), minmode
by padron: egen sdciiu = sd(ciiu) 
by padron: egen sumreentry = sum(reentry)
keep if mciiu==3220
*keep if sdciiu==0|sdciiu==.
*keep if sumreentry==0

replace L_bl = 0 if L_bl==.
replace L_wh = 0 if L_wh==.
replace w_bl = 0 if w_bl==.
replace w_wh = 0 if w_wh==.
gen kk = rkm+rkt+rkb if L>0
gen ln_L = log(L) if L>0
gen ln_Y = log(rY) if L>0
gen ln_va_L = log(rVA/L) if rVA>0&L>0
gen ww = (w_wh*L_wh + w_bl*L_bl)/L if L>0
gen ln_ww = log(ww) if L>0
gen ln_bl_L = log(L_bl/L) if L>0&L_bl~=.
gen bl_L = L_bl/L  if L>0&L_bl~=.
gen ln_wh_L = log(L_wh/L) if L>0&L_wh~=.
gen wh_L = L_wh/L  if L>0&L_wh~=.

gen ln_kk_L = log(kk/L) if L>0&kk>0
tab year, gen(dyear)

gen d_x1 = (1-d_m)*d_e if L>0
gen d_m1 = (1-d_e)*d_m if L>0
gen d_xm = d_e*d_m if L>0

* Pooled OLS (1990-1996)
reg ln_L d_x1 d_m1 d_xm dyear2-dyear7 if L>0, robust 
reg ln_Y d_x1 d_m1 d_xm dyear2-dyear7 if L>0, robust 
reg ln_va_L d_x1 d_m1 d_xm ln_L dyear2-dyear7 if L>0, robust 
reg ln_ww d_x1 d_m1 d_xm ln_L dyear2-dyear7 if L>0, robust 
reg ln_wh_L d_x1 d_m1 d_xm ln_L dyear2-dyear7 if L>0, robust 
reg ln_kk_L d_x1 d_m1 d_xm ln_L dyear2-dyear7 if L>0, robust  

* Fixed Effects (1990-1996)
xtreg ln_L d_x1 d_m1 d_xm dyear2-dyear7 if L>0, fe i(padron) robust  
xtreg ln_Y d_x1 d_m1 d_xm ln_L dyear2-dyear7 if L>0, fe i(padron) robust  
xtreg ln_va_L d_x1 d_m1 d_xm ln_L dyear2-dyear7 if L>0, fe i(padron) robust  
xtreg ln_ww d_x1 d_m1 d_xm ln_L dyear2-dyear7 if L>0, fe i(padron) robust   
xtreg ln_wh_L d_x1 d_m1 d_xm ln_L dyear2-dyear7 if L>0, fe i(padron) robust  
xtreg ln_kk_L d_x1 d_m1 d_xm ln_L dyear2-dyear7 if L>0, fe i(padron) robust  

gen d_L = 0 
replace d_L = 1 if L>0
sort padron
by padron: egen sum_d_L = sum(d_L)
* Note Stata reports the `wrong' number of observations for fixed effects: the following number should be subtracted
summarize sum_d_L if year==90&sum_d_L==1



* Old Table 1
 
sort year
by year: summarize d_x1 d_m1 d_xm
summarize d_x1 d_m1 d_xm

sort year
by year: egen sum_E = sum(export) if L>0
by year: egen sum_E_x1 = sum(export) if d_x1>0
by year: egen sum_E_xm = sum(export) if d_xm>0
by year: egen sum_Mi = sum(Mi) if L>0
by year: egen sum_M_m1 = sum(Mi) if d_m1>0
by year: egen sum_M_xm = sum(Mi) if d_xm>0

gen E_x1 = sum_E_x1/sum_E
gen E_xm = sum_E_xm/sum_E
gen M_m1 = sum_M_m1/sum_Mi
gen M_xm = sum_M_xm/sum_Mi

by year: summarize E_x1 E_xm M_m1 M_xm

sort year
by year: egen sum_Y = sum(Y) if L>0
by year: egen sum_Y_x1 = sum(Y) if d_x1>0&L>0
by year: egen sum_Y_m1 = sum(Y) if d_m1>0&L>0
by year: egen sum_Y_xm = sum(Y) if d_xm>0&L>0

gen Y_x1 = sum_Y_x1/sum_Y 
gen Y_m1 = sum_Y_m1/sum_Y 
gen Y_xm = sum_Y_xm/sum_Y 
by year: summarize Y_x1 Y_m1 Y_xm




replace entry = 0
replace exit = 0

keep entry exit rY L rexport rM rMi padron year
reshape wide entry exit rY L rexport rM rMi, i(padron) j(year)

replace entry91 =1 if L91>0&L90==0
replace entry92 =1 if L92>0&L91==0&L90==0
replace entry93 =1 if L93>0&L92==0&L91==0&L90==0
replace entry94 =1 if L94>0&L93==0&L92==0&L91==0&L90==0
replace entry95 =1 if L95>0&L94==0&L93==0&L92==0&L91==0&L90==0
replace entry96 =1 if L96>0&L95==0&L94==0&L93==0&L92==0&L91==0&L90==0

replace exit91 = 1 if L90>0&L91==0
replace L92 = 0 if exit91==1
replace exit92 = 1 if L91>0&L92==0
replace L93 = 0 if exit92==1|exit91==1
replace exit93 = 1 if L92>0&L93==0
replace L94 = 0 if exit93==1|exit92==1|exit91==1
replace exit94 = 1 if L93>0&L94==0 
replace L95 = 0 if exit94==1|exit93==1|exit92==1|exit91==1
replace exit95 = 1 if L94>0&L95==0  
replace L96 = 0 if exit95==1|exit94==1|exit93==1|exit92==1|exit91==1
replace exit96 = 1 if L95>0&L96==0  

reshape long entry exit rY L rexport rM rMi, i(padron) j(year) 
 
gen d_L = 0 
replace d_L = 1 if L>0

sort year
by year: egen sum_entry = sum(entry)
by year: egen sum_exit = sum(exit)
by year: egen nplant = sum(d_L)
gen entry_nplant = sum_entry/nplant
gen exit_nplant = sum_exit/nplant

by year: summarize nplant sum_entry sum_exit 
by year: summarize entry_nplant exit_nplant 

gen rf_r = rexport/rY if L>0&rexport>0
gen mi_m = rMi/rM if L>0&rMi>0

*sort year
*by year: summarize rY rM L entry_nplant exit_nplant  if L>0
*by year: summarize rexport rf_r if rexport>0
*by year: summarize rMi mi_m if rMi>0

summarize rY rM L if L>0&year==90
summarize rexport rf_r if L>0&rexport>0&year==90
summarize rMi mi_m if rMi>0&year==90

summarize rY rM L if L>0&year==96
summarize rexport rf_r if rexport>0&year==96
summarize rMi mi_m if rMi>0&year==96

sort year
by year: summarize rf_r if L>0&rexport>0
by year: summarize mi_m if rMi>0
 

****** (2) Industry 3560 ********

* Creat Summary Statistics
clear 
set memory 500m
set matsize 800
use "C:\Documents and Settings\hkasahar\My Documents\work\project\exportimport\data\DATA_KL_EXIM5.dta"
drop if year<90 
replace L = 0 if L ==.
sort padron
by padron: egen min_M = min(M) if L>0
drop if min_M == 0
by padron: egen min_Y = min(Y) if L>0
drop if min_Y == 0
by padron: egen sum_L = sum(L)
drop if sum_L == 0

sort padron
by padron: egen mciiu = mode(ciiu), minmode
by padron: egen sdciiu = sd(ciiu) 
by padron: egen sumreentry = sum(reentry)
keep if mciiu==3560
*keep if sdciiu==0|sdciiu==.
*keep if sumreentry==0 


replace L_bl = 0 if L_bl==.
replace L_wh = 0 if L_wh==.
replace w_bl = 0 if w_bl==.
replace w_wh = 0 if w_wh==.
gen kk = rkm+rkt+rkb if L>0
gen ln_L = log(L) if L>0
gen ln_Y = log(rY) if L>0
gen ln_va_L = log(rVA/L) if rVA>0&L>0
gen ww = (w_wh*L_wh + w_bl*L_bl)/L if L>0
gen ln_ww = log(ww) if L>0
gen ln_bl_L = log(L_bl/L) if L>0&L_bl~=.
gen bl_L = L_bl/L  if L>0&L_bl~=.
gen ln_wh_L = log(L_wh/L) if L>0&L_wh~=.
gen wh_L = L_wh/L  if L>0&L_wh~=.

gen ln_kk_L = log(kk/L) if L>0&kk>0
tab year, gen(dyear)

gen d_x1 = (1-d_m)*d_e if L>0
gen d_m1 = (1-d_e)*d_m if L>0
gen d_xm = d_e*d_m if L>0

* Pooled OLS (1990-1996)
reg ln_L d_x1 d_m1 d_xm dyear2-dyear7 if L>0, robust 
reg ln_Y d_x1 d_m1 d_xm dyear2-dyear7 if L>0, robust 
reg ln_va_L d_x1 d_m1 d_xm ln_L dyear2-dyear7 if L>0, robust 
reg ln_ww d_x1 d_m1 d_xm ln_L dyear2-dyear7 if L>0, robust 
reg ln_wh_L d_x1 d_m1 d_xm ln_L dyear2-dyear7 if L>0, robust 
reg ln_kk_L d_x1 d_m1 d_xm ln_L dyear2-dyear7 if L>0, robust  

* Fixed Effects (1990-1996)
xtreg ln_L d_x1 d_m1 d_xm dyear2-dyear7 if L>0, fe i(padron) robust  
xtreg ln_Y d_x1 d_m1 d_xm ln_L dyear2-dyear7 if L>0, fe i(padron) robust  
xtreg ln_va_L d_x1 d_m1 d_xm ln_L dyear2-dyear7 if L>0, fe i(padron) robust  
xtreg ln_ww d_x1 d_m1 d_xm ln_L dyear2-dyear7 if L>0, fe i(padron) robust   
xtreg ln_wh_L d_x1 d_m1 d_xm ln_L dyear2-dyear7 if L>0, fe i(padron) robust  
xtreg ln_kk_L d_x1 d_m1 d_xm ln_L dyear2-dyear7 if L>0, fe i(padron) robust  

gen d_L = 0 
replace d_L = 1 if L>0
sort padron
by padron: egen sum_d_L = sum(d_L)
* Note Stata reports the `wrong' number of observations for fixed effects: the following number should be subtracted
summarize sum_d_L if year==90&sum_d_L==1



* Old Table 1
 
sort year
by year: summarize d_x1 d_m1 d_xm
summarize d_x1 d_m1 d_xm

sort year
by year: egen sum_E = sum(export) if L>0
by year: egen sum_E_x1 = sum(export) if d_x1>0
by year: egen sum_E_xm = sum(export) if d_xm>0
by year: egen sum_Mi = sum(Mi) if L>0
by year: egen sum_M_m1 = sum(Mi) if d_m1>0
by year: egen sum_M_xm = sum(Mi) if d_xm>0

gen E_x1 = sum_E_x1/sum_E
gen E_xm = sum_E_xm/sum_E
gen M_m1 = sum_M_m1/sum_Mi
gen M_xm = sum_M_xm/sum_Mi

by year: summarize E_x1 E_xm M_m1 M_xm

sort year
by year: egen sum_Y = sum(Y) if L>0
by year: egen sum_Y_x1 = sum(Y) if d_x1>0&L>0
by year: egen sum_Y_m1 = sum(Y) if d_m1>0&L>0
by year: egen sum_Y_xm = sum(Y) if d_xm>0&L>0

gen Y_x1 = sum_Y_x1/sum_Y 
gen Y_m1 = sum_Y_m1/sum_Y 
gen Y_xm = sum_Y_xm/sum_Y 
by year: summarize Y_x1 Y_m1 Y_xm





replace entry = 0
replace exit = 0

keep entry exit rY L rexport rM rMi padron year
reshape wide entry exit rY L rexport rM rMi, i(padron) j(year)

replace entry91 =1 if L91>0&L90==0
replace entry92 =1 if L92>0&L91==0&L90==0
replace entry93 =1 if L93>0&L92==0&L91==0&L90==0
replace entry94 =1 if L94>0&L93==0&L92==0&L91==0&L90==0
replace entry95 =1 if L95>0&L94==0&L93==0&L92==0&L91==0&L90==0
replace entry96 =1 if L96>0&L95==0&L94==0&L93==0&L92==0&L91==0&L90==0

replace exit91 = 1 if L90>0&L91==0
replace L92 = 0 if exit91==1
replace exit92 = 1 if L91>0&L92==0
replace L93 = 0 if exit92==1|exit91==1
replace exit93 = 1 if L92>0&L93==0
replace L94 = 0 if exit93==1|exit92==1|exit91==1
replace exit94 = 1 if L93>0&L94==0 
replace L95 = 0 if exit94==1|exit93==1|exit92==1|exit91==1
replace exit95 = 1 if L94>0&L95==0  
replace L96 = 0 if exit95==1|exit94==1|exit93==1|exit92==1|exit91==1
replace exit96 = 1 if L95>0&L96==0  

reshape long entry exit rY L rexport rM rMi, i(padron) j(year) 
 
gen d_L = 0 
replace d_L = 1 if L>0

sort year
by year: egen sum_entry = sum(entry)
by year: egen sum_exit = sum(exit)
by year: egen nplant = sum(d_L)
gen entry_nplant = sum_entry/nplant
gen exit_nplant = sum_exit/nplant

by year: summarize nplant sum_entry sum_exit 
by year: summarize entry_nplant exit_nplant 

gen rf_r = rexport/rY if L>0&rexport>0
gen mi_m = rMi/rM if L>0&rMi>0

*sort year
*by year: summarize rY rM L entry_nplant exit_nplant  if L>0
*by year: summarize rexport rf_r if rexport>0
*by year: summarize rMi mi_m if rMi>0


summarize rY rM L if L>0&year==90
summarize rexport rf_r if rexport>0&year==90
summarize rMi mi_m if rMi>0&year==90

sort year
by year: summarize rf_r if rexport>0
by year: summarize mi_m if rMi>0


****** (3) Industry 3813 ********

* Creat Summary Statistics
clear
set memory 500m
set matsize 800
use "C:\Documents and Settings\hkasahar\My Documents\work\project\exportimport\data\DATA_KL_EXIM5.dta"
drop if year<90 
replace L = 0 if L ==.
sort padron
by padron: egen min_M = min(M) if L>0
drop if min_M == 0
by padron: egen min_Y = min(Y) if L>0
drop if min_Y == 0
by padron: egen sum_L = sum(L)
drop if sum_L == 0


replace ciiu = 3813 if ciiu==3814|ciiu==3815

sort padron
by padron: egen mciiu = mode(ciiu), minmode
by padron: egen sdciiu = sd(ciiu) 
by padron: egen sumreentry = sum(reentry)
keep if mciiu==3813 
*keep if sdciiu==0|sdciiu==.
*keep if sumreentry==0

*drop if dplnt==1
replace L = 0 if L ==.
sort padron
by padron: egen min_M = min(M) if L>0
drop if min_M == 0
by padron: egen min_Y = min(Y) if L>0
drop if min_Y == 0
by padron: egen sum_L = sum(L)
drop if sum_L == 0


replace L_bl = 0 if L_bl==.
replace L_wh = 0 if L_wh==.
replace w_bl = 0 if w_bl==.
replace w_wh = 0 if w_wh==.
gen kk = rkm+rkt+rkb if L>0
gen ln_L = log(L) if L>0
gen ln_Y = log(rY) if L>0
gen ln_va_L = log(rVA/L) if rVA>0&L>0
gen ww = (w_wh*L_wh + w_bl*L_bl)/L if L>0
gen ln_ww = log(ww) if L>0
gen ln_bl_L = log(L_bl/L) if L>0&L_bl~=.
gen bl_L = L_bl/L  if L>0&L_bl~=.
gen ln_wh_L = log(L_wh/L) if L>0&L_wh~=.
gen wh_L = L_wh/L  if L>0&L_wh~=.

gen ln_kk_L = log(kk/L) if L>0&kk>0
tab year, gen(dyear)

gen d_x1 = (1-d_m)*d_e if L>0
gen d_m1 = (1-d_e)*d_m if L>0
gen d_xm = d_e*d_m if L>0

* Pooled OLS (1990-1996)
reg ln_L d_x1 d_m1 d_xm dyear2-dyear7 if L>0, robust 
reg ln_Y d_x1 d_m1 d_xm dyear2-dyear7 if L>0, robust 
reg ln_va_L d_x1 d_m1 d_xm ln_L dyear2-dyear7 if L>0, robust 
reg ln_ww d_x1 d_m1 d_xm ln_L dyear2-dyear7 if L>0, robust 
reg ln_wh_L d_x1 d_m1 d_xm ln_L dyear2-dyear7 if L>0, robust 
reg ln_kk_L d_x1 d_m1 d_xm ln_L dyear2-dyear7 if L>0, robust  

* Fixed Effects (1990-1996)
xtreg ln_L d_x1 d_m1 d_xm dyear2-dyear7 if L>0, fe i(padron) robust  
xtreg ln_Y d_x1 d_m1 d_xm ln_L dyear2-dyear7 if L>0, fe i(padron) robust  
xtreg ln_va_L d_x1 d_m1 d_xm ln_L dyear2-dyear7 if L>0, fe i(padron) robust  
xtreg ln_ww d_x1 d_m1 d_xm ln_L dyear2-dyear7 if L>0, fe i(padron) robust   
xtreg ln_wh_L d_x1 d_m1 d_xm ln_L dyear2-dyear7 if L>0, fe i(padron) robust  
xtreg ln_kk_L d_x1 d_m1 d_xm ln_L dyear2-dyear7 if L>0, fe i(padron) robust  

gen d_L = 0 
replace d_L = 1 if L>0
sort padron
by padron: egen sum_d_L = sum(d_L)
* Note Stata reports the `wrong' number of observations for fixed effects: the following number should be subtracted
summarize sum_d_L if year==90&sum_d_L==1




replace entry = 0
replace exit = 0

keep entry exit rY L rexport rM rMi padron year
reshape wide entry exit rY L rexport rM rMi, i(padron) j(year)

replace entry91 =1 if L91>0&L90==0
replace entry92 =1 if L92>0&L91==0&L90==0
replace entry93 =1 if L93>0&L92==0&L91==0&L90==0
replace entry94 =1 if L94>0&L93==0&L92==0&L91==0&L90==0
replace entry95 =1 if L95>0&L94==0&L93==0&L92==0&L91==0&L90==0
replace entry96 =1 if L96>0&L95==0&L94==0&L93==0&L92==0&L91==0&L90==0

replace exit91 = 1 if L90>0&L91==0
replace L92 = 0 if exit91==1
replace exit92 = 1 if L91>0&L92==0
replace L93 = 0 if exit92==1|exit91==1
replace exit93 = 1 if L92>0&L93==0
replace L94 = 0 if exit93==1|exit92==1|exit91==1
replace exit94 = 1 if L93>0&L94==0 
replace L95 = 0 if exit94==1|exit93==1|exit92==1|exit91==1
replace exit95 = 1 if L94>0&L95==0  
replace L96 = 0 if exit95==1|exit94==1|exit93==1|exit92==1|exit91==1
replace exit96 = 1 if L95>0&L96==0  

reshape long entry exit rY L rexport rM rMi, i(padron) j(year) 
 
gen d_L = 0 
replace d_L = 1 if L>0

sort year
by year: egen sum_entry = sum(entry)
by year: egen sum_exit = sum(exit)
by year: egen nplant = sum(d_L)
gen entry_nplant = sum_entry/nplant
gen exit_nplant = sum_exit/nplant

by year: summarize nplant sum_entry sum_exit 
by year: summarize entry_nplant exit_nplant 

gen rf_r = rexport/rY if L>0&rexport>0
gen mi_m = rMi/rM if L>0&rMi>0

*sort year
*by year: summarize rY rM L entry_nplant exit_nplant  if L>0
*by year: summarize rexport rf_r if rexport>0
*by year: summarize rMi mi_m if rMi>0


summarize rY rM L if L>0&year==90
summarize rexport rf_r if rexport>0&year==90
summarize rMi mi_m if rMi>0&year==90


sort year
by year: summarize rf_r if rexport>0
by year: summarize mi_m if rMi>0



****** (4) Industry 381 ********

* Creat Summary Statistics
clear
set memory 500m
set matsize 800
use "C:\Documents and Settings\hkasahar\My Documents\work\project\exportimport\data\DATA_KL_EXIM5.dta"
drop if year<90 
replace L = 0 if L ==.
sort padron
by padron: egen min_M = min(M) if L>0
drop if min_M == 0
by padron: egen min_Y = min(Y) if L>0
drop if min_Y == 0
by padron: egen sum_L = sum(L)
drop if sum_L == 0

sort padron
by padron: egen sumreentry = sum(reentry)
keep if mcu==381 
 

replace L_bl = 0 if L_bl==.
replace L_wh = 0 if L_wh==.
replace w_bl = 0 if w_bl==.
replace w_wh = 0 if w_wh==.
gen kk = rkm+rkt+rkb if L>0
gen ln_L = log(L) if L>0
gen ln_Y = log(rY) if L>0
gen ln_va_L = log(rVA/L) if rVA>0&L>0
gen ww = (w_wh*L_wh + w_bl*L_bl)/L if L>0
gen ln_ww = log(ww) if L>0
gen ln_bl_L = log(L_bl/L) if L>0&L_bl~=.
gen bl_L = L_bl/L  if L>0&L_bl~=.
gen ln_wh_L = log(L_wh/L) if L>0&L_wh~=.
gen wh_L = L_wh/L  if L>0&L_wh~=.

gen ln_kk_L = log(kk/L) if L>0&kk>0
tab year, gen(dyear)

gen d_x1 = (1-d_m)*d_e if L>0
gen d_m1 = (1-d_e)*d_m if L>0
gen d_xm = d_e*d_m if L>0

* Pooled OLS (1990-1996)
reg ln_L d_x1 d_m1 d_xm dyear2-dyear7 if L>0, robust 
reg ln_Y d_x1 d_m1 d_xm dyear2-dyear7 if L>0, robust 
reg ln_va_L d_x1 d_m1 d_xm ln_L dyear2-dyear7 if L>0, robust 
reg ln_ww d_x1 d_m1 d_xm ln_L dyear2-dyear7 if L>0, robust 
reg ln_wh_L d_x1 d_m1 d_xm ln_L dyear2-dyear7 if L>0, robust 
reg ln_kk_L d_x1 d_m1 d_xm ln_L dyear2-dyear7 if L>0, robust  

* Fixed Effects (1990-1996)
xtreg ln_L d_x1 d_m1 d_xm dyear2-dyear7 if L>0, fe i(padron) robust  
xtreg ln_Y d_x1 d_m1 d_xm ln_L dyear2-dyear7 if L>0, fe i(padron) robust  
xtreg ln_va_L d_x1 d_m1 d_xm ln_L dyear2-dyear7 if L>0, fe i(padron) robust  
xtreg ln_ww d_x1 d_m1 d_xm ln_L dyear2-dyear7 if L>0, fe i(padron) robust   
xtreg ln_wh_L d_x1 d_m1 d_xm ln_L dyear2-dyear7 if L>0, fe i(padron) robust  
xtreg ln_kk_L d_x1 d_m1 d_xm ln_L dyear2-dyear7 if L>0, fe i(padron) robust  

gen d_L = 0 
replace d_L = 1 if L>0
sort padron
by padron: egen sum_d_L = sum(d_L)
* Note Stata reports the `wrong' number of observations for fixed effects: the following number should be subtracted
summarize sum_d_L if year==90&sum_d_L==1




replace entry = 0
replace exit = 0

keep entry exit rY L rexport rM rMi padron year
reshape wide entry exit rY L rexport rM rMi, i(padron) j(year)

replace entry91 =1 if L91>0&L90==0
replace entry92 =1 if L92>0&L91==0&L90==0
replace entry93 =1 if L93>0&L92==0&L91==0&L90==0
replace entry94 =1 if L94>0&L93==0&L92==0&L91==0&L90==0
replace entry95 =1 if L95>0&L94==0&L93==0&L92==0&L91==0&L90==0
replace entry96 =1 if L96>0&L95==0&L94==0&L93==0&L92==0&L91==0&L90==0

replace exit91 = 1 if L90>0&L91==0
replace L92 = 0 if exit91==1
replace exit92 = 1 if L91>0&L92==0
replace L93 = 0 if exit92==1|exit91==1
replace exit93 = 1 if L92>0&L93==0
replace L94 = 0 if exit93==1|exit92==1|exit91==1
replace exit94 = 1 if L93>0&L94==0 
replace L95 = 0 if exit94==1|exit93==1|exit92==1|exit91==1
replace exit95 = 1 if L94>0&L95==0  
replace L96 = 0 if exit95==1|exit94==1|exit93==1|exit92==1|exit91==1
replace exit96 = 1 if L95>0&L96==0  

reshape long entry exit rY L rexport rM rMi, i(padron) j(year) 
 
gen d_L = 0 
replace d_L = 1 if L>0

sort year
by year: egen sum_entry = sum(entry)
by year: egen sum_exit = sum(exit)
by year: egen nplant = sum(d_L)
gen entry_nplant = sum_entry/nplant
gen exit_nplant = sum_exit/nplant

by year: summarize nplant sum_entry sum_exit 
by year: summarize entry_nplant exit_nplant 

gen rf_r = rexport/rY if L>0&rexport>0
gen mi_m = rMi/rM if L>0&rMi>0

*sort year
*by year: summarize rY rM L entry_nplant exit_nplant  if L>0
*by year: summarize rexport rf_r if rexport>0
*by year: summarize rMi mi_m if rMi>0


summarize rY rM L if L>0&year==90
summarize rexport rf_r if rexport>0&year==90
summarize rMi mi_m if rMi>0&year==90


sort year
by year: summarize rf_r if rexport>0
by year: summarize mi_m if rMi>0



****** (5) Industry 311 ********

* Creat Summary Statistics
clear
set memory 500m
set matsize 800
use "C:\Documents and Settings\hkasahar\My Documents\work\project\exportimport\data\DATA_KL_EXIM5.dta"
drop if year<90 
replace L = 0 if L ==.
sort padron
by padron: egen min_M = min(M) if L>0
drop if min_M == 0
by padron: egen min_Y = min(Y) if L>0
drop if min_Y == 0
by padron: egen sum_L = sum(L)
drop if sum_L == 0

sort padron
by padron: egen mciiu = mode(ciiu), minmode
by padron: egen sumreentry = sum(reentry)
keep if mciiu==3111|mciiu==3112|mciiu==3113|mciiu==3114|mciiu==3115|mciiu==3116|mciiu==3118|mciiu==3119
 

replace L_bl = 0 if L_bl==.
replace L_wh = 0 if L_wh==.
replace w_bl = 0 if w_bl==.
replace w_wh = 0 if w_wh==.
gen kk = rkm+rkt+rkb if L>0
gen ln_L = log(L) if L>0
gen ln_Y = log(rY) if L>0
gen ln_va_L = log(rVA/L) if rVA>0&L>0
gen ww = (w_wh*L_wh + w_bl*L_bl)/L if L>0
gen ln_ww = log(ww) if L>0
gen ln_bl_L = log(L_bl/L) if L>0&L_bl~=.
gen bl_L = L_bl/L  if L>0&L_bl~=.
gen ln_wh_L = log(L_wh/L) if L>0&L_wh~=.
gen wh_L = L_wh/L  if L>0&L_wh~=.

gen ln_kk_L = log(kk/L) if L>0&kk>0
tab year, gen(dyear)

gen d_x1 = (1-d_m)*d_e if L>0
gen d_m1 = (1-d_e)*d_m if L>0
gen d_xm = d_e*d_m if L>0

* Pooled OLS (1990-1996)
reg ln_L d_x1 d_m1 d_xm dyear2-dyear7 if L>0, robust 
reg ln_Y d_x1 d_m1 d_xm dyear2-dyear7 if L>0, robust 
reg ln_va_L d_x1 d_m1 d_xm ln_L dyear2-dyear7 if L>0, robust 
reg ln_ww d_x1 d_m1 d_xm ln_L dyear2-dyear7 if L>0, robust 
reg ln_wh_L d_x1 d_m1 d_xm ln_L dyear2-dyear7 if L>0, robust 
reg ln_kk_L d_x1 d_m1 d_xm ln_L dyear2-dyear7 if L>0, robust  

* Fixed Effects (1990-1996)
xtreg ln_L d_x1 d_m1 d_xm dyear2-dyear7 if L>0, fe i(padron) robust  
xtreg ln_Y d_x1 d_m1 d_xm ln_L dyear2-dyear7 if L>0, fe i(padron) robust  
xtreg ln_va_L d_x1 d_m1 d_xm ln_L dyear2-dyear7 if L>0, fe i(padron) robust  
xtreg ln_ww d_x1 d_m1 d_xm ln_L dyear2-dyear7 if L>0, fe i(padron) robust   
xtreg ln_wh_L d_x1 d_m1 d_xm ln_L dyear2-dyear7 if L>0, fe i(padron) robust  
xtreg ln_kk_L d_x1 d_m1 d_xm ln_L dyear2-dyear7 if L>0, fe i(padron) robust  

gen d_L = 0 
replace d_L = 1 if L>0
sort padron
by padron: egen sum_d_L = sum(d_L)
* Note Stata reports the `wrong' number of observations for fixed effects: the following number should be subtracted
summarize sum_d_L if year==90&sum_d_L==1





replace entry = 0
replace exit = 0

keep entry exit rY L rexport rM rMi padron year
reshape wide entry exit rY L rexport rM rMi, i(padron) j(year)

replace entry91 =1 if L91>0&L90==0
replace entry92 =1 if L92>0&L91==0&L90==0
replace entry93 =1 if L93>0&L92==0&L91==0&L90==0
replace entry94 =1 if L94>0&L93==0&L92==0&L91==0&L90==0
replace entry95 =1 if L95>0&L94==0&L93==0&L92==0&L91==0&L90==0
replace entry96 =1 if L96>0&L95==0&L94==0&L93==0&L92==0&L91==0&L90==0

replace exit91 = 1 if L90>0&L91==0
replace L92 = 0 if exit91==1
replace exit92 = 1 if L91>0&L92==0
replace L93 = 0 if exit92==1|exit91==1
replace exit93 = 1 if L92>0&L93==0
replace L94 = 0 if exit93==1|exit92==1|exit91==1
replace exit94 = 1 if L93>0&L94==0 
replace L95 = 0 if exit94==1|exit93==1|exit92==1|exit91==1
replace exit95 = 1 if L94>0&L95==0  
replace L96 = 0 if exit95==1|exit94==1|exit93==1|exit92==1|exit91==1
replace exit96 = 1 if L95>0&L96==0  

reshape long entry exit rY L rexport rM rMi, i(padron) j(year) 
 
gen d_L = 0 
replace d_L = 1 if L>0

sort year
by year: egen sum_entry = sum(entry)
by year: egen sum_exit = sum(exit)
by year: egen nplant = sum(d_L)
gen entry_nplant = sum_entry/nplant
gen exit_nplant = sum_exit/nplant

by year: summarize nplant sum_entry sum_exit 
by year: summarize entry_nplant exit_nplant 




****** (6) Industry 321 ********

* Creat Summary Statistics
clear
set memory 500m
set matsize 800
use "C:\Documents and Settings\hkasahar\My Documents\work\project\exportimport\data\DATA_KL_EXIM5.dta"
drop if year<90 
replace L = 0 if L ==.
sort padron
by padron: egen min_M = min(M) if L>0
drop if min_M == 0
by padron: egen min_Y = min(Y) if L>0
drop if min_Y == 0
by padron: egen sum_L = sum(L)
drop if sum_L == 0

sort padron
by padron: egen sumreentry = sum(reentry)
keep if mcu==321 
 

replace L_bl = 0 if L_bl==.
replace L_wh = 0 if L_wh==.
replace w_bl = 0 if w_bl==.
replace w_wh = 0 if w_wh==.
gen kk = rkm+rkt+rkb if L>0
gen ln_L = log(L) if L>0
gen ln_Y = log(rY) if L>0
gen ln_va_L = log(rVA/L) if rVA>0&L>0
gen ww = (w_wh*L_wh + w_bl*L_bl)/L if L>0
gen ln_ww = log(ww) if L>0
gen ln_bl_L = log(L_bl/L) if L>0&L_bl~=.
gen bl_L = L_bl/L  if L>0&L_bl~=.
gen ln_wh_L = log(L_wh/L) if L>0&L_wh~=.
gen wh_L = L_wh/L  if L>0&L_wh~=.

gen ln_kk_L = log(kk/L) if L>0&kk>0
tab year, gen(dyear)

gen d_x1 = (1-d_m)*d_e if L>0
gen d_m1 = (1-d_e)*d_m if L>0
gen d_xm = d_e*d_m if L>0

* Pooled OLS (1990-1996)
reg ln_L d_x1 d_m1 d_xm dyear2-dyear7 if L>0, robust 
reg ln_Y d_x1 d_m1 d_xm dyear2-dyear7 if L>0, robust 
reg ln_va_L d_x1 d_m1 d_xm ln_L dyear2-dyear7 if L>0, robust 
reg ln_ww d_x1 d_m1 d_xm ln_L dyear2-dyear7 if L>0, robust 
reg ln_wh_L d_x1 d_m1 d_xm ln_L dyear2-dyear7 if L>0, robust 
reg ln_kk_L d_x1 d_m1 d_xm ln_L dyear2-dyear7 if L>0, robust  

* Fixed Effects (1990-1996)
xtreg ln_L d_x1 d_m1 d_xm dyear2-dyear7 if L>0, fe i(padron) robust  
xtreg ln_Y d_x1 d_m1 d_xm ln_L dyear2-dyear7 if L>0, fe i(padron) robust  
xtreg ln_va_L d_x1 d_m1 d_xm ln_L dyear2-dyear7 if L>0, fe i(padron) robust  
xtreg ln_ww d_x1 d_m1 d_xm ln_L dyear2-dyear7 if L>0, fe i(padron) robust   
xtreg ln_wh_L d_x1 d_m1 d_xm ln_L dyear2-dyear7 if L>0, fe i(padron) robust  
xtreg ln_kk_L d_x1 d_m1 d_xm ln_L dyear2-dyear7 if L>0, fe i(padron) robust  

gen d_L = 0 
replace d_L = 1 if L>0
sort padron
by padron: egen sum_d_L = sum(d_L)
* Note Stata reports the `wrong' number of observations for fixed effects: the following number should be subtracted
summarize sum_d_L if year==90&sum_d_L==1



replace entry = 0
replace exit = 0

keep entry exit rY L rexport rM rMi padron year
reshape wide entry exit rY L rexport rM rMi, i(padron) j(year)

replace entry91 =1 if L91>0&L90==0
replace entry92 =1 if L92>0&L91==0&L90==0
replace entry93 =1 if L93>0&L92==0&L91==0&L90==0
replace entry94 =1 if L94>0&L93==0&L92==0&L91==0&L90==0
replace entry95 =1 if L95>0&L94==0&L93==0&L92==0&L91==0&L90==0
replace entry96 =1 if L96>0&L95==0&L94==0&L93==0&L92==0&L91==0&L90==0

replace exit91 = 1 if L90>0&L91==0
replace L92 = 0 if exit91==1
replace exit92 = 1 if L91>0&L92==0
replace L93 = 0 if exit92==1|exit91==1
replace exit93 = 1 if L92>0&L93==0
replace L94 = 0 if exit93==1|exit92==1|exit91==1
replace exit94 = 1 if L93>0&L94==0 
replace L95 = 0 if exit94==1|exit93==1|exit92==1|exit91==1
replace exit95 = 1 if L94>0&L95==0  
replace L96 = 0 if exit95==1|exit94==1|exit93==1|exit92==1|exit91==1
replace exit96 = 1 if L95>0&L96==0  

reshape long entry exit rY L rexport rM rMi, i(padron) j(year) 
 
gen d_L = 0 
replace d_L = 1 if L>0

sort year
by year: egen sum_entry = sum(entry)
by year: egen sum_exit = sum(exit)
by year: egen nplant = sum(d_L)
gen entry_nplant = sum_entry/nplant
gen exit_nplant = sum_exit/nplant

by year: summarize nplant sum_entry sum_exit 
by year: summarize entry_nplant exit_nplant 




****** (4) Industry 331 ********

* Creat Summary Statistics
clear
set memory 500m
set matsize 800
use "C:\Documents and Settings\hkasahar\My Documents\work\project\exportimport\data\DATA_KL_EXIM5.dta"
drop if year<90 
replace L = 0 if L ==.
sort padron
by padron: egen min_M = min(M) if L>0
drop if min_M == 0
by padron: egen min_Y = min(Y) if L>0
drop if min_Y == 0
by padron: egen sum_L = sum(L)
drop if sum_L == 0

sort padron
by padron: egen sumreentry = sum(reentry)
keep if mcu==331 


replace L_bl = 0 if L_bl==.
replace L_wh = 0 if L_wh==.
replace w_bl = 0 if w_bl==.
replace w_wh = 0 if w_wh==.
gen kk = rkm+rkt+rkb if L>0
gen ln_L = log(L) if L>0
gen ln_Y = log(rY) if L>0
gen ln_va_L = log(rVA/L) if rVA>0&L>0
gen ww = (w_wh*L_wh + w_bl*L_bl)/L if L>0
gen ln_ww = log(ww) if L>0
gen ln_bl_L = log(L_bl/L) if L>0&L_bl~=.
gen bl_L = L_bl/L  if L>0&L_bl~=.
gen ln_wh_L = log(L_wh/L) if L>0&L_wh~=.
gen wh_L = L_wh/L  if L>0&L_wh~=.

gen ln_kk_L = log(kk/L) if L>0&kk>0
tab year, gen(dyear)

gen d_x1 = (1-d_m)*d_e if L>0
gen d_m1 = (1-d_e)*d_m if L>0
gen d_xm = d_e*d_m if L>0

* Pooled OLS (1990-1996)
reg ln_L d_x1 d_m1 d_xm dyear2-dyear7 if L>0, robust 
reg ln_Y d_x1 d_m1 d_xm dyear2-dyear7 if L>0, robust 
reg ln_va_L d_x1 d_m1 d_xm ln_L dyear2-dyear7 if L>0, robust 
reg ln_ww d_x1 d_m1 d_xm ln_L dyear2-dyear7 if L>0, robust 
reg ln_wh_L d_x1 d_m1 d_xm ln_L dyear2-dyear7 if L>0, robust 
reg ln_kk_L d_x1 d_m1 d_xm ln_L dyear2-dyear7 if L>0, robust  

* Fixed Effects (1990-1996)
xtreg ln_L d_x1 d_m1 d_xm dyear2-dyear7 if L>0, fe i(padron) robust  
xtreg ln_Y d_x1 d_m1 d_xm ln_L dyear2-dyear7 if L>0, fe i(padron) robust  
xtreg ln_va_L d_x1 d_m1 d_xm ln_L dyear2-dyear7 if L>0, fe i(padron) robust  
xtreg ln_ww d_x1 d_m1 d_xm ln_L dyear2-dyear7 if L>0, fe i(padron) robust   
xtreg ln_wh_L d_x1 d_m1 d_xm ln_L dyear2-dyear7 if L>0, fe i(padron) robust  
xtreg ln_kk_L d_x1 d_m1 d_xm ln_L dyear2-dyear7 if L>0, fe i(padron) robust  

gen d_L = 0 
replace d_L = 1 if L>0
sort padron
by padron: egen sum_d_L = sum(d_L)
* Note Stata reports the `wrong' number of observations for fixed effects: the following number should be subtracted
summarize sum_d_L if year==90&sum_d_L==1





 
replace entry = 0
replace exit = 0

keep entry exit rY L rexport rM rMi padron year
reshape wide entry exit rY L rexport rM rMi, i(padron) j(year)

replace entry91 =1 if L91>0&L90==0
replace entry92 =1 if L92>0&L91==0&L90==0
replace entry93 =1 if L93>0&L92==0&L91==0&L90==0
replace entry94 =1 if L94>0&L93==0&L92==0&L91==0&L90==0
replace entry95 =1 if L95>0&L94==0&L93==0&L92==0&L91==0&L90==0
replace entry96 =1 if L96>0&L95==0&L94==0&L93==0&L92==0&L91==0&L90==0

replace exit91 = 1 if L90>0&L91==0
replace L92 = 0 if exit91==1
replace exit92 = 1 if L91>0&L92==0
replace L93 = 0 if exit92==1|exit91==1
replace exit93 = 1 if L92>0&L93==0
replace L94 = 0 if exit93==1|exit92==1|exit91==1
replace exit94 = 1 if L93>0&L94==0 
replace L95 = 0 if exit94==1|exit93==1|exit92==1|exit91==1
replace exit95 = 1 if L94>0&L95==0  
replace L96 = 0 if exit95==1|exit94==1|exit93==1|exit92==1|exit91==1
replace exit96 = 1 if L95>0&L96==0  

reshape long entry exit rY L rexport rM rMi, i(padron) j(year) 
 
gen d_L = 0 
replace d_L = 1 if L>0

sort year
by year: egen sum_entry = sum(entry)
by year: egen sum_exit = sum(exit)
by year: egen nplant = sum(d_L)
gen entry_nplant = sum_entry/nplant
gen exit_nplant = sum_exit/nplant

by year: summarize nplant sum_entry sum_exit 
by year: summarize entry_nplant exit_nplant 


 

***********************
* Industry-level Distribution of Export and Import
***********************

****** (1) Industry 3220 ********

* Creat Summary Statistics
clear
set memory 500m
set matsize 800
use "C:\Documents and Settings\hkasahar\My Documents\work\project\exportimport\data\DATA_KL_EXIM5.dta"
drop if year<90
sort padron
by padron: egen mciiu = mode(ciiu), minmode
by padron: egen sdciiu = sd(ciiu) 
by padron: egen sumreentry = sum(reentry)
keep if mciiu==3220 

*drop if dplnt==1
replace L = 0 if L ==.
sort padron
by padron: egen min_M = min(M) if L>0
drop if min_M == 0
by padron: egen min_Y = min(Y) if L>0
drop if min_Y == 0
by padron: egen sum_L = sum(L)
drop if sum_L == 0

gen d_x1 = (1-d_m)*d_e if L>0
gen d_m1 = (1-d_e)*d_m if L>0
gen d_xm = d_e*d_m if L>0

* EXPORTS

sort year

by year: egen pct_ex_01 = pctile(rexport) if rexport>0&L>0, p(99) 
by year: egen pct_ex_05 = pctile(rexport) if rexport>0&L>0, p(95)
by year: egen pct_ex_10 = pctile(rexport) if rexport>0&L>0, p(90)

by year: egen ex_01 = total(rexport) if rexport>=pct_ex_01&L>0
by year: egen ex_05 = total(rexport) if rexport>=pct_ex_05&L>0
by year: egen ex_10 = total(rexport) if rexport>=pct_ex_10&L>0  
by year: egen ex_ag = total(rexport) if rexport>0&L>0
gen ex_01_ex_ag = ex_01/ex_ag if rexport>=pct_ex_01&L>0
gen ex_05_ex_ag = ex_05/ex_ag if rexport>=pct_ex_05&L>0
gen ex_10_ex_ag = ex_10/ex_ag if rexport>=pct_ex_10&L>0

summarize ex_01_ex_ag ex_05_ex_ag ex_10_ex_ag if rexport>0&L>0 
 
* IMPORTS

sort year

by year: egen pct_im_01 = pctile(rMi) if rMi>0&L>0, p(99) 
by year: egen pct_im_05 = pctile(rMi) if rMi>0&L>0, p(95)
by year: egen pct_im_10 = pctile(rMi) if rMi>0&L>0, p(90)

by year: egen im_01 = total(rMi) if rMi>=pct_im_01&L>0
by year: egen im_05 = total(rMi) if rMi>=pct_im_05&L>0
by year: egen im_10 = total(rMi) if rMi>=pct_im_10&L>0  
by year: egen im_ag = total(rMi) if rMi>0&L>0
gen im_01_im_ag = im_01/im_ag if rMi>=pct_im_01&L>0
gen im_05_im_ag = im_05/im_ag if rMi>=pct_im_05&L>0
gen im_10_im_ag = im_10/im_ag if rMi>=pct_im_10&L>0

summarize im_01_im_ag im_05_im_ag im_10_im_ag if rMi>0&L>0 
 

 

****** (2) Industry 3560 ********

* Creat Summary Statistics
clear
set memory 500m
set matsize 800
use "C:\Documents and Settings\hkasahar\My Documents\work\project\exportimport\data\DATA_KL_EXIM5.dta"
drop if year<90
sort padron
by padron: egen mciiu = mode(ciiu), minmode
by padron: egen sdciiu = sd(ciiu) 
by padron: egen sumreentry = sum(reentry)
keep if mciiu==3560
keep if sdciiu==0|sdciiu==.
keep if sumreentry==0

*drop if dplnt==1
replace L = 0 if L ==.
sort padron
by padron: egen min_M = min(M) if L>0
drop if min_M == 0
by padron: egen min_Y = min(Y) if L>0
drop if min_Y == 0
by padron: egen sum_L = sum(L)
drop if sum_L == 0

gen d_x1 = (1-d_m)*d_e if L>0
gen d_m1 = (1-d_e)*d_m if L>0
gen d_xm = d_e*d_m if L>0

* EXPORTS

sort year

by year: egen pct_ex_01 = pctile(rexport) if rexport>0&L>0, p(99) 
by year: egen pct_ex_05 = pctile(rexport) if rexport>0&L>0, p(95)
by year: egen pct_ex_10 = pctile(rexport) if rexport>0&L>0, p(90)

by year: egen ex_01 = total(rexport) if rexport>=pct_ex_01&L>0
by year: egen ex_05 = total(rexport) if rexport>=pct_ex_05&L>0
by year: egen ex_10 = total(rexport) if rexport>=pct_ex_10&L>0  
by year: egen ex_ag = total(rexport) if rexport>0&L>0
gen ex_01_ex_ag = ex_01/ex_ag if rexport>=pct_ex_01&L>0
gen ex_05_ex_ag = ex_05/ex_ag if rexport>=pct_ex_05&L>0
gen ex_10_ex_ag = ex_10/ex_ag if rexport>=pct_ex_10&L>0

summarize ex_01_ex_ag ex_05_ex_ag ex_10_ex_ag if rexport>0&L>0 
 
* IMPORTS

sort year

by year: egen pct_im_01 = pctile(rMi) if rMi>0&L>0, p(99) 
by year: egen pct_im_05 = pctile(rMi) if rMi>0&L>0, p(95)
by year: egen pct_im_10 = pctile(rMi) if rMi>0&L>0, p(90)

by year: egen im_01 = total(rMi) if rMi>=pct_im_01&L>0
by year: egen im_05 = total(rMi) if rMi>=pct_im_05&L>0
by year: egen im_10 = total(rMi) if rMi>=pct_im_10&L>0  
by year: egen im_ag = total(rMi) if rMi>0&L>0
gen im_01_im_ag = im_01/im_ag if rMi>=pct_im_01&L>0
gen im_05_im_ag = im_05/im_ag if rMi>=pct_im_05&L>0
gen im_10_im_ag = im_10/im_ag if rMi>=pct_im_10&L>0

summarize im_01_im_ag im_05_im_ag im_10_im_ag if rMi>0&L>0 
 

****** (3) Industry 3813 ********

* Creat Summary Statistics
clear
set memory 500m
set matsize 800
use "C:\Documents and Settings\hkasahar\My Documents\work\project\exportimport\data\DATA_KL_EXIM5.dta"
drop if year<90


replace ciiu = 3813 if ciiu==3814|ciiu==3815

sort padron
by padron: egen mciiu = mode(ciiu), minmode
by padron: egen sdciiu = sd(ciiu) 
by padron: egen sumreentry = sum(reentry)
keep if mciiu==3813
*keep if sdciiu==0|sdciiu==.
*keep if sumreentry==0

*drop if dplnt==1
replace L = 0 if L ==.
sort padron
by padron: egen min_M = min(M) if L>0
drop if min_M == 0
by padron: egen min_Y = min(Y) if L>0
drop if min_Y == 0
by padron: egen sum_L = sum(L)
drop if sum_L == 0

gen d_x1 = (1-d_m)*d_e if L>0
gen d_m1 = (1-d_e)*d_m if L>0
gen d_xm = d_e*d_m if L>0

* EXPORTS

sort year

by year: egen pct_ex_01 = pctile(rexport) if rexport>0&L>0, p(99) 
by year: egen pct_ex_05 = pctile(rexport) if rexport>0&L>0, p(95)
by year: egen pct_ex_10 = pctile(rexport) if rexport>0&L>0, p(90)

by year: egen ex_01 = total(rexport) if rexport>=pct_ex_01&L>0
by year: egen ex_05 = total(rexport) if rexport>=pct_ex_05&L>0
by year: egen ex_10 = total(rexport) if rexport>=pct_ex_10&L>0  
by year: egen ex_ag = total(rexport) if rexport>0&L>0
gen ex_01_ex_ag = ex_01/ex_ag if rexport>=pct_ex_01&L>0
gen ex_05_ex_ag = ex_05/ex_ag if rexport>=pct_ex_05&L>0
gen ex_10_ex_ag = ex_10/ex_ag if rexport>=pct_ex_10&L>0

summarize ex_01_ex_ag ex_05_ex_ag ex_10_ex_ag if rexport>0&L>0 
 
* IMPORTS

sort year

by year: egen pct_im_01 = pctile(rMi) if rMi>0&L>0, p(99) 
by year: egen pct_im_05 = pctile(rMi) if rMi>0&L>0, p(95)
by year: egen pct_im_10 = pctile(rMi) if rMi>0&L>0, p(90)

by year: egen im_01 = total(rMi) if rMi>=pct_im_01&L>0
by year: egen im_05 = total(rMi) if rMi>=pct_im_05&L>0
by year: egen im_10 = total(rMi) if rMi>=pct_im_10&L>0  
by year: egen im_ag = total(rMi) if rMi>0&L>0
gen im_01_im_ag = im_01/im_ag if rMi>=pct_im_01&L>0
gen im_05_im_ag = im_05/im_ag if rMi>=pct_im_05&L>0
gen im_10_im_ag = im_10/im_ag if rMi>=pct_im_10&L>0

summarize im_01_im_ag im_05_im_ag im_10_im_ag if rMi>0&L>0 
 


 
****** (4) Industry 381 ********

* Creat Summary Statistics
clear
set memory 500m
set matsize 800
use "C:\Documents and Settings\hkasahar\My Documents\work\project\exportimport\data\DATA_KL_EXIM5.dta"
drop if year<90 
replace L = 0 if L ==.
sort padron
by padron: egen min_M = min(M) if L>0
drop if min_M == 0
by padron: egen min_Y = min(Y) if L>0
drop if min_Y == 0
by padron: egen sum_L = sum(L)
drop if sum_L == 0

sort padron
by padron: egen sumreentry = sum(reentry)
keep if mcu==381 


* EXPORTS

sort year

by year: egen pct_ex_01 = pctile(rexport) if rexport>0&L>0, p(99) 
by year: egen pct_ex_05 = pctile(rexport) if rexport>0&L>0, p(95)
by year: egen pct_ex_10 = pctile(rexport) if rexport>0&L>0, p(90)

by year: egen ex_01 = total(rexport) if rexport>=pct_ex_01&L>0
by year: egen ex_05 = total(rexport) if rexport>=pct_ex_05&L>0
by year: egen ex_10 = total(rexport) if rexport>=pct_ex_10&L>0  
by year: egen ex_ag = total(rexport) if rexport>0&L>0
gen ex_01_ex_ag = ex_01/ex_ag if rexport>=pct_ex_01&L>0
gen ex_05_ex_ag = ex_05/ex_ag if rexport>=pct_ex_05&L>0
gen ex_10_ex_ag = ex_10/ex_ag if rexport>=pct_ex_10&L>0


 
* IMPORTS

sort year

by year: egen pct_im_01 = pctile(rMi) if rMi>0&L>0, p(99) 
by year: egen pct_im_05 = pctile(rMi) if rMi>0&L>0, p(95)
by year: egen pct_im_10 = pctile(rMi) if rMi>0&L>0, p(90)

by year: egen im_01 = total(rMi) if rMi>=pct_im_01&L>0
by year: egen im_05 = total(rMi) if rMi>=pct_im_05&L>0
by year: egen im_10 = total(rMi) if rMi>=pct_im_10&L>0  
by year: egen im_ag = total(rMi) if rMi>0&L>0
gen im_01_im_ag = im_01/im_ag if rMi>=pct_im_01&L>0
gen im_05_im_ag = im_05/im_ag if rMi>=pct_im_05&L>0
gen im_10_im_ag = im_10/im_ag if rMi>=pct_im_10&L>0

summarize ex_01_ex_ag ex_05_ex_ag ex_10_ex_ag if rexport>0&L>0 
summarize im_01_im_ag im_05_im_ag im_10_im_ag if rMi>0&L>0 

 

****** (5) Industry 311 ********

* Creat Summary Statistics
clear
set memory 500m
set matsize 800
use "C:\Documents and Settings\hkasahar\My Documents\work\project\exportimport\data\DATA_KL_EXIM5.dta"
drop if year<90 
replace L = 0 if L ==.
sort padron
by padron: egen min_M = min(M) if L>0
drop if min_M == 0
by padron: egen min_Y = min(Y) if L>0
drop if min_Y == 0
by padron: egen sum_L = sum(L)
drop if sum_L == 0

sort padron
by padron: egen mciiu = mode(ciiu), minmode
by padron: egen sumreentry = sum(reentry)
keep if mciiu==3111|mciiu==3112|mciiu==3113|mciiu==3114|mciiu==3115|mciiu==3116|mciiu==3118|mciiu==3119
 

* EXPORTS

sort year

by year: egen pct_ex_01 = pctile(rexport) if rexport>0&L>0, p(99) 
by year: egen pct_ex_05 = pctile(rexport) if rexport>0&L>0, p(95)
by year: egen pct_ex_10 = pctile(rexport) if rexport>0&L>0, p(90)

by year: egen ex_01 = total(rexport) if rexport>=pct_ex_01&L>0
by year: egen ex_05 = total(rexport) if rexport>=pct_ex_05&L>0
by year: egen ex_10 = total(rexport) if rexport>=pct_ex_10&L>0  
by year: egen ex_ag = total(rexport) if rexport>0&L>0
gen ex_01_ex_ag = ex_01/ex_ag if rexport>=pct_ex_01&L>0
gen ex_05_ex_ag = ex_05/ex_ag if rexport>=pct_ex_05&L>0
gen ex_10_ex_ag = ex_10/ex_ag if rexport>=pct_ex_10&L>0


 
* IMPORTS

sort year

by year: egen pct_im_01 = pctile(rMi) if rMi>0&L>0, p(99) 
by year: egen pct_im_05 = pctile(rMi) if rMi>0&L>0, p(95)
by year: egen pct_im_10 = pctile(rMi) if rMi>0&L>0, p(90)

by year: egen im_01 = total(rMi) if rMi>=pct_im_01&L>0
by year: egen im_05 = total(rMi) if rMi>=pct_im_05&L>0
by year: egen im_10 = total(rMi) if rMi>=pct_im_10&L>0  
by year: egen im_ag = total(rMi) if rMi>0&L>0
gen im_01_im_ag = im_01/im_ag if rMi>=pct_im_01&L>0
gen im_05_im_ag = im_05/im_ag if rMi>=pct_im_05&L>0
gen im_10_im_ag = im_10/im_ag if rMi>=pct_im_10&L>0

summarize ex_01_ex_ag ex_05_ex_ag ex_10_ex_ag if rexport>0&L>0 
summarize im_01_im_ag im_05_im_ag im_10_im_ag if rMi>0&L>0 



****** (6) Industry 321 ********

* Creat Summary Statistics
clear
set memory 500m
set matsize 800
use "C:\Documents and Settings\hkasahar\My Documents\work\project\exportimport\data\DATA_KL_EXIM5.dta"
drop if year<90 
replace L = 0 if L ==.
sort padron
by padron: egen min_M = min(M) if L>0
drop if min_M == 0
by padron: egen min_Y = min(Y) if L>0
drop if min_Y == 0
by padron: egen sum_L = sum(L)
drop if sum_L == 0

sort padron
by padron: egen sumreentry = sum(reentry)
keep if mcu==321 
 

* EXPORTS

sort year

by year: egen pct_ex_01 = pctile(rexport) if rexport>0&L>0, p(99) 
by year: egen pct_ex_05 = pctile(rexport) if rexport>0&L>0, p(95)
by year: egen pct_ex_10 = pctile(rexport) if rexport>0&L>0, p(90)

by year: egen ex_01 = total(rexport) if rexport>=pct_ex_01&L>0
by year: egen ex_05 = total(rexport) if rexport>=pct_ex_05&L>0
by year: egen ex_10 = total(rexport) if rexport>=pct_ex_10&L>0  
by year: egen ex_ag = total(rexport) if rexport>0&L>0
gen ex_01_ex_ag = ex_01/ex_ag if rexport>=pct_ex_01&L>0
gen ex_05_ex_ag = ex_05/ex_ag if rexport>=pct_ex_05&L>0
gen ex_10_ex_ag = ex_10/ex_ag if rexport>=pct_ex_10&L>0
 
* IMPORTS

sort year

by year: egen pct_im_01 = pctile(rMi) if rMi>0&L>0, p(99) 
by year: egen pct_im_05 = pctile(rMi) if rMi>0&L>0, p(95)
by year: egen pct_im_10 = pctile(rMi) if rMi>0&L>0, p(90)

by year: egen im_01 = total(rMi) if rMi>=pct_im_01&L>0
by year: egen im_05 = total(rMi) if rMi>=pct_im_05&L>0
by year: egen im_10 = total(rMi) if rMi>=pct_im_10&L>0  
by year: egen im_ag = total(rMi) if rMi>0&L>0
gen im_01_im_ag = im_01/im_ag if rMi>=pct_im_01&L>0
gen im_05_im_ag = im_05/im_ag if rMi>=pct_im_05&L>0
gen im_10_im_ag = im_10/im_ag if rMi>=pct_im_10&L>0


summarize ex_01_ex_ag ex_05_ex_ag ex_10_ex_ag if rexport>0&L>0 
summarize im_01_im_ag im_05_im_ag im_10_im_ag if rMi>0&L>0 



****** (4) Industry 331 ********

* Creat Summary Statistics
clear
set memory 500m
set matsize 800
use "C:\Documents and Settings\hkasahar\My Documents\work\project\exportimport\data\DATA_KL_EXIM5.dta"
drop if year<90 
replace L = 0 if L ==.
sort padron
by padron: egen min_M = min(M) if L>0
drop if min_M == 0
by padron: egen min_Y = min(Y) if L>0
drop if min_Y == 0
by padron: egen sum_L = sum(L)
drop if sum_L == 0

sort padron
by padron: egen sumreentry = sum(reentry)
keep if mcu==331 


* EXPORTS

sort year

by year: egen pct_ex_01 = pctile(rexport) if rexport>0&L>0, p(99) 
by year: egen pct_ex_05 = pctile(rexport) if rexport>0&L>0, p(95)
by year: egen pct_ex_10 = pctile(rexport) if rexport>0&L>0, p(90)

by year: egen ex_01 = total(rexport) if rexport>=pct_ex_01&L>0
by year: egen ex_05 = total(rexport) if rexport>=pct_ex_05&L>0
by year: egen ex_10 = total(rexport) if rexport>=pct_ex_10&L>0  
by year: egen ex_ag = total(rexport) if rexport>0&L>0
gen ex_01_ex_ag = ex_01/ex_ag if rexport>=pct_ex_01&L>0
gen ex_05_ex_ag = ex_05/ex_ag if rexport>=pct_ex_05&L>0
gen ex_10_ex_ag = ex_10/ex_ag if rexport>=pct_ex_10&L>0

 
 
* IMPORTS

sort year

by year: egen pct_im_01 = pctile(rMi) if rMi>0&L>0, p(99) 
by year: egen pct_im_05 = pctile(rMi) if rMi>0&L>0, p(95)
by year: egen pct_im_10 = pctile(rMi) if rMi>0&L>0, p(90)

by year: egen im_01 = total(rMi) if rMi>=pct_im_01&L>0
by year: egen im_05 = total(rMi) if rMi>=pct_im_05&L>0
by year: egen im_10 = total(rMi) if rMi>=pct_im_10&L>0  
by year: egen im_ag = total(rMi) if rMi>0&L>0
gen im_01_im_ag = im_01/im_ag if rMi>=pct_im_01&L>0
gen im_05_im_ag = im_05/im_ag if rMi>=pct_im_05&L>0
gen im_10_im_ag = im_10/im_ag if rMi>=pct_im_10&L>0

summarize ex_01_ex_ag ex_05_ex_ag ex_10_ex_ag if rexport>0&L>0 
summarize im_01_im_ag im_05_im_ag im_10_im_ag if rMi>0&L>0 




***********************
* Industry-level Switchers etc.  
***********************

 


****** (1) Industry 3220 ********

clear
set memory 500m
set matsize 800
use "C:\Documents and Settings\hkasahar\My Documents\work\project\exportimport\data\DATA_KL_EXIM5.dta"

drop if year<90
replace L = 0 if L ==.
sort padron
by padron: egen min_M = min(M) if L>0
drop if min_M == 0
by padron: egen min_Y = min(Y) if L>0
drop if min_Y == 0
by padron: egen sum_L = sum(L)
drop if sum_L == 0

sort padron
by padron: egen mciiu = mode(ciiu), minmode
by padron: egen sdciiu = sd(ciiu) 
by padron: egen sumreentry = sum(reentry)
keep if mciiu==3220
 

gen d_x1 = (1-d_m)*d_e if L~=0
gen d_m1 = (1-d_e)*d_m if L~=0
gen d_xm = d_e*d_m if L~=0
gen dxm = 1 if L~=0
replace dxm = 2 if L~=0&d_x1==1
replace dxm = 3 if L~=0&d_m1==1
replace dxm = 4 if L~=0&d_xm == 1

gen switch = 0
gen switch_e = 0
gen switch_m = 0

* define ``switch"
keep switch switch_e switch_m rk dxm d_m d_e export rexport lnkm lnkt lnkb rent fta L_wh L_bl w_wh w_bl rVA rY rM rMi rMd rE Y M Mi Md E rim rit rib im it rkm rkt rkb ciiu cu  exit entry defmach_f defcons_f L padron year
reshape wide switch switch_e switch_m rk dxm d_m d_e export rexport lnkm lnkt lnkb rent fta L_wh L_bl w_wh w_bl rVA rY rM rMi rMd rE Y M Mi Md E rim rit rib im it rkm rkt rkb ciiu cu  exit entry defmach_f defcons_f L, i(padron) j(year)

replace switch91 = 1 if dxm90~=dxm91&L90~=0&L91~=0
replace switch92 = 1 if dxm91~=dxm92&L91~=0&L92~=0
replace switch93 = 1 if dxm92~=dxm93&L92~=0&L93~=0
replace switch94 = 1 if dxm93~=dxm94&L93~=0&L94~=0
replace switch95 = 1 if dxm94~=dxm95&L94~=0&L95~=0
replace switch96 = 1 if dxm95~=dxm96&L95~=0&L96~=0

replace switch_e91 = 1 if d_e90~=d_e91&L90~=0&L91~=0
replace switch_e92 = 1 if d_e91~=d_e92&L91~=0&L92~=0
replace switch_e93 = 1 if d_e92~=d_e93&L92~=0&L93~=0
replace switch_e94 = 1 if d_e93~=d_e94&L93~=0&L94~=0
replace switch_e95 = 1 if d_e94~=d_e95&L94~=0&L95~=0
replace switch_e96 = 1 if d_e95~=d_e96&L95~=0&L96~=0

replace switch_m91 = 1 if d_m90~=d_m91&L90~=0&L91~=0
replace switch_m92 = 1 if d_m91~=d_m92&L91~=0&L92~=0
replace switch_m93 = 1 if d_m92~=d_m93&L92~=0&L93~=0
replace switch_m94 = 1 if d_m93~=d_m94&L93~=0&L94~=0
replace switch_m95 = 1 if d_m94~=d_m95&L94~=0&L95~=0
replace switch_m96 = 1 if d_m95~=d_m96&L95~=0&L96~=0

gen dplant = 0
replace dplant = 1 if L90~=0&L91~=0&L92~=0&L93~=0&L94~=0&L95~=0&L96~=0

reshape long switch switch_e switch_m rk dxm d_m d_e export rexport lnkm lnkt lnkb rent fta L_wh L_bl w_wh w_bl rVA rY rM rMi rMd rE Y M Mi Md E rim rit rib im it rkm rkt rkb ciiu cu exit entry defmach_f defcons_f L, i(padron) j(year)

summarize switch if dplant==1, detail

by padron: egen sum_switch = total(switch) 
by padron: egen sum_switch_e = total(switch_e) 
by padron: egen sum_switch_m = total(switch_m)  
 
summarize sum_switch if year==90&sum_switch>=0
summarize sum_switch if year==90&sum_switch>=1
summarize sum_switch if year==90&sum_switch>=2 
summarize sum_switch if year==90&sum_switch>=3 

summarize sum_switch_e if year==90&sum_switch_e>=0
summarize sum_switch_e if year==90&sum_switch_e>=1
summarize sum_switch_e if year==90&sum_switch_e>=2 
summarize sum_switch_e if year==90&sum_switch_e>=3 
 
summarize sum_switch_m if year==90&sum_switch_m>=0
summarize sum_switch_m if year==90&sum_switch_m>=1
summarize sum_switch_m if year==90&sum_switch_m>=2 
summarize sum_switch_m if year==90&sum_switch_m>=3
 
 
 
****** (2) Industry 3560 ********

clear
set memory 500m
set matsize 800
use "C:\Documents and Settings\hkasahar\My Documents\work\project\exportimport\data\DATA_KL_EXIM5.dta"

drop if year<90
replace L = 0 if L ==.
sort padron
by padron: egen min_M = min(M) if L>0
drop if min_M == 0
by padron: egen min_Y = min(Y) if L>0
drop if min_Y == 0
by padron: egen sum_L = sum(L)
drop if sum_L == 0

sort padron
by padron: egen mciiu = mode(ciiu), minmode
by padron: egen sdciiu = sd(ciiu) 
by padron: egen sumreentry = sum(reentry)
keep if mciiu==3560 

gen d_x1 = (1-d_m)*d_e if L~=0
gen d_m1 = (1-d_e)*d_m if L~=0
gen d_xm = d_e*d_m if L~=0
gen dxm = 1 if L~=0
replace dxm = 2 if L~=0&d_x1==1
replace dxm = 3 if L~=0&d_m1==1
replace dxm = 4 if L~=0&d_xm == 1

gen switch = 0
gen switch_e = 0
gen switch_m = 0

* define ``switch"
keep switch switch_e switch_m rk dxm d_m d_e export rexport lnkm lnkt lnkb rent fta L_wh L_bl w_wh w_bl rVA rY rM rMi rMd rE Y M Mi Md E rim rit rib im it rkm rkt rkb ciiu cu  exit entry defmach_f defcons_f L padron year
reshape wide switch switch_e switch_m rk dxm d_m d_e export rexport lnkm lnkt lnkb rent fta L_wh L_bl w_wh w_bl rVA rY rM rMi rMd rE Y M Mi Md E rim rit rib im it rkm rkt rkb ciiu cu  exit entry defmach_f defcons_f L, i(padron) j(year)

replace switch91 = 1 if dxm90~=dxm91&L90~=0&L91~=0
replace switch92 = 1 if dxm91~=dxm92&L91~=0&L92~=0
replace switch93 = 1 if dxm92~=dxm93&L92~=0&L93~=0
replace switch94 = 1 if dxm93~=dxm94&L93~=0&L94~=0
replace switch95 = 1 if dxm94~=dxm95&L94~=0&L95~=0
replace switch96 = 1 if dxm95~=dxm96&L95~=0&L96~=0

replace switch_e91 = 1 if d_e90~=d_e91&L90~=0&L91~=0
replace switch_e92 = 1 if d_e91~=d_e92&L91~=0&L92~=0
replace switch_e93 = 1 if d_e92~=d_e93&L92~=0&L93~=0
replace switch_e94 = 1 if d_e93~=d_e94&L93~=0&L94~=0
replace switch_e95 = 1 if d_e94~=d_e95&L94~=0&L95~=0
replace switch_e96 = 1 if d_e95~=d_e96&L95~=0&L96~=0

replace switch_m91 = 1 if d_m90~=d_m91&L90~=0&L91~=0
replace switch_m92 = 1 if d_m91~=d_m92&L91~=0&L92~=0
replace switch_m93 = 1 if d_m92~=d_m93&L92~=0&L93~=0
replace switch_m94 = 1 if d_m93~=d_m94&L93~=0&L94~=0
replace switch_m95 = 1 if d_m94~=d_m95&L94~=0&L95~=0
replace switch_m96 = 1 if d_m95~=d_m96&L95~=0&L96~=0

gen dplant = 0
replace dplant = 1 if L90~=0&L91~=0&L92~=0&L93~=0&L94~=0&L95~=0&L96~=0

reshape long switch switch_e switch_m rk dxm d_m d_e export rexport lnkm lnkt lnkb rent fta L_wh L_bl w_wh w_bl rVA rY rM rMi rMd rE Y M Mi Md E rim rit rib im it rkm rkt rkb ciiu cu exit entry defmach_f defcons_f L, i(padron) j(year)

summarize switch if dplant==1, detail

by padron: egen sum_switch = total(switch) 
by padron: egen sum_switch_e = total(switch_e) 
by padron: egen sum_switch_m = total(switch_m)  
 
summarize sum_switch if year==90&sum_switch>=0
summarize sum_switch if year==90&sum_switch>=1
summarize sum_switch if year==90&sum_switch>=2 
summarize sum_switch if year==90&sum_switch>=3 

summarize sum_switch_e if year==90&sum_switch_e>=0
summarize sum_switch_e if year==90&sum_switch_e>=1
summarize sum_switch_e if year==90&sum_switch_e>=2 
summarize sum_switch_e if year==90&sum_switch_e>=3 
 
summarize sum_switch_m if year==90&sum_switch_m>=0
summarize sum_switch_m if year==90&sum_switch_m>=1
summarize sum_switch_m if year==90&sum_switch_m>=2 
summarize sum_switch_m if year==90&sum_switch_m>=3
 




 
****** (3) Industry 3813 ********

clear
set memory 500m
set matsize 800

use "C:\Documents and Settings\hkasahar\My Documents\work\project\exportimport\data\DATA_KL_EXIM5.dta"
drop if year<90
replace L = 0 if L ==.
sort padron
by padron: egen min_M = min(M) if L>0
drop if min_M == 0
by padron: egen min_Y = min(Y) if L>0
drop if min_Y == 0
by padron: egen sum_L = sum(L)
drop if sum_L == 0


replace ciiu = 3813 if ciiu==3814|ciiu==3815

sort padron
by padron: egen mciiu = mode(ciiu), minmode
by padron: egen sdciiu = sd(ciiu) 
by padron: egen sumreentry = sum(reentry)
keep if mciiu==3813
 

gen d_x1 = (1-d_m)*d_e if L~=0
gen d_m1 = (1-d_e)*d_m if L~=0
gen d_xm = d_e*d_m if L~=0
gen dxm = 1 if L~=0
replace dxm = 2 if L~=0&d_x1==1
replace dxm = 3 if L~=0&d_m1==1
replace dxm = 4 if L~=0&d_xm == 1

gen switch = 0
gen switch_e = 0
gen switch_m = 0

* define ``switch"
keep switch switch_e switch_m rk dxm d_m d_e export rexport lnkm lnkt lnkb rent fta L_wh L_bl w_wh w_bl rVA rY rM rMi rMd rE Y M Mi Md E rim rit rib im it rkm rkt rkb ciiu cu  exit entry defmach_f defcons_f L padron year
reshape wide switch switch_e switch_m rk dxm d_m d_e export rexport lnkm lnkt lnkb rent fta L_wh L_bl w_wh w_bl rVA rY rM rMi rMd rE Y M Mi Md E rim rit rib im it rkm rkt rkb ciiu cu  exit entry defmach_f defcons_f L, i(padron) j(year)

replace switch91 = 1 if dxm90~=dxm91&L90~=0&L91~=0
replace switch92 = 1 if dxm91~=dxm92&L91~=0&L92~=0
replace switch93 = 1 if dxm92~=dxm93&L92~=0&L93~=0
replace switch94 = 1 if dxm93~=dxm94&L93~=0&L94~=0
replace switch95 = 1 if dxm94~=dxm95&L94~=0&L95~=0
replace switch96 = 1 if dxm95~=dxm96&L95~=0&L96~=0

replace switch_e91 = 1 if d_e90~=d_e91&L90~=0&L91~=0
replace switch_e92 = 1 if d_e91~=d_e92&L91~=0&L92~=0
replace switch_e93 = 1 if d_e92~=d_e93&L92~=0&L93~=0
replace switch_e94 = 1 if d_e93~=d_e94&L93~=0&L94~=0
replace switch_e95 = 1 if d_e94~=d_e95&L94~=0&L95~=0
replace switch_e96 = 1 if d_e95~=d_e96&L95~=0&L96~=0

replace switch_m91 = 1 if d_m90~=d_m91&L90~=0&L91~=0
replace switch_m92 = 1 if d_m91~=d_m92&L91~=0&L92~=0
replace switch_m93 = 1 if d_m92~=d_m93&L92~=0&L93~=0
replace switch_m94 = 1 if d_m93~=d_m94&L93~=0&L94~=0
replace switch_m95 = 1 if d_m94~=d_m95&L94~=0&L95~=0
replace switch_m96 = 1 if d_m95~=d_m96&L95~=0&L96~=0

gen dplant = 0
replace dplant = 1 if L90~=0&L91~=0&L92~=0&L93~=0&L94~=0&L95~=0&L96~=0

reshape long switch switch_e switch_m rk dxm d_m d_e export rexport lnkm lnkt lnkb rent fta L_wh L_bl w_wh w_bl rVA rY rM rMi rMd rE Y M Mi Md E rim rit rib im it rkm rkt rkb ciiu cu exit entry defmach_f defcons_f L, i(padron) j(year)

summarize switch if dplant==1, detail

by padron: egen sum_switch = total(switch) 
by padron: egen sum_switch_e = total(switch_e) 
by padron: egen sum_switch_m = total(switch_m)  
 
summarize sum_switch if year==90&sum_switch>=0
summarize sum_switch if year==90&sum_switch>=1
summarize sum_switch if year==90&sum_switch>=2 
summarize sum_switch if year==90&sum_switch>=3 

summarize sum_switch_e if year==90&sum_switch_e>=0
summarize sum_switch_e if year==90&sum_switch_e>=1
summarize sum_switch_e if year==90&sum_switch_e>=2 
summarize sum_switch_e if year==90&sum_switch_e>=3 
 
summarize sum_switch_m if year==90&sum_switch_m>=0
summarize sum_switch_m if year==90&sum_switch_m>=1
summarize sum_switch_m if year==90&sum_switch_m>=2 
summarize sum_switch_m if year==90&sum_switch_m>=3
 


****** (4) Industry 381 ********

clear
set memory 500m
set matsize 800

use "C:\Documents and Settings\hkasahar\My Documents\work\project\exportimport\data\DATA_KL_EXIM5.dta"
drop if year<90
replace L = 0 if L ==.
sort padron
by padron: egen min_M = min(M) if L>0
drop if min_M == 0
by padron: egen min_Y = min(Y) if L>0
drop if min_Y == 0
by padron: egen sum_L = sum(L)
drop if sum_L == 0

sort padron
by padron: egen sumreentry = sum(reentry)
keep if mcu==381
 
gen d_x1 = (1-d_m)*d_e if L~=0
gen d_m1 = (1-d_e)*d_m if L~=0
gen d_xm = d_e*d_m if L~=0
gen dxm = 1 if L~=0
replace dxm = 2 if L~=0&d_x1==1
replace dxm = 3 if L~=0&d_m1==1
replace dxm = 4 if L~=0&d_xm == 1

gen switch = 0
gen switch_e = 0
gen switch_m = 0

* define ``switch"
keep switch switch_e switch_m rk dxm d_m d_e export rexport lnkm lnkt lnkb rent fta L_wh L_bl w_wh w_bl rVA rY rM rMi rMd rE Y M Mi Md E rim rit rib im it rkm rkt rkb ciiu cu  exit entry defmach_f defcons_f L padron year
reshape wide switch switch_e switch_m rk dxm d_m d_e export rexport lnkm lnkt lnkb rent fta L_wh L_bl w_wh w_bl rVA rY rM rMi rMd rE Y M Mi Md E rim rit rib im it rkm rkt rkb ciiu cu  exit entry defmach_f defcons_f L, i(padron) j(year)

replace switch91 = 1 if dxm90~=dxm91&L90~=0&L91~=0
replace switch92 = 1 if dxm91~=dxm92&L91~=0&L92~=0
replace switch93 = 1 if dxm92~=dxm93&L92~=0&L93~=0
replace switch94 = 1 if dxm93~=dxm94&L93~=0&L94~=0
replace switch95 = 1 if dxm94~=dxm95&L94~=0&L95~=0
replace switch96 = 1 if dxm95~=dxm96&L95~=0&L96~=0

replace switch_e91 = 1 if d_e90~=d_e91&L90~=0&L91~=0
replace switch_e92 = 1 if d_e91~=d_e92&L91~=0&L92~=0
replace switch_e93 = 1 if d_e92~=d_e93&L92~=0&L93~=0
replace switch_e94 = 1 if d_e93~=d_e94&L93~=0&L94~=0
replace switch_e95 = 1 if d_e94~=d_e95&L94~=0&L95~=0
replace switch_e96 = 1 if d_e95~=d_e96&L95~=0&L96~=0

replace switch_m91 = 1 if d_m90~=d_m91&L90~=0&L91~=0
replace switch_m92 = 1 if d_m91~=d_m92&L91~=0&L92~=0
replace switch_m93 = 1 if d_m92~=d_m93&L92~=0&L93~=0
replace switch_m94 = 1 if d_m93~=d_m94&L93~=0&L94~=0
replace switch_m95 = 1 if d_m94~=d_m95&L94~=0&L95~=0
replace switch_m96 = 1 if d_m95~=d_m96&L95~=0&L96~=0

gen dplant = 0
replace dplant = 1 if L90~=0&L91~=0&L92~=0&L93~=0&L94~=0&L95~=0&L96~=0

reshape long switch switch_e switch_m rk dxm d_m d_e export rexport lnkm lnkt lnkb rent fta L_wh L_bl w_wh w_bl rVA rY rM rMi rMd rE Y M Mi Md E rim rit rib im it rkm rkt rkb ciiu cu exit entry defmach_f defcons_f L, i(padron) j(year)

summarize switch if dplant==1, detail

by padron: egen sum_switch = total(switch) 
by padron: egen sum_switch_e = total(switch_e) 
by padron: egen sum_switch_m = total(switch_m)  
 
summarize sum_switch if year==90&sum_switch>=0
summarize sum_switch if year==90&sum_switch>=1
summarize sum_switch if year==90&sum_switch>=2 
summarize sum_switch if year==90&sum_switch>=3 

summarize sum_switch_e if year==90&sum_switch_e>=0
summarize sum_switch_e if year==90&sum_switch_e>=1
summarize sum_switch_e if year==90&sum_switch_e>=2 
summarize sum_switch_e if year==90&sum_switch_e>=3 
 
summarize sum_switch_m if year==90&sum_switch_m>=0
summarize sum_switch_m if year==90&sum_switch_m>=1
summarize sum_switch_m if year==90&sum_switch_m>=2 
summarize sum_switch_m if year==90&sum_switch_m>=3
 


****** (5) Industry 311 ********

clear
set memory 500m
set matsize 800

use "C:\Documents and Settings\hkasahar\My Documents\work\project\exportimport\data\DATA_KL_EXIM5.dta"
drop if year<90
replace L = 0 if L ==.
sort padron
by padron: egen min_M = min(M) if L>0
drop if min_M == 0
by padron: egen min_Y = min(Y) if L>0
drop if min_Y == 0
by padron: egen sum_L = sum(L)
drop if sum_L == 0

sort padron
by padron: egen mciiu = mode(ciiu), minmode
by padron: egen sumreentry = sum(reentry)
keep if mciiu==3111|mciiu==3112|mciiu==3113|mciiu==3114|mciiu==3115|mciiu==3116|mciiu==3118|mciiu==3119
 
gen d_x1 = (1-d_m)*d_e if L~=0
gen d_m1 = (1-d_e)*d_m if L~=0
gen d_xm = d_e*d_m if L~=0
gen dxm = 1 if L~=0
replace dxm = 2 if L~=0&d_x1==1
replace dxm = 3 if L~=0&d_m1==1
replace dxm = 4 if L~=0&d_xm == 1

gen switch = 0
gen switch_e = 0
gen switch_m = 0

* define ``switch"
keep switch switch_e switch_m rk dxm d_m d_e export rexport lnkm lnkt lnkb rent fta L_wh L_bl w_wh w_bl rVA rY rM rMi rMd rE Y M Mi Md E rim rit rib im it rkm rkt rkb ciiu cu  exit entry defmach_f defcons_f L padron year
reshape wide switch switch_e switch_m rk dxm d_m d_e export rexport lnkm lnkt lnkb rent fta L_wh L_bl w_wh w_bl rVA rY rM rMi rMd rE Y M Mi Md E rim rit rib im it rkm rkt rkb ciiu cu  exit entry defmach_f defcons_f L, i(padron) j(year)

replace switch91 = 1 if dxm90~=dxm91&L90~=0&L91~=0
replace switch92 = 1 if dxm91~=dxm92&L91~=0&L92~=0
replace switch93 = 1 if dxm92~=dxm93&L92~=0&L93~=0
replace switch94 = 1 if dxm93~=dxm94&L93~=0&L94~=0
replace switch95 = 1 if dxm94~=dxm95&L94~=0&L95~=0
replace switch96 = 1 if dxm95~=dxm96&L95~=0&L96~=0

replace switch_e91 = 1 if d_e90~=d_e91&L90~=0&L91~=0
replace switch_e92 = 1 if d_e91~=d_e92&L91~=0&L92~=0
replace switch_e93 = 1 if d_e92~=d_e93&L92~=0&L93~=0
replace switch_e94 = 1 if d_e93~=d_e94&L93~=0&L94~=0
replace switch_e95 = 1 if d_e94~=d_e95&L94~=0&L95~=0
replace switch_e96 = 1 if d_e95~=d_e96&L95~=0&L96~=0

replace switch_m91 = 1 if d_m90~=d_m91&L90~=0&L91~=0
replace switch_m92 = 1 if d_m91~=d_m92&L91~=0&L92~=0
replace switch_m93 = 1 if d_m92~=d_m93&L92~=0&L93~=0
replace switch_m94 = 1 if d_m93~=d_m94&L93~=0&L94~=0
replace switch_m95 = 1 if d_m94~=d_m95&L94~=0&L95~=0
replace switch_m96 = 1 if d_m95~=d_m96&L95~=0&L96~=0

gen dplant = 0
replace dplant = 1 if L90~=0&L91~=0&L92~=0&L93~=0&L94~=0&L95~=0&L96~=0

reshape long switch switch_e switch_m rk dxm d_m d_e export rexport lnkm lnkt lnkb rent fta L_wh L_bl w_wh w_bl rVA rY rM rMi rMd rE Y M Mi Md E rim rit rib im it rkm rkt rkb ciiu cu exit entry defmach_f defcons_f L, i(padron) j(year)

summarize switch if dplant==1, detail

by padron: egen sum_switch = total(switch) 
by padron: egen sum_switch_e = total(switch_e) 
by padron: egen sum_switch_m = total(switch_m)  
 
summarize sum_switch if year==90&sum_switch>=0
summarize sum_switch if year==90&sum_switch>=1
summarize sum_switch if year==90&sum_switch>=2 
summarize sum_switch if year==90&sum_switch>=3 

summarize sum_switch_e if year==90&sum_switch_e>=0
summarize sum_switch_e if year==90&sum_switch_e>=1
summarize sum_switch_e if year==90&sum_switch_e>=2 
summarize sum_switch_e if year==90&sum_switch_e>=3 
 
summarize sum_switch_m if year==90&sum_switch_m>=0
summarize sum_switch_m if year==90&sum_switch_m>=1
summarize sum_switch_m if year==90&sum_switch_m>=2 
summarize sum_switch_m if year==90&sum_switch_m>=3
 



****** (6) Industry 321 ********

clear
set memory 500m
set matsize 800

use "C:\Documents and Settings\hkasahar\My Documents\work\project\exportimport\data\DATA_KL_EXIM5.dta"
drop if year<90
replace L = 0 if L ==.
sort padron
by padron: egen min_M = min(M) if L>0
drop if min_M == 0
by padron: egen min_Y = min(Y) if L>0
drop if min_Y == 0
by padron: egen sum_L = sum(L)
drop if sum_L == 0

sort padron
by padron: egen sumreentry = sum(reentry)
keep if mcu==321
 
gen d_x1 = (1-d_m)*d_e if L~=0
gen d_m1 = (1-d_e)*d_m if L~=0
gen d_xm = d_e*d_m if L~=0
gen dxm = 1 if L~=0
replace dxm = 2 if L~=0&d_x1==1
replace dxm = 3 if L~=0&d_m1==1
replace dxm = 4 if L~=0&d_xm == 1

gen switch = 0
gen switch_e = 0
gen switch_m = 0

* define ``switch"
keep switch switch_e switch_m rk dxm d_m d_e export rexport lnkm lnkt lnkb rent fta L_wh L_bl w_wh w_bl rVA rY rM rMi rMd rE Y M Mi Md E rim rit rib im it rkm rkt rkb ciiu cu  exit entry defmach_f defcons_f L padron year
reshape wide switch switch_e switch_m rk dxm d_m d_e export rexport lnkm lnkt lnkb rent fta L_wh L_bl w_wh w_bl rVA rY rM rMi rMd rE Y M Mi Md E rim rit rib im it rkm rkt rkb ciiu cu  exit entry defmach_f defcons_f L, i(padron) j(year)

replace switch91 = 1 if dxm90~=dxm91&L90~=0&L91~=0
replace switch92 = 1 if dxm91~=dxm92&L91~=0&L92~=0
replace switch93 = 1 if dxm92~=dxm93&L92~=0&L93~=0
replace switch94 = 1 if dxm93~=dxm94&L93~=0&L94~=0
replace switch95 = 1 if dxm94~=dxm95&L94~=0&L95~=0
replace switch96 = 1 if dxm95~=dxm96&L95~=0&L96~=0

replace switch_e91 = 1 if d_e90~=d_e91&L90~=0&L91~=0
replace switch_e92 = 1 if d_e91~=d_e92&L91~=0&L92~=0
replace switch_e93 = 1 if d_e92~=d_e93&L92~=0&L93~=0
replace switch_e94 = 1 if d_e93~=d_e94&L93~=0&L94~=0
replace switch_e95 = 1 if d_e94~=d_e95&L94~=0&L95~=0
replace switch_e96 = 1 if d_e95~=d_e96&L95~=0&L96~=0

replace switch_m91 = 1 if d_m90~=d_m91&L90~=0&L91~=0
replace switch_m92 = 1 if d_m91~=d_m92&L91~=0&L92~=0
replace switch_m93 = 1 if d_m92~=d_m93&L92~=0&L93~=0
replace switch_m94 = 1 if d_m93~=d_m94&L93~=0&L94~=0
replace switch_m95 = 1 if d_m94~=d_m95&L94~=0&L95~=0
replace switch_m96 = 1 if d_m95~=d_m96&L95~=0&L96~=0

gen dplant = 0
replace dplant = 1 if L90~=0&L91~=0&L92~=0&L93~=0&L94~=0&L95~=0&L96~=0

reshape long switch switch_e switch_m rk dxm d_m d_e export rexport lnkm lnkt lnkb rent fta L_wh L_bl w_wh w_bl rVA rY rM rMi rMd rE Y M Mi Md E rim rit rib im it rkm rkt rkb ciiu cu exit entry defmach_f defcons_f L, i(padron) j(year)

summarize switch if dplant==1, detail

by padron: egen sum_switch = total(switch) 
by padron: egen sum_switch_e = total(switch_e) 
by padron: egen sum_switch_m = total(switch_m)  
 
summarize sum_switch if year==90&sum_switch>=0
summarize sum_switch if year==90&sum_switch>=1
summarize sum_switch if year==90&sum_switch>=2 
summarize sum_switch if year==90&sum_switch>=3 

summarize sum_switch_e if year==90&sum_switch_e>=0
summarize sum_switch_e if year==90&sum_switch_e>=1
summarize sum_switch_e if year==90&sum_switch_e>=2 
summarize sum_switch_e if year==90&sum_switch_e>=3 
 
summarize sum_switch_m if year==90&sum_switch_m>=0
summarize sum_switch_m if year==90&sum_switch_m>=1
summarize sum_switch_m if year==90&sum_switch_m>=2 
summarize sum_switch_m if year==90&sum_switch_m>=3
 




****** (7) Industry 331 ********

clear
set memory 500m
set matsize 800

use "C:\Documents and Settings\hkasahar\My Documents\work\project\exportimport\data\DATA_KL_EXIM5.dta"
drop if year<90
replace L = 0 if L ==.
sort padron
by padron: egen min_M = min(M) if L>0
drop if min_M == 0
by padron: egen min_Y = min(Y) if L>0
drop if min_Y == 0
by padron: egen sum_L = sum(L)
drop if sum_L == 0

sort padron
by padron: egen sumreentry = sum(reentry)
keep if mcu==331
 
gen d_x1 = (1-d_m)*d_e if L~=0
gen d_m1 = (1-d_e)*d_m if L~=0
gen d_xm = d_e*d_m if L~=0
gen dxm = 1 if L~=0
replace dxm = 2 if L~=0&d_x1==1
replace dxm = 3 if L~=0&d_m1==1
replace dxm = 4 if L~=0&d_xm == 1

gen switch = 0
gen switch_e = 0
gen switch_m = 0

* define ``switch"
keep switch switch_e switch_m rk dxm d_m d_e export rexport lnkm lnkt lnkb rent fta L_wh L_bl w_wh w_bl rVA rY rM rMi rMd rE Y M Mi Md E rim rit rib im it rkm rkt rkb ciiu cu  exit entry defmach_f defcons_f L padron year
reshape wide switch switch_e switch_m rk dxm d_m d_e export rexport lnkm lnkt lnkb rent fta L_wh L_bl w_wh w_bl rVA rY rM rMi rMd rE Y M Mi Md E rim rit rib im it rkm rkt rkb ciiu cu  exit entry defmach_f defcons_f L, i(padron) j(year)

replace switch91 = 1 if dxm90~=dxm91&L90~=0&L91~=0
replace switch92 = 1 if dxm91~=dxm92&L91~=0&L92~=0
replace switch93 = 1 if dxm92~=dxm93&L92~=0&L93~=0
replace switch94 = 1 if dxm93~=dxm94&L93~=0&L94~=0
replace switch95 = 1 if dxm94~=dxm95&L94~=0&L95~=0
replace switch96 = 1 if dxm95~=dxm96&L95~=0&L96~=0

replace switch_e91 = 1 if d_e90~=d_e91&L90~=0&L91~=0
replace switch_e92 = 1 if d_e91~=d_e92&L91~=0&L92~=0
replace switch_e93 = 1 if d_e92~=d_e93&L92~=0&L93~=0
replace switch_e94 = 1 if d_e93~=d_e94&L93~=0&L94~=0
replace switch_e95 = 1 if d_e94~=d_e95&L94~=0&L95~=0
replace switch_e96 = 1 if d_e95~=d_e96&L95~=0&L96~=0

replace switch_m91 = 1 if d_m90~=d_m91&L90~=0&L91~=0
replace switch_m92 = 1 if d_m91~=d_m92&L91~=0&L92~=0
replace switch_m93 = 1 if d_m92~=d_m93&L92~=0&L93~=0
replace switch_m94 = 1 if d_m93~=d_m94&L93~=0&L94~=0
replace switch_m95 = 1 if d_m94~=d_m95&L94~=0&L95~=0
replace switch_m96 = 1 if d_m95~=d_m96&L95~=0&L96~=0

gen dplant = 0
replace dplant = 1 if L90~=0&L91~=0&L92~=0&L93~=0&L94~=0&L95~=0&L96~=0

reshape long switch switch_e switch_m rk dxm d_m d_e export rexport lnkm lnkt lnkb rent fta L_wh L_bl w_wh w_bl rVA rY rM rMi rMd rE Y M Mi Md E rim rit rib im it rkm rkt rkb ciiu cu exit entry defmach_f defcons_f L, i(padron) j(year)

summarize switch if dplant==1, detail

by padron: egen sum_switch = total(switch) 
by padron: egen sum_switch_e = total(switch_e) 
by padron: egen sum_switch_m = total(switch_m)  
 
summarize sum_switch if year==90&sum_switch>=0
summarize sum_switch if year==90&sum_switch>=1
summarize sum_switch if year==90&sum_switch>=2 
summarize sum_switch if year==90&sum_switch>=3 

summarize sum_switch_e if year==90&sum_switch_e>=0
summarize sum_switch_e if year==90&sum_switch_e>=1
summarize sum_switch_e if year==90&sum_switch_e>=2 
summarize sum_switch_e if year==90&sum_switch_e>=3 
 
summarize sum_switch_m if year==90&sum_switch_m>=0
summarize sum_switch_m if year==90&sum_switch_m>=1
summarize sum_switch_m if year==90&sum_switch_m>=2 
summarize sum_switch_m if year==90&sum_switch_m>=3
 





 
