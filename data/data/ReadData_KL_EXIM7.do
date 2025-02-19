* ReadData_KL_EXIM7.do
* Last Modified: Sept, 2010
* Notes: This do-file creates the data set for estimating an empirical model of firm's import/export
*	   ``Productivity and the Decision to Import and Export: Theory and Evidence'' by Hiro Kasahara and Bev Lapham
* 
* This is for the third revision at JIE
*
* Written by Hiro Kasahara, UBC
 


clear
set memory 500m
use "C:\Users\Hiro\Documents\Work\Project\exportimport\data\DATA_KL_EXIM5.dta"
drop if year<90
*drop if dplnt==1
replace L = 0 if L ==.
sort padron
by padron: egen min_M = min(M) if L>0
drop if min_M == 0
by padron: egen min_Y = min(Y) if L>0
drop if min_Y == 0
by padron: egen sum_L = sum(L)
drop if sum_L == 0
replace rY = 0 if rY ==.
replace year = 0 if year ==.
replace ciiu = 0 if ciiu ==.
replace mcu = 0 if mcu==. 
replace export = 0 if export==.
replace rexport = 0 if rexport==.
replace L = 0 if L==.
replace rM = 0 if rM==.
replace rMi = 0 if rMi==.
replace d_m_ind = 0 if d_m_ind ==.
replace Mi_M_ind = 0 if Mi_M_ind ==.
replace d_e_ind = 0 if d_e_ind ==.
replace E_Y_ind = 0 if E_Y_ind ==.
replace d_e = -1 if d_e==.
replace d_m = -1 if d_m==.
replace Mi = 0 if Mi==.
replace M = 0 if M==.
replace markup = 0 if markup==.
sort padron year
keep padron year ciiu rY d_e d_m rexport rMi L rM markup
order padron year ciiu rY d_e d_m rexport rMi L rM markup 
outsheet using "C:\Users\Hiro\Documents\Work\Project\exportimport\data\DATA_KL_EXIM5.xl", replace 


by year, sort: summarize markup if L>0
by year, sort: summarize d_e d_m if L>0
gen Mi_M = rMi/rM
gen e_Y = rexport/rY
by year, sort: egen med_e_Y = median(e_Y) if d_e>0&L>0
by year, sort: egen med_Mi_M = median(Mi_M) if d_m>0&L>0
by year, sort: egen ave_e_Y = mean(e_Y) if d_e>0&L>0
by year, sort: egen ave_Mi_M = mean(Mi_M) if d_m>0&L>0
by year: summarize ave_e_Y med_e_Y if d_e>0&L>0
by year: summarize ave_Mi_M med_Mi_M if d_m>0&L>0



************************* Industry 3220 (Wearing Apparel) *************************************

clear
set memory 500m
use "C:\Users\Hiro\Documents\Work\Project\exportimport\data\DATA_KL_EXIM5.dta"
drop if year<90
*drop if dplnt==1
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

sort year
merge year using "C:\Users\Hiro\Documents\Work\Project\exportimport\data\deflate9096.dta"
drop _merge
sort year
merge year using "C:\Users\Hiro\Documents\Work\Project\exportimport\data\deflate.dta"
replace rY = defmanuf*rY/y3220
replace rM = defmanuf*rM/y3220
replace rMi = defmanuf*rMi/y3220
replace rexport = defmanuf*rexport/y3220
drop if year<90

replace rY = 0 if rY ==.
replace year = 0 if year ==.
replace ciiu = 0 if ciiu ==.
replace mcu = 0 if mcu==. 
replace export = 0 if export==.
replace rexport = 0 if rexport==.
replace L = 0 if L==.
replace rM = 0 if rM==.
replace rMi = 0 if rMi==.
replace d_m_ind = 0 if d_m_ind ==.
replace Mi_M_ind = 0 if Mi_M_ind ==.
replace d_e_ind = 0 if d_e_ind ==.
replace E_Y_ind = 0 if E_Y_ind ==.
replace d_e = -1 if d_e==.
replace d_m = -1 if d_m==.
replace Mi = 0 if Mi==.
replace M = 0 if M==.
replace markup = 0 if markup==.
sort padron year
keep padron year ciiu rY d_e d_m rexport rMi L rM markup
order padron year ciiu rY d_e d_m rexport rMi L rM markup 
outsheet using "C:\Users\Hiro\Documents\Work\Project\exportimport\data\DATA_KL_3220_ciiu.xl", replace


summarize markup if L>0&year==90
summarize markup if L>0&year==96




summarize padron year ciiu rY d_e d_m rexport rMi L rM markup
gen lny = log(rY)
reg markup lny d_e d_m if L>0
sort year
by year: summarize markup if L>0
by year: summarize d_e d_m if L>0

gen Mi_M = rMi/rM
gen e_Y = rexport/rY

summarize rY rM L markup if year == 90&L>0
summarize rexport e_Y if year == 90&d_e>0&L>0
summarize rMi Mi_M if year == 90&d_m>0&L>0
summarize rY rM L markup if year == 96&L>0
summarize rexport e_Y if year == 96&d_e>0&L>0
summarize rMi Mi_M if year == 96&d_m>0&L>0




by year, sort: egen med_e_Y = median(e_Y) if d_e>0&L>0
by year, sort: egen med_Mi_M = median(Mi_M) if d_m>0&L>0
by year, sort: egen ave_e_Y = mean(e_Y) if d_e>0&L>0
by year, sort: egen ave_Mi_M = mean(Mi_M) if d_m>0&L>0
by year: summarize ave_e_Y med_e_Y if d_e>0&L>0
by year: summarize ave_Mi_M med_Mi_M if d_m>0&L>0



************************* Industry 3560 (Plastic Products) *************************************

clear
set memory 500m
use "C:\Users\Hiro\Documents\Work\Project\exportimport\data\DATA_KL_EXIM5.dta"
drop if year<90
*drop if dplnt==1
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
keep if mciiu==3560
*keep if sdciiu==0|sdciiu==.
*keep if sumreentry==0

*gen markup = (rY-( w_wh*L_wh + w_bl*L_bl + rM + rE))/rY if L~=.
*gen sigma = 1/markup if L~=.
*gen M_share = rM/(w_wh*L_wh + w_bl*L_bl + rM + rE) if L>0
*summarize markup M_share if L>0

sort year
merge year using "C:\Users\Hiro\Documents\Work\Project\exportimport\data\deflate9096.dta"
drop _merge
sort year
merge year using "C:\Users\Hiro\Documents\Work\Project\exportimport\data\deflate.dta"
replace rY = defmanuf*rY/y3560
replace rM = defmanuf*rM/y3560
replace rMi = defmanuf*rMi/y3560
replace rexport = defmanuf*rexport/y3560
drop if year<90

replace rY = 0 if rY ==.
replace year = 0 if year ==.
replace ciiu = 0 if ciiu ==.
replace mcu = 0 if mcu==. 
replace export = 0 if export==.
replace rexport = 0 if rexport==.
replace L = 0 if L==.
replace rM = 0 if rM==.
replace rMi = 0 if rMi==.
replace d_m_ind = 0 if d_m_ind ==.
replace Mi_M_ind = 0 if Mi_M_ind ==.
replace d_e_ind = 0 if d_e_ind ==.
replace E_Y_ind = 0 if E_Y_ind ==.
replace d_e = -1 if d_e==.
replace d_m = -1 if d_m==.
replace Mi = 0 if Mi==.
replace M = 0 if M==.
replace markup = 0 if markup==.
sort padron year
keep padron year ciiu rY d_e d_m rexport rMi L rM markup
order padron year ciiu rY d_e d_m rexport rMi L rM markup 
outsheet using "C:\Users\Hiro\Documents\Work\Project\exportimport\data\DATA_KL_3560_ciiu.xl", replace

summarize  padron year ciiu rY d_e d_m rexport rMi L rM markup

gen lny = log(rY)
reg markup lny d_e d_m if L>0
sort year
by year: summarize markup if L>0
by year: summarize d_e d_m if L>0

gen Mi_M = rMi/rM
gen e_Y = rexport/rY
summarize rY rM L markup if year == 90&L>0
summarize rexport e_Y if year == 90&d_e>0&L>0
summarize rMi Mi_M if year == 90&d_m>0&L>0
summarize rY rM L markup if year == 96&L>0
summarize rexport e_Y if year == 96&d_e>0&L>0
summarize rMi Mi_M if year == 96&d_m>0&L>0



summarize rY rM L rexport rMi e_Y Mi_M markup if year == 90
summarize rY rM L rexport rMi e_Y Mi_M markup if year == 96

sort year
by year, sort: egen med_e_Y = median(e_Y) if d_e>0&L>0
by year, sort: egen med_Mi_M = median(Mi_M) if d_m>0&L>0
by year, sort: egen ave_e_Y = mean(e_Y) if d_e>0&L>0
by year, sort: egen ave_Mi_M = mean(Mi_M) if d_m>0&L>0
by year: summarize ave_e_Y med_e_Y if d_e>0&L>0
by year: summarize ave_Mi_M med_Mi_M if d_m>0&L>0


************************* Industry 3813 (Structural Metal) *************************************

clear
set memory 500m
use "C:\Users\Hiro\Documents\Work\Project\exportimport\data\DATA_KL_EXIM5.dta"
drop if year<90
*drop if dplnt==1
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
replace mciiu = 3813 if mciiu==3814|mciiu==3815
keep if mciiu==3813
*keep if sdciiu==0|sdciiu==.
*keep if sumreentry==0

*gen markup = (rY-( w_wh*L_wh + w_bl*L_bl + rM + rE))/rY if L~=.
*gen sigma = 1/markup if L~=.
*gen M_share = rM/(w_wh*L_wh + w_bl*L_bl + rM + rE) if L>0
*summarize markup M_share if L>0

sort year
merge year using "C:\Users\Hiro\Documents\Work\Project\exportimport\data\deflate9096.dta"
drop _merge
sort year
merge year using "C:\Users\Hiro\Documents\Work\Project\exportimport\data\deflate.dta"
replace rY = defmanuf*rY/y3813
replace rM = defmanuf*rM/y3813
replace rMi = defmanuf*rMi/y3813
replace rexport = defmanuf*rexport/y3813
drop if year<90

replace rY = 0 if rY ==.
replace year = 0 if year ==.
replace ciiu = 0 if ciiu ==.
replace mcu = 0 if mcu==. 
replace export = 0 if export==.
replace rexport = 0 if rexport==.
replace L = 0 if L==.
replace rM = 0 if rM==.
replace rMi = 0 if rMi==.
replace d_m_ind = 0 if d_m_ind ==.
replace Mi_M_ind = 0 if Mi_M_ind ==.
replace d_e_ind = 0 if d_e_ind ==.
replace E_Y_ind = 0 if E_Y_ind ==.
replace d_e = -1 if d_e==.
replace d_m = -1 if d_m==.
replace Mi = 0 if Mi==.
replace M = 0 if M==.
replace markup = 0 if markup==.
sort padron year
keep padron year ciiu rY d_e d_m rexport rMi L rM markup
order padron year ciiu rY d_e d_m rexport rMi L rM markup 
outsheet using "C:\Users\Hiro\Documents\Work\Project\exportimport\data\DATA_KL_3813_ciiu.xl", replace

summarize  padron year ciiu rY d_e d_m rexport rMi L rM markup

gen lny = log(rY)
reg markup lny d_e d_m if L>0
sort year
by year: summarize markup if L>0
by year: summarize d_e d_m if L>0

gen Mi_M = rMi/rM
gen e_Y = rexport/rY
summarize rY rM L markup if year == 90&L>0
summarize rexport e_Y if year == 90&d_e>0&L>0
summarize rMi Mi_M if year == 90&d_m>0&L>0
summarize rY rM L markup if year == 96&L>0
summarize rexport e_Y if year == 96&d_e>0&L>0
summarize rMi Mi_M if year == 96&d_m>0&L>0


summarize rY rM L rexport rMi e_Y Mi_M markup if year == 90
summarize rY rM L rexport rMi e_Y Mi_M markup if year == 96

sort year
by year, sort: egen med_e_Y = median(e_Y) if d_e>0&L>0
by year, sort: egen med_Mi_M = median(Mi_M) if d_m>0&L>0
by year, sort: egen ave_e_Y = mean(e_Y) if d_e>0&L>0
by year, sort: egen ave_Mi_M = mean(Mi_M) if d_m>0&L>0
by year: summarize ave_e_Y med_e_Y if d_e>0&L>0
by year: summarize ave_Mi_M med_Mi_M if d_m>0&L>0




************************* Industry 381 (Fabricated Metal Products) *************************************
clear
set memory 500m
use "C:\Users\Hiro\Documents\Work\Project\exportimport\data\DATA_KL_EXIM5.dta"
drop if year<90
*drop if dplnt==1
replace L = 0 if L ==.
sort padron
by padron: egen min_M = min(M) if L>0
drop if min_M == 0
by padron: egen min_Y = min(Y) if L>0
drop if min_Y == 0
by padron: egen sum_L = sum(L)
drop if sum_L == 0


by padron: egen sumreentry = sum(reentry)
keep if mcu==381


*gen markup = (rY-( w_wh*L_wh + w_bl*L_bl + rM + rE))/rY if L~=.
*gen sigma = 1/markup if L~=.
gen M_share = rM/(w_wh*L_wh + w_bl*L_bl + rM + rE) if L>0
summarize markup M_share if L>0

sort year
merge year using "C:\Users\Hiro\Documents\Work\Project\exportimport\data\deflate9096.dta"
drop _merge
sort year
merge year using "C:\Users\Hiro\Documents\Work\Project\exportimport\data\deflate.dta"
replace rY = defmanuf*rY/y381
replace rM = defmanuf*rM/y381
replace rMi = defmanuf*rMi/y381
replace rexport = defmanuf*rexport/y381
drop if year<90

replace rY = 0 if rY ==.
replace year = 0 if year ==.
replace ciiu = 0 if ciiu ==.
replace mcu = 0 if mcu==. 
replace export = 0 if export==.
replace rexport = 0 if rexport==.
replace L = 0 if L==.
replace rM = 0 if rM==.
replace rMi = 0 if rMi==.
replace d_m_ind = 0 if d_m_ind ==.
replace Mi_M_ind = 0 if Mi_M_ind ==.
replace d_e_ind = 0 if d_e_ind ==.
replace E_Y_ind = 0 if E_Y_ind ==.
replace d_e = -1 if d_e==.
replace d_m = -1 if d_m==.
replace Mi = 0 if Mi==.
replace M = 0 if M==.
replace markup = 0 if markup==.
sort padron year
keep padron year ciiu rY d_e d_m rexport rMi L rM markup
order padron year ciiu rY d_e d_m rexport rMi L rM markup 
outsheet using "C:\Users\Hiro\Documents\Work\Project\exportimport\data\DATA_KL_381_ciiu.xl", replace

summarize  padron year ciiu rY d_e d_m rexport rMi L rM markup

gen lny = log(rY)
reg markup lny d_e d_m if L>0
sort year
by year: summarize markup if L>0
by year: summarize d_e d_m if L>0

gen Mi_M = rMi/rM
gen e_Y = rexport/rY
summarize rY rM L markup if year == 90&L>0
summarize rexport e_Y if year == 90&d_e>0&L>0
summarize rMi Mi_M if year == 90&d_m>0&L>0
summarize rY rM L markup if year == 96&L>0
summarize rexport e_Y if year == 96&d_e>0&L>0
summarize rMi Mi_M if year == 96&d_m>0&L>0


summarize rY rM L rexport rMi e_Y Mi_M markup if year == 90
summarize rY rM L rexport rMi e_Y Mi_M markup if year == 96

sort year
by year, sort: egen med_e_Y = median(e_Y) if d_e>0&L>0
by year, sort: egen med_Mi_M = median(Mi_M) if d_m>0&L>0
by year, sort: egen ave_e_Y = mean(e_Y) if d_e>0&L>0
by year, sort: egen ave_Mi_M = mean(Mi_M) if d_m>0&L>0
by year: summarize ave_e_Y med_e_Y if d_e>0&L>0
by year: summarize ave_Mi_M med_Mi_M if d_m>0&L>0



****************** Industry 311 except for 3117 (Food products except for Bakery products) *******************


clear
set memory 500m
use "C:\Users\Hiro\Documents\Work\Project\exportimport\data\DATA_KL_EXIM5.dta"
drop if year<90
*sort ciiu
*by ciiu: egen sumY = sum(rY) if year==90
*summarize sumY if year==90&ciiu==3111
*summarize sumY if year==90&ciiu==3112
*summarize sumY if year==90&ciiu==3113
*summarize sumY if year==90&ciiu==3114
*summarize sumY if year==90&ciiu==3115
*summarize sumY if year==90&ciiu==3116
*summarize sumY if year==90&ciiu==3118
*summarize sumY if year==90&ciiu==3119

*drop if dplnt==1
replace L = 0 if L ==.
sort padron
by padron: egen min_M = min(M) if L>0
drop if min_M == 0
by padron: egen min_Y = min(Y) if L>0
drop if min_Y == 0
by padron: egen sum_L = sum(L)
drop if sum_L == 0

by padron: egen mciiu = mode(ciiu), minmode
by padron: egen sumreentry = sum(reentry)
keep if mciiu==3111|mciiu==3112|mciiu==3113|mciiu==3114|mciiu==3115|mciiu==3116|mciiu==3118|mciiu==3119

gen M_share = rM/(w_wh*L_wh + w_bl*L_bl + rM + rE) if L>0
summarize markup M_share if L>0


sort year
merge year using "C:\Users\Hiro\Documents\Work\Project\exportimport\data\deflate9096.dta"
drop _merge
sort year
merge year using "C:\Users\Hiro\Documents\Work\Project\exportimport\data\deflate.dta"
*gen y3110 = 0.180*y3111 + 0.166*y3112 + 0.096*y3113 + 0.105*y3114 + 0.223*y3115 + 0.117*y3117 + 0.073*y3118 + 0.039*y3119 
replace rY = defmanuf*rY/y3110
replace rM = defmanuf*rM/y3110
replace rMi = defmanuf*rMi/y3110
replace rexport = defmanuf*rexport/y3110
drop if year<90

replace rY = 0 if rY ==.
replace year = 0 if year ==.
replace ciiu = 0 if ciiu ==.
replace mcu = 0 if mcu==. 
replace export = 0 if export==.
replace rexport = 0 if rexport==.
replace L = 0 if L==.
replace rM = 0 if rM==.
replace rMi = 0 if rMi==.
replace d_m_ind = 0 if d_m_ind ==.
replace Mi_M_ind = 0 if Mi_M_ind ==.
replace d_e_ind = 0 if d_e_ind ==.
replace E_Y_ind = 0 if E_Y_ind ==.
replace d_e = -1 if d_e==.
replace d_m = -1 if d_m==.
replace Mi = 0 if Mi==.
replace M = 0 if M==.
replace markup = 0 if markup==.
sort padron year
keep padron year ciiu rY d_e d_m rexport rMi L rM markup
order padron year ciiu rY d_e d_m rexport rMi L rM markup 
outsheet using "C:\Users\Hiro\Documents\Work\Project\exportimport\data\DATA_KL_311_ciiu.xl", replace

gen Mi_M = rMi/rM
gen e_Y = rexport/rY
summarize rY rM L d_e d_m markup if year == 90&L>0
summarize rexport e_Y if year == 90&d_e>0&L>0
summarize rMi Mi_M if year == 90&d_m>0&L>0
summarize rY rM L d_e d_m markup if year == 96&L>0
summarize rexport e_Y if year == 96&d_e>0&L>0
summarize rMi Mi_M if year == 96&d_m>0&L>0







************************* Industry 321 (Textile) *************************************


clear
set memory 500m
use "C:\Users\Hiro\Documents\Work\Project\exportimport\data\DATA_KL_EXIM5.dta"
drop if year<90
*drop if dplnt==1
replace L = 0 if L ==.
sort padron
by padron: egen min_M = min(M) if L>0
drop if min_M == 0
by padron: egen min_Y = min(Y) if L>0
drop if min_Y == 0
by padron: egen sum_L = sum(L)
drop if sum_L == 0


by padron: egen sumreentry = sum(reentry)
keep if mcu==321

gen M_share = rM/(w_wh*L_wh + w_bl*L_bl + rM + rE) if L>0
summarize markup M_share if L>0


sort year
merge year using "C:\Users\Hiro\Documents\Work\Project\exportimport\data\deflate9096.dta"
drop _merge
sort year
merge year using "C:\Users\Hiro\Documents\Work\Project\exportimport\data\deflate.dta"
replace rY = defmanuf*rY/y321
replace rM = defmanuf*rM/y321
replace rMi = defmanuf*rMi/y321
replace rexport = defmanuf*rexport/y321
drop if year<90

replace rY = 0 if rY ==.
replace year = 0 if year ==.
replace ciiu = 0 if ciiu ==.
replace mcu = 0 if mcu==. 
replace export = 0 if export==.
replace rexport = 0 if rexport==.
replace L = 0 if L==.
replace rM = 0 if rM==.
replace rMi = 0 if rMi==.
replace d_m_ind = 0 if d_m_ind ==.
replace Mi_M_ind = 0 if Mi_M_ind ==.
replace d_e_ind = 0 if d_e_ind ==.
replace E_Y_ind = 0 if E_Y_ind ==.
replace d_e = -1 if d_e==.
replace d_m = -1 if d_m==.
replace Mi = 0 if Mi==.
replace M = 0 if M==.
replace markup = 0 if markup==.
sort padron year
keep padron year ciiu rY d_e d_m rexport rMi L rM markup
order padron year ciiu rY d_e d_m rexport rMi L rM markup 
outsheet using "C:\Users\Hiro\Documents\Work\Project\exportimport\data\DATA_KL_321_ciiu.xl", replace


gen Mi_M = rMi/rM
gen e_Y = rexport/rY
summarize rY rM L d_e d_m markup if year == 90&L>0
summarize rexport e_Y if year == 90&d_e>0&L>0
summarize rMi Mi_M if year == 90&d_m>0&L>0
summarize rY rM L d_e d_m markup if year == 96&L>0
summarize rexport e_Y if year == 96&d_e>0&L>0
summarize rMi Mi_M if year == 96&d_m>0&L>0




************************* Industry 331 (Wood Products) *************************************


clear
set memory 500m
use "C:\Users\Hiro\Documents\Work\Project\exportimport\data\DATA_KL_EXIM5.dta"
drop if year<90
*drop if dplnt==1
replace L = 0 if L ==.
sort padron
by padron: egen min_M = min(M) if L>0
drop if min_M == 0
by padron: egen min_Y = min(Y) if L>0
drop if min_Y == 0
by padron: egen sum_L = sum(L)
drop if sum_L == 0


by padron: egen sumreentry = sum(reentry)
keep if mcu==331

gen M_share = rM/(w_wh*L_wh + w_bl*L_bl + rM + rE) if L>0
summarize markup M_share if L>0


sort year
merge year using "C:\Users\Hiro\Documents\Work\Project\exportimport\data\deflate9096.dta"
drop _merge
sort year
merge year using "C:\Users\Hiro\Documents\Work\Project\exportimport\data\deflate.dta"
replace rY = defmanuf*rY/y331
replace rM = defmanuf*rM/y331
replace rMi = defmanuf*rMi/y331
replace rexport = defmanuf*rexport/y331
drop if year<90

replace rY = 0 if rY ==.
replace year = 0 if year ==.
replace ciiu = 0 if ciiu ==.
replace mcu = 0 if mcu==. 
replace export = 0 if export==.
replace rexport = 0 if rexport==.
replace L = 0 if L==.
replace rM = 0 if rM==.
replace rMi = 0 if rMi==.
replace d_m_ind = 0 if d_m_ind ==.
replace Mi_M_ind = 0 if Mi_M_ind ==.
replace d_e_ind = 0 if d_e_ind ==.
replace E_Y_ind = 0 if E_Y_ind ==.
replace d_e = -1 if d_e==.
replace d_m = -1 if d_m==.
replace Mi = 0 if Mi==.
replace M = 0 if M==.
replace markup = 0 if markup==.
sort padron year
keep padron year ciiu rY d_e d_m rexport rMi L rM markup
order padron year ciiu rY d_e d_m rexport rMi L rM markup 
outsheet using "C:\Users\Hiro\Documents\Work\Project\exportimport\data\DATA_KL_331_ciiu.xl", replace


gen Mi_M = rMi/rM
gen e_Y = rexport/rY
summarize rY rM L d_e d_m markup if year == 90&L>0
summarize rexport e_Y if year == 90&d_e>0&L>0
summarize rMi Mi_M if year == 90&d_m>0&L>0
summarize rY rM L d_e d_m markup if year == 96&L>0
summarize rexport e_Y if year == 96&d_e>0&L>0
summarize rMi Mi_M if year == 96&d_m>0&L>0












