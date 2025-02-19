* ReadData_KL_EXIM5.do
* Last Modified: April, 2008
* Notes: This do-file creates the data set for estimating an empirical model of firm's import/export
*	   ``Productivity and the Decision to Import and Export: Theory and Evidence'' by Hiro Kasahara and Bev Lapham
* 
* This is for the first revision at JIE
*
* Written by Hiro Kasahara, UWO

clear 
set memory 500m
set matsize 800
*******
use "C:\Documents and Settings\hkasahar\My Documents\work\project\exportimport\data\Enia86.dta"
keep year xx1 xx2 xx4 xx49 xx52 xx109-xx129 xx131-xx138 xx157 xx159 xx160 xx161 xx165-xx167 xx37-xx43 xx45 xx75 xx77 xx80 xx82 xx52 xx54 xx56 xx58 xx60 xx62 xx64 xx66 xx68 xx70 xx72 xx73 xx13-xx32
append using "C:\Documents and Settings\hkasahar\My Documents\work\project\exportimport\data\Enia87.dta"
keep year xx1 xx2 xx4 xx49 xx52 xx109-xx129 xx131-xx138 xx157 xx159 xx160 xx161 xx165-xx167 xx37-xx43 xx45 xx75 xx77 xx80 xx82 xx52 xx54 xx56 xx58 xx60 xx62 xx64 xx66 xx68 xx70 xx72 xx73 xx13-xx32
append using "C:\Documents and Settings\hkasahar\My Documents\work\project\exportimport\data\Enia88.dta"
keep year xx1 xx2 xx4 xx49 xx52 xx109-xx129 xx131-xx138 xx157 xx159 xx160 xx161 xx165-xx167 xx37-xx43 xx45 xx75 xx77 xx80 xx82 xx52 xx54 xx56 xx58 xx60 xx62 xx64 xx66 xx68 xx70 xx72 xx73 xx13-xx32
append using "C:\Documents and Settings\hkasahar\My Documents\work\project\exportimport\data\Enia89.dta"
keep year xx1 xx2 xx4 xx49 xx52 xx109-xx129 xx131-xx138 xx157 xx159 xx160 xx161 xx165-xx167 xx37-xx43 xx45 xx75 xx77 xx80 xx82 xx52 xx54 xx56 xx58 xx60 xx62 xx64 xx66 xx68 xx70 xx72 xx73 xx13-xx32 xx152-xx156 
append using "C:\Documents and Settings\hkasahar\My Documents\work\project\exportimport\data\Enia90.dta"
keep year xx1 xx2 xx4 xx49 xx52 xx109-xx129 xx131-xx138 xx157 xx159 xx160 xx161 xx165-xx167 xx37-xx43 xx45 xx75 xx77 xx80 xx82 xx52 xx54 xx56 xx58 xx60 xx62 xx64 xx66 xx68 xx70 xx72 xx73 xx13-xx32 xx152-xx156
append using "C:\Documents and Settings\hkasahar\My Documents\work\project\exportimport\data\Enia91.dta"
keep year xx1 xx2 xx4 xx49 xx52 xx109-xx129 xx131-xx138 xx157 xx159 xx160 xx161 xx165-xx167 xx37-xx43 xx45 xx75 xx77 xx80 xx82 xx52 xx54 xx56 xx58 xx60 xx62 xx64 xx66 xx68 xx70 xx72 xx73 xx13-xx32 xx152-xx156
append using "C:\Documents and Settings\hkasahar\My Documents\work\project\exportimport\data\Enia92.dta"
keep year xx1 xx2 xx4 xx49 xx52 xx109-xx129 xx131-xx138 xx157 xx159 xx160 xx161 xx165-xx167 xx37-xx43 xx45 xx75 xx77 xx80 xx82 xx52 xx54 xx56 xx58 xx60 xx62 xx64 xx66 xx68 xx70 xx72 xx73 xx13-xx32 xx152-xx156 xx236-xx239
append using "C:\Documents and Settings\hkasahar\My Documents\work\project\exportimport\data\Enia93.dta"
keep year xx1 xx2 xx4 xx49 xx52 xx109-xx129 xx131-xx138 xx157 xx159 xx160 xx161 xx165-xx167 xx37-xx43 xx45 xx75 xx77 xx80 xx82 xx52 xx54 xx56 xx58 xx60 xx62 xx64 xx66 xx68 xx70 xx72 xx73 xx13-xx32 xx152-xx156 xx236-xx239
append using "C:\Documents and Settings\hkasahar\My Documents\work\project\exportimport\data\Enia94.dta"
keep year xx1 xx2 xx4 xx49 xx52 xx109-xx129 xx131-xx138 xx157 xx159 xx160 xx161 xx165-xx167 xx37-xx43 xx45 xx75 xx77 xx80 xx82 xx52 xx54 xx56 xx58 xx60 xx62 xx64 xx66 xx68 xx70 xx72 xx73 xx13-xx32 xx152-xx156 xx236-xx239
append using "C:\Documents and Settings\hkasahar\My Documents\work\project\exportimport\data\Enia95.dta"
keep year xx1 xx2 xx4 xx49 xx52 xx109-xx129 xx131-xx138 xx157 xx159 xx160 xx161 xx165-xx167 xx37-xx43 xx45 xx75 xx77 xx80 xx82 xx52 xx54 xx56 xx58 xx60 xx62 xx64 xx66 xx68 xx70 xx72 xx73 xx13-xx32 xx152-xx156 xx236-xx239
append using "C:\Documents and Settings\hkasahar\My Documents\work\project\exportimport\data\Enia96.dta"
keep year xx1 xx2 xx4 xx49 xx52 xx109-xx129 xx131-xx138 xx157 xx159 xx160 xx161 xx165-xx167 xx37-xx43 xx45 xx75 xx77 xx80 xx82 xx52 xx54 xx56 xx58 xx60 xx62 xx64 xx66 xx68 xx70 xx72 xx73 xx13-xx32 xx152-xx156 xx236-xx239

******** Nominal yriables *******
gen padron = xx1
gen ciiu = xx2
gen L = xx165
gen ib = xx109+xx112-xx115+xx118+xx121 + xx152+xx153-xx154+xx155+xx156 if year>=89
replace ib = xx109+xx112-xx115+xx118+xx121 if year<89
gen im = xx110+xx113-xx116+xx119+xx122 
gen it = xx111+xx114-xx117+xx120+xx123

*gen ib = xx109+xx112+xx118+xx121+xx152+xx153+xx155+xx156 if year>=89
*replace ib = xx109+xx112-xx115+xx118+xx121 if year<89
*gen im = xx110+xx113+xx119+xx122 
*gen it = xx111+xx114+xx120+xx123

gen nim = xx110+xx119+xx122
gen nit = xx111+xx120+xx123 
gen kb = xx236+xx237
gen km = xx238
gen kt = xx239
gen FUEL = xx54+xx56+xx58+xx60+xx62+xx64+xx66+xx68+xx70+xx72+xx73
replace FUEL = xx167 if FUEL==.|FUEL==0
gen E = xx49 + xx45 + xx75 + xx77 + FUEL
gen MATT = xx40 
gen Mi = xx41
gen Md = (MATT - Mi) + xx42 + xx43
*gen Y = xx124 + xx125 + xx126 + xx127 + xx128 + xx129 + xx52 + xx121 + xx122 + xx123
*gen mky = xx122 + xx123
*gen bky = xx121
*gen ey = xx52
*replace Y = Y + mky + bky + ey
gen rent = xx80
gen L_wh = xx13+xx14+xx15+xx16+xx17+xx18+xx19+xx20
gen L_bl = xx21+xx22+xx23+xx24
gen w_wh = (xx29+xx31)/L_wh
gen w_bl = (xx30+xx32)/L_bl
gen fta = xx82

gen export = xx157

gen ivy_i = xx135 + xx137
gen ivy_f = xx136 + xx138
gen ivm_i = xx131 + xx133
gen ivm_f = xx132 + xx134

gen Y = xx161 + xx52 + xx121 + xx122 + xx123 
gen M = xx166 + xx49 if year<88
replace M = xx166 if year>=88


***** Construct missing book values of capital stock *****
gen x1 = ln(L_wh) if L~=.
gen x2 = ln(L_bl) if L~=.
gen x3 = ln(M) if L~=.
gen x4 = ln(E) if L~=.
gen x5 = 0 if L~=.
replace x5 = 1 if Mi>0&L~=.
gen x6 = ln(Y) if L~=.
replace x1 = 1 if L_wh<=0&L~=.
replace x2 = 1 if L_bl<=0&L~=.
replace x3 = 1 if M<=0&L~=.
replace x4 = 1 if E<=0&L~=.
replace x6 = 1 if Y<=0&L~=.
tab ciiu, gen(dind)
tab year, gen(dyear)
gen y = ln(km)
reg y x1 x2 x3 x4 x5 x6 dind1-dind86 if year==92&km>0&L~=.
predict yhat if year==92&L~=.
gen lnkm = yhat if year==92&L~=. 
drop yhat
reg y x1 x2 x3 x4 x5 x6 dind1-dind86 if year==93&km>0&L~=.
predict yhat if year==93&L~=.
replace lnkm = yhat if year==93&L~=. 
drop yhat
reg y x1 x2 x3 x4 x5 x6 dind1-dind86 if year==94&km>0&L~=.
predict yhat if year==94&L~=.
replace lnkm = yhat if year==94&L~=. 
drop yhat
reg y x1 x2 x3 x4 x5 x6 dind1-dind86 if year==95&km>0&L~=.
predict yhat if year==95&L~=.
replace lnkm = yhat if year==95&L~=. 
drop yhat
reg y x1 x2 x3 x4 x5 x6 dind1-dind86 if year==96&km>0&L~=.
predict yhat if year==96&L~=.
replace lnkm = yhat if year==96&L~=. 
drop yhat
replace y = ln(kt)
reg y x1 x2 x3 x4 x5 x6 dind1-dind86 if year==92&kt>0&L~=.
predict yhat if year==92&L~=.
gen lnkt = yhat if year==92&L~=. 
drop yhat
reg y x1 x2 x3 x4 x5 x6 dind1-dind86 if year==93&kt>0&L~=.
predict yhat if year==93&L~=.
replace lnkt = yhat if year==93&L~=. 
drop yhat
reg y x1 x2 x3 x4 x5 x6 dind1-dind86 if year==94&kt>0&L~=.
predict yhat if year==94&L~=.
replace lnkt = yhat if year==94&L~=. 
drop yhat
reg y x1 x2 x3 x4 x5 x6 dind1-dind86 if year==95&kt>0&L~=.
predict yhat if year==95&L~=.
replace lnkt = yhat if year==95&L~=. 
drop yhat
reg y x1 x2 x3 x4 x5 x6 dind1-dind86 if year==96&kt>0&L~=.
predict yhat if year==96&L~=.
replace lnkt = yhat if year==96&L~=. 
drop yhat
replace y = ln(kb)
reg y x1 x2 x3 x4 x5 x6 dind1-dind86 if year==92&kb>0&L~=.
predict yhat if year==92&L~=.
gen lnkb = yhat if year==92&L~=. 
drop yhat
reg y x1 x2 x3 x4 x5 x6 dind1-dind86 if year==93&kb>0&L~=.
predict yhat if year==93&L~=.
replace lnkb = yhat if year==93&L~=. 
drop yhat
reg y x1 x2 x3 x4 x5 x6 dind1-dind86 if year==94&kb>0&L~=.
predict yhat if year==94&L~=.
replace lnkb = yhat if year==94&L~=. 
drop yhat
reg y x1 x2 x3 x4 x5 x6 dind1-dind86 if year==95&kb>0&L~=.
predict yhat if year==95&L~=.
replace lnkb = yhat if year==95&L~=. 
drop yhat
reg y x1 x2 x3 x4 x5 x6 dind1-dind86 if year==96&kb>0&L~=.
predict yhat if year==96&L~=.
replace lnkb = yhat if year==96&L~=. 
drop yhat

drop xx1 xx2 xx4 xx49 xx52 xx109-xx129 xx131-xx138 xx157 xx159 xx160 xx161 xx165-xx167 xx37-xx43 xx45 xx75 xx77 xx80 xx82 xx52 xx54 xx56 xx58 xx60 xx62 xx64 xx66 xx68 xx70 xx72 xx73 xx13-xx32 xx80

gen cu = 312
replace cu = 313 if ciiu>=3130
replace cu = 314 if ciiu>=3140
replace cu = 321 if ciiu>=3210
replace cu = 322 if ciiu>=3220
replace cu = 323 if ciiu>=3230
replace cu = 324 if ciiu>=3240
replace cu = 331 if ciiu>=3310
replace cu = 332 if ciiu>=3320
replace cu = 341 if ciiu>=3410
replace cu = 342 if ciiu>=3420
replace cu = 351 if ciiu>=3510
replace cu = 352 if ciiu>=3520
replace cu = 353 if ciiu>=3530
replace cu = 354 if ciiu>=3540
replace cu = 355 if ciiu>=3550
replace cu = 356 if ciiu>=3560
replace cu = 361 if ciiu>=3610
replace cu = 362 if ciiu>=3620
replace cu = 369 if ciiu>=3690
replace cu = 371 if ciiu>=3710
replace cu = 372 if ciiu>=3720
replace cu = 381 if ciiu>=3810
replace cu = 382 if ciiu>=3820
replace cu = 383 if ciiu>=3830
replace cu = 384 if ciiu>=3840
replace cu = 385 if ciiu>=3850
replace cu = 390 if ciiu>=3900

 

save "C:\Documents and Settings\hkasahar\My Documents\work\project\exportimport\data\DATA8696.dta", replace
sort year
merge year using "C:\Documents and Settings\hkasahar\My Documents\work\project\exportimport\data\deflate.dta"
drop _merge

replace rent = rent/defcons
gen rim = im/defmach
gen rit = it/defmach
gen rib = ib/defcons
replace rim = 0 if L==.
replace rit = 0 if L==.
replace rib = 0 if L==.

gen padrony=year*0.01+padron
sort padrony
save "C:\Documents and Settings\hkasahar\My Documents\work\project\exportimport\data\DATA8696.dta", replace
* Construction of Capital Stock Variables -- the beginning of the periods
* Notes: depreciation rate == 0.1 for machinery & 0.2 for vehicle

gen exit = 0
gen entry = 0


keep export lnkm lnkt lnkb rent fta L_wh L_bl w_wh w_bl Y M Mi Md E rim rit rib im it km kt kb ciiu cu padron year exit entry defmach_f defcons_f L ivy_i ivy_f ivm_i ivm_f
reshape wide export lnkm lnkt lnkb rent fta L_wh L_bl w_wh w_bl Y M Mi Md E rim rit rib im it km kt kb ciiu cu exit entry defmach_f defcons_f L ivy_i ivy_f ivm_i ivm_f, i(padron) j(year)
 
replace exit87 =1 if L86~=.&L87==.
replace exit88 =1 if L87~=.&L88==.
replace exit89 =1 if L88~=.&L89==.
replace exit90 =1 if L89~=.&L90==.
replace exit91 =1 if L90~=.&L91==.
replace exit92 =1 if L91~=.&L92==.
replace exit93 =1 if L92~=.&L93==.
replace exit94 =1 if L93~=.&L94==.
replace exit95 =1 if L94~=.&L95==.
replace exit96 =1 if L95~=.&L96==.

replace entry87 =1 if L86==.&L87~=.
replace entry88 =1 if L86==.&L87==.&L88~=.
replace entry89 =1 if L86==.&L87==.&L88==.&L89~=.
replace entry90 =1 if L86==.&L87==.&L88==.&L89==.&L90~=.
replace entry91 =1 if L86==.&L87==.&L88==.&L89==.&L90==.&L91~=.
replace entry92 =1 if L86==.&L87==.&L88==.&L89==.&L90==.&L91==.&L92~=.
replace entry93 =1 if L86==.&L87==.&L88==.&L89==.&L90==.&L91==.&L92==.&L93~=.
replace entry94 =1 if L86==.&L87==.&L88==.&L89==.&L90==.&L91==.&L92==.&L93==.&L94~=.
replace entry95 =1 if L86==.&L87==.&L88==.&L89==.&L90==.&L91==.&L92==.&L93==.&L94==.&L95~=.
replace entry96 =1 if L86==.&L87==.&L88==.&L89==.&L90==.&L91==.&L92==.&L93==.&L94==.&L95==.&L96~=.

replace km92 = exp(lnkm92) if (km92 == 0)&L92~=.
replace kt92 = exp(lnkt92) if (kt92 == 0)&L92~=.
replace kb92 = exp(lnkb92) if (kb92 == 0)&(rent92==0)&L92~=.
replace km93 = exp(lnkm93) if (km93 == 0)&L93~=.
replace kt93 = exp(lnkt93) if (kt93 == 0)&L93~=.
replace kb93 = exp(lnkb93) if (kb93 == 0)&(rent93==0)&L93~=.
replace km94 = exp(lnkm94) if (km94 == 0)&L94~=.
replace kt94 = exp(lnkt94) if (kt94 == 0)&L94~=.
replace kb94 = exp(lnkb94) if (kb94 == 0)&(rent94==0)&L94~=.
replace km95 = exp(lnkm95) if (km95 == 0)&L95~=.
replace kt95 = exp(lnkt95) if (kt95 == 0)&L95~=.
replace kb95 = exp(lnkb95) if (kb95 == 0)&(rent95==0)&L95~=.
replace km96 = exp(lnkm96) if (km96 == 0)&L96~=.
replace kt96 = exp(lnkt96) if (kt96 == 0)&L96~=.
replace kb96 = exp(lnkb96) if (kb96 == 0)&(rent96==0)&L96~=.

gen rkm92 = km92/defmach_f92 if L92~=.
gen rkt92 = kt92/defmach_f92 if L92~=. 
gen rkb92 = kb92/defcons_f92 if L92~=.
gen rkm91 = (rkm92-rim91)/0.9 if L92~=.
gen rkt91 = (rkt92-rit91)/0.8 if L92~=.
gen rkb91 = (rkb92-rib91)/0.95 if L92~=.
replace rkm91 = 0 if rkm91<0&L92~=.
replace rkt91 = 0 if rkt91<0&L92~=.
replace rkb91 = 0 if rkb91<0&L92~=.
gen rkm90 = (rkm91-rim90)/0.9 if L92~=.
gen rkt90 = (rkt91-rit90)/0.8 if L92~=.
gen rkb90 = (rkb91-rib90)/0.95 if L92~=.
replace rkm90 = 0 if rkm90<0&L92~=.
replace rkt90 = 0 if rkt90<0&L92~=.
replace rkb90 = 0 if rkb90<0&L92~=.
gen rkm89 = (rkm90-rim89)/0.9 if L92~=.
gen rkt89 = (rkt90-rit89)/0.8 if L92~=.
gen rkb89 = (rkb90-rib89)/0.95 if L92~=.
replace rkm89 = 0 if rkm89<0&L92~=.
replace rkt89 = 0 if rkt89<0&L92~=.
replace rkb89 = 0 if rkb89<0&L92~=.
gen rkm93 = 0.9*rkm92+rim92 if L92~=.
gen rkt93 = 0.8*rkt92+rit92 if L92~=.
gen rkb93 = 0.95*rkb92+rib92 if L92~=.
replace rkm93 = km93/defmach_f93 if L92==.&L93~=.
replace rkt93 = kt93/defmach_f93 if L92==.&L93~=. 
replace rkb93 = kb93/defcons_f93 if L92==.&L93~=.
replace rkm93 = 0 if rkm93<0
replace rkt93 = 0 if rkt93<0
replace rkb93 = 0 if rkb93<0
gen rkm94 = 0.9*rkm93+rim93 if L92~=.|L93~=.
gen rkt94 = 0.8*rkt93+rit93 if L92~=.|L93~=.
gen rkb94 = 0.95*rkb93+rib93 if L92~=.|L93~=.
replace rkm94 = km94/defmach_f94 if L92==.&L93==.&L94~=.
replace rkt94 = kt94/defmach_f94 if L92==.&L93==.&L94~=. 
replace rkb94 = kb94/defcons_f94 if L92==.&L93==.&L94~=.
replace rkm94 = 0 if rkm94<0
replace rkt94 = 0 if rkt94<0
replace rkb94 = 0 if rkb94<0
gen rkm95 = 0.9*rkm94+rim95 if L92~=.|L93~=.|L94~=.
gen rkt95 = 0.8*rkt94+rit95 if L92~=.|L93~=.|L94~=.
gen rkb95 = 0.95*rkb94+rib95 if L92~=.|L93~=.|L94~=.
replace rkm95 = km95/defmach_f95 if L92==.&L93==.&L94==.&L95~=.
replace rkt95 = kt95/defmach_f95 if L92==.&L93==.&L94==.&L95~=.
replace rkb95 = kb95/defcons_f95 if L92==.&L93==.&L94==.&L95~=.
replace rkm95 = 0 if rkm95<0
replace rkt95 = 0 if rkt95<0
replace rkb95 = 0 if rkb95<0
gen rkm96 = 0.9*rkm95+rim96 if L92~=.|L93~=.|L94~=.|L95~=.
gen rkt96 = 0.8*rkt95+rit96 if L92~=.|L93~=.|L94~=.|L95~=.
gen rkb96 = 0.95*rkb95+rib96 if L92~=.|L93~=.|L94~=.|L95~=.
replace rkm96 = km96/defmach_f96 if L92==.&L93==.&L94==.&L95==.&L96~=.
replace rkt96 = kt96/defmach_f96 if L92==.&L93==.&L94==.&L95==.&L96~=.
replace rkb96 = kb96/defcons_f96 if L92==.&L93==.&L94==.&L95==.&L96~=.
replace rkm96 = 0 if rkm96<0
replace rkt96 = 0 if rkt96<0
replace rkb96 = 0 if rkb96<0

reshape long export rkm rkt rkb lnkm lnkt lnkb rent fta L_wh L_bl w_wh w_bl Y M Mi Md E rim rit rib im it km kt kb ciiu cu exit entry defmach_f defcons_f L ivy_i ivy_f ivm_i ivm_f, i(padron) j(year)

drop lnkm lnkt lnkt
drop if year<89|L==.

gen padrony=padron+year*0.01
drop if padron==.
sort padrony
save "C:\Documents and Settings\hkasahar\My Documents\work\project\exportimport\data\DATA8996_temp.dta", replace

use "C:\Documents and Settings\hkasahar\My Documents\work\project\exportimport\data\DATA8996_temp.dta"
sort year
merge year using "C:\Documents and Settings\hkasahar\My Documents\work\project\exportimport\data\deflate.dta"
drop _merge
sort year

drop if year<89

sort padron
by padron: egen mcu = mode(cu), minmode


gen VA = Y-M-E

gen rVA = VA/y312 + (ivy_f/y312_f) - (ivy_i/y312_i)
replace rVA =  VA/y313 + (ivy_f/y313_f) - (ivy_i/y313_i) if mcu==313
replace rVA =  VA/y314 + (ivy_f/y314_f) - (ivy_i/y314_i) if mcu==314
replace rVA =  VA/y321 + (ivy_f/y321_f) - (ivy_i/y321_i) if mcu==321
replace rVA =  VA/y322 + (ivy_f/y322_f) - (ivy_i/y322_i) if mcu==322
replace rVA =  VA/y323 + (ivy_f/y323_f) - (ivy_i/y323_i) if mcu==323
replace rVA =  VA/y324 + (ivy_f/y324_f) - (ivy_i/y324_i) if mcu==324
replace rVA =  VA/y331 + (ivy_f/y331_f) - (ivy_i/y331_i) if mcu==331
replace rVA =  VA/y332 + (ivy_f/y332_f) - (ivy_i/y332_i) if mcu==332
replace rVA =  VA/y341 + (ivy_f/y341_f) - (ivy_i/y341_i) if mcu==341
replace rVA =  VA/y342 + (ivy_f/y342_f) - (ivy_i/y342_i) if mcu==342
replace rVA =  VA/y351 + (ivy_f/y351_f) - (ivy_i/y351_i) if mcu==351
replace rVA =  VA/y352 + (ivy_f/y352_f) - (ivy_i/y352_i) if mcu==352
replace rVA =  VA/y353 + (ivy_f/y353_f) - (ivy_i/y353_i) if mcu==353
replace rVA =  VA/y354 + (ivy_f/y354_f) - (ivy_i/y354_i) if mcu==354
replace rVA =  VA/y355 + (ivy_f/y355_f) - (ivy_i/y355_i) if mcu==355
replace rVA =  VA/y356 + (ivy_f/y356_f) - (ivy_i/y356_i) if mcu==356
replace rVA =  VA/y361 + (ivy_f/y361_f) - (ivy_i/y361_i) if mcu==361
replace rVA =  VA/y362 + (ivy_f/y362_f) - (ivy_i/y362_i) if mcu==362
replace rVA =  VA/y369 + (ivy_f/y369_f) - (ivy_i/y369_i) if mcu==369
replace rVA =  VA/y371 + (ivy_f/y371_f) - (ivy_i/y371_i) if mcu==371
replace rVA =  VA/y372 + (ivy_f/y372_f) - (ivy_i/y372_i) if mcu==372
replace rVA =  VA/y381 + (ivy_f/y381_f) - (ivy_i/y381_i) if mcu==381
replace rVA =  VA/y382 + (ivy_f/y382_f) - (ivy_i/y382_i) if mcu==382
replace rVA =  VA/y383 + (ivy_f/y383_f) - (ivy_i/y383_i) if mcu==383
replace rVA =  VA/y384 + (ivy_f/y384_f) - (ivy_i/y384_i) if mcu==384
replace rVA =  VA/y385 + (ivy_f/y385_f) - (ivy_i/y385_i) if mcu==385
replace rVA =  VA/y390 + (ivy_f/y390_f) - (ivy_i/y390_i) if mcu==390

gen rY = Y/y312 + (ivy_f/y312_f) - (ivy_i/y312_i)
replace rY =  Y/y313 + (ivy_f/y313_f) - (ivy_i/y313_i) if mcu==313
replace rY =  Y/y314 + (ivy_f/y314_f) - (ivy_i/y314_i) if mcu==314
replace rY =  Y/y321 + (ivy_f/y321_f) - (ivy_i/y321_i) if mcu==321
replace rY =  Y/y322 + (ivy_f/y322_f) - (ivy_i/y322_i) if mcu==322
replace rY =  Y/y323 + (ivy_f/y323_f) - (ivy_i/y323_i) if mcu==323
replace rY =  Y/y324 + (ivy_f/y324_f) - (ivy_i/y324_i) if mcu==324
replace rY =  Y/y331 + (ivy_f/y331_f) - (ivy_i/y331_i) if mcu==331
replace rY =  Y/y332 + (ivy_f/y332_f) - (ivy_i/y332_i) if mcu==332
replace rY =  Y/y341 + (ivy_f/y341_f) - (ivy_i/y341_i) if mcu==341
replace rY =  Y/y342 + (ivy_f/y342_f) - (ivy_i/y342_i) if mcu==342
replace rY =  Y/y351 + (ivy_f/y351_f) - (ivy_i/y351_i) if mcu==351
replace rY =  Y/y352 + (ivy_f/y352_f) - (ivy_i/y352_i) if mcu==352
replace rY =  Y/y353 + (ivy_f/y353_f) - (ivy_i/y353_i) if mcu==353
replace rY =  Y/y354 + (ivy_f/y354_f) - (ivy_i/y354_i) if mcu==354
replace rY =  Y/y355 + (ivy_f/y355_f) - (ivy_i/y355_i) if mcu==355
replace rY =  Y/y356 + (ivy_f/y356_f) - (ivy_i/y356_i) if mcu==356
replace rY =  Y/y361 + (ivy_f/y361_f) - (ivy_i/y361_i) if mcu==361
replace rY =  Y/y362 + (ivy_f/y362_f) - (ivy_i/y362_i) if mcu==362
replace rY =  Y/y369 + (ivy_f/y369_f) - (ivy_i/y369_i) if mcu==369
replace rY =  Y/y371 + (ivy_f/y371_f) - (ivy_i/y371_i) if mcu==371
replace rY =  Y/y372 + (ivy_f/y372_f) - (ivy_i/y372_i) if mcu==372
replace rY =  Y/y381 + (ivy_f/y381_f) - (ivy_i/y381_i) if mcu==381
replace rY =  Y/y382 + (ivy_f/y382_f) - (ivy_i/y382_i) if mcu==382
replace rY =  Y/y383 + (ivy_f/y383_f) - (ivy_i/y383_i) if mcu==383
replace rY =  Y/y384 + (ivy_f/y384_f) - (ivy_i/y384_i) if mcu==384
replace rY =  Y/y385 + (ivy_f/y385_f) - (ivy_i/y385_i) if mcu==385
replace rY =  Y/y390 + (ivy_f/y390_f) - (ivy_i/y390_i) if mcu==390

replace rY = Y/defmanuf

gen rE = E/e312 
replace rE =  E/e313  if mcu==313
replace rE =  E/e314  if mcu==314
replace rE =  E/e321  if mcu==321
replace rE =  E/e322  if mcu==322
replace rE =  E/e323  if mcu==323
replace rE =  E/e324  if mcu==324
replace rE =  E/e331  if mcu==331
replace rE =  E/e332  if mcu==332
replace rE =  E/e341  if mcu==341
replace rE =  E/e342  if mcu==342
replace rE =  E/e351  if mcu==351
replace rE =  E/e352  if mcu==352
replace rE =  E/e353  if mcu==353
replace rE =  E/e354  if mcu==354
replace rE =  E/e355  if mcu==355
replace rE =  E/e356  if mcu==356
replace rE =  E/e361  if mcu==361
replace rE =  E/e362  if mcu==362
replace rE =  E/e369  if mcu==369
replace rE =  E/e371  if mcu==371
replace rE =  E/e372  if mcu==372
replace rE =  E/e381  if mcu==381
replace rE =  E/e382  if mcu==382
replace rE =  E/e383  if mcu==383
replace rE =  E/e384  if mcu==384
replace rE =  E/e385  if mcu==385
replace rE =  E/e390  if mcu==390

gen rMi = Mi/m312 
replace rMi =  Mi/m313  if mcu==313
replace rMi =  Mi/m314  if mcu==314
replace rMi =  Mi/m321  if mcu==321
replace rMi =  Mi/m322  if mcu==322
replace rMi =  Mi/m323  if mcu==323
replace rMi =  Mi/m324  if mcu==324
replace rMi =  Mi/m331  if mcu==331
replace rMi =  Mi/m332  if mcu==332
replace rMi =  Mi/m341  if mcu==341
replace rMi =  Mi/m342  if mcu==342
replace rMi =  Mi/m351  if mcu==351
replace rMi =  Mi/m352  if mcu==352
replace rMi =  Mi/m353  if mcu==353
replace rMi =  Mi/m354  if mcu==354
replace rMi =  Mi/m355  if mcu==355
replace rMi =  Mi/m356  if mcu==356
replace rMi =  Mi/m361  if mcu==361
replace rMi =  Mi/m362  if mcu==362
replace rMi =  Mi/m369  if mcu==369
replace rMi =  Mi/m371  if mcu==371
replace rMi =  Mi/m372  if mcu==372
replace rMi =  Mi/m381  if mcu==381
replace rMi =  Mi/m382  if mcu==382
replace rMi =  Mi/m383  if mcu==383
replace rMi =  Mi/m384  if mcu==384
replace rMi =  Mi/m385  if mcu==385
replace rMi =  Mi/m390  if mcu==390

replace rMi = Mi/defimport

gen rMd = Md/m312
replace rMd =  Md/m313 if mcu==313
replace rMd =  Md/m314 if mcu==314
replace rMd =  Md/m321 if mcu==321
replace rMd =  Md/m322 if mcu==322
replace rMd =  Md/m323 if mcu==323
replace rMd =  Md/m324 if mcu==324
replace rMd =  Md/m331 if mcu==331
replace rMd =  Md/m332 if mcu==332
replace rMd =  Md/m341 if mcu==341
replace rMd =  Md/m342 if mcu==342
replace rMd =  Md/m351 if mcu==351
replace rMd =  Md/m352 if mcu==352
replace rMd =  Md/m353 if mcu==353
replace rMd =  Md/m354 if mcu==354
replace rMd =  Md/m355 if mcu==355
replace rMd =  Md/m356 if mcu==356
replace rMd =  Md/m361 if mcu==361
replace rMd =  Md/m362 if mcu==362
replace rMd =  Md/m369 if mcu==369
replace rMd =  Md/m371 if mcu==371
replace rMd =  Md/m372 if mcu==372
replace rMd =  Md/m381 if mcu==381
replace rMd =  Md/m382 if mcu==382
replace rMd =  Md/m383 if mcu==383
replace rMd =  Md/m384 if mcu==384
replace rMd =  Md/m385 if mcu==385
replace rMd =  Md/m390 if mcu==390

*replace w_wh = w_wh/defmanuf
*replace w_bl = w_bl/defmanuf

gen rexport = export/defmanuf
replace fta = fta/defmanuf

gen rM = M/defmanuf
replace rMi = Mi/defmanuf

drop if year<89

***** Construct missing book values of capital stock *****
gen x1 = ln(L_wh) if L~=.
gen x2 = ln(L_bl) if L~=.
gen x3 = ln(rM) if L~=.
gen x4 = ln(rE) if L~=.
gen x5 = 0 if L~=.
replace x5 = 1 if Mi>0&L~=.
gen x6 = ln(rY) if L~=.
replace x1 = 1 if L_wh<=0&L~=.
replace x2 = 1 if L_bl<=0&L~=.
replace x3 = 1 if M<=0&L~=.
replace x4 = 1 if E<=0&L~=.
replace x6 = 1 if Y<=0&L~=.
tab ciiu, gen(dind)
gen y = ln(rkm)
reg y x1 x2 x3 x4 x5 x6 dind1-dind86 if year==91&rkm>0&L~=.
predict yhat if year==91&L~=.
gen lnkm = yhat if year==91&L~=. 
drop yhat
reg y x1 x2 x3 x4 x5 x6 dind1-dind86 if year==90&rkm>0&L~=.
predict yhat if year==90&L~=.
replace lnkm = yhat if year==90&L~=. 
drop yhat
reg y x1 x2 x3 x4 x5 x6 dind1-dind86 if year==89&rkm>0&L~=.
predict yhat if year==89&L~=.
replace lnkm = yhat if year==89&L~=. 
drop yhat
replace y = ln(rkt)
reg y x1 x2 x3 x4 x5 x6 dind1-dind86 if year==91&rkt>0&L~=.
predict yhat if year==91&L~=.
gen lnkt = yhat if year==91&L~=. 
drop yhat
reg y x1 x2 x3 x4 x5 x6 dind1-dind86 if year==90&rkt>0&L~=.
predict yhat if year==90&L~=.
replace lnkt = yhat if year==90&L~=. 
drop yhat
reg y x1 x2 x3 x4 x5 x6 dind1-dind86 if year==89&rkt>0&L~=.
predict yhat if year==89&L~=.
replace lnkt = yhat if year==89&L~=. 
drop yhat
replace y = ln(rkb)
reg y x1 x2 x3 x4 x5 x6 dind1-dind86 if year==91&rkb>0&L~=.
predict yhat if year==91&L~=.
replace lnkb = yhat if year==91&L~=. 
drop yhat
reg y x1 x2 x3 x4 x5 x6 dind1-dind86 if year==90&rkb>0&L~=.
predict yhat if year==90&L~=.
replace lnkb = yhat if year==90&L~=. 
drop yhat
reg y x1 x2 x3 x4 x5 x6 dind1-dind86 if year==89&rkb>0&L~=.
predict yhat if year==89&L~=.
replace lnkb = yhat if year==89&L~=. 
drop yhat

gen reentry = 0

keep mcu reentry export rexport lnkm lnkt lnkb rent fta L_wh L_bl w_wh w_bl rVA rY rM rMi rMd rE Y M Mi Md E rim rit rib im it rkm rkt rkb ciiu cu padron year exit entry defmach_f defcons_f L
reshape wide reentry export rexport lnkm lnkt lnkb rent fta L_wh L_bl w_wh w_bl rVA rY rM rMi rMd rE Y M Mi Md E rim rit rib im it rkm rkt rkb ciiu cu  exit entry defmach_f defcons_f L, i(padron) j(year)

replace rkm91 = exp(lnkm91) if L92==.&L91~=.
replace rkt91 = exp(lnkt91) if L92==.&L91~=.
replace rkb91 = exp(lnkb91) if L92==.&L91~=.
replace rkm91 = 0 if rkm91<0
replace rkt91 = 0 if rkt91<0
replace rkb91 = 0 if rkb91<0
replace rkm90 = (rkm91-rim90)/0.9
replace rkt90 = (rkt91-rit90)/0.8
replace rkb90 = (rkb91-rib90)/0.95
replace rkm90 = 0 if rkm90<0
replace rkt90 = 0 if rkt90<0
replace rkb90 = 0 if rkb90<0
replace rkm90 = exp(lnkm90) if L92==.&L91==.&L90~=.
replace rkt90 = exp(lnkt90) if L92==.&L91==.&L90~=.
replace rkb90 = exp(lnkb90) if L92==.&L91==.&L90~=.
replace rkm89 = (rkm90-rim89)/0.9
replace rkt89 = (rkt90-rit89)/0.8
replace rkb89 = (rkb90-rib89)/0.95
replace rkm89 = 0 if rkm89<0
replace rkt89 = 0 if rkt89<0
replace rkb89 = 0 if rkb89<0
replace rkm89 = exp(lnkm89) if L92==.&L91==.&L90==.&L89~=.
replace rkt89 = exp(lnkt89) if L92==.&L91==.&L90==.&L89~=.
replace rkb89 = exp(lnkb89) if L92==.&L91==.&L90==.&L89~=.

gen dplnt=0
replace dplnt=1 if L90==.&L91==.&L92==.&L93==.&L94==.&L95==.&L96==.


replace exit91 =1 if L90~=.&L91==.
replace exit92 =1 if L91~=.&L92==.
replace exit93 =1 if L92~=.&L93==.
replace exit94 =1 if L93~=.&L94==.
replace exit95 =1 if L94~=.&L95==.
replace exit96 =1 if L95~=.&L96==.


replace reentry91 = 1 if exit91==1&L92~=.
replace reentry91 = 1 if exit91==1&L93~=.
replace reentry91 = 1 if exit91==1&L94~=.
replace reentry91 = 1 if exit91==1&L95~=.
replace reentry91 = 1 if exit91==1&L96~=. 
replace reentry92 = 1 if exit92==1&L93~=.
replace reentry92 = 1 if exit92==1&L94~=.
replace reentry92 = 1 if exit92==1&L95~=.
replace reentry92 = 1 if exit92==1&L96~=.
replace reentry93 = 1 if exit93==1&L94~=.
replace reentry93 = 1 if exit93==1&L95~=.
replace reentry93 = 1 if exit93==1&L96~=. 
replace reentry94 = 1 if exit94==1&L95~=.
replace reentry94 = 1 if exit94==1&L96~=. 
replace reentry95 = 1 if exit95==1&L96~=. 

reshape long reentry export rexport lnkm lnkt lnkb rent fta L_wh L_bl w_wh w_bl rVA rY rM rMi rMd rE Y M Mi Md E rim rit rib im it rkm rkt rkb ciiu cu exit entry defmach_f defcons_f L, i(padron) j(year)


*drop if L==.
gen rk = rkm + rkt
gen d_m = 0
replace d_m = 1 if (rMi~=.)&(rMi>0)
replace d_m = . if L==.|L==0
gen d_e = 0
replace d_e = 1 if export>0
replace d_e = . if L==.|L==0
replace d_e = . if year<90
replace M = . if L==.|L==0
gen cuy = year*0.01+mcu
replace cuy = . if L==.|L==0 
replace export = . if year<90
sort cuy
by cuy: egen d_m_ind = mean(d_m)
by cuy: egen d_e_ind = mean(d_e)
by cuy: egen sum_M_ind = sum(M)
by cuy: egen sum_Mi_ind = sum(Mi) 
by cuy: egen sum_Y_ind = sum(Y)
by cuy: egen sum_E_ind = sum(export) 
gen Mi_M_ind = sum_Mi_ind/sum_M_ind
gen E_Y_ind = sum_E_ind/sum_Y_ind


sort year
by year: egen wage_wh = mean(w_wh)
by year: egen wage_bl = mean(w_bl)
by year: egen wage_wh_md = median(w_wh)
by year: egen wage_bl_md = median(w_bl)

 
replace w_wh = 0 if w_wh==.
replace w_bl = 0 if w_bl==.
gen vc = w_wh*L_wh + w_bl*L_bl + M + E if L~=.
gen profit = Y - vc  & L~=. 
gen sigma = Y/profit if L~=.
gen markup = profit/Y if L~=.
 
 

sort padron year

save "C:\Documents and Settings\hkasahar\My Documents\work\project\exportimport\data\DATA_KL_EXIM5.dta", replace

stop




clear 
insheet using "C:\Documents and Settings\hkasahar\My Documents\work\project\exportimport\data\rer.txt"
sort year
save "C:\Documents and Settings\hkasahar\My Documents\work\project\exportimport\data\rer.dta",replace

clear 
insheet using "C:\Documents and Settings\hkasahar\My Documents\work\project\exportimport\data\tfp.txt"
sort year
save "C:\Documents and Settings\hkasahar\My Documents\work\project\exportimport\data\tfp.dta",replace


clear 

use "C:\Documents and Settings\hkasahar\My Documents\work\project\exportimport\data\rer.dta"
gen lnrer = ln(rer)
gen lnrer_1 = ln(rer_1)
reg lnrer lnrer_1
predict r, resid
summarize r


clear
use "C:\Documents and Settings\hkasahar\My Documents\work\project\exportimport\data\tfp.dta"
gen time = year-1960
reg tfp tfp_1, noconst
predict r, resid
summarize r



sort year
by year: summarize rk rkb if L>0, detail

drop rM cuy statey d_m d_e d_m_ind d_e_ind sum_M_ind sum_Mi_ind sum_Y_ind sum_E_ind Mi_M_ind E_Y_ind 
drop d_m_st d_e_st sum_M_st sum_Mi_st sum_Y_st sum_E_st Mi_M_st E_Y_st

sort year
by year: summarize d_m_ind, detail
by year: summarize d_e_ind, detail
by year: summarize Mi_M_ind, detail
by year: summarize E_Y_ind, detail
by year: summarize d_m_st, detail
by year: summarize d_e_st, detail
by year: summarize Mi_M_st, detail
by year: summarize E_Y_st, detail

gen Y = ln(rY)
gen l_s = ln(L_wh)
gen l_u = ln(L_bl)
gen k = ln(rk+rkb)
gen e = ln(rE)
gen m = ln(rMi+rMd)
gen d = 0
replace d = 1 if (rMi>0)&(rMi~=.)

reg Y l_s l_u k e m d if dplnt==0

drop if year<90
sort year
by year: summarize export 
by year: summarize export if Mi>0 
by year: summarize rMi 
by year: summarize rMi if export>0


clear
set memory 500m
use "C:\Documents and Settings\hkasahar\My Documents\work\project\exportimport\data\DATA_KL_EXIM5.dta"
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
outsheet using "C:\Documents and Settings\hkasahar\My Documents\work\project\exportimport\data\DATA_KL_EXIM5.xl", replace 




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



clear
set memory 500m
use "C:\Documents and Settings\hkasahar\My Documents\work\project\exportimport\data\DATA_KL_EXIM5.dta"
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
merge year using "C:\Documents and Settings\hkasahar\My Documents\work\project\exportimport\data\deflate9096.dta"
drop _merge
sort year
merge year using "C:\Documents and Settings\hkasahar\My Documents\work\project\exportimport\data\deflate.dta"
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
outsheet using "C:\Documents and Settings\hkasahar\My Documents\work\project\exportimport\data\DATA_KL_3220_ciiu.xl", replace





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


clear
set memory 500m
use "C:\Documents and Settings\hkasahar\My Documents\work\project\exportimport\data\DATA_KL_EXIM5.dta"
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
merge year using "C:\Documents and Settings\hkasahar\My Documents\work\project\exportimport\data\deflate9096.dta"
drop _merge
sort year
merge year using "C:\Documents and Settings\hkasahar\My Documents\work\project\exportimport\data\deflate.dta"
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
outsheet using "C:\Documents and Settings\hkasahar\My Documents\work\project\exportimport\data\DATA_KL_3560_ciiu.xl", replace

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




clear
set memory 500m
use "C:\Documents and Settings\hkasahar\My Documents\work\project\exportimport\data\DATA_KL_EXIM5.dta"
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
merge year using "C:\Documents and Settings\hkasahar\My Documents\work\project\exportimport\data\deflate9096.dta"
drop _merge
sort year
merge year using "C:\Documents and Settings\hkasahar\My Documents\work\project\exportimport\data\deflate.dta"
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
outsheet using "C:\Documents and Settings\hkasahar\My Documents\work\project\exportimport\data\DATA_KL_3813_ciiu.xl", replace

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














