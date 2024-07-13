cap log close
 log using "~/Desktop/Stata/Metrics/HW1/HW.smcl", replace
clear
set more off

cd "~/Desktop/Stata/Metrics/HW1"

use "~/Desktop/Stata/Metrics/HW1/cps08.dta",clear

reg ahe age, r
est store a1

reg ahe age if bachelor == 0
est store a2

reg ahe age if bachelor == 1
est store a3

outreg2 [a1 a2 a3] using HW1.tex, replace
