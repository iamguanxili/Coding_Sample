*** Econometrics HW2
*** June 1, 2024
*** Guanxi Li

cd "~/Desktop/Stata/Metrics/HW3"

*** P292, 2

use "~/Desktop/Stata/Metrics/HW3/SeatBelts.dta",clear

** Data cleaning

egen state_dum = group(state)

gen lnincome = ln(income)

** Main regressions

reg fatalityrate sb_useage speed65 speed70 ba08 drinkage21 lnincome age, r
est store a1

reghdfe fatalityrate sb_useage speed65 speed70 ba08 drinkage21 lnincome age, absorb(state_dum)
est store a2

reghdfe fatalityrate sb_useage speed65 speed70 ba08 drinkage21 lnincome age, absorb state i.year)

est store a3

outreg2 [a1 a2 a3] using HW3-2-1.tex, replace

reghdfe sb_useage primary secondary speed65 speed70 ba08 drinkage21 lnincome age, absorb(state year)
est store a4

outreg2 [a4] using HW3-2-2.tex, replace

** F-test

reghdfe fatalityrate sb_useage speed65 speed70 ba08 drinkage21 lnincome age i.state_dum i.year

** Mean of millions of traffic miles per year. (Note: Number of fatalities = fatalityrate \times vmt)
mean(vmt)

testparm i.year

*** P356, 2

cd "~/Desktop/Stata/Metrics/HW3"

use "~/Desktop/Stata/Metrics/HW3/fertility.dta",clear

reg weeksm1 morekids, r
est store z1

ivreg weeksm1 (morekids = samesex), first
est store z2

global X1 "agem1 black hispan othrace"

ivreg weeksm1 $X1 (morekids = samesex), first
est store z3

outreg2 [z1 z2 z3] using HW3-5.tex, replace
