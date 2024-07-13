*** Econometrics HW2
*** May 3, 2024
*** Guanxi Li

cd "~/Desktop/Stata/Metrics/HW2"

use "~/Desktop/Stata/Metrics/HW2/Growth.dta",clear

* drop Malta
drop if country_name == "Malta"

** main regressions
reg growth tradeshare yearsschool, r
est store a1

gen lnyearsschool = ln(yearsschool)
reg growth tradeshare lnyearsschool, r
est store a2

gen lnrgdp60 = ln(rgdp60)
reg growth tradeshare lnyearsschool rev_coups assasinations lnrgdp60, r
est store a3

gen trdslnys = tradeshare*ln(yearsschool)
reg growth tradeshare lnyearsschool rev_coups assasinations lnrgdp60 trdslnys, r
est store a4

gen tradeshare2 = tradeshare^2
gen tradeshare3 = tradeshare^3
reg growth tradeshare tradeshare2 tradeshare3 lnyearsschool rev_coups assasinations lnrgdp60, r
est store a5

* F test
test tradeshare2 = tradeshare3 = 0

outreg2 [a1 a2 a3 a4 a5] using HW2.tex, replace keep(tradeshare lnyearsschool rev_coups assasinations lnrgdp60 yearsschool trdslnys tradeshare2 tradeshare3)

** scatter
graph twoway ///
    (lfit growth yearsschool, lpattern(dash) lcolor(blue)) ///
    (scatter growth yearsschool, mcolor(red) msize(small)) ///
    , ///
    title(" " " ") ///
    xtitle("Years of Schooling") ytitle("Growth") ///
    legend(on)
	
graph twoway ///
    (lfit growth lnyearsschool, lpattern(dash) lcolor(blue)) ///
    (scatter growth lnyearsschool, mcolor(red) msize(small)) ///
    , ///
    title(" " " ") ///
    xtitle("ln(Years of Schooling)") ytitle("Growth") ///
    legend(on)




