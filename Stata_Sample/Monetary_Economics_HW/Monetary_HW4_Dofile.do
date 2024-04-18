cd "/Users/liguanxi/Dropbox/00 大学/大三上/课程（大三上）/货币经济学/货币经济学 作业/货币经济学HW4/货币经济学HW4_Data"
clear
insheet using 4-1-部分实证数据-KHY.csv
browse


tsset year
gen lnm2=ln(m2)
gen lngdp=ln(gdp)
gen lninf=ln(inf)
gen M = m


* Specify the max n of lag
di 12 * 0.35^0.25



* check the constant and term of time trend
twoway(line m2 year)
twoway(line gdp year)
twoway(line i year)
twoway(line inf year)
twoway(line M year)
twoway(line lnm2 year)
twoway(line lngdp year)
twoway(line lninf year)


* ADF
dfuller lnm2, lags(9) trend
dfuller lnm2, lags(8) trend
dfuller lnm2, lags(7) trend
dfuller lnm2, lags(6) trend
dfuller lnm2, lags(5) trend
dfuller lnm2, lags(4) trend
dfuller lnm2, lags(3) trend
dfuller lnm2, lags(2) trend
dfuller lnm2, lags(1) trend
dfuller lnm2, trend



dfuller lngdp, lags(9) trend 
dfuller lngdp, lags(8) trend
dfuller lngdp, lags(7) trend
dfuller lngdp, lags(6) trend
dfuller lngdp, lags(5) trend
dfuller lngdp, lags(4) trend
dfuller lngdp, lags(3) trend
dfuller lngdp, lags(2) trend
dfuller lngdp, lags(1) trend
dfuller lngdp, trend


dfuller i, lags(9) reg
dfuller i, lags(8) reg
dfuller i, lags(7) reg
dfuller i, lags(6) reg
dfuller i, lags(5) reg
dfuller i, lags(4) reg
dfuller i, lags(3) reg
dfuller i, lags(2) reg
dfuller i, lags(1) reg
adf i, reg


dfuller lninf, lags(9) reg
dfuller lninf, lags(8) reg
dfuller lninf, lags(7) reg
dfuller lninf, lags(6) reg
dfuller lninf, lags(5) reg
dfuller lninf, lags(4) reg
dfuller lninf, lags(3) reg
dfuller lninf, lags(2) reg
dfuller lninf, lags(1) reg
dfuller lninf, reg


dfuller M, lags(9) trend reg
dfuller M, lags(8) trend reg
dfuller M, lags(7) trend reg
dfuller M, lags(6) trend reg
dfuller M, lags(5) trend reg
dfuller M, lags(4) trend reg
dfuller M, lags(3) trend reg
dfuller M, lags(2) trend reg
dfuller M, lags(1) trend reg
dfuller M, reg


dfuller d.lnm2, lags(9) trend reg
dfuller d.lnm2, lags(8) trend reg
dfuller d.lnm2, lags(7) trend reg
dfuller d.lnm2, lags(6) trend reg
dfuller d.lnm2, lags(5) trend reg
dfuller d.lnm2, lags(4) trend reg
dfuller d.lnm2, lags(3) trend reg
dfuller d.lnm2, lags(2) trend reg
dfuller d.lnm2, lags(1) trend reg
dfuller d.lnm2, trend reg


dfuller d.lngdp, lags(9) trend reg
dfuller d.lngdp, lags(8) trend reg
dfuller d.lngdp, lags(7) trend reg
dfuller d.lngdp, lags(6) trend reg
dfuller d.lngdp, lags(5) trend reg
dfuller d.lngdp, lags(4) trend reg
dfuller d.lngdp, lags(3) trend reg
dfuller d.lngdp, lags(2) trend reg
dfuller d.lngdp, lags(1) trend reg
dfuller d.lngdp, trend reg


dfuller d.i, lags(9) reg
dfuller d.i, lags(8) reg
dfuller d.i, lags(7) reg
dfuller d.i, lags(6) reg
dfuller d.i, lags(5) reg
dfuller d.i, lags(4) reg
dfuller d.i, lags(3) reg
dfuller d.i, lags(2) reg
dfuller d.i, lags(1) reg
dfuller d.i, reg


dfuller d.lninf, lags(9) reg
dfuller d.lninf, lags(8) reg
dfuller d.lninf, lags(7) reg
dfuller d.lninf, lags(6) reg
dfuller d.lninf, lags(5) reg
dfuller d.lninf, lags(4) reg
dfuller d.lninf, lags(3) reg
dfuller d.lninf, lags(2) reg
dfuller d.lninf, lags(1) reg
dfuller d.lninf, reg


dfuller d.M, lags(9) trend reg
dfuller d.M, lags(8) trend reg
dfuller d.M, lags(7) trend reg
dfuller d.M, lags(6) trend reg
dfuller d.M, lags(5) trend reg
dfuller d.M, lags(4) trend reg
dfuller d.M, lags(3) trend reg
dfuller d.M, lags(2) trend reg
dfuller d.M, lags(1) trend reg
dfuller d.M, trend reg


* Long-term broad money demand function and residual stationarity test
reg lnm2 lngdp i, r
est store z1
reg lnm2 lngdp i lninf, r
est store z2
reg lnm2 lngdp i lninf M, r
est store z3
esttab z1 z2 z3 using regn.tex, r2 ar2 se star(* 0.1 ** 0.05 *** 0.01) replace nogap


predict e,resid
twoway(line e year)

dfuller e, lags(9) nocon
dfuller e, lags(8) nocon
dfuller e, lags(7) nocon
dfuller e, lags(6) nocon
dfuller e, lags(5) nocon
dfuller e, lags(4) nocon
dfuller e, lags(3) nocon
dfuller e, lags(2) nocon
dfuller e, lags(1) nocon
dfuller e, nocon

gen estlnm2 = -3.427 + 1.254 * lngdp - 0.0173 * i + 0.0163 * lninf + 0.348 * M
twoway(line estlnm2 year) (line lnm2 year)



* Solution to short-term broad money demand function
gen ecm=e
reg d.lnm2 d.lngdp d.i d.lninf d.M, nocon
est store z4
esttab z4 using regn.tex, r2 ar2 se star(* 0.1 ** 0.05 *** 0.01) replace nogap

predict e,resid
dfuller e, lags(9) nocon
dfuller e, lags(8) nocon
dfuller e, lags(7) nocon
dfuller e, lags(6) nocon
dfuller e, lags(5) nocon
dfuller e, lags(4) nocon
dfuller e, lags(3) nocon
dfuller e, lags(2) nocon
dfuller e, lags(1) nocon
dfuller e, nocon

gen estdlnm2 = 1.028 * d.lngdp - 0.00684 * d.i + 0.00407 * d.lninf + 1.070 * d.M
twoway(line estdlnm2 year) (line d.lnm2 year)
