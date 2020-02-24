Converting Brandon Kelly's IDL code supou.pro (described in general terms in Kelly+11) into Python

Input file: lightcurve, of format [day, flux, flux err]

2018-07-30
The code works, and produces the same output as the IDL version for data on the quasar 0238+166. It does however take about 20 minutes to run, which I think is longer than the IDL version. Can we optimize the emcee part?
