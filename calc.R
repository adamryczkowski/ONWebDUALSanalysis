#devtools::install_github('adamryczkowski/ONWebDUALSanalysis')
library(ONWebDUALSanalysis)
library(doMC)
registerDoMC(2)
#debugonce(make_rap)
make_rap(dv_nr = 5, rap_path = 'reports')

for(i in 1:6) do_calc(i)
