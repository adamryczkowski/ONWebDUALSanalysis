gen_geom_series<-function(n, start, end, steepness=0.7) {
  offset<-(-1/steepness + 2)*start
  s<-offset + exp(seq(from = log(start-offset), by = log((end-offset)/(start-offset))/(n-1), length.out = n))
  return(s)
}

gen_int_geom_series<-function(n, count, steepness=0.7) {
  unique(c(round(gen_geom_series(ncol(ads)-1, 1, ncol(ads)-1))), count)
}

