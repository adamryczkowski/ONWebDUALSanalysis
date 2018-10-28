process_q48<-function(dt) {
  #groupings
  #1. Oro-facial
  #2. Limb
  #3. Upper limb
  #4. Lower limb
  #5. Respiratory
  #6. Neck
  #7. Thoracic/abdomen
  #8. Dyscognition
  groupings=list(
    "Oro-facial"=1,
    "Limb"=c(2:13),
    "Upper limb"=c(2:7),
    "Lower limb"=c(8:13),
    "Respiratory symptoms"=14,
    "Neck"=15,
    "Thoracic/abdomen"=16,
    "Dyscognition"=17
  )


#  dt<-as.data.frame(joined_df)[stringr::str_detect(colnames(joined_df), stringr::regex('^q_48._[rt]$'))]
  n_regions<-colnames(dt)[stringr::str_detect(colnames(dt), stringr::regex('^q_48._r$'))]
  n_timings<-colnames(dt)[stringr::str_detect(colnames(dt), stringr::regex('^q_48._t$'))]
  get_set<-function(rownr) {
    regions<-as.integer((dt %>% select(!!n_regions))[(rownr),])
    timings<-as.integer((dt %>% select(!!n_timings))[(rownr),])
    maxidx_r<-which.max(is.na(regions))-1
    maxidx_t<-which.max(is.na(timings))-1
    if(maxidx_r==0 && maxidx_t==0) {
      return(list(regions=integer(0),
             timings=integer(0)))
    }
    if(maxidx_r != maxidx_t +1 ) {
      browser()
    }
    timings<-timings[seq(1, maxidx_t)]
    timings[timings==2]<-0
    timings<-cumsum(c(1, timings)) #"After" increases counter, "Same time" does not. That's why we replace it with "0" in cumsum.
    return(list(
      regions=regions[seq(1, maxidx_r)],
      timings=timings
      ))
    #Function, that for a given row returns a two vectors: places and timings
  }
  sets<-purrr::map(seq(1, nrow(dt)), get_set)

  locate_grouping<-function(a_set, grouping) {
    #Returns index of the first particular grouping or NA if absent in the given set "a_set"
    ans<-which(a_set$regions %in% grouping)
    if(length(ans)>0) {
      r<-ans[[1]]
      return(as.integer(a_set$timings[[r]]))
    } else {
      return(NA)
    }
  }
  for(i in seq_along(groupings)) {
    gr_name<-names(groupings)[[i]]
    gr_list<-groupings[[i]]
    var_name<-paste0('q_48_', i)
    var<-purrr::map_int(sets, locate_grouping, grouping=gr_list)
    dt[,(var_name):=var]
    data.table::setattr(dt[[var_name]], 'label', paste0("Ordinal number of ", gr_name, " symptoms"))
  }

  return(dt)
}

