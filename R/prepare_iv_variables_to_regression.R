select_variables_sep2018<-function(joined_df) {
  db<-nodupdb<-joined_df[!duplicated(joined_df$q_0),]

#1.	Bulbar (pyt. 16), upper limb, lower limb onset
#   (pyt. 16: limb onset oraz pyt. 48a: cervical=upper limb, lumbar=lower limb);
#   Jeśli pacjent w pytaniu 36a ma co innego niż 1, 2, 3 albo 4 – brak danych.
#   Jeśli w 16b jest tak, to przyjmujemy Bulbar i koniec.
#   Jeśli w 16a jest tak, to patrzymy na pytanie 49a.
#   Jeśli 49a = Cervical -> Upper Limb.
#   Jeśli 49a = Lumbo-sacral-> Lower Limb. Brak danych w innym razie.



  flags2=(db$q_16b=='Yes' & !is.na(db$q_36a)) #bulbat
  flags1=(!flags2 & db$q_16a=='Yes' & db$q_49a=='Cervical' & !is.na(db$q_36a)) #lower limb
  flags3=(!flags2 & !flags2 & db$q_16a=='Yes' & db$q_49a=='Lumbo-sacral' & !is.na(db$q_36a)) #upper limb
  ans=flags1*1 + flags2*2 + flags3*3
  db[,iv56:=factor(ifelse(ans==0, NA, ans), labels=c("Bulbar", 'Upper limb', 'Lower limb'),  levels=1:3)]
#  setattr(db$iv56, 'label', 'Fenotype classification 1')
  setattr(db$iv56, 'label', 'Fenotype classification Bulbar/Limb')

  # 2. UMN or LMN involvement on onset (pyt. 17).
  #    Jeśli pacjent w pytaniu 36a ma co innego niż 1, 2, 3 albo 4 – brak danych.
  #    Dla każdej grupy pyt 17 po jednej grupie fenotypu: UMN, LMN, UMN+LMN

  flags1=(!is.na(db$q_36a) & db$q_17=="UMN")
  flags2=(!is.na(db$q_36a) & db$q_17=="LMN")
  flags3=(!is.na(db$q_36a) & db$q_17=="UMN+LMN")
  ans=flags1*1 + flags2*2 + flags3*3
  db[,iv57:=factor(ifelse(ans==0, NA, ans), labels=c("UMN", 'LMN', 'UMN+LMN'),  levels=1:3)]
  #setattr(db$iv57, 'label', 'Fenotype classification 2')
#  db[,iv57:=factor( (db$iv57=='UMN')*1+(db$iv57=='LMN')*2, labels=c('UMN', 'LMN'), levels=c(1, 2))]
  setattr(db$iv57, 'label', 'Fenotype classification UMN/LMN')

  db[,dv1:=db$q_51_1.score]
  db[,dv2:=db$q_51_1.1 + db$q_51_1.2 + db$q_51_1.3]
  setattr(db$dv2, 'label', "ALSFRS Speech, Salivation and Swallowing")
  db[,dv3:=db$q_51_1.4 + db$q_51_1.5 + db$q_51_1.6]
  setattr(db$dv3, 'label', "ALSFRS Handwriting, Cutting Food and Hygiene")
  db[,dv4:=db$q_51_1.7 + db$q_51_1.8 + db$q_51_1.9]
  setattr(db$dv4, 'label', "ALSFRS Turning in Bed, Walking and Climbing Stairs")
  db[,dv5:=db$q_51_1.10 + db$q_51_1.11 + db$q_51_1.12]
  setattr(db$dv5, 'label', "ALSFRS Dyspnea, Orthopnea and Respiratory insufficiency")
  db[,dv6:=db$q_51_1.rate]

  db[,iv7:=(as.numeric(db$q_11 - db$q_5)/365.25)]
  setattr(db$iv7, 'label', "Age at first symptoms")

  src_vars=c('q_73', 'q_74', 'q_75', 'q_76', 'q_77', 'q_78', 'q_79', 'q_80', 'q_81', 'q_82', 'q_83', 'q_84')
  dst_vars=paste0('iv', seq(8, 7 + length(src_vars)))

  for(i in seq_along(src_vars)) {
    src_var=src_vars[[i]]
    dst_var=dst_vars[[i]]
    db[,(dst_var):=db[[src_var]]]
  }

  flags = dplyr::coalesce(db$q_85a_1!='No', FALSE) |
    dplyr::coalesce(db$q_85a_2!='No', FALSE) |
    dplyr::coalesce(db$q_85a_3!='No', FALSE)

  db[,iv20:=factor(flags, labels=c("No cancer", "Any cancer"), levels=c(FALSE,TRUE))]
  setattr(db$iv20, 'label', "Cancer")

  db[,iv21:=db$q_138]
  db[,iv22:=db$q_139]
  db[,iv23:=factor(db$q_139=='Yes' | db$q_138 == 'Yes', labels=c('No', 'Yes'), levels=c(FALSE, TRUE))]
  setattr(db$iv23, 'label', 'Emotional stress')


  db[,iv24:=db$q_20]
  db[,iv25:=db$q_21]
  db[,iv26:=db$q_22]
  db[,iv27:=db$q_23]
  db[,iv28:=db$q_37]
  db[,iv29:=db$q_38.1a]
  db[,iv30:=db$q_38.1b]
  db[,iv31:=factor(db$q_38.1a=='Yes' | db$q_38.1b == 'Yes', labels=c('No', 'Yes'), levels=c(FALSE, TRUE))]
  setattr(db$iv31, 'label', 'Resting tremor')
  db[,iv32:=db$q_38.2a]
  db[,iv33:=db$q_38.2b]
  db[,iv34:=factor(db$q_38.2a=='Yes' | db$q_38.2b == 'Yes', labels=c('No', 'Yes'), levels=c(FALSE, TRUE))]
  setattr(db$iv34, 'label', 'Bradykinesia')
  db[,iv35:=db$q_39]
  db[,iv36:=db$q_40]
  db[,iv37:=db$q_41]
  db[,iv38:=db$q_42]
  db[,iv39:=db$q_43]
  db[,iv40:=db$q_44]

  db[,iv41:=db$q_89a]
  db[,iv42:=db$q_90a]
  db[,iv43:=db$q_91a]
  db[,iv44:=db$q_92a]
  db[,iv45:=db$q_93a]
  db[,iv46:=db$q_94a]
  db[,iv47:=db$q_95a]
  db[,iv48:=db$q_96a]

  db[,iv49:=db$q_100a_1]
  setattr(db$iv49, 'f.o.b', 1)
  db[,iv50:=db$q_101a_1]
  setattr(db$iv50, 'f.o.b', 1)
  db[,iv51:=db$q_102a_1]
  setattr(db$iv51, 'f.o.b', 1)
  db[,iv52:=db$q_103a_1]
  setattr(db$iv52, 'f.o.b', 1)
  db[,iv53:=db$q_104a_1]
  setattr(db$iv53, 'f.o.b', 1)

  db[,iv54:=factor(db$q_115=='Yes' |
                     db$q_116 == 'Yes' |
                     db$q_117a == 'Yes' |
                     db$q_117b == 'Yes' |
                     db$q_117c == 'Yes' |
                     db$q_118a == 'Yes' |
                     db$q_118b == 'Yes' |
                     db$q_118c == 'Yes' , labels=c('No', 'Yes'), levels=c(FALSE, TRUE))]
  setattr(db$iv54, 'label', 'Trauma')


  db[,iv55:=factor(db$q_119a == 'Yes' |
                     db$q_119b == 'Yes' |
                     db$q_119c == 'Yes' |
                     db$q_120a == 'Yes' |
                     db$q_120b == 'Yes' |
                     db$q_120c == 'Yes' |
                     db$q_121a == 'Yes' |
                     db$q_121b == 'Yes' |
                     db$q_122a == 'Yes' |
                     db$q_122b == 'Yes' |
                     db$q_123 == 'Yes' |
                     db$q_124 == 'Yes' |
                     db$q_125 == 'Yes' |
                     db$q_126 == 'Yes' |
                     db$q_127 == 'Yes' |
                     db$q_128 == 'Yes' |
                     db$q_129 == 'Yes' |
                     db$q_130 == 'Yes', labels=c('No', 'Yes'), levels=c(FALSE, TRUE))]

  setattr(db$iv55, 'label', 'Surgery')




  # db[,dv1:=factor( (db$depvar1=='Bulbar'), labels=c('Not bulbar', 'Bulbar'), levels=c(FALSE, TRUE))]
  # setattr(db$dv1, 'label', 'Is fenotype bulbar')
  # db[,dv2:=factor( (db$depvar1=='Upper limb'), labels=c('Not upper limb', 'Upper limb'), levels=c(FALSE, TRUE))]
  # setattr(db$dv2, 'label', 'Is fenotype upper limb')
  # db[,dv3:=factor( (db$depvar1=='Lower limb'), labels=c('Not lower limb', 'Lower limb'), levels=c(FALSE, TRUE))]
  # setattr(db$dv3, 'label', 'Is fenotype lower limb')




  args<-list(
    zn=names(db)[stringr::str_detect(names(db), '^iv.*')],
    zz=names(db)[stringr::str_detect(names(db), '^dv.*')])
  return(list(db=as_tibble(db)[,c(args$zz, args$zn)], iv_names=args$zn, dv_names=args$zz))
}


regression_recode_ordinal<-function(varnames, depvar, mydt) {
  cat('Binning the variables...\n')
  woedt <- mydt %>% select_(.dots = varnames)
  mybins <- binnr::bin(data = woedt,
                       y = depvar,
                       mono = 1, max.bin=20)
  woes<-mybins$predict(type = 'woe')

  ans<-as.data.table(woes)
  danesurowe::copy_dt_attributes(mydt, ans)
  return(list(df=ans, woe_info=mybins$variables))
}

padded_count<-function(n, ntotal) {
  ndigits<-floor(log10(ntotal))+1
  return(stringr::str_pad(n, ndigits, pad='0'))
}

regression_recode_factor_1<-function(varname, dt, keep_nominal=character(0)) {

  # tmp_df<-data.frame(x=factor(c(1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,3,NA,NA,NA,NA),
  #                             levels=c(1,2,3), labels=c('A', 'B', 'C')))
  #as.integer(dt[[varname]])
  tmp_fac<-addNA(as.integer(dt[[varname]]), ifany=TRUE)
  if(length(unique(tmp_fac))>1) {
    tmp_df<-setNames(tibble::tibble(tmp_fac), varname)
    tmp_df<-predict(caret::dummyVars(~ ., tmp_df, fullRank=FALSE, levelsOnly = FALSE, drop2nd=TRUE), tmp_df)
    to_remove<-caret::nearZeroVar(as.data.table(tmp_df), saveMetrics = FALSE)
    if(varname %in% keep_nominal) {
      #    browser()
      for(to_rm in to_remove) {
        tmp_fac[as.integer(tmp_fac)==to_rm]<-NA
      }
      class(tmp_fac)<-class(dt[[varname]])
      tmp_fac<-addNA(tmp_fac, ifany = TRUE)
      if(any(is.na(levels(tmp_fac)))) {
        setattr(tmp_fac, 'levels', c(attr(dt[[varname]], 'levels'), as.character(NA)))
      } else {
        setattr(tmp_fac, 'levels', attr(dt[[varname]], 'levels'))
      }
      setNames(data.table(x=tmp_fac), varname)
      #We still need to replace small groups with NA.
    } else {
      to_keep<-setdiff(seq(ncol(tmp_df)), to_remove)
      data.table::as.data.table(tmp_df[,to_keep, drop=FALSE])
    }
  } else {
    data.table()
  }
}

# regression_recode_factor_1<-function(varname, dt, count_threshold=1) {
#   var<-dt[[varname]]
#   mylabel<-danesurowe::GetVarLabel(dt, varname)
#   levels<-as_tibble(table(var)) %>% filter(n>count_threshold)
#
#   if(nrow(levels)>2) {
#     base_level<-arrange(levels, -n)$var[[1]]
#     levels<-levels[2:nrow(levels),]
#     newname<-paste0(varname,'_', padded_count(1, nrow(levels)))
#     newvar<-as.integer(var==levels$var[[1]])
#     newvar[is.na(newvar)]<-0
#     ans<-data.table(var=newvar)
#     data.table::setnames(ans, 'var',  newname)
#     setattr(ans[[newname]], 'label', paste0(mylabel,' ',levels$var[[1]] ))
#     if(nrow(levels)>1) {
#       for(i in seq(2, nrow(levels))) {
#         newname<-paste0(varname,'_', padded_count(i, nrow(levels)))
#         newvar<-as.integer(var==levels$var[[i]])
#         newvar[is.na(newvar)]<-0
#         ans[,(newname):=newvar]
#         setattr(ans[[newname]], 'label', paste0(mylabel,' ',levels$var[[i]] ))
#       }
#     }
#   } else if (nrow(levels)==2) {
#     ans<-data.table(.tmp=as.integer(var))
#     data.table::setnames(ans, '.tmp', varname)
#   } else {
#     ans<-data.table(.tmp=rep(NA, nrow(dt)))
#   }
#   if(sum(is.na(var))>count_threshold) {
#     newname<-paste0(varname,'_NA')
#     ans[,(newname):=as.integer(is.na(var))]
#     setattr(ans[[newname]], 'label', paste0('Braki danych w ',mylabel))
#   }
#   if('.tmp' %in% names(ans)) {
#     ans[,.tmp:=NULL]
#   }
#   return(ans)
# }

regression_recode_binary_1<-function(varname, dt, count_threshold=1) {
  # cat(varname)
  # if(varname=='iv53') {
  #   browser()
  # }
  var<-as.numeric(dt[[varname]])
  levels<-unique(var)
  levels<-levels[!is.na(levels)]
  mylabel<-danesurowe::GetVarLabel(dt, varname)
  if(length(levels)>1) {
    newvar<-as.integer(var==levels[[2]])
    ans<-data.table(var=newvar)
    danesurowe::copy_var_attributes(dt[[varname]], var_dest_name = 'var', dt_dest = ans)
    data.table::setnames(ans, 'var',  varname)
  } else {
    ans<-data.table(.tmp=rep(NA, nrow(dt)))
  }

  if(sum(is.na(var))>count_threshold) {
    newname<-paste0(varname,'_NA')
    ans[,(newname):=as.integer(is.na(var))]
    setattr(ans[[newname]], 'label', paste0('Braki danych w ',mylabel))
  }
  if('.tmp' %in% names(ans)) {
    ans[,.tmp:=NULL]
  }

  return(ans)
}

prepare_variables_hybrid<-function(dt, iv_names, dv_name, keep_nominal=character(0)) {
  dv<-dt[[dv_name]]
  nas<-is.na(dv)
  mydt<-as.data.table(dt[!nas,iv_names])
  dv<-dv[!nas]

  fobs<-purrr::map_dbl(iv_names, ~danesurowe::GetFOB(mydt[[.]], flag_recalculate_uniques = TRUE))

  ordinals<-iv_names[fobs==2]
  ans<-list()
  if(length(ordinals)>0) {
    bins<-regression_recode_ordinal(ordinals, depvar, mydt)
    ans<-c(ans, bins$df)
  }
  nominals<-iv_names[fobs==1]
  if(length(nominals)>0) {
    a<-flatten(purrr::map(nominals, regression_recode_factor_1, dt=mydt, keep_nominal=keep_nominal))
    if ('V1' %in% names(a)) {
      browser()
    }
    ans<-c(ans,a)
  }
  binaries<-iv_names[fobs==3]
  if(length(binaries)>0) {
    a<-flatten(purrr::map(binaries, regression_recode_factor_1, dt=mydt, keep_nominal=keep_nominal))
#    a<-flatten(purrr::map(binaries, regression_recode_binary_1, dt=mydt))
    ans<-c(ans,a)
  }
  numerics<-iv_names[fobs==0]
  if(length(numerics)>0) {
    a<-mydt %>% select_(.dots=numerics) %>% as_tibble() %>% mutate_all(funs(as.numeric))
    #    a<-flatten(purrr::map(numerics, ~scale(mydt[[.]])))
    ans<-c(ans,a)
    missing_pattern<-mice::md.pattern(a)
    ans2<-as.data.table(ans)

    imputes<-NULL
    if(nrow(missing_pattern)>2) {
      if(nrow(mydt)>40000) {
        m<-missing_pattern[nrow(missing_pattern),1:ncol(missing_pattern)-1]
        vars_miss<-names(m[m>0])
        ms<-list()
        for(colname in vars_miss) {
          var<-ans2[[colname]]
          m<-mean(var, na.rm = TRUE)
          imputes<-c(imputes, list(sumNA=sum(is.na(var))))
          var[is.na(var)]<-m
          ms<-c(ms, list(m))
          ans2[,(colname):=as.numeric(var)]
        }
        imputes<-data.frame(varname=vars_miss, missing_count=unlist(imputes), mean=unlist(ms))

      } else {
        library(mice)
        imp<-mice::mice(ans2)
        ans2<-complete(imp)
      }
    }
  } else{
    ans2<-mydt
  }

  mydt<-as.data.table(ans2)
  return(list(db=mydt, dv=dv))
}

make_ads<-function(dt, iv_names, dv_name, keep_nominal=character(0)) {
  rmvar<-caret::nearZeroVar(dt)
  groupvar<-dt$iv56
  nearZeroVars<-colnames(dt)[ seq_along(colnames(dt)) %in% rmvar]
  vars_to_keep<-colnames(dt)[ ! seq_along(colnames(dt)) %in% rmvar]
  ads<-dt[c(vars_to_keep, dv_name)]
  iv_names=intersect(vars_to_keep, iv_names)

  ans<-prepare_variables_hybrid(ads, iv_names, dv_name, keep_nominal=keep_nominal)
  depvar<-ans$dv
  nas<-is.na(depvar)
  # if(length(keep)>0) {
  #   browser()
  #   to_keep<-setdiff(colnames(ans$db), keep)
  #   ans$db<-cbind(ans$db[to_keep], setNames(data.frame(x=dt[keep]), keep))
  # }
  cbind(data.frame(dv=depvar[!nas]), ans$db)
}
