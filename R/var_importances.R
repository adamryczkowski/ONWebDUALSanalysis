#Gathers variable importances across all models
gather_variable_importances<-function(models, adf, df) {
  var_imp<-function(model, adf) {
    dfexp=tryCatch({
#      cat(paste0(model$method, '\n'))
      explainer<-DALEX::explain(model, data=adf %>% select(-dv), y=adf$dv)
      dfexp<-dplyr::select(DALEX::variable_importance(explainer, loss_function = DALEX::loss_root_mean_square, type = "ratio"),-label)
      names(dfexp)<-c('variable', model$method)
      return(dfexp)},
      error=function(e) {
        return(paste0('error ', e$message, ' in model ', model$method))
      }
    )
    dfexp
  }


  var_imps<-foreach(model = models) %dopar%
    var_imp(model=model,  adf=adf)

  names(var_imps) <- names(models)

  var_imps<-var_imps[!purrr::map_lgl(var_imps, is_character)]

  comb_imps<-reduce(var_imps, function(imps1, imps2) {
    dplyr::full_join(imps1, imps2, by = 'variable')
  })

  for(i in seq(2, ncol(comb_imps))) {
    lab_idx<-which(df$model == colnames(comb_imps)[[i]])
    lab<-df$name[[lab_idx]]
    data.table::setattr(comb_imps[[i]], 'label', lab)
  }

  da<-data.matrix(comb_imps[seq(2, ncol(comb_imps))])
  zero_models<-c(1, 1+which(abs(plyr::aaply(da, 2, function(x) sum(x^2) )-0)>1.0E-10))
  comb_imps<-comb_imps[zero_models]

  comb_imps$variable<-as.character(comb_imps$variable)

  comb_imps<-dplyr::right_join(tibble(variable=colnames(adf), varlabel=Hmisc::label(adf) ), comb_imps, by='variable')

  return(comb_imps)
}

mds_from_varimps<-function(comb_imps, df, flag_variables=TRUE) {
  tmpdf1<-comb_imps %>% filter(variable != '_full_model_' & variable != '_baseline_')
  tmpdf<-data.matrix(tmpdf1 %>% dplyr::select(-variable))
  colnames(tmpdf)<-colnames(comb_imps %>% dplyr::select(-variable))
  rownames(tmpdf)<-tmpdf1$variable

  if(!flag_variables) {
    tmpdf<-t(tmpdf)
  }
  comb_dists<-as.matrix(dist(tmpdf))
  comb_weights<-plyr::aaply(tmpdf, 1, mean)
  comb_weights_sd<-plyr::aaply(tmpdf, 1, sd)
  comb_weights_cnt<-plyr::aaply(tmpdf, 1, function(x) sum(x != 0))
  min_comb_weights_cnt<-min(comb_weights_cnt)
  comb_weights_cnt<-comb_weights_cnt - min(comb_weights_cnt)


  library(ggrepel)

  mds_comb<-cbind(as.data.frame(cmdscale(comb_dists)),
                  weight=comb_weights, weight_sd=comb_weights_sd, weight_cnt=comb_weights_cnt,
                  label=Hmisc::label(comb_imps)[2:length(comb_imps)] ) %>% arrange(-comb_weights_cnt)
  plot1 <- ggplot(mds_comb, aes(V1, V2, label=label, size=weight_cnt)) +
    geom_point(colour="blue", alpha=0.2) +
    scale_size('weight_cnt')+
    geom_text_repel(colour="black", size=2.5,
                    hjust = "center", vjust = "bottom", nudge_x = 0, nudge_y = 0.025) +
    labs(x="", y="") + theme_bw()

  plot1
}

gather_detailed_variable_importances<-function(model, adf, df, how_many_resamples=500) {
  var_imp<-function(model, adf, nr) {
    explainer<-DALEX::explain(model, data=adf %>% select(-dv), y=adf$dv)
    dfexp<-dplyr::select(DALEX::variable_importance(explainer, loss_function = DALEX::loss_root_mean_square, type = "raw"),-label)
    #    dfexp<-dplyr::select(DALEX::variable_importance(explainer, loss_function = DALEX::loss_root_mean_square, type = "difference"),-label)
    names(dfexp)<-c('variable', paste0('var_', nr))
    dfexp
  }

  var_imps<-foreach(nr = seq(how_many_resamples)) %dopar%
    var_imp(model=model,  adf=adf, nr=nr)

  comb_imps<-reduce(var_imps, function(imps1, imps2) {
    dplyr::full_join(imps1, imps2, by = 'variable')
  })

  a<-data.matrix(comb_imps[,seq(2, ncol(comb_imps))])
  fullmodels<-a[1,]
  baselines<-a[nrow(a),]

  a2<-((a - fullmodels)/(baselines - fullmodels))
#  a2<-a
  ans<-tibble(variable=as.character(comb_imps[[1]]),
              # raw_importances=plyr::aaply(a2,1,mean),
              # sd_raw_importances=plyr::aaply(a2,1,sd),
              # md_raw_importances=plyr::aaply(a2,1,median),
              # low_raw_importances=plyr::aaply(a2,1,function(x) quantile(x, probs = 0.25)),
              # high_raw_importances=plyr::aaply(a2,1,function(x) quantile(x, probs = 0.75))
              importances=plyr::aaply(a2,1,function(x) mean(x, na.rm=TRUE)),
              sd_importances=plyr::aaply(a2,1,function(x) sd(x, na.rm=TRUE)),
              md_importances=plyr::aaply(a2,1,function(x) median(x, na.rm=TRUE)),
              low_importances=plyr::aaply(a2,1,function(x) quantile(x, probs = 0.25, na.rm=TRUE)),
              high_importances=plyr::aaply(a2,1,function(x) quantile(x, probs = 0.75, na.rm=TRUE))
  )

  ans2<-dplyr::arrange(dplyr::right_join(tibble(variable=colnames(adf), varlabel=Hmisc::label(adf)), ans, by='variable'), -importances)
  return(ans2)
}
