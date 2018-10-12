#Gathers variable importances across all models
gather_variable_importances<-function(models, adf, df) {
  var_imp<-function(model, adf) {
    cat(paste0(model$method, '\n'))
    explainer<-DALEX::explain(model, data=adf %>% select(-dv), y=adf$dv)
    dfexp<-dplyr::select(DALEX::variable_importance(explainer, loss_function = DALEX::loss_root_mean_square, type = "difference"),-label)
    names(dfexp)<-c('variable', model$method)
    dfexp
  }

  var_imps<-foreach(model = models) %dopar%
    var_imp(model=model,  adf=adf)
  names(var_imps) <- names(models)

  comb_imps<-reduce(var_imps, function(imps1, imps2) {
    dplyr::full_join(imps1, imps2, by = 'variable')
  })


  for(i in seq(2, ncol(comb_imps))) {
    lab_idx<-which(df$model == colnames(comb_imps)[[i]])
    lab<-df$name[[lab_idx]]
    data.table::setattr(comb_imps[[i]], 'label', lab)
  }
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
