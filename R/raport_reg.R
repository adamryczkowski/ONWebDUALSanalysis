#Function that generates raport for one dv
make_rap<-function(dv_nr, rap_path) {
  time_consuming_models<-c('ANFIS', 'DENFIS', 'FIR.DM', 'FS.HGD', 'GFS.FR.MOGUL', 'GFS.LT.RS', 'HYFIS', 'Rborist', 'xgbDART', 'xgbLinear', 'xgbTree')
  rather_long<-c('brnn', 'nodeHarvest', 'qrnn', 'rfRules')
  #empty_models<-c('avNNet', 'ANFIS')
  really_long<-c('DENFIS', 'FIR.DM', 'FS.HGD', 'leapSeq')

  empty_models<-character(0)
  tensor_flow<-c('mlpKerasDecay', 'mlpKerasDropout')
  not_parallel<-c('M5', 'M5Rules')
  mem_insufficient<-c('randomGLM')
  models_with_broken_packages<-c('elm', 'mxnet', 'mxnetAdam', 'mlpSGD')
  model_blacklist<-c('bag', 'bam', 'bartMachine', 'blackboost','bstSm', 'bstTree', 'elm', 'extraTrees',
                     'gam', 'gamboost', 'gamLoess', 'gamSpline', 'gbm_h2o', 'GFS.THRIFT',
                     'glmboost', 'glmnet_h2o', 'HYFIS', 'krlsRadial',
                     'logicBag', 'logreg', 'mlpSGD', 'mxnet', 'mxnetAdam',
                     'neuralnet', 'rlm', 'svmBoundrangeString', 'svmExpoString', 'svmSpectrumString',
                     #   'xgbLinear',
                     time_consuming_models, empty_models, not_parallel, tensor_flow, mem_insufficient)

  all_models<-unique((caret::modelLookup() %>% filter(forReg==TRUE & ! model %in% model_blacklist))$model)
  models_that_hangs<-c('bartMachine', 'extraTrees')
  model_names<-setdiff(all_models, c(time_consuming_models, really_long, rather_long,  models_that_hangs,
                                        tensor_flow,not_parallel, mem_insufficient, models_with_broken_packages)  )
  #Reads all models for a given dv_nr
  ans<-calc_models(model_names, dv_nr=dv_nr, adaptive = NA, assume_calculated = TRUE)
  models<-ans$models
  adf<-ans$ads


  #Gathers all metadata for models
  df<-model_perfs(ans)

  #Filters models that are much worse
  models<-models[df$model]

  #Imports all the possibly missing packages
  packages<-unlist(purrr::map(models, function(m) m$modelInfo$library ))
  new.packages <- packages[!(packages %in% installed.packages()[,"Package"])]
  if(length(new.packages)) install.packages(new.packages)

  #Show all models' performance
  comp_models<-function(models, df) {
    r<-caret::resamples(models)
    o<-order(names(models))

    model_family_colors<-RColorBrewer::brewer.pal(name="Dark2", n = length(unique(df$model_family)))
    a<-ggplot(r,metric="RMSE",main="GBM vs xgboost")
    a$data<-cbind(a$data, model_family = df$model_family, model_family_colors = model_family_colors[as.integer(df$model_family)])
    a + aes(color=model_family) + scale_color_manual(values = model_family_colors)+ labs(color = "Model family") +
      ylab(paste0("RMSE of predicting ", Hmisc::label(adf$dv), " (less is better)")) +
      theme(axis.text.y = element_text(colour=model_family_colors[as.integer(df$model_family)][order(a$data$Estimate)]))

  }
  plot<-comp_models(models, df)
  ggsave(filename=paste0('all_models_perf_', dv_nr, '.png'), plot=plot, device='png', path=rap_path)
  ggsave(filename=paste0('all_models_perf_', dv_nr, '.svg'), plot=plot, device='svg', path=rap_path)

  #Gather variable importance statistics
  if(file.exists(paste0('tmp_comp_imps_', dv_nr, '.rds'))) {
    comb_imps<-readRDS(paste0('tmp_comp_imps_', dv_nr, '.rds'))
  } else {
    comb_imps<-gather_variable_importances(models=models, adf=adf, df=df)
    saveRDS(comb_imps, paste0('tmp_comp_imps_', dv_nr, '.rds'))
  }

  mds_from_varimps<-function(df_comb, df, adf, flag_variables=TRUE) {
    tmpdf1<-df_comb %>% filter(variable != '_full_model_' & variable != '_baseline_')
    tmpdf<-data.matrix(tmpdf1 %>% dplyr::select(-variable))
    colnames(tmpdf)<-colnames(df_comb %>% dplyr::select(-variable))
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

    if(flag_variables) {
      labels<-Hmisc::label(df_comb)[2:length(df_comb)]
    } else {
      pos<-which(purrr::map_int(as.character(tmpdf1$variable), function(x) which(colnames(adf)  %in% x )))
      labels<-Hmisc::label(adf)[pos]
    }

    mds_comb<-cbind(as.data.frame(cmdscale(comb_dists)),
                    weight=comb_weights, weight_sd=comb_weights_sd, weight_cnt=comb_weights_cnt,
                    label=labels ) %>% arrange(-comb_weights_cnt)
    plot1 <- ggplot(mds_comb, aes(V1, V2, label=label, size=weight_cnt)) +
      geom_point(colour="blue", alpha=0.2) +
      scale_size('weight_cnt')+
      ggrepel::geom_text_repel(colour="black", size=2.5,
                      hjust = "center", vjust = "bottom", nudge_x = 0, nudge_y = 0.025) +
      labs(x="", y="") + theme_bw()

    plot1
  }

  a<-which(colnames(comb_imps) %in% (df %>% dplyr::filter(model_family=='Linear Regression'))$model)
  df_comb<-comb_imps[,c(1,a)]
  mds_from_varimps(df_comb, df = df, flag_variables = TRUE)

#  dplyr::full_join(comb_imps  colnames(comb_imps)

  summary(res, metric='RMSE', decreasing = FALSE)
  bwplot(res, metric='RMSE', decreasing = TRUE)
  dotplot(res)

  best_model_name<-df$model[[1]]
  best_model<-models[[ best_model_name ]]
  glmnet_model<-models[['glmnet']]

  ads_pred<-cbind(ads, dv_predict=predict(best_model))

}