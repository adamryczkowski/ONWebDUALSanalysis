do_calc<-function(dv_nr, cache_path='models') {
  registerDoMC(12)
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
  ans<-calc_models(model_names, dv_nr=dv_nr, adaptive = NA, assume_calculated = FALSE, path_prefix=paste0(cache_path, '/'))

}

calc_models<-function(model_names, dv_nr, path_prefix='models/', adaptive=NA, assume_calculated=FALSE, keep_nominal=character(0)) {
  library(caret)
  joined_df<-readRDS(system.file('db_no_duplicates.rds', package='ONWebDUALSanalysis'))
  ans<-select_variables_sep2018(joined_df)
  dt<-ans$db
  iv_names<-ans$iv_names
  dv_names<-ans$dv_names
  dv_name<-dv_names[[dv_nr]]

  do_model_inner<-function(dv_name, model_name, ads, tc) {
    plik_name<-paste0(path_prefix, 'model_', dv_name, '_', model_name, '.rds')
    if(file.exists(plik_name)) {
      if (assume_calculated) {
        if(file.size(plik_name)>0) {
          cat(paste0("Reading in already calculated model ", model_name, "...\n"))
          return(readRDS(plik_name))
        } else {
          msg<-paste0("File ", plik_name, " is still being calculated")
          cat(paste0(msg, "\n"))
          return(msg)
        }
      } else {
        msg<-paste0("Skipping already calculated model ", model_name)
        cat(paste0(msg, "\n"))
        return(msg)
      }
    } else {
      if(assume_calculated) {
        msg<-paste0("Model ", model_name, " is not computed, skipping it.")
        cat(paste0(msg,'\n'))
        return(msg)
      } else {
        write.table(data.frame(), file=plik_name, col.names=FALSE)
        return(tryCatch(
          {
            if(is.na(adaptive)) {
              msg<-paste0("Trying adaptive train of model ", model_name, "...")
              cat(paste0(msg, "\n"))
              model<-caret::train(dv ~ ., data = ads, method = model_name,
                                  trControl = tc_adaptive, tuneLength=15)
            } else if(adaptive) {
              msg<-paste0("Calculating adaptive train of model ", model_name, "..")
              cat(paste0(msg, "\n"))
              model<-caret::train(dv ~ ., data = ads, method = model_name,
                                  trControl = tc_adaptive, tuneLength=15)
            } else {
              msg<-paste0("Calculating non-adaptive train of model ", model_name, "...")
              cat(paste0(msg, "\n"))
              model<-caret::train(dv ~ ., data = ads, method = model_name, trControl = tc)
            }
            saveRDS(model, plik_name)
            return(msg)
          },
          error=function(e) {
            if(stringr::str_detect(e$message, stringr::fixed('For adaptive resampling, there needs to be more than one tuning parameter'))) {
              msg<-paste0("Adaptive train failed. Calculating non-adaptive train of model ", model_name, "...")
              cat(paste0(msg, "\n"))
              return(tryCatch(
                {
                  model<-caret::train(dv ~ ., data = ads, method = model_name, trControl = tc, tuneLength=15)
                  saveRDS(model, plik_name)
                  return(msg)
                },
                error=function(e) {
                  msg=paste0("Non-adaptive run of model ", model_name, " returned error: ", e$message)
                  cat(paste0(msg, '\n'))
                  saveRDS(msg, plik_name)
                  return(msg)
                }
              ))
            }
            msg=paste0("Adaptive run of model ", model_name, " returned error: ", e$message)
            cat(paste0(msg, '\n'))
            saveRDS(msg, plik_name)
            return(msg)
          }
        ))
      }
    }
  }
  #ads<-make_ads(dt, iv_names, dv_name, keep_nominal ='iv56')
  ads<-make_ads(dt, c(iv_names, 'dv3'), dv_name, keep_nominal=keep_nominal)

  #  which(is.na(data.matrix(ads)),arr.ind = TRUE)

  if(length(keep_nominal)>0) {
    groupvar<-ads[[keep_nominal]]
    cvIndex<-caret::createMultiFolds(groupvar, k = 10, times=10)
  } else {
    cvIndex<-caret::createMultiFolds(ads$dv, times=10,  k = 10)
  }
  selFun<-function(x, metric,  maximize) caret::oneSE(x=x, metric = metric, num=10, maximize = maximize)
  tc_adaptive <- caret::trainControl(index = cvIndex,
                                     method = 'adaptive_cv',
                                     adaptive = list(min = 5, alpha = 0.05,
                                                     method = "gls", complete = TRUE),
                                     search = "random",
                                     selectionFunction = selFun
  )
  tc <- caret::trainControl(index = cvIndex,
                            method = 'cv',
                            number = 10, repeats = 10)
  models=setNames(purrr::map(model_names, do_model_inner, ads=ads, tc=tc, dv_name=dv_name), model_names)

  list(ads=ads, models=models)
}

model_perfs<- function(ans) {
  valid_models<-which(!purrr::map_lgl(ans$models, is.character))
  models<-ans$models[valid_models]
  model_names<-names(models)
  ads<-ans$ads

  a1<-purrr::map_dbl(models, function(x) {caret::getTrainPerf(x)$TrainRMSE} )
  a3<-purrr::map_dbl(models, function(x) {caret::getTrainPerf(x)$TrainRsquared} )
  a2<-purrr::map_dbl(models, function(x) {caret::getTrainPerf(x)$TrainMAE} )
  n1<-purrr::map_chr(models, function(x) {x$modelInfo$label} )
  b1<-purrr::map_dbl(models, function(x) {as.numeric(x$times$everything['elapsed'])} )
  b2<-purrr::map_dbl(models, function(x) {as.numeric(x$times$everything['user.self'])} )
  b3<-purrr::map_dbl(models, function(x) {as.numeric(x$times$everything['user.child'])} )
  b4<-purrr::map_chr(models, function(x) {paste(x$modelInfo$tags, collapse = ', ')} )
  b4_1<-purrr::map_lgl(models, function(x) {'Neural Network' %in% x$modelInfo$tags} )
  b4_2<-purrr::map_lgl(models, function(x) {'Bagging' %in% x$modelInfo$tags} )
  b4_3<-purrr::map_lgl(models, function(x) {'Random Forest' %in% x$modelInfo$tags} )
  b4_4<-purrr::map_lgl(models, function(x) {'Linear Regression' %in% x$modelInfo$tags} )
  b4_5<-purrr::map_lgl(models, function(x) {'Bayesian Model' %in% x$modelInfo$tags} )
  b4_6<-purrr::map_lgl(models, function(x) {'Implicit Feature Selection' %in% x$modelInfo$tags} )
  b4_7<-purrr::map_lgl(models, function(x) {'Boosting' %in% x$modelInfo$tags} )
  b4_8<-purrr::map_chr(models, function(x) {x$modelType} )
  b5<-purrr::map(models, function(x) {list(setdiff(x$modelInfo$tags,
                                                       c('Bagging', 'Implicit Feature Selection', 'Boosting')))} )

  b5_grp<-rep(NA, length(b5))

  fn<-function(b5_grp, tag, tag_level){
    lgls<-which(purrr::map_lgl(b5, function(x) tag %in% x[[1]]))

    tmp<-lgls[which(!is.na(b5_grp[lgls]))]
    if(length(tmp)>0) {
      for(i in tmp) {
        cat(paste0("Model ", n1[[i]], "(",model_names[[i]],") is also of type ", b5_grp[[i]], ".\n"))
      }
    }
    cat(paste0(length(lgls), " ", tag, "s\n"))
    b5_grp[lgls]<-tag_level
    return(b5_grp)
  }

  b5_grp<-fn(b5_grp, 'Linear Regression', 1)
  b5_grp<-fn(b5_grp, 'Relevance Vector Machines', 2)
  b5_grp<-fn(b5_grp, 'Support Vector Machines', 2)
  b5_grp<-fn(b5_grp, 'Gaussian Process',3)
  b5_grp<-fn(b5_grp, 'Random Forest', 4)
  b5_grp<-fn(b5_grp, 'Multivariate Adaptive Regression Splines', 5)
  b5_grp<-fn(b5_grp, 'Tree-Based Model', 4)
  b5_grp<-fn(b5_grp, 'Neural Network', 6)
  b5_grp[is.na(b5_grp)]<-7



  df<-dplyr::arrange(tibble(model=model_names, name=n1,  rmse=a1, rsq=a2, mae=a3, elapsed_time=b1, user_time=b2+b3,
                            is_nn=b4_1, is_bagging=b4_2, is_boost=b4_7,
                            is_rf=b4_3, is_lm=b4_4, is_bayes=b4_5, is_feature_sel=b4_6, modelType = b4_8, tags = b5,
                            model_family=factor(b5_grp, levels=sort(unique(b5_grp)), labels=c('Linear Regression', 'Support Vector Machines',
                                                                             'Gaussian Process', 'Random Forest and trees',
                                                                             'Multivariate Adaptive Regression Splines',
                                                                             'Neural Network', 'Other'))), rmse)

  sub_df <- df %>% filter(b5_grp != "Linear Regression")
  models_to_get_metric<-c('glmnet', sub_df$model)
  models_hard_to_compare<-c('lasso', 'pls', 'rpart', 'rpart2', 'simpls')
  models_hard_to_compare<-c('bagEarth', 'blasso', 'earth', 'gcvEarth', 'kernelpls', 'lars2', 'lasso', 'pls', 'pcr', 'rpart', 'rpart2', 'simpls', 'spikeslab')
  #models_hard_to_compare<-c('bagEarth', 'blasso', 'earth', 'gcvEarth', 'kernelpls', 'lars2', 'lasso', 'pls', 'rpart', 'rpart2', 'simpls', 'spikeslab')
  models_to_get_metric<-setdiff(sample(df$model), models_hard_to_compare)

  df_to_remove<- df %>% filter(model %in% models_hard_to_compare)

#  browser()
  # res_big<-summary(resamples(models[
  #   c('rpart', 'rpart2', 'blasso', 'kernelpls')]), metric='RMSE', decreasing=TRUE)
  # c('gcvEarth', 'earth', 'rpart', 'rpart2', 'simpls', 'pls', 'lasso', 'blasso', 'kernelpls')
  # res_big<-summary(resamples(models[models_hard_to_compare[1:4]]), metric='RMSE', decreasing=TRUE)


  # old_timing<-0
  # old_size<-0
  # for(i in seq(2, length(models_to_get_metric))) {
  #   subsecik<-models_to_get_metric[seq(1, i)]
  #   cat(paste0("Trying adding ", length(subsecik), "th model ", subsecik[[i]], " of size ", utils:::format.object_size(object.size(models[[i]]), "auto"), "\n"))
  #   ts<-system.time(res<-summary(resamples(models[subsecik]), metric='RMSE', decreasing=TRUE))
  #   new_timing<-ts[['elapsed']]
  #   new_size<-object.size(res)
  #   cat(paste0("Elapsed time: ", round(new_timing), " which is ", round(100*(new_timing-old_timing)/(length(subsecik)-1))/100, " sec per sample.\n"))
  #   cat(paste0("Extra object size time: ", gdata::humanReadable(new_size - old_size),
  #              " which is ", gdata::humanReadable(round((new_size-old_size)/(length(subsecik)-1))), " per sample.\n\n"))
  #   old_timing<-new_timing
  #   old_size<-new_size
  # }


  res<-summary(resamples(models[models_to_get_metric]), metric='RMSE', decreasing=TRUE)



#  res<-summary(resamples(models[df$model]), metric='RMSE', decreasing=TRUE)
  rmse_3rd<-min(res$statistics$RMSE[,5]) #3rd quantile of the best model's RMSE
  idx_ok<-rownames(res$statistics$RMSE)[(res$statistics$RMSE[,2]<rmse_3rd)] #Which models are not statistically worse then the best

  # best_model<-df$model[[1]]
  # mem_size<-0
  # dftmp<-NULL
  # for(i in seq(2, length(df$model))) {
  #   ms<-models[c(best_model, df$model[[i]])]
  #   res<-summary(resamples(ms), metric='RMSE', decreasing=TRUE)
  #   if(is.null(dftmp)) {
  #     dftmp<-res$statistics$RMSE
  #   } else {
  #     dftmp<-rbind(dftmp, res$statistics$RMSE[2,])
  #     rownames(dftmp)<-c(rownames(dftmp)[seq(1, nrow(dftmp)-1)], rownames(res$statistics$RMSE)[[2]])
  #   }
  # }
  # res<-dftmp
  # #res<-summary(resamples(models[df$model]), metric='RMSE', decreasing=TRUE)
  # rmse_3rd<-res[1,5] #3rd quantile of the best model's RMSE
  # idx_ok<-which(res[,2]<rmse_3rd) #Which models are not statistically worse then the best
  cat(paste0("Discarded ", nrow(df)-length(idx_ok), " models that are not as good as the best model\n"))
#  df<-dplyr::arrange(df[idx_ok,], rmse)
  df<-dplyr::inner_join(df, tibble::tibble(model=idx_ok), by='model')
  res<-caret::resamples(models[df$model])


  return(list(df=df, resamples=res))
}

