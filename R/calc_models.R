calc_models<-function(model_names, dv_nr, path_prefix='models/', adaptive=NA, assume_calculated=FALSE) {
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
      cat(paste0("Reading in already calculated model ", model_name, "...\n"))
      return(readRDS(plik_name))
    } else {
      if(assume_calculated) {
        msg<-paste0("Model ", model_name, " is not computed, skipping it.")
        cat(paste0(msg,'\n'))
        saveRDS(msg, plik_name)
        return(msg)
      } else {
        return(tryCatch(
          {
            if(is.na(adaptive)) {
              cat(paste0("Trying adaptive train of model ", model_name, "...\n"))
              model<-caret::train(dv ~ ., data = ads, method = model_name,
                                  trControl = tc_adaptive, tuneLength=15)
            } else if(adaptive) {
              cat(paste0("Calculating adaptive train of model ", model_name, "...\n"))
              model<-caret::train(dv ~ ., data = ads, method = model_name,
                                  trControl = tc_adaptive, tuneLength=15)
            } else {
              cat(paste0("Calculating non-adaptive train of model ", model_name, "...\n"))
              model<-caret::train(dv ~ ., data = ads, method = model_name, trControl = tc)
            }
            saveRDS(model, plik_name)
            return(model)
          },
          error=function(e) {
            if(stringr::str_detect(e$message, stringr::fixed('For adaptive resampling, there needs to be more than one tuning parameter'))) {
              cat(paste0("Adaptive train failed. Calculating non-adaptive train of model ", model_name, "...\n"))
              return(tryCatch(
                {
                  model<-caret::train(dv ~ ., data = ads, method = model_name, trControl = tc, tuneLength=15)
                  saveRDS(model, plik_name)
                  return(model)
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

  ads<-make_ads(dt, iv_names, dv_name, keep_nominal ='iv56')

  #  which(is.na(data.matrix(ads)),arr.ind = TRUE)

  groupvar<-ads$iv56
  cvIndex<-caret::createFolds(groupvar, 10, returnTrain = T)
  selFun<-function(x, metric,  maximize) caret::oneSE(x=x, metric = metric, num=10, maximize = maximize)
  tc_adaptive <- caret::trainControl(index = cvIndex,
                                     method = 'adaptive_cv',
                                     number = 10, repeats = 10,
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

  df<-dplyr::arrange(tibble(model=model_names, name=n1,  rmse=a1, rsq=a2, mae=a3, elapsed_time=b1, user_time=b2+b3,
                            is_nn=b4_1, is_bagging=b4_2, is_boost=b4_7,
                            is_rf=b4_3, is_lm=b4_4, is_bayes=b4_5, is_feature_sel=b4_6, modelType = b4_8), rmse)

  res<-summary(resamples(models[df$model]), metric='RMSE', decreasing=TRUE)
  rmse_3rd<-res$statistics$RMSE[1,5] #3rd quantile of the best model's RMSE
  idx_ok<-which(res$statistics$RMSE[,2]<rmse_3rd) #Which models are not statistically worse then the best
  cat(paste0("Discarded ", nrow(df)-length(idx_ok), " models that are not as good as the best model\n"))
  return(dplyr::arrange(df[idx_ok,], rmse))
}
