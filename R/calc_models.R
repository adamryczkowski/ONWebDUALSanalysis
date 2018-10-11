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
        return(msg)
      } else {
        return(tryCatch(
          {
            if(is.na(adaptive)) {
              cat(paste0("Trying adaptive train of model ", model_name, "...\n"))
              model<-caret::train(dv ~ ., data = ads, method = model_name, trControl = tc_adaptive, tuneLength=15)
            } else if(adaptive) {
              cat(paste0("Calculating adaptive train of model ", model_name, "...\n"))
              model<-caret::train(dv ~ ., data = ads, method = model_name, trControl = tc_adaptive, tuneLength=15)
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
                  return(msg)
                }
              ))
            }
            msg=paste0("Adaptive run of model ", model_name, " returned error: ", e$message)
            cat(paste0(msg, '\n'))
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

  tc_adaptive <- caret::trainControl(index = cvIndex,
                                     method = 'adaptive_cv',
                                     number = 10, repeats = 10,
                                     adaptive = list(min = 5, alpha = 0.05,
                                                     method = "gls", complete = TRUE),
                                     search = "random"
  )
  tc <- caret::trainControl(index = cvIndex,
                            method = 'cv',
                            number = 10, repeats = 10)
  models=setNames(purrr::map(model_names, do_model_inner, ads=ads, tc=tc, dv_name=dv_name), model_names)

  list(ads=ads, models=models)
}
