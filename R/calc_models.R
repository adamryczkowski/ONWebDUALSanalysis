calc_models<-function(model_names, dv_nr, path_prefix='models/', adaptive=TRUE) {
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
      readRDS(plik_name)
    } else {
      cat(paste0("Calculating model ", model_name, "...\n"))
      tryCatch(
        {
          if(adaptive) {
            model<-caret::train(dv ~ ., data = ads, method = model_name, trControl = tc, tuneLength=15)
          } else {
            model<-caret::train(dv ~ ., data = ads, method = model_name, trControl = tc)
          }
          saveRDS(model, plik_name)
          model
        },
        error=function(e) {
          cat(paste0("Model ", model_name, " returned error: ", e$message, '\n'))
        }
      )
    }
  }

  ads<-make_ads(dt, iv_names, dv_name, keep_nominal ='iv56')

  #  which(is.na(data.matrix(ads)),arr.ind = TRUE)

  groupvar<-ads$iv56
  cvIndex<-caret::createFolds(groupvar, 10, returnTrain = T)

  if(adaptive) {
    tc <- caret::trainControl(index = cvIndex,
                              method = 'adaptive_cv',
                              number = 10, repeats = 10,
                              adaptive = list(min = 5, alpha = 0.05,
                                              method = "gls", complete = TRUE),
                              search = "random"
    )
  } else {
    tc <- caret::trainControl(index = cvIndex,
                              method = 'cv',
                              number = 10, repeats = 10)
  }


  list(ads=ads, models=purrr::map(model_names, do_model_inner, ads=ads, tc=tc, dv_name=dv_name))
}
