#devtools::install_github('adamryczkowski/ONWebDUALSanalysis')
library(ONWebDUALSanalysis)
library(doMC)
registerDoMC()


time_consuming_models<-c('ANFIS', 'DENFIS', 'FIR.DM', 'FS.HGD', 'GFS.FR.MOGUL', 'GFS.LT.RS', 'HYFIS', 'Rborist', 'xgbDART', 'xgbLinear', 'xgbTree')
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
dubious_models<-setdiff(all_models, c(time_consuming_models, really_long,  models_that_hangs,
                                      tensor_flow,not_parallel, mem_insufficient, models_with_broken_packages)  )
model_names=time_consuming_models
model_names=c(all_models, c('mlpKerasDropout','M5', 'M5Rules'))
model_names=all_models
model_names=not_parallel
model_names=c("nodeHarvest", "earth", "enet", "BstLm", "glmnet", "RRF", "ctree2", "ranger", "RRFglobal", "rf", "cforest",
              "evtree", "gbm", "spikeslab", "ctree", "lars", "lasso", "rpart2", "parRF",
              "cubist", "lars2", "rqnc", "rpart", 'krlsRadial',
              "penalized", "msaenet", "rfRules", "qrf", "relaxo")
ans<-calc_models(model_names, dv_nr=5, adaptive = NA)
ans<-calc_models('gbm_h2o', dv_nr=5, adaptive = NA)
ans<-calc_models(dubious_models, dv_nr=1, adaptive = NA)

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

  length(m$finalModel$nodes)

  df<-dplyr::arrange(tibble(model=model_names, name=n1,  rmse=a1, rsq=a2, mae=a3, elapsed_time=b1, user_time=b2+b3,
             is_nn=b4_1, is_bagging=b4_2, is_boost=b4_7,
             is_rf=b4_3, is_lm=b4_4, is_bayes=b4_5, is_feature_sel=b4_6, modelType = b4_8), rmse)

  res<-summary(resamples(models[df$model]), metric='RMSE', decreasing=TRUE)
  rmse_3rd<-res$statistics$RMSE[1,5] #3rd quantile of the best model's RMSE
  idx_ok<-which(res$statistics$RMSE[,2]<rmse_3rd) #Which models are not statistically worse then the best
  return(dplyr::arrange(df[idx_ok,], rmse))
}

df<-model_perfs(ans)
m<-models[[df$model[[3]]]]
descr_model<-function(m) {
  tuned_pars<-setNames(m$finalModel$tuneValue, as.character(m$modelInfo$parameters$label))
  return(tuned_pars)
}

descr_regr_model<-function(m) {
  if(m$modelType!='Regression') {
    return(NA)
  }
  coef(m$finalModel)
}

a<-unclass(m$finalModel)


summary(res, metric='RMSE', decreasing = FALSE)
bwplot(res, metric='RMSE', decreasing = TRUE)
dotplot(res)

best_model_name<-df$model[[1]]
best_model<-models[[ best_model_name ]]
glmnet_model<-models[['glmnet']]

ads_pred<-cbind(ads, dv_predict=predict(best_model))

#ggplot2::qplot(dv, dv_predict, data = ads)


min_point<-c(min(ads$dv), min(ads_pred$dv_predict))
max_point<-c(max(ads$dv), max(ads_pred$dv_predict))

min_dist<-((ads$dv - min_point[[1]])^2 + (ads_pred$dv_predict - min_point[[2]])^2)^0.5
max_dist<-((ads$dv - max_point[[1]])^2 + (ads_pred$dv_predict - max_point[[2]])^2)^0.5

ktory_quantil=0.05

dvquantiles<-c(quantile(min_dist, probs = ktory_quantil),quantile(max_dist, probs = ktory_quantil))

min_records<-which(min_dist<=dvquantiles[[1]])
max_records<-which(max_dist<=dvquantiles[[2]])


explainer<-DALEX::explain(best_model, data=ads, label=best_model_name)

library(doParallel)
registerDoParallel(cores=4)

make_dalex_single_pred_df<-function(ads, records) {
  rec_dalex<-foreach(rec_nr=records) %dopar%
    DALEX::single_prediction(explainer, observation = ads[rec_nr,])

  a<-unlist(purrr::map(rec_dalex, function(x) x$contribution))
  dim(a)<-c(nrow(rec_dalex[[1]]), length(rec_dalex))
  cum_min_dalex<-tibble(name=rec_dalex[[1]]$variable,
                        mean_contr = plyr::aaply(a, 1, mean),
                        sd_contr = plyr::aaply(a, 1, sd),
                        min_contr = plyr::aaply(a, 1, min),
                        max_contr = plyr::aaply(a, 1, max),
                        med_contr = plyr::aaply(a, 1, median),
                        q05_contr = plyr::aaply(a, 1, function(x) quantile(x, probs=0.05)),
                        q95_contr = plyr::aaply(a, 1, function(x) quantile(x, probs=0.95))
  )
}

min_dalex<-make_dalex_single_pred_df(ads, min_records)
max_dalex<-make_dalex_single_pred_df(ads, max_records)

class(min_dalex[[1]])
m<-models$glmnet

importance <- caret::varImp(best_model, scale=FALSE)

