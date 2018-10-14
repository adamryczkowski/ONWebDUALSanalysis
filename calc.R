#devtools::install_github('adamryczkowski/ONWebDUALSanalysis')
library(ONWebDUALSanalysis)
library(doMC)
registerDoMC(8)
dv_nr<-1
#debugonce(make_rap)
make_rap(dv_nr = 1, rap_path = 'reports')

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
dubious_models<-setdiff(all_models, c(time_consuming_models, really_long, rather_long,  models_that_hangs,
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
ans<-calc_models(dubious_models, dv_nr=1, adaptive = NA, assume_calculated = TRUE)
ans<-calc_models(dubious_models, dv_nr=1, adaptive = NA, assume_calculated = FALSE)

models<-ans$models
ads<-ans$ads

df<-model_perfs(ans)
models<-models[df$model]

packages<-unlist(purrr::map(models, function(m) m$modelInfo$library ))
new.packages <- packages[!(packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)


m<-models[[df$model[[1]]]]
m<-models$glmnet
descr_model<-function(m) {
  tuned_pars<-setNames(m$finalModel$tuneValue, as.character(m$modelInfo$parameters$label))
  return(tuned_pars)
}

get_coefs<-function(m) {
  if(m$method=='glmnet') {
    ans<-list()
    a<-coef(m$finalModel, m$bestTune$lambda)
    ans$coefs<-setNames(attr(a, 'x'), attr(a, 'Dimnames')[[1]][1+attr(a, 'i')])
    ans$rmse<-caret::RMSE(predict(m), ads$dv)
    ans$train_rmse<-caret::getTrainPerf(m)$TrainRMSE
    ans$r2<-caret::R2(predict(m), ads$dv)
    ans$train_rsq<-caret::getTrainPerf(m)$TrainRsquared
    ans$rmse<-caret::MAE(predict(m), ads$dv)
    ans$train_mae<-caret::getTrainPerf(m)$TrainMAE
    ans$cputime<-as.numeric(m$times$everything['user.self']) + as.numeric(m$times$everything['user.child'])
    pred<-predict(m)
    ans
  }
}

var_imp<-function(model, ads) {
  cat(paste0(model$method, '\n'))
  explainer<-DALEX::explain(model, data=ads %>% select(-dv), y=ads$dv)
  dfexp<-dplyr::select(DALEX::variable_importance(explainer, loss_function = DALEX::loss_root_mean_square, type = "difference"),-label)
  names(dfexp)<-c('variable', model$method)
  dfexp
}







comp_models<-function(models) {
  r<-caret::resamples(models)
  bwplot(r,metric="RMSE",main="GBM vs xgboost")
}

summary(res, metric='RMSE', decreasing = FALSE)
bwplot(res, metric='RMSE', decreasing = TRUE)
dotplot(res)

best_model_name<-df$model[[1]]
best_model<-models[[ best_model_name ]]
glmnet_model<-models[['glmnet']]

ads_pred<-cbind(ads, dv_predict=predict(best_model))

#ggplot2::qplot(dv, dv_predict, data = ads)







# #Explanations
#
# min_point<-c(min(ads$dv), min(ads_pred$dv_predict))
# max_point<-c(max(ads$dv), max(ads_pred$dv_predict))
#
# min_dist<-((ads$dv - min_point[[1]])^2 + (ads_pred$dv_predict - min_point[[2]])^2)^0.5
# max_dist<-((ads$dv - max_point[[1]])^2 + (ads_pred$dv_predict - max_point[[2]])^2)^0.5
#
# ktory_quantil=0.05
#
# dvquantiles<-c(quantile(min_dist, probs = ktory_quantil),quantile(max_dist, probs = ktory_quantil))
#
# min_records<-which(min_dist<=dvquantiles[[1]])
# max_records<-which(max_dist<=dvquantiles[[2]])
#
#
# explainer<-DALEX::explain(best_model, data=ads, label=best_model_name)
#
# library(doParallel)
# registerDoParallel(cores=4)
#
# make_dalex_single_pred_df<-function(ads, records) {
#   rec_dalex<-foreach(rec_nr=records) %dopar%
#     DALEX::single_prediction(explainer, observation = ads[rec_nr,])
#
#   a<-unlist(purrr::map(rec_dalex, function(x) x$contribution))
#   dim(a)<-c(nrow(rec_dalex[[1]]), length(rec_dalex))
#   cum_min_dalex<-tibble(name=rec_dalex[[1]]$variable,
#                         mean_contr = plyr::aaply(a, 1, mean),
#                         sd_contr = plyr::aaply(a, 1, sd),
#                         min_contr = plyr::aaply(a, 1, min),
#                         max_contr = plyr::aaply(a, 1, max),
#                         med_contr = plyr::aaply(a, 1, median),
#                         q05_contr = plyr::aaply(a, 1, function(x) quantile(x, probs=0.05)),
#                         q95_contr = plyr::aaply(a, 1, function(x) quantile(x, probs=0.95))
#   )
# }
#
# min_dalex<-make_dalex_single_pred_df(ads, min_records)
# max_dalex<-make_dalex_single_pred_df(ads, max_records)
#
# class(min_dalex[[1]])
# m<-models$glmnet
#
# importance <- caret::varImp(best_model, scale=FALSE)
#
