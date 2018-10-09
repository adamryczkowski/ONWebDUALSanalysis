#devtools::install_github('adamryczkowski/ONWebDUALSimport')
library(ONWebDUALSanalysis)
library(data.table)
library(caret)
#contr.ltfr<-caret::contr.ltfr
joined_df<-readRDS(system.file('db_no_duplicates.rds', package='ONWebDUALSanalysis'))
ans<-select_variables_sep2018(joined_df)
dt<-ans$db
iv_names<-ans$iv_names
dv_names<-ans$dv_names
dv_name<-dv_names[[1]]
library(doMC)
registerDoMC(4)


time_consuming_models<-c('ANFIS', 'DENFIS', 'FIR.DM', 'FS.HGD', 'GFS.FR.MOGUL', 'GFS.LT.RS', 'Rborist', 'xgbDART', 'xgbLinear', 'xgbTree')
#empty_models<-c('avNNet', 'ANFIS')
really_long<-c('DENFIS', 'FIR.DM', 'FS.HGD')

empty_models<-character(0)
tensor_flow<-c('mlpKerasDecay', 'mlpKerasDropout')
not_parallel<-c('M5', 'M5Rules')
mem_insufficient<-c('randomGLM')
model_blacklist<-c('bag', 'bam', 'bartMachine', 'blackboost','bstSm', 'bstTree', 'elm', 'extraTrees',
                   'gam', 'gamboost', 'gamLoess', 'gamSpline', 'gbm_h2o', 'GFS.THRIFT',
                   'glmboost', 'glmnet_h2o', 'HYFIS', 'krlsRadial',
                   'logicBag', 'logreg', 'mlpSGD', 'mxnet', 'mxnetAdam',
                   'neuralnet', 'rlm', 'svmBoundrangeString', 'svmExpoString', 'svmSpectrumString',
                   #   'xgbLinear',
                   time_consuming_models, empty_models, not_parallel, tensor_flow, mem_insufficient)

all_models<-unique((caret::modelLookup() %>% filter(forReg==TRUE & ! model %in% model_blacklist))$model)


model_names=time_consuming_models
model_names=c(all_models, c('mlpKerasDropout','M5', 'M5Rules'))
model_names=all_models
model_names=not_parallel
model_names=c("earth", "enet", "blasso", "BstLm", "glmnet", "RRF", "ctree2", "ranger", "RRFglobal", "rf", "cforest", "evtree", "gbm", "bagEarth", "spikeslab", "ctree", "nodeHarvest", "lars", "lasso", "rpart2", "parRF", "cubist", "lars2", "rqlasso", "bagEarthGCV", "rqnc", "rpart", "gcvEarth", "glmStepAIC", "penalized", "msaenet", "rpart1SE", "rfRules", "qrf", "relaxo")

do_model_inner<-function(dv_name, model_name, ads, tc) {
  plik_name<-paste0('model_', dv_name, '_', model_name, '.rds')
  if(file.exists(plik_name)) {
    cat(paste0("Reading in already calculated model ", model_name, "...\n"))
    readRDS(plik_name)
  } else {
    cat(paste0("Calculating model ", model_name, "...\n"))
    tryCatch(
      {
        model<-caret::train(dv ~ ., data = ads, method = model_name, trControl = tc, tuneLength=15)
#        model<-caret::rfe(x = ads, y = ads$dv, method = model_name, rfeControl = tc, sizes=gen_int_geom_series(52, 20))
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

tc <- caret::trainControl(index = cvIndex,
                          method = 'adaptive_cv',
                          number = 10, repeats = 10,
                          adaptive = list(min = 5, alpha = 0.05,
                                          method = "gls", complete = TRUE),
                          search = "random"
                          )


#caret::checkInstall(model_names)
#caret::checkInstall('mlp')
#do_model_inner(dv_name, 'mlpKerasDropout', ads, tc)
#do_model_inner(dv_name, 'xgbTree', ads, tc)
m<-do_model_inner(dv_name, 'glmnet', ads, tc)
models<-
  purrr::map(model_names, do_model_inner, ads=ads, tc=tc, dv_name=dv_name)
names(models)<-model_names
a<-purrr::map_dbl(models, function(x) {caret::getTrainPerf(x)$TrainRMSE} )
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



df<-tibble(model=model_names, rmsa=a, elapsed=b1, user=b2+b3,
           is_nn=b4_1, is_bagging=b4_2, is_boost=b4_7,
           is_rf=b4_3, is_lm=b4_4, is_bayes=b4_5, is_feature_sel=b4_6,
           ) %>% arrange(rmsa)

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

