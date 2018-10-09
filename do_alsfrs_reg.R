#devtools::install_github('adamryczkowski/ONWebDUALSimport')
library(ONWebDUALSanalysis)
library(data.table)
contr.ltfr<-caret::contr.ltfr
joined_df<-readRDS(system.file('db_no_duplicates.rds', package='ONWebDUALSanalysis'))
ans<-select_variables_sep2018(joined_df)
dt<-ans$db
iv_names<-ans$iv_names
dv_names<-ans$dv_names
dv_name<-dv_names[[1]]
library(doMC)
registerDoMC(2)


time_consuming_models<-c('ANFIS', 'DENFIS', 'FIR.DM', 'FS.HGD', 'GFS.FR.MOGUL', 'GFS.LT.RS', 'Rborist', 'xgbDART', 'xgbLinear', 'xgbTree')
#empty_models<-c('avNNet', 'ANFIS')
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

do_model_inner<-function(dv_name, model_name, ads, tc) {
  plik_name<-paste0('model_', dv_name, '_', model_name, '.rds')
  if(file.exists(plik_name)) {
    cat(paste0("Reading in already calculated model ", model_name, "...\n"))
    readRDS(plik_name)
  } else {
    cat(paste0("Calculating model ", model_name, "...\n"))
    model<-caret::train(dv ~ ., data = ads, method = model_name, trControl = tc)
    saveRDS(model, plik_name)
    model
  }
}

ads<-make_ads(dt, iv_names, dv_name, keep_nominal ='iv56')

#  which(is.na(data.matrix(ads)),arr.ind = TRUE)

groupvar<-ads$iv56
cvIndex<-caret::createFolds(groupvar, 10, returnTrain = T)

tc <- caret::trainControl(index = cvIndex,
                          method = 'cv',
                          number = 10)


caret::checkInstall(model_names)
do_model_inner(dv_name, 'mlpKerasDropout', ads, tc)
models<-
  purrr::map(model_names, do_model_inner, ads=ads, tc=tc, dv_name=dv_name)

dal<-DALEX::explain(ans$mod, data=dt, y=depvar)
