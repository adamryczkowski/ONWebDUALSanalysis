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
