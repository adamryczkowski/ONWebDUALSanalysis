get_coefs<-function(m, dv) {
  if(m$method=='glmnet') {
    ans<-list()
    a<-coef(m$finalModel, m$bestTune$lambda)
    ans$coefs<-setNames(attr(a, 'x'), attr(a, 'Dimnames')[[1]][1+attr(a, 'i')])
    ans$rmse<-caret::RMSE(predict(m), dv)
    ans$train_rmse<-caret::getTrainPerf(m)$TrainRMSE
    ans$r2<-caret::R2(predict(m), dv)
    ans$train_rsq<-caret::getTrainPerf(m)$TrainRsquared
    ans$rmse<-caret::MAE(predict(m), dv)
    ans$train_mae<-caret::getTrainPerf(m)$TrainMAE
    ans$cputime<-as.numeric(m$times$everything['user.self']) + as.numeric(m$times$everything['user.child'])
    pred<-predict(m)
    ans
  }
}

plot_coefs<-function(m, adf) {
  ans<-get_coefs(m = m, dv = adf$dv)
  glmnet_imps<-gather_detailed_variable_importances(model = m, adf = adf)

  # glmnet_baseline<-glmnet_imps$raw_importances[[which(glmnet_imps$variable=='_baseline_')]]
  # glmnet_fullmodel<-glmnet_imps$raw_importances[[which(glmnet_imps$variable=='_full_model_')]]
  # glmnet_imps$importances<-(glmnet_imps$raw_importances-glmnet_fullmodel)/(glmnet_baseline-glmnet_fullmodel)
  # glmnet_imps$md_importances<-(glmnet_imps$md_raw_importances-glmnet_fullmodel)/(glmnet_baseline-glmnet_fullmodel)
  # glmnet_imps$low_importances<-(glmnet_imps$low_raw_importances-glmnet_fullmodel)/(glmnet_baseline-glmnet_fullmodel)
  # glmnet_imps$high_importances<-(glmnet_imps$high_raw_importances-glmnet_fullmodel)/(glmnet_baseline-glmnet_fullmodel)
  # glmnet_imps$sd_importances<-glmnet_imps$sd_raw_importances/(glmnet_baseline-glmnet_fullmodel)

  coefs_db<-dplyr::arrange(dplyr::left_join(tibble(variable=names(ans$coefs), coef=as.numeric(ans$coefs)), glmnet_imps, by='variable'), -importances)
  coefs_db<-dplyr::left_join(coefs_db, tibble(variable=colnames(adf), varlabel=Hmisc::label(adf)), by=c('variable', 'varlabel'))
  coefs_db$varlabel <- factor(coefs_db$varlabel, levels=coefs_db$varlabel[order(coefs_db$importances)])

  plotdf<-coefs_db %>% filter(variable !='(Intercept)')
  browser()

  glmnet_plot<-ggplot2::ggplot(data = plotdf, aes(x=varlabel, y=importances)) +
    ggplot2::geom_bar(stat="identity") + coord_flip() +
    geom_text(aes(label=scales::comma_format(accuracy=0.001)(plotdf$coef),
                  hjust=ifelse(importances>max(plotdf$importances, plotdf$high_importances)*0.3,1.2,ifelse(importances>0, -0.2, -1)),
                  colour = ifelse(importances>max(plotdf$importances, plotdf$high_importances)*0.3,"black","white")
    )
    ) +
    guides(colour=FALSE)+
    scale_colour_manual(values=c("#FFFFFF", "#000000")) +
#    scale_colour_manual(values=c("#000000", "#000000")) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    ylab(paste0("Regression coefficients sorted by importances\n",
                "for model LASSO (glmnet)\n",
                Hmisc::label(adf$dv))) +
    xlab(NULL) +
    geom_errorbar(aes(x=varlabel, ymin=low_importances, ymax=high_importances), colour='dark grey')

  return(list(plot=glmnet_plot, df=plotdf, stats=ans))
}
