
find_roc_range=function(data,k,seed,data_name){
require(caret)
require(pROC)
require(dplyr)
require(data.table)
set.seed(seed)
#data.names.list=c('Metabolomics 230','Metabolimcs_elastic_net','Covariates','Covariates and Metabolimcs_elastic_net')
#data.names.list_id=0
splits <- createFolds(data[,ncol(data)], returnTrain = TRUE,k = k)
results <- lapply(splits, 
                  function(x, dat) {
                    holdout <- (1:nrow(dat))[-unique(x)]
                    data.frame(index = holdout, 
                               obs = dat$lable[holdout])
                  },
                  dat = data)
mods <- vector(mode = "list", length = length(splits))
roc <- vector(mode = "list", length = length(splits))
## foreach or lapply would do this faster
for(i in seq(along = splits)) {
  in_train <- unique(splits[[i]])
  set.seed(seed)
  mod <- train(lable ~ ., data = data[in_train, ],
               method = "regLogistic", #cost (Cost) loss (Loss Function)  epsilon (Tolerance)
               preProc = c("center", "scale"),trControl=trainControl(method="LOOCV", 
                     summaryFunction=twoClassSummary, 
                     classProbs=T,
                     savePredictions = T,allowParallel = TRUE)
               )
  results[[i]]$pred <- predict(mod, data[-in_train,!(colnames(data) %in% c('lable')) ],type="prob")
  roc[[i]]<- roc(predictor=results[[i]]$pred$Normal,response=data[-in_train,'lable' ],levels=rev(levels(data[-in_train,'lable' ])))#,smooth=TRUE)
  mods[[i]] <- mod
  
}

###############How to plot concatante variables importance from 5 models
important_featurs=(lapply(mods,function(xx) (varImp(xx,scale=T))))
important_featurs_list=(lapply(important_featurs,function(x) data.frame(x$importance)))
important_featurs_list2=data.frame((lapply(important_featurs_list,function(x) (x[1]))))
#print(important_featurs_list2)
important_featurs2=data.frame(apply(important_featurs_list2,1,mean ))
#print(important_featurs2)

colnames(important_featurs2)='importance'
zz=important_featurs2[order(-important_featurs2$importance),,drop=F]

zz=(zz[apply(zz, 1, function(row) all(row !=0 )),,drop=F])
zz=(zz[apply(zz, 1, function(row) all(row >0 )),,drop=F])

#pdf("metabolites_230_importance.pdf")
#pdf("covariates_importance.pdf")
#pdf("metabolites_lasso_importance.pdf")
#pdf("covariates_and_lasso_importance.pdf")
p=ggplot(data=zz,aes(x=reorder(rownames(zz),importance), y=(importance))) +
geom_bar(stat="identity",fill="steelblue")+
xlab("")+theme(axis.text=element_text(size=15),axis.title=element_text(size=10,face="bold"))
p1=p+coord_flip()
print(p1)
#dev.off()

predicted_57=(rbindlist(lapply(results,function(x) data.frame(index=x$index,obs=x$obs,pre=x$pred))))
#pdf("metabolites_230_ROC.pdf")
#pdf("covariates_ROC.pdf")
#pdf("metabolites_lasso_ROC.pdf")
#pdf("covariates_and_lasso_ROC.pdf")

plot(roc(predictor=predicted_57$pre.Normal,response=predicted_57$obs,
levels=rev(levels(predicted_57$obs)),smooth=TRUE,ci=TRUE),print.auc = TRUE,main = paste0(data_name),col='blue',cex.lab=2, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
#dev.off()

#print(predicted_57)  
#print(results)
#print(roc)
#lapply(roc, function(x)plot(x,print.auc = TRUE))
#cat(sprintf("The average AUC of %d fold changes is %f ", k, mean(as.numeric(lapply(roc,function(x) (x$auc))))))
clab = 0.1
cmain = 0.1
caxis = 0.1
#lapply(mods,function(x) plot(varImp(x,scale=T),cex.lab=clab,cex.main =cmain,cex.axis=caxis))
#lapply(mods,function(x) plot(varImp(x,scale=T),20,cex.lab=clab,cex.main =cmain,cex.axis=caxis))
#:lapply(mods,function(x) plot(x))	
    }
   