
########### Fig5a
filtet_data = TCGA_exp[c(gene,gene2),]
library(caret)
highCorr = findCorrelation(cor(data.frame(t(filtet_data))), 0.8) 
filtet_data2 = filtet_data[filter_gene[-highCorr],]
filtet_data2 = data.frame(t(filtet_data2))
findLinearCombos(filtet_data2)  
filtet_data2_scale = predict(preProcess(filtet_data2), filtet_data2) #标准化训练集
inTrain = createDataPartition(TCGA_sur2$class, p = 2/3, list = FALSE)

trainx = filtet_data2_scale[inTrain,]
testx = filtet_data2_scale[-inTrain,]
trainy = as.character(TCGA_sur$class)[inTrain]
testy = as.character(TCGA_sur$class)[-inTrain]

ctrl= rfeControl(functions = rfFuncs, method = "cv",verbose = FALSE, returnResamp = "final")

Profile = rfe(filtet_data2_scale, TCGA_sur2$class, sizes = 1:40, rfeControl = ctrl) 

gene_filter = gene[gene %in% Profile$optVariables]
gene2_filter = gene2[gene2 %in% Profile$optVariables]

filtet_data3 = filtet_data2[,Profile$optVariables]

require(xgboost)
require(Matrix)
require(data.table)

filtet_data3 = filtet_data3[rownames(TCGA_sur),]
filtet_data3$class = TCGA_sur$class

sparse_matrix <- sparse.model.matrix(class~ ., data = filtet_data3)[,-1]
filtet_data3$class = as.character(filtet_data3$class)
filtet_data3[filtet_data3$class == "L-ITH","class"] = 0
filtet_data3[filtet_data3$class == "H-ITH","class"] = 1
filtet_data3$class = as.numeric(filtet_data3$class)
bst <- xgboost(data = sparse_matrix, label = filtet_data3$class, max_depth = 10,eta = 0.3,  nrounds =1000)

importance <- xgb.importance(feature_names = colnames(sparse_matrix), model = bst)
importance = data.frame(importance)
xgb.plot.importance(importance, rel_to_first = TRUE, xlab = "Relative importance")


########### Fig5b
validation  = exGSE22138[importance$Feature,]%>% t()%>% data.frame()
validation = validation[rownames(meta_exGSE22138),]
validation$group = meta_exGSE22138$class

validation[1:9] = decostand(validation[,1:9] , 'standardize') 

validation_matrix <- sparse.model.matrix(group ~., data = validation)
validation_label <-  as.numeric(validation$group)

validation_fin <- list(data=validation_matrix,label=validation_label) 
dvalidation <- xgb.DMatrix(data = validation_fin$data, label = validation_fin$label)


pre_xgb = round(predict(xgb,newdata = dvalidation))
meta_exGSE22138$pre = pre_xgb
table(meta_exGSE22138$pre)

ggsurvplot(survfit(Surv(meta_exGSE22138$os.time,event = meta_exGSE22138$os  )~pre ,
                   data= meta_exGSE22138)  ,pval = T)


validation  = GSE84976[importance$Feature,]%>% t()%>% data.frame()
validation[,1:9] = decostand(validation[,1:9] , 'standardize') 

validation_matrix <- sparse.model.matrix(group ~., data = validation)
validation_label <-  as.numeric(validation$group)

validation_fin <- list(data=validation_matrix,label=validation_label) 
dvalidation <- xgb.DMatrix(data = validation_fin$data, label = validation_fin$label)

pre_xgb = round(predict(xgb,newdata = dvalidation))
meta_GSE84976$pre = pre_xgb
table(meta_GSE84976$pre)

ggsurvplot(survfit(Surv(meta_GSE84976$os.time,event = meta_GSE84976$os  )~pre ,
                   data= meta_GSE84976)  ,pval = T)


validation  = GSE27831_log[importance$Feature,]%>% t()%>% data.frame()
validation[,1:9] = decostand(validation[,1:9] , 'standardize') 

validation_matrix <- sparse.model.matrix(group ~., data = validation)
validation_label <-  as.numeric(validation$group)

validation_fin <- list(data=validation_matrix,label=validation_label) 
dvalidation <- xgb.DMatrix(data = validation_fin$data, label = validation_fin$label)

pre_xgb = round(predict(xgb,newdata = dvalidation))
meta_GSE27831$pre = pre_xgb
table(meta_GSE27831$pre)

ggsurvplot(survfit(Surv(meta_GSE27831$DFS_months,event = meta_GSE27831$metastasis  )~pre ,
                   data= meta_GSE27831)  ,pval = T,palette=c("#00468B","#ED0000"))








































