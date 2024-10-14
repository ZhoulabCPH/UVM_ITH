########### Fig2a
marker = read.csv('DEGs.csv')
ggplot(data=marker,aes(x=avg_log2FC ,
                          y= abs(pct.1-pct.2)))+
  geom_point(alpha=0.8,size=3)


########### Fig2c
library(vegan)
load('TCGA-UVM2.Rdata')
gene = filter(marker,avg_log2FC > 0.58 & p_val_adj < 0.05 & pct.1 > 0.6)$Genes 
gene2 = filter(marker,avg_log2FC < (-0.58) & p_val_adj < 0.05& pct.2 > 0.6)$Genes

t_test = data.frame(test = rep(10,length(colnames(TCGA_exp))),row.names = colnames(TCGA_exp))
t_test$test = as.numeric(t_test$test)

for (i in colnames(TCGA_exp)) {
  a = t.test(TCGA_exp[gene,i],TCGA_exp[gene2,i],alternative = c("two.sided"))
  t_test[i,1] = a$statistic
}
t_test = arrange(t_test,test)
t_test$test = decostand(t_test$test,method = "standardize",MARGIN = 2)
TCGA_meta = TCGA_meta[rownames(t_test),]
TCGA_meta$t_test = t_test$test

bk <- c(seq(-2,2,by=0.01))
color = c(colorRampPalette(colors = c("#0000FF","white"))(length(bk)/2),colorRampPalette(colors = c("white","#FF0000"))(length(bk)/2))
pheatmap::pheatmap(TCGA_exp[c(gene,gene2),rownames(t_test)],border = F,
                   show_rownames = F,show_colnames = F,cluster_rows = F,cluster_cols = F,
                   scale = "row",breaks = bk,color = color,
                   annotation_col = TCGA_meta[,c("t_test","stage","diagnosis","age",
                                                 "meta_all.tumor_basal_diameter",
                                                 "meta_all.gender.demographic")])

########### Fig2d
library(survival)
library(survminer)
ggsurvplot(survfit(Surv(TCGA_sur$OS.time/12,event = TCGA_sur$OS  )~class ,
                   data= TCGA_sur)  ,pval = T,
           conf.int = F,
           conf.int.style="step",
           conf.int.alpha =0.3,
           pval.method=TRUE,
           risk.table =T ,
           risk.table.col = "class",risk.table.y.text.col = T)



########### Fig2g
load('exGSE22138.Rdata')
t_test = data.frame(test = rep(10,length(colnames(exGSE22138))),row.names = colnames(exGSE22138))
t_test$test = as.numeric(t_test$test)

gene = gene[gene %in% rownames(exGSE22138)]
gene2 = gene2[gene2 %in% rownames(exGSE22138)]
for (i in colnames(exGSE22138)) {
  a = t.test(exGSE22138[gene,i],exGSE22138[gene2,i],alternative = c("two.sided"))
  t_test[i,1] = a$statistic
}
t_test = arrange(t_test,test)

bk <- c(seq(-2,2,by=0.01))
color = c(colorRampPalette(colors = c("#0000FF","white"))(length(bk)/2),colorRampPalette(colors = c("white","#FF0000"))(length(bk)/2))
pheatmap::pheatmap(exGSE22138[c(gene,gene2),rownames(t_test)],border = F,
                   show_rownames = F,show_colnames = F,cluster_rows = F,cluster_cols = F,
                   scale = "row",breaks = bk,color = color)


########### Fig2f
library(survival)
library(survminer)
ggsurvplot(survfit(Surv(meta_exGSE22138$os.time,event = meta_exGSE22138$os  )~class ,
                   data= meta_exGSE22138)  ,pval = T,
           conf.int = F,
           conf.int.style="step",
           conf.int.alpha =0.3,
           pval.method=TRUE,
           risk.table =T ,
           risk.table.col = "class",risk.table.y.text.col = T)


########### Fig2f
cox_data = filter(TCGA_meta,meta_all.tumor_basal_diameter != "NA", diagnosis != "Malignant melanoma, NOS" & meta_all.tumor_stage.diagnoses != "not reported")

cox_data$sample = rownames(cox_data)
TCGA_sur$sample = rownames(TCGA_sur)
cox_data = inner_join(TCGA_sur,cox_data)
rownames(cox_data) = cox_data$sample

cox_data[cox_data$diagnosis == "Spindle","diagnosis"] = "1_Spindle"
cox_data[cox_data$diagnosis == "Mixed","diagnosis"] = "2_Mixed"
cox_data[cox_data$diagnosis == "Epithelioid","diagnosis"] = "3_Epithelioid"

cox_data$class = as.factor(cox_data$class)
covariates = names(cox_data)[c(6,10,11,12,17,19,20)]

univ_formulas  <- sapply(covariates,                
                         function(x) as.formula(paste('Surv(OS.time2, OS)~', x)))
univ_models <- lapply( univ_formulas, function(x){coxph(x, data = cox_data)}) 

res = list()
for (i in noquote(names(univ_models))) {
  out_multi <- cbind(
    coef=summary(univ_models[[i]])$coefficients[,"coef"],
    HR= summary(univ_models[[i]])$conf.int[,"exp(coef)"],
    HR.95L=summary(univ_models[[i]])$conf.int[,"lower .95"],
    HR.95H=summary(univ_models[[i]])$conf.int[,"upper .95"],
    pvalue=summary(univ_models[[i]])$coefficients[,"Pr(>|z|)"])
  res[[i]] = out_multi
}

rownames(res$meta_all.gender.demographic) = "gender"
rownames(res$meta_all.age_at_diagnosis.diagnoses) = "Age"
rownames(res$meta_all.tumor_basal_diameter) = "tumor_basal_diameter"
res = rbind(res$meta_all.age_at_diagnosis.diagnoses,res$meta_all.tumor_basal_diameter,res$meta_all.gender.demographic,res$diagnosis,res$stage,res$class)

cox.mod = coxph(Surv(OS.time,OS) ~meta_all.age_at_diagnosis.diagnoses+meta_all.tumor_basal_diameter+meta_all.gender.demographic+
                  stage+diagnosis+class,data = cox_data)
multiCoxSum = summary(cox.mod)

res2 <- cbind(as.data.frame(signif(multiCoxSum$coefficients[,"coef"], digits=2)),
              as.data.frame(signif(multiCoxSum$conf.int[,"exp(coef)"], digits=2)),
              as.data.frame(signif(multiCoxSum$conf.int[,"lower .95"], digits=2)),
              as.data.frame(signif(multiCoxSum$conf.int[,"upper .95"], digits=2)),
              as.data.frame(signif(multiCoxSum$coefficients[,"Pr(>|z|)"], digits=2)))
names(res2) = c("coef","HR","HR.95L","HR.95H","pvalue")



cox_DATA = meta_exGSE22138
cox_DATA = filter(cox_DATA,tumor_dia != "NA")

cox_DATA$tumor_dia = as.numeric(cox_DATA$tumor_dia)
cox_DATA$age = as.numeric(cox_DATA$age)
cox_DATA$tumor_thickness = as.numeric(cox_DATA$tumor_thickness)
cox_DATA$class = as.factor(cox_DATA$class)

covariates = names(cox_DATA)[c(5,1,2,6,7,8,9)]

univ_formulas  <- sapply(covariates,                   # 单COX
                         function(x) as.formula(paste('Surv(os.time, os)~', x)))
univ_models <- lapply( univ_formulas, function(x){coxph(x, data = cox_DATA)})  

res = list()
for (i in noquote(names(univ_models))) {
  out_multi <- cbind(
    coef=summary(univ_models[[i]])$coefficients[,"coef"],
    HR= summary(univ_models[[i]])$conf.int[,"exp(coef)"],
    HR.95L=summary(univ_models[[i]])$conf.int[,"lower .95"],
    HR.95H=summary(univ_models[[i]])$conf.int[,"upper .95"],
    pvalue=summary(univ_models[[i]])$coefficients[,"Pr(>|z|)"])
  res[[i]] = out_multi
}
rownames(res$metastasis) = "metastasis"
rownames(res$tumor_dia) = "tumor_basal_diameter"
rownames(res$age) = "age"
rownames(res$gender) = c("gender")
rownames(res$tumor_thickness) = "tumor_thickness"

res = rbind(res$class,res$metastasis,res$tumor_dia,res$age,res$gender,res$tumor_thickness)

res[,1:4] = round(res[,1:4],2)
res[,5]= round(res[,5],3)

res = data.frame(res)
rownames(res) = c("class","metastasis","tumor_basal_diameter","gender","tumor_thickness")

coxph(Surv(os.time,os) ~class+tumor_dia +age +gender + tumor_thickness,data = cox_DATA)

cox.mod = coxph(Surv(os.time,os) ~class+tumor_dia +age +gender + tumor_thickness,data = cox_DATA)
multiCoxSum = summary(cox.mod)
res3 <- cbind(as.data.frame(signif(multiCoxSum$coefficients[,"coef"], digits=2)),
              as.data.frame(signif(multiCoxSum$conf.int[,"exp(coef)"], digits=2)),
              as.data.frame(signif(multiCoxSum$conf.int[,"lower .95"], digits=2)),
              as.data.frame(signif(multiCoxSum$conf.int[,"upper .95"], digits=2)),
              as.data.frame(signif(multiCoxSum$coefficients[,"Pr(>|z|)"], digits=2)))
names(res3) = c("coef","HR","HR.95L","HR.95H","pvalue")


cox_DATA = meta_GSE84976
cox_DATA$class = as.factor(cox_DATA$class)
cox_DATA$age = as.numeric(cox_DATA$age)

covariates = names(cox_DATA)[c(5,6)]

univ_formulas  <- sapply(covariates,                   # 单COX
                         function(x) as.formula(paste('Surv(os.time, os)~', x)))
univ_models <- lapply( univ_formulas, function(x){coxph(x, data = cox_DATA)})  
res = list()
for (i in noquote(names(univ_models))) {
  out_multi <- cbind(
    coef=summary(univ_models[[i]])$coefficients[,"coef"],
    HR= summary(univ_models[[i]])$conf.int[,"exp(coef)"],
    HR.95L=summary(univ_models[[i]])$conf.int[,"lower .95"],
    HR.95H=summary(univ_models[[i]])$conf.int[,"upper .95"],
    pvalue=summary(univ_models[[i]])$coefficients[,"Pr(>|z|)"])
  res[[i]] = out_multi
}

rownames(res$age) = "age"
rownames(res$class) = c("class")

res = rbind(res$class,res$age)

res[,1:4] = round(res[,1:4],2)
res[,5]= round(res[,5],3)
res = data.frame(res)

res_list[[5]] = res

cox.mod = coxph(Surv(os.time,os) ~class  +age ,data = cox_DATA)
multiCoxSum = summary(cox.mod)
res4 <- cbind(as.data.frame(signif(multiCoxSum$coefficients[,"coef"], digits=2)),
              as.data.frame(signif(multiCoxSum$conf.int[,"exp(coef)"], digits=2)),
              as.data.frame(signif(multiCoxSum$conf.int[,"lower .95"], digits=2)),
              as.data.frame(signif(multiCoxSum$conf.int[,"upper .95"], digits=2)),
              as.data.frame(signif(multiCoxSum$coefficients[,"Pr(>|z|)"], digits=2)))
names(res4) = c("coef","HR","HR.95L","HR.95H","pvalue")



rbind(res2,res3,res4) %>% write.csv('E:/FOREST2.csv')
a = read.csv('E:/FOREST2.csv')
a$Subgroups <- ifelse(is.na(a$HR),
                      a$Subgroup,
                      paste0("      ", a$Subgroup))

a$` ` <- paste(rep(" ", 24), collapse = " ")
a$`HR (95% CI)` <- ifelse(is.na(a$HR), "",
                          sprintf("%.2f (%.2f to %.2f)",
                                  a$HR, a$HR.95L, a$HR.95H))

a$`HR (95% CI)`
forest(a[,c(1,8,6)],
       est = a$HR,
       lower = a$HR.95L,
       upper = a$HR.95H,
       ref_line  = 1,
       xlim = c(0, 30),
       ci_column = 2)












