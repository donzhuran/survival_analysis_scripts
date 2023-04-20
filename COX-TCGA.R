library("survival")
library("survminer")

setwd("C:/Users/shiva/Desktop/Survival-workshop")
TCGAMergedata <- read.csv(file ="TCGA-3.csv", sep = "\t",check.names = F,row.names = 1)

dim(TCGAMergedata)
View(TCGAMergedata)

getwd()

table(TCGAMergedata$vital_status)
table(TCGAMergedata$`Overall survival`)
table(TCGAMergedata$race)

##Clinical stage HR analysis
res.cox1 <- coxph(Surv(dataStage$`Overall survival`, dataStage$vital_status) ~ dataStage$ajcc_pathologic_stage, data = dataStage)
res.cox1
summary(res.cox1)
a <- summary(res.cox1)
b <- coef(a)
print (paste("Group=",rownames(b)))
print (paste("HR=",b[2]))
##the patients with late stage of BRCA are at high risk of dying

##Gene expression
dataGene=mutate(TCGAMergedata, `BRCA1|672`=ifelse((TCGAMergedata$`BRCA1|672`>median(TCGAMergedata$`BRCA1|672`)),"High","Low"))
res.cox2 <- coxph(Surv(dataGene$`Overall survival`, dataGene$vital_status) ~ dataGene$`BRCA1|672`, data = dataGene)
res.cox2
summary(res.cox2)
a1 <- summary(res.cox2)
b1 <- coef(a1)
print (paste("Group=",rownames(b1)))
print (paste("HR=",b1[2]))
###the patients with low expression of BRCA1 gene are at low risk of dying, not making sense
surv_obj <- Surv(dataGene$`Overall survival`, dataGene$vital_status)
fit4 <- survfit(surv_obj~dataGene$`BRCA1|672`)
ggsurvplot(fit4,data=dataGene,pval = T)

##########Multivariate COX regression, clinical stage+gene expression
dataMV=mutate(TCGAMergedata, ajcc_pathologic_stage=ifelse(((ajcc_pathologic_stage=="Stage I")|(ajcc_pathologic_stage=="Stage IA")|(ajcc_pathologic_stage=="Stage IB")|(ajcc_pathologic_stage=="Stage II")|(ajcc_pathologic_stage=="Stage IIA")|(ajcc_pathologic_stage=="Stage IIA")|(ajcc_pathologic_stage=="Stage IIB")),"Early-Stage","Late-Stage"))
dataMV=mutate(dataMV, `BRCA1|672`=ifelse((dataMV$`BRCA1|672`>median(dataMV$`BRCA1|672`)),"High","Low"))

res.cox2 <- coxph(Surv(dataMV$`Overall survival`, dataMV$vital_status) ~ dataMV$ajcc_pathologic_stage + dataMV$`BRCA1|672`, data =  dataMV)
summary(res.cox2)
ggforest(res.cox2,data=dataMV)


--------------------
###### Nomograms for prognostic biomarker discovery - combined HR with Prognostic Index

## PI = b1g1 + b2g2

# Calculate HR and b for gene 1
BRAF <- coxph(Surv(TCGAMergedata$`Overall survival`,TCGAMergedata$vital_status)~TCGAMergedata$`BRAF|673`>median(TCGAMergedata$`BRAF|673`))
#BRCA1 <- coxph(Surv(dataGene$`Overall survival`, dataGene$vital_status) ~ dataGene$`BRCA1|672`, data = dataGene)
A <- summary(BRAF)
b1 = coef(A)[2]
b1

# Calculate HR and b for gene 2
BRCA1 <- coxph(Surv(TCGAMergedata$`Overall survival`,TCGAMergedata$vital_status)~TCGAMergedata$`BRCA1|672`>median(TCGAMergedata$`BRCA1|672`))
A1 <- summary(BRCA1)
b2 = coef(A1)[2]
b2
TCGAMergedata2 <- TCGAMergedata
TCGAMergedata2$PI = (b1*TCGAMergedata$`BRAF|673`+b2*TCGAMergedata$`BRCA1|672`)

# Calculate PI (additive effect)
PI <- coxph(Surv(TCGAMergedata2$`Overall survival`,TCGAMergedata2$vital_status)~TCGAMergedata2$PI>median(TCGAMergedata2$PI))
coef(summary(PI))[2]
summary(PI)
median(TCGAMergedata2$PI)
##Patients who have a PI value more than the median of PI (722.3185) are at 1.5 times more risk of dying as compared to patients who have PI value less than 722.3185

